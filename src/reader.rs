use std::io::{self, Read};

use crate::{selection::AtomSelection, BoxVec};

struct DecodeState {
    count: usize,
    lastbits: usize,
    lastbyte: u8,
}

// TODO: I have a constexpr laying around for this somewhere.
#[rustfmt::skip]
pub const MAGICINTS: [i32; 73] = [
    0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216
];
pub const FIRSTIDX: usize = 9; // Note that MAGICINTS[FIRSTIDX-1] == 0.

// TODO: Amortize the read_opaque call such that not all data is read in at once if that's wasteful
// given the atom_selection.
#[inline]
pub(crate) fn read_compressed_positions<R: Read>(
    file: &mut R,
    positions: &mut [f32],
    precision: f32,
    scratch: &mut Vec<u8>,
    atom_selection: &AtomSelection,
) -> io::Result<()> {
    let n = positions.len();
    assert_eq!(n % 3, 0, "the length of `positions` must be divisible by 3");
    let natoms = n / 3;

    let invprecision = precision.recip();

    let minint = [0; 3].try_map(|_| read_i32(file))?;
    let maxint = [0; 3].try_map(|_| read_i32(file))?;
    let mut smallidx = read_u32(file)? as usize;
    assert!(smallidx < MAGICINTS.len());

    let mut sizeint = [0u32; 3];
    let mut bitsizeint = [0u32; 3];
    let bitsize = calc_sizeint(minint, maxint, &mut sizeint, &mut bitsizeint);

    let tmpidx = smallidx - 1;
    let tmpidx = if FIRSTIDX > tmpidx { FIRSTIDX } else { tmpidx };

    let mut smaller = MAGICINTS[tmpidx] / 2;
    let mut smallnum = MAGICINTS[smallidx] / 2;
    let mut sizesmall = [MAGICINTS[smallidx] as u32; 3];

    let compressed_data = scratch;
    read_opaque(file, compressed_data)?;

    let mut state = DecodeState {
        count: 0,
        lastbits: 0,
        lastbyte: 0,
    };
    let mut run: i32 = 0;
    let mut prevcoord;
    let mut write_idx = 0;
    let mut read_idx = 0;
    while read_idx < natoms {
        let mut coord = [0i32; 3];
        let mut position: &mut [f32; 3] = positions.array_chunks_mut().nth(write_idx).unwrap();
        if bitsize == 0 {
            coord[0] = decodebits(compressed_data, &mut state, bitsizeint[0] as usize);
            coord[1] = decodebits(compressed_data, &mut state, bitsizeint[1] as usize);
            coord[2] = decodebits(compressed_data, &mut state, bitsizeint[2] as usize);
        } else {
            decodeints(compressed_data, &mut state, bitsize, sizeint, &mut coord);
        }

        coord[0] += minint[0];
        coord[1] += minint[1];
        coord[2] += minint[2];
        prevcoord = coord;

        macro_rules! write_position {
            ($position:ident, $write_idx:ident, $coord:ident  ) => {
                match atom_selection.is_included($write_idx) {
                    None => return Ok(()),
                    Some(false) => {}
                    Some(true) => {
                        *$position = $coord.map(|v| v as f32 * invprecision);
                        $write_idx += 1;
                    }
                };
            };
        }

        let flag: bool = decodebits::<u8>(compressed_data, &mut state, 1) > 0;
        let mut is_smaller = 0;
        if flag {
            run = decodebits(compressed_data, &mut state, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller -= 1;
        }
        if run > 0 {
            if write_idx * 3 + run as usize > n {
                // eprintln!("may attempt to write a run beyond the positions buffer");
                // dbg!(write_idx, run, n, write_idx * 3 + run as usize);
            }

            // Let's read the next coordinate.
            coord.fill(0);

            for k in (0..run).step_by(3) {
                decodeints(
                    compressed_data,
                    &mut state,
                    smallidx as u32,
                    sizesmall,
                    &mut coord,
                );
                read_idx += 1;
                coord[0] += prevcoord[0] - smallnum;
                coord[1] += prevcoord[1] - smallnum;
                coord[2] += prevcoord[2] - smallnum;
                if k == 0 {
                    // Swap the first and second atom. This is done to achieve better compression
                    // for water atoms. Waters are stored as OHH, but right now we want to swap the
                    // atoms such that e.g., water will become HOH again.
                    std::mem::swap(&mut coord[0], &mut prevcoord[0]);
                    std::mem::swap(&mut coord[1], &mut prevcoord[1]);
                    std::mem::swap(&mut coord[2], &mut prevcoord[2]);
                    write_position!(position, write_idx, prevcoord);
                    position = match positions.array_chunks_mut().nth(write_idx) {
                        Some(c) => c,
                        None => break,
                    };
                } else {
                    prevcoord = coord;
                }
                write_position!(position, write_idx, coord);
                position = match positions.array_chunks_mut().nth(write_idx) {
                    Some(c) => c,
                    None => break,
                };
            }
        } else {
            write_position!(position, write_idx, coord);
        }

        match is_smaller.cmp(&0) {
            std::cmp::Ordering::Less => {
                smallidx -= 1;
                smallnum = smaller;
                if smallidx > FIRSTIDX {
                    smaller = MAGICINTS[smallidx - 1] / 2;
                } else {
                    smaller = 0;
                }
            }
            std::cmp::Ordering::Greater => {
                smallidx += 1;
                smaller = smallnum;
                smallnum = MAGICINTS[smallidx] / 2;
            }
            std::cmp::Ordering::Equal => {}
        }

        assert_ne!(MAGICINTS[smallidx], 0, "found an invalid size");
        sizesmall.fill(MAGICINTS[smallidx] as u32);
        read_idx += 1;
    }

    Ok(())
}

#[inline]
pub(crate) fn read_boxvec<R: Read>(file: &mut R) -> io::Result<BoxVec> {
    let mut boxvec = [0.0; 9];
    read_f32s(file, &mut boxvec)?;
    let cols = [
        [boxvec[0], boxvec[1], boxvec[2]],
        [boxvec[3], boxvec[4], boxvec[5]],
        [boxvec[6], boxvec[7], boxvec[8]],
    ];
    Ok(BoxVec::from_cols_array_2d(&cols))
}

fn read_opaque<R: Read>(file: &mut R, data: &mut Vec<u8>) -> io::Result<()> {
    let count = read_u32(file)? as usize;
    let padding = (4 - (count % 4)) % 4; // FIXME: Why, and also, can we do this better?
    data.resize(count + padding, 0);
    file.read_exact(data)
}

pub(crate) fn read_f32s<R: Read>(file: &mut R, buf: &mut [f32]) -> io::Result<()> {
    for value in buf {
        *value = read_f32(file)?
    }
    Ok(())
}

// FIXME: These read_* functions are prime targets for a macro tbh.
pub(crate) fn read_f32<R: Read>(file: &mut R) -> io::Result<f32> {
    let mut buf: [u8; 4] = Default::default();
    file.read_exact(&mut buf)?;
    Ok(f32::from_be_bytes(buf))
}

pub(crate) fn read_i32<R: Read>(file: &mut R) -> io::Result<i32> {
    let mut buf: [u8; 4] = Default::default();
    file.read_exact(&mut buf)?;
    Ok(i32::from_be_bytes(buf))
}

fn read_u32<R: Read>(file: &mut R) -> io::Result<u32> {
    let mut buf: [u8; 4] = Default::default();
    file.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

// CHECKED(2024-03-07 11:51): Looks good.
fn calc_sizeint(
    minint: [i32; 3],
    maxint: [i32; 3],
    sizeint: &mut [u32; 3],
    bitsizeint: &mut [u32; 3],
) -> u32 {
    sizeint[0] = (maxint[0] - minint[0]) as u32 + 1;
    sizeint[1] = (maxint[1] - minint[1]) as u32 + 1;
    sizeint[2] = (maxint[2] - minint[2]) as u32 + 1;

    bitsizeint.fill(0);

    // Check if one of the sizes is too big to be multiplied.
    if (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        return 0; // This flags the use of large sizes. // FIXME: This can become an enum to be more explicit?
    }

    sizeofints(*sizeint)
}

#[inline]
const fn sizeofint(size: u32) -> u32 {
    let mut n = 1;
    let mut nbits = 0;

    while size >= n && nbits < 32 {
        nbits += 1;
        n <<= 1;
    }

    nbits
}

fn sizeofints(sizes: [u32; 3]) -> u32 {
    let mut nbytes = 1;
    let mut bytes = [0u8; 32];
    bytes[0] = 1;
    let mut nbits = 0;

    for size in sizes {
        let mut tmp = 0;
        let mut bytecount = 0;
        while bytecount < nbytes {
            tmp += bytes[bytecount] as u32 * size;
            bytes[bytecount] = (tmp & 0xff) as u8;
            tmp >>= 8;
            bytecount += 1;
        }
        while tmp != 0 {
            bytes[{
                let t = bytecount;
                bytecount += 1;
                t
            }] = (tmp & 0xff) as u8;
            tmp >>= 8;
        }
        nbytes = bytecount;
    }

    nbytes -= 1;
    let mut num = 1;
    while bytes[nbytes] as u32 >= num {
        nbits += 1;
        num *= 2;
    }

    nbytes as u32 * 8 + nbits // FIXME: Check whether it is okay for nbytes to have the type of usize not u32
}

fn decodebyte(buf: &[u8], state: &mut DecodeState) -> u8 {
    let mask = 0xff;

    let DecodeState {
        mut count,
        mut lastbits,
        lastbyte,
    } = *state;
    let mut lastbyte = lastbyte as u32;

    let mut num = 0;
    let mut nbits = 8;
    while nbits >= 8 {
        lastbyte = (lastbyte << 8) | buf[count] as u32;
        count += 1;
        num |= (lastbyte >> lastbits) << (nbits - 8);
        nbits -= 8;
    }

    if nbits > 0 {
        if lastbits < nbits {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | buf[count] as u32;
            count += 1;
        }
        lastbits -= nbits;
        num |= (lastbyte >> lastbits) & mask;
    }

    num &= mask;
    *state = DecodeState {
        count,
        lastbits,
        lastbyte: (lastbyte & 0xff) as u8, // We don't care about anything but the last byte.
    };

    debug_assert_eq!(num & 0xff, num);
    num as u8
}

fn decodebits<T: TryFrom<u32>>(buf: &[u8], state: &mut DecodeState, mut nbits: usize) -> T {
    let mask = (1 << nbits) - 1; // A string of ones that is nbits long.

    let DecodeState {
        mut count,
        mut lastbits,
        lastbyte,
    } = *state;
    let mut lastbyte = lastbyte as u32;

    let mut num = 0;
    while nbits >= 8 {
        lastbyte = (lastbyte << 8) | buf[count] as u32;
        count += 1;
        num |= (lastbyte >> lastbits) << (nbits - 8);
        nbits -= 8;
    }

    if nbits > 0 {
        if lastbits < nbits {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | buf[count] as u32;
            count += 1;
        }
        lastbits -= nbits;
        num |= (lastbyte >> lastbits) & mask;
    }

    num &= mask;
    *state = DecodeState {
        count,
        lastbits,
        lastbyte: (lastbyte & 0xff) as u8, // We don't care about anything but the last byte.
    };

    match num.try_into() {
        Ok(n) => n,
        Err(_) => unreachable!(), // We just checked for that!
    }
}

fn decodeints(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: u32,
    sizes: [u32; 3],
    nums: &mut [i32; 3],
) {
    if nbits <= 32 {
        unpack_from_int_into_u32(buf, state, nbits, sizes, nums);
        return;
    }
    if nbits <= 64 {
        unpack_from_int_into_u64(buf, state, nbits, sizes, nums);
        return;
    }

    let mut bytes = [0u8; 32];
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        bytes[nbytes] = decodebyte(buf, state);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        bytes[nbytes] = decodebits(buf, state, nbits as usize);
        nbytes += 1;
    }

    for i in (0..2).rev() {
        let mut num: u32 = 0;
        for j in 0..nbytes {
            let k = nbytes - 1 - j;
            num = (num << 8) | bytes[k] as u32;
            let p = num / sizes[i];
            bytes[k] = p as u8;
            num -= p * sizes[i];
        }
        nums[i] = num as i32;
    }

    nums[0] = i32::from_le_bytes(bytes[..4].try_into().unwrap());
}

fn unpack_from_int_into_u32(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: u32,
    sizes: [u32; 3],
    nums: &mut [i32; 3],
) {
    type T = u32;
    let mut v: T = 0;
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        let byte: T = decodebyte(buf, state) as T;
        v |= byte << (8 * nbytes as u32);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        let byte: T = decodebits(buf, state, nbits as usize);
        v |= byte << (8 * nbytes as u32);
    }

    // FIXME: What's up with the whole FastType stuff here?
    let sz: T = sizes[2];
    let sy: T = sizes[1];
    let szy: T = sz * sy;
    let x1 = v / szy;
    let q1 = v - x1 * szy;
    let y1 = q1 / sz;
    let z1 = q1 - y1 * sz;

    *nums = [x1, y1, z1].map(|v| v as i32);
}

fn unpack_from_int_into_u64(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: u32,
    sizes: [u32; 3],
    nums: &mut [i32; 3],
) {
    type T = u64;
    let mut v: T = 0;
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        let byte: T = decodebyte(buf, state) as T;
        v |= byte << (8 * nbytes as u32);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        let byte: T = decodebits(buf, state, nbits as usize);
        v |= byte << (8 * nbytes as u32);
    }

    // FIXME: What's up with the whole FastType stuff here?
    let sz: T = sizes[2] as u64;
    let sy: T = sizes[1] as u64;
    let szy: T = sz * sy;
    let x1 = v / szy;
    let q1 = v - x1 * szy;
    let y1 = q1 / sz;
    let z1 = q1 - y1 * sz;

    *nums = [x1, y1, z1].map(|v| v as i32);
}
