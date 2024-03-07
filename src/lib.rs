#![feature(array_chunks, iter_array_chunks, array_try_map)]

use glam::{Mat3, Vec3};

pub type BoxVec = Mat3;

#[derive(Debug, Default, Clone)]
pub struct Frame {
    pub step: u32,
    /// Time in picoseconds.
    pub time: f32,
    pub boxvec: BoxVec,
    pub positions: Vec<Vec3>,
}

#[derive(Debug, Clone)]
pub struct XTCReader<R> {
    pub file: R,
    pub frames: Vec<Frame>,
    step: usize,
}

impl<R: std::io::Read> XTCReader<R> {
    pub const MAGIC: i32 = 1995;

    pub fn new(reader: R) -> Self {
        Self {
            file: reader,
            frames: Vec::new(),
            step: 0,
        }
    }

    /// Reads and returns a [`Frame`] and advances one step.
    pub fn read_frame(&mut self) -> std::io::Result<Frame> {
        let file = &mut self.file;

        // Start of by reading the header.
        let magic = read_i32(file)?;
        assert_eq!(magic, Self::MAGIC);
        let natoms = read_i32(file)?.try_into().unwrap();
        let step: u32 = read_i32(file)?.try_into().unwrap();
        let time = read_f32(file)?;

        // Read the frame data.
        let boxvec = read_boxvec(file)?;
        let natoms_repeated = read_i32(file)?.try_into().unwrap();
        assert_eq!(natoms, natoms_repeated);

        // Now, we read the atoms.
        let mut atoms = Vec::with_capacity(natoms);
        if natoms <= 9 {
            // In case the number of atoms is very small, just read their uncompressed positions.
            let positions = read_f32s(file, natoms * 3)?;
            atoms.extend(positions.array_chunks().map(|pos| Vec3::from_array(pos)));
        } else {
            let precision = read_f32(file)?;
            let positions = read_compressed_floats(file, natoms * 3, precision)?;
            atoms.extend(positions.array_chunks().map(|&pos| Vec3::from_array(pos)));
        }

        self.step += 1;

        Ok(Frame {
            step,
            time,
            boxvec,
            positions: atoms,
        })
    }
}

fn read_boxvec(file: &mut impl std::io::Read) -> std::io::Result<BoxVec> {
    let boxvec: Vec<_> = read_f32s(file, 9)?.collect();
    let cols = [
        [boxvec[0], boxvec[3], boxvec[6]],
        [boxvec[1], boxvec[4], boxvec[7]],
        [boxvec[2], boxvec[5], boxvec[8]],
    ];
    Ok(BoxVec::from_cols_array_2d(&cols))
}

fn read_opaque(file: &mut impl std::io::Read, data: &mut Vec<u8>) -> std::io::Result<()> {
    let count = read_u32(file)? as usize;
    let padding = (4 - (count % 4)) % 4; // FIXME: Why, and also, can we do this better?
    data.resize(count + padding, 0);
    file.read_exact(data)
}

fn read_f32s(
    file: &mut impl std::io::Read,
    n: usize,
) -> std::io::Result<impl Iterator<Item = f32>> {
    let mut buf = vec![0; n * 4]; // TODO: Gotta remove this.
    file.read_exact(&mut buf)?;
    Ok(buf
        .into_iter()
        .array_chunks()
        .map(|chunk| f32::from_be_bytes(chunk)))
}

// FIXME: These read_* functions are prime targets for a macro tbh.
fn read_f32(file: &mut impl std::io::Read) -> std::io::Result<f32> {
    let mut buf: [u8; 4] = Default::default();
    file.read_exact(&mut buf)?;
    Ok(f32::from_be_bytes(buf))
}

fn read_i32(file: &mut impl std::io::Read) -> std::io::Result<i32> {
    let mut buf: [u8; 4] = Default::default();
    file.read_exact(&mut buf)?;
    Ok(i32::from_be_bytes(buf))
}

fn read_u32(file: &mut impl std::io::Read) -> std::io::Result<u32> {
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

fn sizeofint(size: u32) -> u32 {
    let mut n = 1;
    let mut nbits = 0;

    while size >= n && nbits < 32 {
        nbits += 1;
        n <<= 1;
    }

    eprintln!("nbits = {nbits}");
    nbits
}

// CHECKED(2024-03-07 11:43): Looks good (for at least the first entries).
fn sizeofints(sizes: [u32; 3]) -> u32 {
    let nints = sizes.len(); // FIXME: Inline this.
    let mut nbytes = 1;
    let mut bytes = [0u8; 32];
    bytes[0] = 1;
    let mut nbits = 0;

    for i in 0..nints {
        let mut tmp = 0;
        let mut bytecount = 0;
        while bytecount < nbytes {
            tmp = bytes[bytecount] as u32 * sizes[i] + tmp; // FIXME: +=? or overflows?
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

    // eprintln!("sizeofints -> {}", nbytes as u32 * 8 + nbits);
    nbytes as u32 * 8 + nbits // FIXME: Check whether it is okay for nbytes to have the type of usize not u32
}

// CHECKED(2024-03-07 11:55): Seems fine from a cursory look. May want to look further.
fn decodebits<T: TryFrom<u32> + std::fmt::Debug>(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: usize,
) -> T {
    // eprintln!("decodebits::nbits = {nbits}");
    let mask = (1 << nbits) - 1;

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

    assert!(std::mem::size_of::<T>() * 8 >= nbits);
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
    bytes[..3].fill(0);
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        bytes[nbytes] = decodebits(buf, state, 8);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        bytes[nbytes] = decodebits(buf, state, nbits as usize);
        nbytes += 1;
    }

    for i in 2..0 {
        let mut num: u32 = 0;
        for j in 0..nbytes {
            let k = nbytes - 1 - j;
            num = (num << 8) | bytes[k] as u32;
            let p = num / sizes[i];
            bytes[k] = p as u8;
            num = num - p * sizes[i]; // FIXME: Is there a reason not to use -= here?
        }
        nums[i] = num as i32;
    }

    // FIXME: Isn't this just a little endian cast of [u8; 4] to i32?
    nums[0] = bytes[0] as i32
        | ((bytes[1] as i32) << 8)
        | ((bytes[2] as i32) << 16)
        | ((bytes[3] as i32) << 24);
}

fn unpack_from_int_into_u32(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: u32,
    sizes: [u32; 3],
    nums: &mut [i32; 3],
) 
{
    type T = u32;
    let mut v: T = 0;
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        let byte: T = decodebits(buf, state, 8);
        assert_eq!(nbytes as u32 as usize, nbytes);
        v |= byte << (8 * nbytes as u32);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        let byte: T = decodebits(buf, state, nbits as usize);
        v |= byte << (8 * nbytes as u32);
    }

    // FIXME: What's up with the whole FastType stuff here?
    let sz: T = sizes[2].into();
    let sy: T = sizes[1].into();
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
) 
{
    type T = u64;
    let mut v: T = 0;
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        let byte: T = decodebits(buf, state, 8);
        assert_eq!(nbytes as u32 as usize, nbytes);
        v |= byte << (8 * nbytes as u32);
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        let byte: T = decodebits(buf, state, nbits as usize);
        v |= byte << (8 * nbytes as u32);
    }

    // FIXME: What's up with the whole FastType stuff here?
    let sz: T = sizes[2].into();
    let sy: T = sizes[1].into();
    let szy: T = sz * sy;
    let x1 = v / szy;
    let q1 = v - x1 * szy;
    let y1 = q1 / sz;
    let z1 = q1 - y1 * sz;

    *nums = [x1, y1, z1].map(|v| v as i32);
}

fn unpack_from_int<T>(
    buf: &[u8],
    state: &mut DecodeState,
    mut nbits: u32,
    sizes: [u32; 3],
    nums: &mut [i32; 3],
) where
    T: From<u32>
        + std::ops::Div<T, Output = T>
        + std::ops::Mul<T, Output = T>
        + std::ops::Shl<T, Output = T>
        + std::ops::BitOrAssign<T>
        + std::ops::Shl<T>
        + std::ops::Sub<T, Output = T>
        + Into<u64>
        + Copy
        + std::fmt::Debug,
{
    let mut v: T = unsafe { std::mem::zeroed() };
    let mut nbytes: usize = 0;
    while nbits >= 8 {
        let byte: T = decodebits(buf, state, 8);
        assert_eq!(nbytes as u32 as usize, nbytes);
        v |= byte << (8 * nbytes as u32).into();
        nbytes += 1;
        nbits -= 8;
    }
    if nbits > 0 {
        let byte: T = decodebits(buf, state, nbits as usize);
        v |= byte << (8 * nbytes as u32).into();
    }

    // FIXME: What's up with the whole FastType stuff here?
    let sz: T = sizes[2].into();
    let sy: T = sizes[1].into();
    let szy: T = sz * sy;
    let x1 = v / szy;
    let q1 = v - x1 * szy;
    let y1 = q1 / sz;
    let z1 = q1 - y1 * sz;

    *nums = [x1, y1, z1].map(|v| v.into() as i32);
}

struct DecodeState {
    count: usize,
    lastbits: usize,
    lastbyte: u8,
}

fn read_compressed_floats(
    file: &mut impl std::io::Read,
    n: usize,
    precision: f32,
) -> std::io::Result<Box<[f32]>> {
    // TODO: I have a constexpr laying around for this somewhere.
    #[rustfmt::skip]
    const  MAGICINTS : [i32; 73] = [
        0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
        10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
        101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
        1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
        10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
        104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
        1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
        10568983, 13316085, 16777216
    ];

    const FIRSTIDX: usize = 9; // Note that MAGICINTS[FIRSTIDX-1] == 0.
    const LASTIDX: usize = MAGICINTS.len(); // FIXME: Check if this is actually equivalent to #define LASTIDX (sizeof(MAGICINTS) / sizeof(*MAGICINTS))

    let minint = [0; 3].try_map(|_| read_i32(file))?;
    let maxint = [0; 3].try_map(|_| read_i32(file))?;
    let mut smallidx = read_u32(file)? as usize;
    assert!(smallidx < LASTIDX);

    let mut sizeint = [0u32; 3];
    let mut bitsizeint = [0u32; 3];
    let bitsize = calc_sizeint(minint, maxint, &mut sizeint, &mut bitsizeint);
    // eprintln!("   sizeint = {sizeint:?}");
    // eprintln!("bitsizeint = {bitsizeint:?}");

    let tmpidx = smallidx - 1;
    let tmpidx = if FIRSTIDX > tmpidx { FIRSTIDX } else { tmpidx };

    let mut smaller = MAGICINTS[tmpidx] / 2;
    let mut smallnum = MAGICINTS[smallidx] / 2;
    let mut sizesmall = [MAGICINTS[smallidx] as u32; 3];

    let mut compressed_data = Vec::new();
    read_opaque(file, &mut compressed_data)?;
    let mut data = {
        let mut data = Vec::with_capacity(n);
        data.resize(n, 0.0);
        data.into_boxed_slice()
    };

    assert_eq!(n % 3, 0, "length of data should be divisible by 3");
    let natoms = n / 3;

    let mut state = DecodeState {
        count: 0,
        lastbits: 0,
        lastbyte: 0,
    };
    let mut run: i32 = 0;
    let mut prevcoord = [0i32; 3];
    let inv_precision = 1.0 / precision;
    let mut write_idx = 0;
    let mut read_idx = 0;
    while read_idx < natoms {
        let mut thiscoord = [0i32; 3];
        let mut thiscoord_fl = data.get_mut(write_idx * 3..(write_idx + 1) * 3).unwrap();

        if bitsize == 0 {
            thiscoord[0] = decodebits(&compressed_data, &mut state, bitsizeint[0] as usize);
            thiscoord[1] = decodebits(&compressed_data, &mut state, bitsizeint[1] as usize);
            thiscoord[2] = decodebits(&compressed_data, &mut state, bitsizeint[2] as usize);
        } else {
            decodeints(
                &compressed_data,
                &mut state,
                bitsize,
                sizeint,
                &mut thiscoord,
            );
        }

        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];

        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];

        let flag: bool = decodebits::<u8>(&compressed_data, &mut state, 1) > 0;
        let mut is_smaller = 0;
        if flag {
            run = decodebits(&compressed_data, &mut state, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller -= 1;
        }
        if run > 0 && write_idx * 3 + run as usize > n {
            panic!("attempt to write a run beyond the positions buffer")
        }
        if run > 0 {
            // Let's read the next coordinate.
            thiscoord.fill(0);

            for k in (0..run).step_by(3) {
                decodeints(
                    &compressed_data,
                    &mut state,
                    smallidx as u32,
                    sizesmall,
                    &mut thiscoord,
                );
                read_idx += 1;
                thiscoord[0] += prevcoord[0] - smallnum;
                thiscoord[1] += prevcoord[1] - smallnum;
                thiscoord[2] += prevcoord[2] - smallnum;
                if k == 0 {
                    // Swap the first and second atom. This is done to achieve better compression
                    // for water atoms. Waters are stored as OHH, but right now we want to swap the
                    // atoms such that e.g., water will become HOH again.
                    std::mem::swap(&mut thiscoord[0], &mut prevcoord[0]);
                    std::mem::swap(&mut thiscoord[1], &mut prevcoord[1]);
                    std::mem::swap(&mut thiscoord[2], &mut prevcoord[2]);
                    thiscoord_fl[0] = prevcoord[0] as f32 * inv_precision;
                    thiscoord_fl[1] = prevcoord[1] as f32 * inv_precision;
                    thiscoord_fl[2] = prevcoord[2] as f32 * inv_precision;
                    write_idx += 1;
                    // thiscoord_fl = &mut data[write_idx * 3..(write_idx + 1) * 3];
                    thiscoord_fl = data.get_mut(write_idx * 3..(write_idx + 1) * 3).unwrap();
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                thiscoord_fl[0] = thiscoord[0] as f32 * inv_precision;
                thiscoord_fl[1] = thiscoord[1] as f32 * inv_precision;
                thiscoord_fl[2] = thiscoord[2] as f32 * inv_precision;
                write_idx += 1;
                thiscoord_fl = match data.get_mut(write_idx * 3..(write_idx + 1) * 3) {
                    Some(c) => c,
                    None => break,
                };
            }
        } else {
            thiscoord_fl[0] = thiscoord[0] as f32 * inv_precision;
            thiscoord_fl[1] = thiscoord[1] as f32 * inv_precision;
            thiscoord_fl[2] = thiscoord[2] as f32 * inv_precision;
            write_idx += 1;
        }

        if is_smaller < 0 {
            smallidx -= 1;
            smallnum = smaller;
            if smallidx > FIRSTIDX {
                smaller = MAGICINTS[smallidx - 1] / 2;
            } else {
                smaller = 0;
            }
        } else if is_smaller > 0 {
            smallidx += 1;
            smaller = smallnum;
            smallnum = MAGICINTS[smallidx] / 2;
        }
        if MAGICINTS[smallidx] == 0 {
            panic!("found an invalid size")
        }
        sizesmall.fill(MAGICINTS[smallidx] as u32);
        read_idx += 1;
    }

    Ok(data)
}
