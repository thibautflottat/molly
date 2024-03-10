#![feature(array_chunks)]

fn main() {
    let bytes= [1u8, 2 , 3 , 4, 5, 6, 7];
    let num= bytes[0] as i32
        | ((bytes[1] as i32) << 8)
        | ((bytes[2] as i32) << 16)

        | ((bytes[3] as i32) << 24);
    dbg!(num);
    let head = bytes.array_chunks().nth(0).unwrap();
    let i = i32::from_le_bytes(*head);
    assert_eq!(num, i);
}
