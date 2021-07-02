use std::cmp::{max, min};

pub const CODON_CODE: [u8; 65] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
    b'X',
];

pub const ANTI_CODON_CODE: [u8; 65] = [
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'M', b'W', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'I', b'*', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'X',
];

pub fn nt2int(nt: u8) -> Option<usize> {
    match nt {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub fn trinucleotide(a: u8, b: u8, c: u8) -> usize {
    if let (Some(a_), Some(b_), Some(c_)) = (nt2int(a), nt2int(b), nt2int(c)) {
        16 * a_ + 4 * b_ + c_
    } else {
        0
    }
}

pub fn trinucleotide_pep(a: u8, b: u8, c: u8) -> usize {
    if let (Some(a_), Some(b_), Some(c_)) = (nt2int(a), nt2int(b), nt2int(c)) {
        16 * a_ + 4 * b_ + c_
    } else {
        64
    }
}

pub fn get_rc_dna(dna: &Vec<u8>) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'N' => b'N',
            b'n' => b'n',
            _ => b'x',
        })
        .collect()
}

pub fn count_cg_content(seq: &[u8]) -> usize {
    let mut count = 0;
    for l in seq.iter() {
        if b"CcGg".contains(l) {
            count += 1;
        }
    }
    min(43, max(26, count * 100 / seq.len()) - 26)
}
