use std::cmp::{max, min};

pub const CG_MIN: usize = 26;
pub const CG_MAX: usize = 70;

pub const ACGT: usize = 4;
pub const BI_ACGT: usize = 4 * 4;
pub const TRI_ACGT: usize = 4 * 4 * 4;

#[derive(PartialEq, Eq, Clone, Copy)]
pub enum Nuc {
    A,
    C,
    G,
    T,
    N,
}

impl Nuc {
    pub fn to_int(&self) -> Option<usize> {
        match self {
            Nuc::A => Some(0),
            Nuc::C => Some(1),
            Nuc::G => Some(2),
            Nuc::T => Some(3),
            Nuc::N => None,
        }
    }

    pub fn to_lower(&self) -> u8 {
        match self {
            Nuc::A => b'a',
            Nuc::C => b'c',
            Nuc::G => b'g',
            Nuc::T => b't',
            Nuc::N => b'n',
        }
    }

    pub fn to_upper(&self) -> u8 {
        match self {
            Nuc::A => b'A',
            Nuc::C => b'C',
            Nuc::G => b'G',
            Nuc::T => b'T',
            Nuc::N => b'N',
        }
    }
}

impl From<u8> for Nuc {
    fn from(nt: u8) -> Nuc {
        match nt {
            b'A' => Nuc::A,
            b'C' => Nuc::C,
            b'G' => Nuc::G,
            b'T' => Nuc::T,
            _ => Nuc::N,
        }
    }
}

pub const CODON_CODE: [u8; TRI_ACGT] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
];

pub const ANTI_CODON_CODE: [u8; TRI_ACGT] = [
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'M', b'W', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'I', b'*', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
];

pub fn trinucleotide(a: Nuc, b: Nuc, c: Nuc) -> Option<usize> {
    if let (Some(a_), Some(b_), Some(c_)) = (a.to_int(), b.to_int(), c.to_int()) {
        Some(16 * a_ + 4 * b_ + c_)
    } else {
        None
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

pub fn count_cg_content(seq: &[Nuc]) -> usize {
    let mut count = 0;
    for &l in seq.iter() {
        if l == Nuc::C || l == Nuc::G {
            count += 1;
        }
    }
    min(
        CG_MAX - CG_MIN - 1,
        max(CG_MIN, count * 100 / seq.len()) - CG_MIN,
    )
}
