use std::cmp::{max, min};

pub const CG_MIN: usize = 26;
pub const CG_MAX: usize = 70;

pub const ACGT: usize = 4;
pub const BI_ACGT: usize = 4 * 4;
pub const TRI_ACGT: usize = 4 * 4 * 4;

// Optimization: #[repr(u8)] and specific ordering (A, Ai, C, Ci, ...)
// allows efficient bitwise calculation in to_int() by grouping variants.
#[repr(u8)]
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Nuc {
    A,
    Ai,
    C,
    Ci,
    G,
    Gi,
    T,
    Ti,
    N,
    Ni,
    Xi,
}

impl Nuc {
    #[inline]
    pub fn to_int(&self) -> Option<usize> {
        let val = *self as u8;
        if val < 8 {
            Some((val >> 1) as usize)
        } else {
            None
        }
    }

    pub fn to_lower(&self) -> Nuc {
        match self {
            Nuc::A | Nuc::Ai => Nuc::Ai,
            Nuc::C | Nuc::Ci => Nuc::Ci,
            Nuc::G | Nuc::Gi => Nuc::Gi,
            Nuc::T | Nuc::Ti => Nuc::Ti,
            Nuc::N | Nuc::Ni => Nuc::Ni,
            Nuc::Xi => Nuc::Xi,
        }
    }

    pub fn is_insertion(&self) -> bool {
        match self {
            Nuc::Ai | Nuc::Ci | Nuc::Gi | Nuc::Ti | Nuc::Ni => true,
            _ => false,
        }
    }

    pub fn rc(&self) -> Nuc {
        match self {
            Nuc::A => Nuc::T,
            Nuc::C => Nuc::G,
            Nuc::G => Nuc::C,
            Nuc::T => Nuc::A,
            Nuc::N => Nuc::N,
            Nuc::Ai => Nuc::Ti,
            Nuc::Ci => Nuc::Gi,
            Nuc::Gi => Nuc::Ci,
            Nuc::Ti => Nuc::Ai,
            Nuc::Ni => Nuc::Ni,
            Nuc::Xi => Nuc::Xi,
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
            b'a' => Nuc::Ai,
            b'c' => Nuc::Ci,
            b'g' => Nuc::Gi,
            b't' => Nuc::Ti,
            b'x' => Nuc::Xi,
            b'n' => Nuc::Ni,
            _ => Nuc::N,
        }
    }
}

impl From<Nuc> for u8 {
    fn from(nt: Nuc) -> u8 {
        match nt {
            Nuc::A => b'A',
            Nuc::C => b'C',
            Nuc::G => b'G',
            Nuc::T => b'T',
            Nuc::N => b'N',
            Nuc::Ai => b'a',
            Nuc::Ci => b'c',
            Nuc::Gi => b'g',
            Nuc::Ti => b't',
            Nuc::Ni => b'n',
            Nuc::Xi => b'N',
        }
    }
}

pub fn dna(seq: &str) -> Vec<Nuc> {
    seq.as_bytes().into_iter().copied().map(Nuc::from).collect()
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

// Optimization: Direct arithmetic calculation of codon index (0-63).
// Treats nucleotides as base-4 digits (A=0, C=1, G=2, T=3).
pub fn trinucleotide(n: &[Nuc]) -> Option<usize> {
    let s = n.get(0..3).unwrap();
    let n0 = s[0].to_int()?;
    let n1 = s[1].to_int()?;
    let n2 = s[2].to_int()?;
    Some(n0 * 16 + n1 * 4 + n2)
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
