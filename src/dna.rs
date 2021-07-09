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
    Ai,
    Ci,
    Gi,
    Ti,
    Ni,
    Xi,
}

impl Nuc {
    pub fn to_int(&self) -> Option<usize> {
        match self {
            Nuc::A | Nuc::Ai => Some(0),
            Nuc::C | Nuc::Ci => Some(1),
            Nuc::G | Nuc::Gi => Some(2),
            Nuc::T | Nuc::Ti => Some(3),
            Nuc::N | Nuc::Ni => None,
            Nuc::Xi => None,
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
            b'A' | b'a' => Nuc::A,
            b'C' | b'c' => Nuc::C,
            b'G' | b'g' => Nuc::G,
            b'T' | b't' => Nuc::T,
            b'x' => Nuc::Xi,
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

pub fn trinucleotide(n: &[Nuc]) -> Option<usize> {
    if let (Some(a), Some(b), Some(c)) = (n[0].to_int(), n[1].to_int(), n[2].to_int()) {
        Some(16 * a + 4 * b + c)
    } else {
        None
    }
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
