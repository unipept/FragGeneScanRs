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

    pub fn to_upper(&self) -> Nuc {
        match self {
            Nuc::A | Nuc::Ai => Nuc::A,
            Nuc::C | Nuc::Ci => Nuc::C,
            Nuc::G | Nuc::Gi => Nuc::G,
            Nuc::T | Nuc::Ti => Nuc::T,
            Nuc::N | Nuc::Ni => Nuc::N,
            Nuc::Xi => Nuc::N,
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
    match [n[0].to_upper(), n[1].to_upper(), n[2].to_upper()] {
        [Nuc::A, Nuc::A, Nuc::A] => Some(0),
        [Nuc::A, Nuc::A, Nuc::C] => Some(1),
        [Nuc::A, Nuc::A, Nuc::G] => Some(2),
        [Nuc::A, Nuc::A, Nuc::T] => Some(3),
        [Nuc::A, Nuc::C, Nuc::A] => Some(4),
        [Nuc::A, Nuc::C, Nuc::C] => Some(5),
        [Nuc::A, Nuc::C, Nuc::G] => Some(6),
        [Nuc::A, Nuc::C, Nuc::T] => Some(7),
        [Nuc::A, Nuc::G, Nuc::A] => Some(8),
        [Nuc::A, Nuc::G, Nuc::C] => Some(9),
        [Nuc::A, Nuc::G, Nuc::G] => Some(10),
        [Nuc::A, Nuc::G, Nuc::T] => Some(11),
        [Nuc::A, Nuc::T, Nuc::A] => Some(12),
        [Nuc::A, Nuc::T, Nuc::C] => Some(13),
        [Nuc::A, Nuc::T, Nuc::G] => Some(14),
        [Nuc::A, Nuc::T, Nuc::T] => Some(15),
        [Nuc::C, Nuc::A, Nuc::A] => Some(16),
        [Nuc::C, Nuc::A, Nuc::C] => Some(17),
        [Nuc::C, Nuc::A, Nuc::G] => Some(18),
        [Nuc::C, Nuc::A, Nuc::T] => Some(19),
        [Nuc::C, Nuc::C, Nuc::A] => Some(20),
        [Nuc::C, Nuc::C, Nuc::C] => Some(21),
        [Nuc::C, Nuc::C, Nuc::G] => Some(22),
        [Nuc::C, Nuc::C, Nuc::T] => Some(23),
        [Nuc::C, Nuc::G, Nuc::A] => Some(24),
        [Nuc::C, Nuc::G, Nuc::C] => Some(25),
        [Nuc::C, Nuc::G, Nuc::G] => Some(26),
        [Nuc::C, Nuc::G, Nuc::T] => Some(27),
        [Nuc::C, Nuc::T, Nuc::A] => Some(28),
        [Nuc::C, Nuc::T, Nuc::C] => Some(29),
        [Nuc::C, Nuc::T, Nuc::G] => Some(30),
        [Nuc::C, Nuc::T, Nuc::T] => Some(31),
        [Nuc::G, Nuc::A, Nuc::A] => Some(32),
        [Nuc::G, Nuc::A, Nuc::C] => Some(33),
        [Nuc::G, Nuc::A, Nuc::G] => Some(34),
        [Nuc::G, Nuc::A, Nuc::T] => Some(35),
        [Nuc::G, Nuc::C, Nuc::A] => Some(36),
        [Nuc::G, Nuc::C, Nuc::C] => Some(37),
        [Nuc::G, Nuc::C, Nuc::G] => Some(38),
        [Nuc::G, Nuc::C, Nuc::T] => Some(39),
        [Nuc::G, Nuc::G, Nuc::A] => Some(40),
        [Nuc::G, Nuc::G, Nuc::C] => Some(41),
        [Nuc::G, Nuc::G, Nuc::G] => Some(42),
        [Nuc::G, Nuc::G, Nuc::T] => Some(43),
        [Nuc::G, Nuc::T, Nuc::A] => Some(44),
        [Nuc::G, Nuc::T, Nuc::C] => Some(45),
        [Nuc::G, Nuc::T, Nuc::G] => Some(46),
        [Nuc::G, Nuc::T, Nuc::T] => Some(47),
        [Nuc::T, Nuc::A, Nuc::A] => Some(48),
        [Nuc::T, Nuc::A, Nuc::C] => Some(49),
        [Nuc::T, Nuc::A, Nuc::G] => Some(50),
        [Nuc::T, Nuc::A, Nuc::T] => Some(51),
        [Nuc::T, Nuc::C, Nuc::A] => Some(52),
        [Nuc::T, Nuc::C, Nuc::C] => Some(53),
        [Nuc::T, Nuc::C, Nuc::G] => Some(54),
        [Nuc::T, Nuc::C, Nuc::T] => Some(55),
        [Nuc::T, Nuc::G, Nuc::A] => Some(56),
        [Nuc::T, Nuc::G, Nuc::C] => Some(57),
        [Nuc::T, Nuc::G, Nuc::G] => Some(58),
        [Nuc::T, Nuc::G, Nuc::T] => Some(59),
        [Nuc::T, Nuc::T, Nuc::A] => Some(60),
        [Nuc::T, Nuc::T, Nuc::C] => Some(61),
        [Nuc::T, Nuc::T, Nuc::G] => Some(62),
        [Nuc::T, Nuc::T, Nuc::T] => Some(63),
        _ => None,
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
