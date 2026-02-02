use std::cmp::{max, min};

pub const CG_MIN: usize = 26;
pub const CG_MAX: usize = 70;

pub const ACGT: usize = 4;
pub const BI_ACGT: usize = 4 * 4;
pub const TRI_ACGT: usize = 4 * 4 * 4;

/// Nucleotide representation optimized for bit operations.
///
/// Layout: bit 0 = insertion flag, bits 1-2 = base encoding
/// - A=0, Ai=1, C=2, Ci=3, G=4, Gi=5, T=6, Ti=7
/// - N, Ni, Xi are special values >= 8
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(u8)]
pub enum Nuc {
    A = 0,
    Ai = 1,
    C = 2,
    Ci = 3,
    G = 4,
    Gi = 5,
    T = 6,
    Ti = 7,
    N = 8,
    Ni = 9,
    Xi = 10,
}

impl Nuc {
    /// Convert to base index (0=A, 1=C, 2=G, 3=T), None for N/Ni/Xi
    #[inline]
    pub fn to_int(&self) -> Option<usize> {
        let v = *self as u8;
        if v < 8 {
            Some((v >> 1) as usize)
        } else {
            None
        }
    }

    /// Convert to lowercase (insertion) variant
    #[inline]
    pub fn to_lower(&self) -> Nuc {
        let v = *self as u8;
        // Set bit 0 to make it an insertion variant
        // SAFETY: v | 1 is always a valid Nuc variant
        unsafe { std::mem::transmute(v | 1) }
    }

    /// Check if this is an insertion (lowercase) variant
    #[inline]
    pub fn is_insertion(&self) -> bool {
        (*self as u8) & 1 == 1
    }

    /// Reverse complement: A<->T, C<->G
    #[inline]
    pub fn rc(&self) -> Nuc {
        let v = *self as u8;
        if v < 8 {
            // XOR with 6 swaps A(0)<->T(6) and C(2)<->G(4)
            // while preserving bit 0 (insertion flag)
            // SAFETY: v ^ 6 for v in 0..8 is always in 0..8
            unsafe { std::mem::transmute(v ^ 6) }
        } else {
            // N, Ni, Xi stay the same
            *self
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

/// Convert a trinucleotide to its codon index (0-63).
/// Returns None if any nucleotide is N, Ni, or Xi.
#[inline]
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
