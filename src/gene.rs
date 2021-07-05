use std::fs::File;
use std::io;
use std::io::Write;

extern crate seq_io;
use seq_io::fasta;
use seq_io::fasta::Record;

extern crate thiserror;
use thiserror::Error;

use crate::dna::Nuc::{A, C, G, T};
use crate::dna::{trinucleotide, Nuc, ANTI_CODON_CODE, CODON_CODE};

pub struct Gene {
    pub head: Vec<u8>,
    pub start: usize,
    pub metastart: usize,
    pub end: usize,
    pub frame: usize,
    pub score: f64,
    pub dna: Vec<Nuc>,
    pub forward_strand: bool,
    pub inserted: Vec<usize>,
    pub deleted: Vec<usize>,
}

impl Gene {
    pub fn print_meta(&self, file: &mut File) -> Result<(), GeneError> {
        fasta::OwnedRecord {
            head: self.head.to_owned(),
            seq: format!(
                "{}\t{}\t{}\t{}\t{}\tI:{}\tD:{}",
                self.metastart,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                self.frame,
                self.score,
                self.inserted
                    .iter()
                    .map(|i: &usize| { format!("{},", i) })
                    .collect::<String>(),
                self.deleted
                    .iter()
                    .map(|i: &usize| { format!("{},", i) })
                    .collect::<String>()
            )
            .into_bytes(),
        }
        .write(file)?;
        Ok(())
    }

    pub fn print_dna(&self, file: &mut File, formatted: bool) -> Result<(), GeneError> {
        let dna: Vec<u8> = match (self.forward_strand, formatted) {
            (true, true) => self.dna.iter().map(|&n| u8::from(n)).collect(),
            (true, false) => self
                .dna
                .iter()
                .filter(|n| !n.is_insertion())
                .map(|&n| u8::from(n))
                .collect(),
            (false, true) => self.dna.iter().rev().map(|&n| u8::from(n.rc())).collect(),
            (false, false) => self
                .dna
                .iter()
                .rev()
                .filter(|n| !n.is_insertion())
                .map(|&n| u8::from(n.rc()))
                .collect(),
        };
        fasta::OwnedRecord {
            head: format!(
                "{}_{}_{}_{}",
                std::str::from_utf8(&self.head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' }
            )
            .into_bytes(),
            seq: dna,
        }
        .write(file)?;
        Ok(())
    }

    pub fn print_protein<W: Write>(
        &self,
        whole_genome: bool,
        file: &mut W,
    ) -> Result<(), GeneError> {
        let dna = self
            .dna
            .iter()
            .filter(|n| !n.is_insertion())
            .map(|&n| n)
            .collect::<Vec<Nuc>>();
        let mut protein: Vec<u8> = if self.forward_strand {
            dna.chunks_exact(3)
                .map(|c| {
                    trinucleotide(c[0], c[1], c[2])
                        .map(|i| CODON_CODE[i])
                        .unwrap_or(b'X')
                })
                .collect()
        } else {
            dna.rchunks_exact(3)
                .map(|c| {
                    trinucleotide(c[0], c[1], c[2])
                        .map(|i| ANTI_CODON_CODE[i])
                        .unwrap_or(b'X')
                })
                .collect()
        };
        if protein.last() == Some(&b'*') {
            protein.pop();
        }

        // alternative start codons still encode for Met
        // E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG and one or two others (e.g., an AUU and possibly a CUG)
        // only consider two major alternative ones, GTG and TTG
        if whole_genome {
            if self.forward_strand {
                let s = trinucleotide(self.dna[0], self.dna[1], self.dna[2]);
                if s == trinucleotide(G, T, G) || s == trinucleotide(T, T, G) {
                    protein[0] = b'M';
                }
            } else {
                let s = trinucleotide(
                    self.dna[self.dna.len() - 3],
                    self.dna[self.dna.len() - 2],
                    self.dna[self.dna.len() - 1],
                );
                if s == trinucleotide(C, A, C) || s == trinucleotide(C, A, A) {
                    protein[0] = b'M';
                }
            }
        }

        fasta::OwnedRecord {
            head: format!(
                "{}_{}_{}_{}",
                std::str::from_utf8(&self.head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' }
            )
            .into_bytes(),
            seq: protein,
        }
        .write(file)?;

        Ok(())
    }
}

#[derive(Error, Debug)]
pub enum GeneError {
    #[error("could not write to file")]
    IoError(#[from] io::Error),
    #[error("could not convert header back to UTF-8")]
    Utf8Error(#[from] std::str::Utf8Error),
}
