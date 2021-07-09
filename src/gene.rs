use std::fs::File;
use std::io;
use std::io::Write;

extern crate thiserror;
use thiserror::Error;

use crate::dna::Nuc::{A, C, G, T};
use crate::dna::{trinucleotide, Nuc, ANTI_CODON_CODE, CODON_CODE};

pub struct ReadPrediction {
    pub head: Vec<u8>,
    pub genes: Vec<Gene>,
}

impl ReadPrediction {
    pub fn new(head: Vec<u8>) -> Self {
        ReadPrediction {
            head: head,
            genes: vec![],
        }
    }

    pub fn print_meta(&self, file: &mut File) -> Result<(), GeneError> {
        if !self.genes.is_empty() {
            file.write_all(&format!(">{}\n", std::str::from_utf8(&self.head)?).into_bytes())?;
        }
        for gene in &self.genes {
            gene.print_meta(file)?;
        }
        Ok(())
    }

    pub fn print_dna(&self, file: &mut File, formatted: bool) -> Result<(), GeneError> {
        for gene in &self.genes {
            gene.print_dna(file, &self.head, formatted)?;
        }
        Ok(())
    }

    pub fn print_protein<W: Write>(
        &self,
        whole_genome: bool,
        file: &mut W,
    ) -> Result<(), GeneError> {
        for gene in &self.genes {
            gene.print_protein(file, &self.head, whole_genome)?;
        }
        Ok(())
    }
}

pub struct Gene {
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
        file.write_all(
            &format!(
                "{}\t{}\t{}\t{}\t{:.6}\tI:{}\tD:{}\n",
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
        )?;
        Ok(())
    }

    pub fn print_dna(
        &self,
        file: &mut File,
        head: &Vec<u8>,
        formatted: bool,
    ) -> Result<(), GeneError> {
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

        file.write_all(
            &format!(
                ">{}_{}_{}_{}\n{}\n",
                std::str::from_utf8(head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                std::str::from_utf8(&dna)?,
            )
            .into_bytes(),
        )?;

        Ok(())
    }

    pub fn print_protein<W: Write>(
        &self,
        file: &mut W,
        head: &Vec<u8>,
        whole_genome: bool,
    ) -> Result<(), GeneError> {
        let dna = self
            .dna
            .iter()
            .filter(|n| !n.is_insertion())
            .map(|&n| n)
            .collect::<Vec<Nuc>>();
        let mut protein: Vec<u8> = if self.forward_strand {
            dna.chunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        } else {
            dna.rchunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| ANTI_CODON_CODE[i]).unwrap_or(b'X'))
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
                let s = trinucleotide(&self.dna);
                if s == trinucleotide(&[G, T, G]) || s == trinucleotide(&[T, T, G]) {
                    protein[0] = b'M';
                }
            } else {
                let s = trinucleotide(self.dna.get(self.dna.len() - 3..).unwrap());
                if s == trinucleotide(&[C, A, C]) || s == trinucleotide(&[C, A, A]) {
                    protein[0] = b'M';
                }
            }
        }

        file.write_all(
            &format!(
                ">{}_{}_{}_{}\n{}\n",
                std::str::from_utf8(head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' },
                std::str::from_utf8(&protein)?,
            )
            .into_bytes(),
        )?;
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
