use std::error::Error;
use std::fs::File;
use std::io::Write;

extern crate seq_io;
use seq_io::fasta;
use seq_io::fasta::Record;

use crate::dna::{trinucleotide_pep, ANTI_CODON_CODE, CODON_CODE};

pub struct Gene {
    pub head: Vec<u8>,
    pub start: usize,
    pub metastart: usize,
    pub end: usize,
    pub frame: usize,
    pub score: f64,
    pub dna: Vec<u8>,
    pub dna_ffn: Vec<u8>,
    pub forward_strand: bool,
    pub inserted: Vec<usize>,
    pub deleted: Vec<usize>,
}

impl Gene {
    pub fn print_meta(&self, file: &mut File) -> Result<(), Box<dyn Error>> {
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

    pub fn print_dna(&self, file: &mut File) -> Result<(), Box<dyn Error>> {
        fasta::OwnedRecord {
            head: format!(
                "{}_{}_{}_{}",
                std::str::from_utf8(&self.head)?,
                self.start,
                self.end,
                if self.forward_strand { '+' } else { '-' }
            )
            .into_bytes(),
            seq: self.dna_ffn.clone(),
        }
        .write(file)?;
        Ok(())
    }

    pub fn print_protein<W: Write>(
        &self,
        whole_genome: bool,
        file: &mut W,
    ) -> Result<(), Box<dyn Error>> {
        let mut protein: Vec<u8> = if self.forward_strand {
            self.dna
                .chunks_exact(3)
                .map(|c| CODON_CODE[trinucleotide_pep(c[0], c[1], c[2])])
                .collect()
        } else {
            self.dna
                .rchunks_exact(3)
                .map(|c| ANTI_CODON_CODE[trinucleotide_pep(c[0], c[1], c[2])])
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
                let s = trinucleotide_pep(self.dna[0], self.dna[1], self.dna[2]);
                if s == trinucleotide_pep(b'G', b'T', b'G')
                    || s == trinucleotide_pep(b'T', b'T', b'G')
                {
                    protein[0] = b'M';
                }
            } else {
                let s = trinucleotide_pep(
                    self.dna[self.dna.len() - 3],
                    self.dna[self.dna.len() - 2],
                    self.dna[self.dna.len() - 1],
                );
                if s == trinucleotide_pep(b'C', b'A', b'C')
                    || s == trinucleotide_pep(b'C', b'A', b'A')
                {
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
