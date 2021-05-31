//! FragGeneScanRs executable
#![allow(non_snake_case)]

use std::cmp::{max, min};
use std::error::Error;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;

extern crate clap;
use clap::{App, Arg};

extern crate seq_io;
use seq_io::fasta;
use seq_io::policy::StdPolicy;

extern crate frag_gene_scan_rs;
use frag_gene_scan_rs::hmm;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("FragGeneScanRs")
        .version("0.0.1")
        .author("Felix Van der Jeugt <felix.vanderjeugt@ugent.be>")
        .about("Scalable high-throughput short-read open reading frame prediction.")
        .arg(Arg::with_name("seq-file-name")
            .short("s")
            .long("seq-file-name")
            .value_name("seq_file_name")
            .takes_value(true)
            .default_value("stdin")
            .help("Sequence file name including the full path."))
        .arg(Arg::with_name("output-file-name")
            .short("o")
            .long("output-file-name")
            .value_name("output_file_name")
            .takes_value(true)
            .default_value("stdout")
            .help("Output file name including the full path."))
        .arg(Arg::with_name("complete")
            .short("w")
            .long("complete")
            .help("The input sequence has complete genomic sequences; not short sequence reads."))
        .arg(Arg::with_name("train-file")
            .short("t")
            .long("training-file")
            .value_name("train_file_name")
            .takes_value(true)
            .required(true)
            .help("File name that contains model parameters; this file should be in the -r directory or one of the following:
[complete] for complete genomic sequences or short sequence reads without sequencing error
[sanger_5] for Sanger sequencing reads with about 0.5% error rate
[sanger_10] for Sanger sequencing reads with about 1% error rate
[454_5] for 454 pyrosequencing reads with about 0.5% error rate
[454_10] for 454 pyrosequencing reads with about 1% error rate
[454_30] for 454 pyrosequencing reads with about 3% error rate
[illumina_1] for Illumina sequencing reads with about 0.1% error rate
[illumina_5] for Illumina sequencing reads with about 0.5% error rate
[illumina_10] for Illumina sequencing reads with about 1% error rate"))
        .arg(Arg::with_name("train-file-dir")
            .short("r")
            .long("train-file-dir")
            .value_name("train_file_dir")
            .takes_value(true)
            .help("Full path of the directory containing the training model files."))
        .arg(Arg::with_name("thread-num")
            .short("p")
            .long("thread-num")
            .value_name("thread_num")
            .takes_value(true)
            .default_value("1")
            .help("The number of threads used by FragGeneScan++."))
        .arg(Arg::with_name("metadata-file")
            .short("e")
            .long("metadata-file")
            .value_name("metadata_file")
            .takes_value(true)
            .help("Output metadata for sequences."))
        .arg(Arg::with_name("dna-file")
            .short("d")
            .long("dna-file")
            .value_name("dna_file")
            .takes_value(true)
            .help("Output predicted DNA reads."))
        .arg(Arg::with_name("chunk-size")
            .short("c")
            .long("chunk-size")
            .value_name("chunk_size")
            .takes_value(true)
            .default_value("1")
            .help("Number of sequences in a chunk, scales speed and memory usage."))
        .arg(Arg::with_name("translation-table")
            .short("x")
            .long("translation-table")
            .value_name("translation_table")
            .takes_value(true)
            .default_value("11")
            .help("Which translation table to use."))
        .get_matches();

    let (global, locals) = hmm::get_train_from_file(
        PathBuf::from(matches.value_of("train-dir").unwrap_or("train")),
        PathBuf::from(matches.value_of("train-file").unwrap()),
    )?;

    // thread_data_init in run_hmm.c
    // writeOutputFiles in run_hmm.c
    let stdin = io::stdin();
    let stdout = io::stdout();
    run(
        global,
        locals,
        &mut fasta::Reader::new(stdin.lock()),
        &mut stdout.lock(),
    )?;

    Ok(())
}

fn run<R: Read, W: Write>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: &mut fasta::Reader<R, StdPolicy>,
    outputstream: &mut W,
) -> Result<(), Box<dyn Error>> {
    while let Some(record) = inputseqs.next() {
        let record = record?;
        viterbi(&global, &locals, record, outputstream)?;
    }
    Ok(())
}

fn count_cg_content(seq: &[u8]) -> usize {
    let mut count = 0;
    for l in seq.iter() {
        if b"CcGg".contains(l) {
            count += 1;
        }
    }
    min(43, max(0, count * 100 / seq.len() - 26))
}

fn viterbi<W: Write>(
    _global: &hmm::Global,
    _locals: &Vec<hmm::Local>,
    record: fasta::RefRecord,
    outputstream: &mut W,
) -> Result<(), Box<dyn Error>> {
    let seq = record.full_seq();
    let _cg = count_cg_content(&seq);
    record.write_unchanged(&mut *outputstream)?; // TODO
    Ok(())
}
