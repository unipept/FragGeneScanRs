//! FragGeneScanRs executable
#![allow(non_snake_case)]

use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Mutex;

extern crate anyhow;
use anyhow::Result;

extern crate clap;
use clap::{App, Arg};

extern crate seq_io;
use seq_io::fasta;

extern crate rayon;
use rayon::iter::{ParallelBridge, ParallelIterator};

extern crate frag_gene_scan_rs;
use frag_gene_scan_rs::dna::{count_cg_content, Nuc};
use frag_gene_scan_rs::hmm;
use frag_gene_scan_rs::viterbi::viterbi;

fn main() -> Result<()> {
    let matches = App::new("FragGeneScanRs")
        .version("0.0.1")
        .author("Felix Van der Jeugt <felix.vanderjeugt@ugent.be>")
        .about("Scalable high-throughput short-read open reading frame prediction.")
        .arg(Arg::with_name("seq-file")
            .short("s")
            .long("seq-file-name")
            .value_name("seq_file_name")
            .takes_value(true)
            .default_value("stdin")
            .help("Sequence file name including the full path. Using 'stdin' (or not suplying this argument) reads from standard input."))
        .arg(Arg::with_name("output-file")
            .short("o")
            .long("output-file-name")
            .value_name("output_file_name")
            .takes_value(true)
            .help("Output metadata, proteins and dna to files with this base name. Using 'stdout' write the predicted AA reads to standard output."))
        .arg(Arg::with_name("complete")
            .short("w")
            .long("complete")
            .help("The input sequence has complete genomic sequences; not short sequence reads."))
        .arg(Arg::with_name("formatted")
            .short("f")
            .long("formatted")
            .help("Format the DNA output."))
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
        .arg(Arg::with_name("meta-file")
            .short("m")
            .long("meta-file")
            .value_name("meta_file")
            .takes_value(true)
            .help("Output metadata to this file (supersedes -o)."))
        .arg(Arg::with_name("aa-file")
            .short("e")
            .long("aa-file")
            .value_name("aa_file")
            .takes_value(true)
            .help("Output predicted AA reads to this file (supersedes -o)."))
        .arg(Arg::with_name("dna-file")
            .short("d")
            .long("dna-file")
            .value_name("dna_file")
            .takes_value(true)
            .help("Output predicted DNA reads to this file (supersedes -o)."))
        // .arg(Arg::with_name("chunk-size")
        //     .short("c")
        //     .long("chunk-size")
        //     .value_name("chunk_size")
        //     .takes_value(true)
        //     .default_value("1")
        //     .help("Number of sequences in a chunk, scales speed and memory usage."))
        // .arg(Arg::with_name("translation-table")
        //     .short("x")
        //     .long("translation-table")
        //     .value_name("translation_table")
        //     .takes_value(true)
        //     .default_value("11")
        //     .help("Which translation table to use."))
        .get_matches();

    // thread_data_init in run_hmm.c
    // writeOutputFiles in run_hmm.c

    let (global, locals) = hmm::get_train_from_file(
        PathBuf::from(matches.value_of("train-file-dir").unwrap_or("train")),
        PathBuf::from(matches.value_of("train-file").unwrap()),
    )?;

    let inputseqs: Box<dyn Read + Send> = match matches.value_of("seq-file").unwrap_or("stdin") {
        "stdin" => Box::new(io::stdin()),
        filename => Box::new(File::open(filename)?),
    };

    let mut aastream: Option<Box<dyn Write + Send>> =
        match (matches.value_of("aa-file"), matches.value_of("output-file")) {
            (Some(filename), _) => Some(Box::new(File::create(filename)?)),
            (None, Some("stdout")) => Some(Box::new(io::stdout())),
            (None, Some(filename)) => Some(Box::new(File::create(filename.to_owned() + ".faa")?)),
            (None, None) => None,
        };

    let metastream: Option<File> = match (
        matches.value_of("meta-file"),
        matches.value_of("output-file"),
    ) {
        (Some(filename), _) => Some(File::create(filename)?),
        (None, Some("stdout")) => None,
        (None, Some(filename)) => Some(File::create(filename.to_owned() + ".out")?),
        (None, None) => None,
    };

    let dnastream: Option<File> = match (
        matches.value_of("dna-file"),
        matches.value_of("output-file"),
    ) {
        (Some(filename), _) => Some(File::create(filename)?),
        (None, Some("stdout")) => None,
        (None, Some(filename)) => Some(File::create(filename.to_owned() + ".ffn")?),
        (None, None) => None,
    };

    if aastream.is_none() && metastream.is_none() && dnastream.is_none() {
        aastream = Some(Box::new(io::stdout()));
    }

    run(
        global,
        locals,
        inputseqs,
        aastream,
        metastream,
        dnastream,
        matches.is_present("complete"),
        matches.is_present("formatted"),
        usize::from_str_radix(matches.value_of("thread-num").unwrap_or("1"), 10)?,
    )?;

    Ok(())
}

fn run<R: Read + Send, W: Write + Send>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: R,
    aastream: Option<W>,
    metastream: Option<File>,
    dnastream: Option<File>,
    whole_genome: bool,
    formatted: bool,
    _thread_num: usize,
) -> Result<()> {
    let aastream = aastream.map(Mutex::new);
    let metastream = metastream.map(Mutex::new);
    let dnastream = dnastream.map(Mutex::new);
    fasta::Reader::new(inputseqs)
        .into_records()
        .par_bridge()
        .map(|record| {
            let fasta::OwnedRecord { mut head, seq } = record?;
            head = head.into_iter().take_while(u8::is_ascii_graphic).collect();
            let nseq: Vec<Nuc> = seq.into_iter().map(Nuc::from).collect();
            let genes = viterbi(
                &global,
                &locals[count_cg_content(&nseq)],
                head,
                nseq,
                whole_genome,
            );
            for gene in genes {
                if let Some(metastream) = &metastream {
                    gene.print_meta(&mut *metastream.lock().unwrap())?;
                }
                if let Some(dnastream) = &dnastream {
                    gene.print_dna(&mut *dnastream.lock().unwrap(), formatted)?;
                }
                if let Some(aastream) = &aastream {
                    gene.print_protein(whole_genome, &mut *aastream.lock().unwrap())?;
                }
            }
            Ok(())
        })
        .collect()
}
