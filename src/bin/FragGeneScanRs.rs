//! FragGeneScanRs executable
#![allow(non_snake_case)]

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::marker::Sync;
use std::ops::DerefMut;
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
use rayon::ThreadPoolBuilder;

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
        .arg(Arg::with_name("output-prefix")
            .short("o")
            .long("output-prefix")
            .value_name("output_prefix")
            .takes_value(true)
            .help("Output metadata (.out), proteins (.faa) and genes (.ffn) to files with this prefix. Use 'stdout' to write the predicted proteins to standard output."))
        .arg(Arg::with_name("complete")
            .short("w")
            .long("complete")
            .value_name("complete")
            .takes_value(true)
            .default_value("0")
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
            .short("a")
            .long("aa-file")
            .value_name("aa_file")
            .takes_value(true)
            .help("Output predicted proteins to this file (supersedes -o)."))
        .arg(Arg::with_name("nucleotide-file")
            .short("n")
            .long("nucleotide-file")
            .value_name("nucleotide_file")
            .takes_value(true)
            .help("Output predicted genes to this file (supersedes -o)."))
        .get_matches();

    let (global, locals) = hmm::get_train_from_file(
        PathBuf::from(matches.value_of("train-file-dir").unwrap_or("train")),
        PathBuf::from(matches.value_of("train-file").unwrap()),
    )?;

    let inputseqs: Box<dyn Read + Send> = match matches.value_of("seq-file").unwrap_or("stdin") {
        "stdin" => Box::new(io::stdin()),
        filename => Box::new(File::open(filename)?),
    };

    let mut aastream: Option<Box<dyn Write + Send + Sync>> = match (
        matches.value_of("aa-file"),
        matches.value_of("output-prefix"),
    ) {
        (Some(filename), _) => Some(Box::new(File::create(filename)?)),
        (None, Some("stdout")) => Some(Box::new(io::stdout())),
        (None, Some(filename)) => Some(Box::new(File::create(filename.to_owned() + ".faa")?)),
        (None, None) => None,
    };

    let metastream: Option<File> = match (
        matches.value_of("meta-file"),
        matches.value_of("output-prefix"),
    ) {
        (Some(filename), _) => Some(File::create(filename)?),
        (None, Some("stdout")) => None,
        (None, Some(filename)) => Some(File::create(filename.to_owned() + ".out")?),
        (None, None) => None,
    };

    let dnastream: Option<File> = match (
        matches.value_of("nucleotide-file"),
        matches.value_of("output-prefix"),
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
        matches.value_of("complete").unwrap() == "1",
        matches.is_present("formatted"),
        usize::from_str_radix(matches.value_of("thread-num").unwrap(), 10)?,
    )?;

    Ok(())
}

fn run<R: Read + Send, W: Write + Send + Sync>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: R,
    aastream: Option<W>,
    metastream: Option<File>,
    dnastream: Option<File>,
    whole_genome: bool,
    formatted: bool,
    thread_num: usize,
) -> Result<()> {
    ThreadPoolBuilder::new()
        .num_threads(thread_num)
        .build_global()?;

    let hasmeta = metastream.is_some();
    let hasdna = metastream.is_some();
    let hasaa = metastream.is_some();
    let output = Mutex::new(OutputBuffer::new(aastream, metastream, dnastream));

    Chunked::new(100, fasta::Reader::new(inputseqs).into_records())
        .enumerate()
        .par_bridge()
        .map(|(index, recordvec)| {
            let mut metabuf = Vec::new();
            let mut dnabuf = Vec::new();
            let mut aabuf = Vec::new();
            for record in recordvec {
                let fasta::OwnedRecord { mut head, seq } = record?;
                head = head.into_iter().take_while(u8::is_ascii_graphic).collect();
                let nseq: Vec<Nuc> = seq.into_iter().map(Nuc::from).collect();
                let read_prediction = viterbi(
                    &global,
                    &locals[count_cg_content(&nseq)],
                    head,
                    nseq,
                    whole_genome,
                );
                if hasmeta {
                    read_prediction.meta(&mut metabuf)?;
                }
                if hasdna {
                    read_prediction.dna(&mut dnabuf, formatted)?;
                }
                if hasaa {
                    read_prediction.protein(&mut aabuf, whole_genome)?;
                }
            }
            let mut buffer = output.lock().unwrap();
            buffer.set(index, (metabuf, dnabuf, aabuf));
            let bufs = buffer.deref_mut().collect::<Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>>();
            for (metabuf, dnabuf, aabuf) in bufs {
                if let Some(metastream) = &mut buffer.metastream {
                    &mut metastream.write_all(&metabuf)?;
                }
                if let Some(dnastream) = &mut buffer.dnastream {
                    &mut dnastream.write_all(&dnabuf)?;
                }
                if let Some(aastream) = &mut buffer.aastream {
                    &mut aastream.write_all(&aabuf)?;
                }
            }
            Ok(())
        })
        .collect()
}

struct Chunked<I: Iterator> {
    size: usize,
    iterator: I,
}

impl<I: Iterator> Chunked<I> {
    fn new(size: usize, iterator: I) -> Self {
        Chunked { size, iterator }
    }
}

impl<I: Iterator> Iterator for Chunked<I> {
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut items = Vec::with_capacity(self.size);
        for _ in 0..self.size {
            if let Some(item) = self.iterator.next() {
                items.push(item);
            } else {
                break;
            }
        }
        if items.is_empty() {
            None
        } else {
            Some(items)
        }
    }
}

struct OutputBuffer<I, W: Write + Send> {
    next: usize,
    queue: VecDeque<Option<I>>,
    aastream: Option<W>,
    metastream: Option<File>,
    dnastream: Option<File>,
}

impl<I, W: Write + Send> OutputBuffer<I, W> {
    fn new(aastream: Option<W>, metastream: Option<File>, dnastream: Option<File>) -> Self {
        OutputBuffer {
            next: 0,
            queue: VecDeque::new(),
            aastream: aastream,
            metastream: metastream,
            dnastream: dnastream,
        }
    }

    fn set(&mut self, index: usize, item: I) {
        while self.next + self.queue.len() <= index {
            self.queue.push_back(None);
        }
        self.queue[index - self.next] = Some(item);
    }
}

impl <I, W: Write + Send> Iterator for OutputBuffer<I, W> {
    type Item = I;

    fn next(&mut self) -> Option<I> {
        if self.queue.front().map(Option::is_some).unwrap_or(false) {
            let item = self.queue.pop_front().unwrap().unwrap();
            self.next += 1;
            Some(item)
        } else {
            None
        }
    }
}
