//! FragGeneScanRs executable
#![allow(non_snake_case)]

use std::cmp::{max, min};
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;

extern crate clap;
use clap::{App, Arg};

extern crate seq_io;
use seq_io::fasta;
use seq_io::fasta::Record;

extern crate frag_gene_scan_rs;
use frag_gene_scan_rs::hmm;

fn main() -> Result<(), Box<dyn Error>> {
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

    let stdin = io::stdin();
    let mut inputseqs: Box<dyn Read> = match matches.value_of("seq-file").unwrap_or("stdin") {
        "stdin" => Box::new(stdin.lock()),
        filename => Box::new(File::open(filename)?),
    };

    let stdout = io::stdout();
    let mut aastream: Option<Box<dyn Write>> =
        match (matches.value_of("aa-file"), matches.value_of("output-file")) {
            (Some(filename), _) => Some(Box::new(File::create(filename)?)),
            (None, Some("stdout")) => Some(Box::new(stdout.lock())),
            (None, Some(filename)) => Some(Box::new(File::create(filename.to_owned() + ".faa")?)),
            (None, None) => None,
        };

    let mut metastream: Option<File> = match (
        matches.value_of("meta-file"),
        matches.value_of("output-file"),
    ) {
        (Some(filename), _) => Some(File::create(filename)?),
        (None, Some("stdout")) => None,
        (None, Some(filename)) => Some(File::create(filename.to_owned() + ".out")?),
        (None, None) => None,
    };

    let mut dnastream: Option<File> = match (
        matches.value_of("dna-file"),
        matches.value_of("output-file"),
    ) {
        (Some(filename), _) => Some(File::create(filename)?),
        (None, Some("stdout")) => None,
        (None, Some(filename)) => Some(File::create(filename.to_owned() + ".ffn")?),
        (None, None) => None,
    };

    if aastream.is_none() && metastream.is_none() && dnastream.is_none() {
        aastream = Some(Box::new(stdout.lock()));
    }

    run(
        global,
        locals,
        &mut inputseqs,
        &mut aastream,
        &mut metastream,
        &mut dnastream,
        matches.is_present("complete"),
        matches.is_present("formatted"),
        usize::from_str_radix(matches.value_of("thread-num").unwrap_or("1"), 10)?,
    )?;

    Ok(())
}

fn run<R: Read, W: Write>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: &mut R,
    aastream: &mut Option<W>,
    metastream: &mut Option<File>,
    dnastream: &mut Option<File>,
    whole_genome: bool,
    formatted: bool,
    _thread_num: usize,
) -> Result<(), Box<dyn Error>> {
    let mut sequences = fasta::Reader::new(inputseqs);
    while let Some(record) = sequences.next() {
        let record = record?;
        viterbi(
            &global,
            &locals,
            record,
            aastream,
            metastream,
            dnastream,
            whole_genome,
            formatted,
        )?;
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
    min(43, max(26, count * 100 / seq.len()) - 26)
}

fn viterbi<W: Write>(
    global: &hmm::Global,
    locals: &Vec<hmm::Local>,
    record: fasta::RefRecord,
    aastream: &mut Option<W>,
    metastream: &mut Option<File>,
    dnastream: &mut Option<File>,
    whole_genome: bool,
    formatted: bool,
) -> Result<(), Box<dyn Error>> {
    let fasta::OwnedRecord { mut head, seq } = record.to_owned_record();
    head.truncate(head.partition_point(u8::is_ascii_whitespace));
    let cg = count_cg_content(&seq);

    let gene_len = if whole_genome { 120 } else { 60 }; // minimum length to be output

    let mut alpha: Vec<[f64; hmm::NUM_STATE]> = vec![];
    let mut path: Vec<[usize; hmm::NUM_STATE]> = vec![];
    let mut temp_i: [usize; 6] = [0; 6];
    let mut temp_i_1: [usize; 6] = [0; 6];

    for _ in 0..seq.len() {
        alpha.push([0.0; hmm::NUM_STATE]);
        path.push([hmm::NOSTATE; hmm::NUM_STATE]);
    }
    alpha[0].copy_from_slice(&global.pi);
    for i in &mut alpha[0] {
        *i *= -1.0
    } // TODO broadcast operation

    // If the sequence starts with a stop codon
    if seq[0] == b'T'
        && ((seq[1] == b'A' && seq[2] == b'A')
            || (seq[1] == b'A' && seq[2] == b'G')
            || (seq[1] == b'G' && seq[2] == b'A'))
    {
        alpha[0][hmm::E_STATE] = f64::INFINITY;
        alpha[1][hmm::E_STATE] = f64::INFINITY;
        path[1][hmm::E_STATE] = hmm::E_STATE;
        path[2][hmm::E_STATE] = hmm::E_STATE;

        alpha[2][hmm::M6_STATE] = f64::INFINITY;
        alpha[1][hmm::M5_STATE] = f64::INFINITY;
        alpha[0][hmm::M4_STATE] = f64::INFINITY;
        alpha[2][hmm::M3_STATE] = f64::INFINITY;
        alpha[1][hmm::M2_STATE] = f64::INFINITY;
        alpha[0][hmm::M1_STATE] = f64::INFINITY;

        alpha[2][hmm::E_STATE] -= if seq[1] == b'A' && seq[2] == b'A' {
            0.53_f64.ln()
        } else if seq[1] == b'A' && seq[2] == b'G' {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    // If the sequence starts with a reverse stop codon
    if seq[2] == b'A'
        && ((seq[1] == b'T' && seq[0] == b'T')
            || (seq[1] == b'T' && seq[0] == b'C')
            || (seq[1] == b'C' && seq[0] == b'T'))
    {
        alpha[0][hmm::S_STATE_1] = f64::INFINITY;
        alpha[1][hmm::S_STATE_1] = f64::INFINITY;
        alpha[2][hmm::S_STATE_1] = alpha[0][hmm::S_STATE];
        path[1][hmm::S_STATE_1] = hmm::S_STATE_1;
        path[2][hmm::S_STATE_1] = hmm::S_STATE_1;

        alpha[2][hmm::M3_STATE_1] = f64::INFINITY;
        alpha[2][hmm::M6_STATE_1] = f64::INFINITY;

        alpha[2][hmm::S_STATE_1] = if seq[1] == b'T' && seq[0] == b'T' {
            0.53_f64.ln()
        } else if seq[1] == b'T' && seq[0] == b'C' {
            0.16_f64.ln()
        } else {
            0.30_f64.ln()
        }
    }

    let mut num_N = 0; // number of invalid nts in sequence
    for t in 1..seq.len() {
        let from = nt2int(seq[t - 1]).unwrap_or(2);
        let from0 = if t > 1 {
            nt2int(seq[t - 2]).unwrap_or(2)
        } else {
            2
        };
        let to = nt2int(seq[t]).unwrap_or_else(|| {
            num_N += 1;
            2
        });
        let from2 = from0 * 4 + from;

        // M state

        for i in hmm::M1_STATE..=hmm::M6_STATE {
            if alpha[t][i].is_finite() {
                if i == hmm::M1_STATE {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::M6_STATE]
                        - global.tr[hmm::TR_GG]
                        - global.tr[hmm::TR_MM]
                        - locals[cg].e_m[0][from2][to];
                    path[t][i] = hmm::M6_STATE;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE..=hmm::M5_STATE).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr[hmm::TR_MD]
                                    - locals[cg].e_m[0][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr[hmm::TR_DD] * (num_d - 2) as f64
                                    - global.tr[hmm::TR_DM];
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }

                    // from Start state
                    let temp_alpha = alpha[t - 1][hmm::S_STATE] - locals[cg].e_m[0][from2][to];
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = hmm::S_STATE;
                    }
                } else {
                    // from M state
                    alpha[t][i] = alpha[t - 1][i - 1]
                        - global.tr[hmm::TR_MM]
                        - locals[cg].e_m[i - hmm::M1_STATE][from2][to];
                    path[t][i] = i - 1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE..=hmm::M6_STATE).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr[hmm::TR_MD]
                                    - locals[cg].e_m[i - hmm::M1_STATE][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr[hmm::TR_DD] * (num_d - 2) as f64
                                    - global.tr[hmm::TR_DM];
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                }

                // from I state (l251)
                let j = if i == hmm::M1_STATE {
                    hmm::I6_STATE
                } else {
                    hmm::I1_STATE + (i - hmm::M1_STATE - 1)
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::M2_STATE || i == hmm::M5_STATE)
                    && t + 1 < seq.len()
                    && seq[temp_i[j - hmm::I1_STATE]] == b'T'
                    && ((seq[t] == b'A' && seq[t + 1] == b'A')
                        || (seq[t] == b'A' && seq[t + 1] == b'G')
                        || (seq[t] == b'G' && seq[t + 1] == b'A'))
                {
                } else if (i == hmm::M3_STATE || i == hmm::M6_STATE)
                    && temp_i[j - hmm::I1_STATE] > 0
                    && (seq[temp_i[j - hmm::I1_STATE] - 1] == b'T')
                    && ((seq[temp_i[j - hmm::I1_STATE]] == b'A' && seq[t] == b'A')
                        || (seq[temp_i[j - hmm::I1_STATE]] == b'A' && seq[t] == b'G')
                        || (seq[temp_i[j - hmm::I1_STATE]] == b'G' && seq[t] == b'A'))
                {
                } else {
                    let temp_alpha = alpha[t - 1][j] - global.tr[hmm::TR_IM] - 0.25_f64.ln();
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = j;
                    }
                }
            }
        }

        // I state
        for i in hmm::I1_STATE..=hmm::I6_STATE {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr[hmm::TR_II] - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            let temp_alpha = alpha[t - 1][i - hmm::I1_STATE + hmm::M1_STATE]
                - global.tr[hmm::TR_MI]
                - global.tr_mi[from][to]
                - if i == hmm::I6_STATE {
                    global.tr[hmm::TR_GG]
                } else {
                    0.0
                };
            if temp_alpha < alpha[t][i] {
                alpha[t][i] = temp_alpha;
                path[t][i] = i - hmm::I1_STATE + hmm::M1_STATE;
                temp_i[i - hmm::I1_STATE] = t - 1;
            }
        }

        // M' state
        for i in hmm::M1_STATE_1..=hmm::M6_STATE_1 {
            if (i == hmm::M1_STATE_1 || i == hmm::M4_STATE_1)
                && t >= 3
                && seq[t - 1] == b'A'
                && ((seq[t - 2] == b'T' && seq[t - 3] == b'T')
                    || (seq[t - 2] == b'T' && seq[t - 3] == b'C')
                    || (seq[t - 2] == b'C' && seq[t - 3] == b'T'))
            {
                // from Start state since this is actually a stop codon in minus strand
                alpha[t][i] =
                    alpha[t - 1][hmm::S_STATE_1] - locals[cg].e_m1[i - hmm::M1_STATE_1][from2][to];
                path[t][i] = hmm::S_STATE_1;
            } else {
                if i == hmm::M1_STATE_1 {
                    // from M state
                    alpha[t][i] = alpha[t - 1][hmm::M6_STATE_1]
                        - global.tr[hmm::TR_GG]
                        - global.tr[hmm::TR_MM]
                        - locals[cg].e_m1[0][from2][to];
                    path[t][i] = hmm::M6_STATE_1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE_1..=hmm::M5_STATE_1).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j] - global.tr[hmm::TR_MD]
                                               - locals[cg].e_m1[0][from2][to] // TODO different from forward but merge?
                                               - 0.25_f64.ln() * (num_d - 1) as f64
                                               - global.tr[hmm::TR_DD] * (num_d - 2) as f64
                                               - global.tr[hmm::TR_DM];
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                } else {
                    // from M state
                    alpha[t][i] = alpha[t - 1][i - 1]
                        - global.tr[hmm::TR_MM]
                        - locals[cg].e_m1[i - hmm::M1_STATE_1][from2][to];
                    path[t][i] = i - 1;

                    // from D state
                    if !whole_genome {
                        for j in (hmm::M1_STATE_1..=hmm::M6_STATE_1).rev() {
                            let num_d = if j >= i {
                                (i + 6 - j) as i32
                            } else if j + 1 < i {
                                (i - j) as i32
                            } else {
                                -10
                            };
                            if num_d > 0 {
                                let temp_alpha = alpha[t - 1][j]
                                    - global.tr[hmm::TR_MD]
                                    - locals[cg].e_m1[i - hmm::M1_STATE_1][from2][to]
                                    - 0.25_f64.ln() * (num_d - 1) as f64
                                    - global.tr[hmm::TR_DD] * (num_d - 2) as f64
                                    - global.tr[hmm::TR_DM];
                                if temp_alpha < alpha[t][i] {
                                    alpha[t][i] = temp_alpha;
                                    path[t][i] = j;
                                }
                            }
                        }
                    }
                }

                // from I state
                let j = if i == hmm::M1_STATE_1 {
                    hmm::I6_STATE_1
                } else {
                    hmm::I1_STATE_1 + i - hmm::M1_STATE_1 - 1
                };

                // to avoid stop codon
                if t < 2 {
                } else if (i == hmm::M2_STATE_1 || i == hmm::M5_STATE_1)
                    && t + 1 < seq.len()
                    && seq[t + 1] == b'A'
                    && ((seq[t] == b'T' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T')
                        || (seq[t] == b'T' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'C')
                        || (seq[t] == b'A' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T'))
                {
                } else if (i == hmm::M3_STATE_1 || i == hmm::M6_STATE_1)
                    && seq[t] == b'A'
                    && temp_i_1[j - hmm::I1_STATE_1] > 1
                    && ((seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T'
                        && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == b'T')
                        || (seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T'
                            && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == b'C')
                        || (seq[temp_i_1[j - hmm::I1_STATE_1]] == b'C'
                            && seq[temp_i_1[j - hmm::I1_STATE_1] - 1] == b'T'))
                {
                } else {
                    let temp_alpha = alpha[t - 1][j] - global.tr[hmm::TR_IM] - 0.25_f64.ln();
                    if temp_alpha < alpha[t][i] {
                        alpha[t][i] = temp_alpha;
                        path[t][i] = j;
                    }
                }
            }
        }

        // I' state
        for i in hmm::I1_STATE_1..=hmm::I6_STATE_1 {
            // from I state
            alpha[t][i] = alpha[t - 1][i] - global.tr[hmm::TR_II] - global.tr_ii[from][to];
            path[t][i] = i;

            // from M state
            if (t >= 3 && path[t - 3][hmm::S_STATE_1] != hmm::R_STATE)
                && (t >= 4 && path[t - 4][hmm::S_STATE_1] != hmm::R_STATE)
                && (t >= 5 && path[t - 5][hmm::S_STATE_1] != hmm::R_STATE)
            {
                let temp_alpha = alpha[t - 1][i - hmm::I1_STATE_1 + hmm::M1_STATE_1]
                    - global.tr[hmm::TR_MI]
                    - global.tr_mi[from][to]
                    - if i == hmm::I6_STATE_1 {
                        global.tr[hmm::TR_GG]
                    } else {
                        0.0
                    };
                if temp_alpha < alpha[t][i] {
                    alpha[t][i] = temp_alpha;
                    path[t][i] = i - hmm::I1_STATE_1 + hmm::M1_STATE_1;
                    temp_i_1[i - hmm::I1_STATE_1] = t - 1;
                }
            }
        }

        // non_coding state
        // TODO just a minimum of three
        alpha[t][hmm::R_STATE] =
            alpha[t - 1][hmm::R_STATE] - locals[cg].tr_rr[from][to] - global.tr[hmm::TR_RR];
        path[t][hmm::R_STATE] = hmm::R_STATE;

        let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr[hmm::TR_ER];
        if temp_alpha < alpha[t][hmm::R_STATE] {
            alpha[t][hmm::R_STATE] = temp_alpha;
            path[t][hmm::R_STATE] = hmm::E_STATE;
        }

        let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr[hmm::TR_ER];
        if temp_alpha < alpha[t][hmm::R_STATE] {
            alpha[t][hmm::R_STATE] = temp_alpha;
            path[t][hmm::R_STATE] = hmm::E_STATE_1;
        }

        alpha[t][hmm::R_STATE] -= 0.95_f64.ln();

        // end state
        if alpha[t][hmm::E_STATE] == 0.0 {
            alpha[t][hmm::E_STATE] = f64::INFINITY;
            path[t][hmm::E_STATE] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t] == b'T'
                && ((seq[t + 1] == b'A' && seq[t + 2] == b'A')
                    || (seq[t + 1] == b'A' && seq[t + 2] == b'G')
                    || (seq[t + 1] == b'G' && seq[t + 2] == b'A'))
            {
                alpha[t + 2][hmm::E_STATE] = f64::INFINITY;

                // transition from frame4, frame5 and frame6
                let temp_alpha = alpha[t - 1][hmm::M6_STATE] - global.tr[hmm::TR_GE];
                if temp_alpha < alpha[t + 2][hmm::E_STATE] {
                    alpha[t + 2][hmm::E_STATE] = temp_alpha;
                    path[t][hmm::E_STATE] = hmm::M6_STATE;
                }

                // transition from frame1, frame2 and frame3
                let temp_alpha = alpha[t - 1][hmm::M3_STATE] - global.tr[hmm::TR_GE];
                if temp_alpha < alpha[t + 2][hmm::E_STATE] {
                    alpha[t + 2][hmm::E_STATE] = temp_alpha;
                    path[t][hmm::E_STATE] = hmm::M3_STATE;
                }

                alpha[t][hmm::E_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::E_STATE] = f64::INFINITY;
                path[t + 1][hmm::E_STATE] = hmm::E_STATE;
                path[t + 2][hmm::E_STATE] = hmm::E_STATE;

                alpha[t + 2][hmm::M6_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::M5_STATE] = f64::INFINITY;
                alpha[t][hmm::M4_STATE] = f64::INFINITY;
                alpha[t + 2][hmm::M3_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::M2_STATE] = f64::INFINITY;
                alpha[t][hmm::M1_STATE] = f64::INFINITY;

                alpha[t + 2][hmm::E_STATE] -= if seq[t + 1] == b'A' && seq[t + 2] == b'A' {
                    0.54_f64.ln()
                } else if seq[t + 1] == b'A' && seq[t + 2] == b'G' {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 60 {
                    // bug reported by Yu-Wei ------ TODO 60 is the incomplete minimum length? can be merged?
                    for i in (t - 60)..=(t - 3) {
                        if i + 2 < seq.len() {
                            start_freq -= locals[cg].tr_e[i + 60 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else if t > 3 {
                    for i in 0..=(t - 3) {
                        if i + 2 < seq.len() {
                            sub_sum += locals[cg].tr_e[i + 60 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                    sub_sum = sub_sum * 58.0 / (t - 3 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = locals[cg].dist_e[2]
                    * (-1.0 * (start_freq - locals[cg].dist_e[1]).powi(2)
                        / (locals[cg].dist_e[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = locals[cg].dist_e[5]
                    * (-1.0 * (start_freq - locals[cg].dist_e[4]).powi(2)
                        / (locals[cg].dist_e[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::E_STATE] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start' state
        // originally stop codon of genes in - strand
        if alpha[t][hmm::S_STATE_1] == 0.0 {
            alpha[t][hmm::S_STATE_1] = f64::INFINITY;
            path[t][hmm::S_STATE_1] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t + 2] == b'A'
                && ((seq[t + 1] == b'T' && seq[t] == b'T')
                    || (seq[t + 1] == b'T' && seq[t] == b'C')
                    || (seq[t + 1] == b'C' && seq[t] == b'T'))
            {
                alpha[t][hmm::S_STATE_1] = f64::INFINITY;
                alpha[t + 1][hmm::S_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::S_STATE_1] = alpha[t - 1][hmm::R_STATE] - global.tr[hmm::TR_RS];
                path[t][hmm::S_STATE_1] = hmm::R_STATE;
                path[t + 1][hmm::S_STATE_1] = hmm::S_STATE_1;
                path[t + 2][hmm::S_STATE_1] = hmm::S_STATE_1;

                let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr[hmm::TR_ES];
                if temp_alpha < alpha[t + 2][hmm::S_STATE_1] {
                    alpha[t + 2][hmm::S_STATE_1] = temp_alpha;
                    path[t][hmm::S_STATE_1] = hmm::E_STATE_1;
                }

                let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr[hmm::TR_ES1];
                if temp_alpha < alpha[t + 2][hmm::S_STATE_1] {
                    alpha[t + 2][hmm::S_STATE_1] = temp_alpha;
                    path[t][hmm::S_STATE_1] = hmm::E_STATE;
                }

                alpha[t + 2][hmm::M3_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::M6_STATE_1] = f64::INFINITY;

                alpha[t + 2][hmm::S_STATE_1] -= if seq[t + 1] == b'T' && seq[t] == b'T' {
                    0.54_f64.ln()
                } else if seq[t + 1] == b'T' && seq[t] == b'C' {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;

                // TODO needs same 60-edgecase as above?
                for i in 3..=60 {
                    if t + i + 2 < seq.len() {
                        start_freq += locals[cg].tr_s1[i - 3]
                            [trinucleotide(seq[t + i], seq[t + i + 1], seq[t + i + 2])];
                    }
                }

                let h_kd = locals[cg].dist_s1[2]
                    * (-1.0 * (start_freq - locals[cg].dist_s1[1]).powi(2)
                        / (locals[cg].dist_s1[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = locals[cg].dist_s1[5]
                    * (-1.0 * (start_freq - locals[cg].dist_s1[4]).powi(2)
                        / (locals[cg].dist_s1[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::S_STATE_1] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // start state
        if alpha[t][hmm::S_STATE] == 0.0 {
            alpha[t][hmm::S_STATE] = f64::INFINITY;
            path[t][hmm::S_STATE] = hmm::NOSTATE;

            if t < seq.len() - 2
                && (seq[t] == b'A' || seq[t] == b'G' || seq[t] == b'T')
                && seq[t + 1] == b'A'
                && seq[t + 2] == b'G'
            {
                alpha[t][hmm::S_STATE] = f64::INFINITY;
                alpha[t + 1][hmm::S_STATE] = f64::INFINITY;
                alpha[t + 2][hmm::S_STATE] = alpha[t - 1][hmm::R_STATE] - global.tr[hmm::TR_RS];
                path[t][hmm::S_STATE] = hmm::R_STATE;
                path[t + 1][hmm::S_STATE] = hmm::S_STATE;
                path[t + 2][hmm::S_STATE] = hmm::S_STATE;

                let temp_alpha = alpha[t - 1][hmm::E_STATE] - global.tr[hmm::TR_ES];
                if temp_alpha < alpha[t + 2][hmm::S_STATE] {
                    alpha[t + 2][hmm::S_STATE] = temp_alpha;
                    path[t][hmm::S_STATE] = hmm::E_STATE;
                }

                let temp_alpha = alpha[t - 1][hmm::E_STATE_1] - global.tr[hmm::TR_ES1];
                if temp_alpha < alpha[t + 2][hmm::S_STATE] {
                    alpha[t + 2][hmm::S_STATE] = temp_alpha;
                    path[t][hmm::S_STATE] = hmm::E_STATE_1;
                }

                alpha[t + 2][hmm::S_STATE] -= if seq[t] == b'A' {
                    0.83_f64.ln()
                } else if seq[t] == b'G' {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 30 {
                    // TODO why 30?
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -= locals[cg].tr_s[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum += locals[cg].tr_s[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                    sub_sum *= 61.0 / (t + 30 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = locals[cg].dist_s[2]
                    * (-1.0 * (start_freq - locals[cg].dist_s[1]).powi(2)
                        / (locals[cg].dist_s[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = locals[cg].dist_s[5]
                    * (-1.0 * (start_freq - locals[cg].dist_s[4]).powi(2)
                        / (locals[cg].dist_s[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::S_STATE] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        // end' state
        // originally start codon of genes in - strand
        if alpha[t][hmm::E_STATE_1] == 0.0 {
            alpha[t][hmm::E_STATE_1] = f64::INFINITY;
            path[t][hmm::E_STATE_1] = hmm::NOSTATE;

            if t < seq.len() - 2
                && seq[t] == b'C'
                && seq[t + 1] == b'A'
                && (seq[t + 2] == b'T' || seq[t + 2] == b'C' || seq[t + 2] == b'A')
            {
                // transition from frame6
                alpha[t][hmm::E_STATE_1] = f64::INFINITY;
                alpha[t + 1][hmm::E_STATE_1] = f64::INFINITY;
                alpha[t + 2][hmm::E_STATE_1] =
                    alpha[t - 1][hmm::M6_STATE_1] - global.tr[hmm::TR_GE];
                path[t][hmm::E_STATE_1] = hmm::M6_STATE_1;
                path[t + 1][hmm::E_STATE_1] = hmm::E_STATE_1;
                path[t + 2][hmm::E_STATE_1] = hmm::E_STATE_1;

                alpha[t + 2][hmm::E_STATE_1] -= if seq[t + 2] == b'T' {
                    0.83_f64.ln()
                } else if seq[t + 2] == b'C' {
                    0.10_f64.ln()
                } else {
                    0.07_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut sub_sum = 0.0;

                if t >= 30 {
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -= locals[cg].tr_e1[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum += locals[cg].tr_e1[i + 30 - t]
                                [trinucleotide(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                    sub_sum *= 61.0 / (t + 30 + 1) as f64;
                    start_freq -= sub_sum;
                }

                let h_kd = locals[cg].dist_e1[2]
                    * (-1.0 * (start_freq - locals[cg].dist_e1[1]).powi(2)
                        / (locals[cg].dist_e1[0]).powi(2)
                        / 2.0)
                        .exp();
                let r_kd = locals[cg].dist_e1[5]
                    * (-1.0 * (start_freq - locals[cg].dist_e1[4]).powi(2)
                        / (locals[cg].dist_e1[3]).powi(2)
                        / 2.0)
                        .exp();
                alpha[t + 2][hmm::E_STATE_1] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln();
            }
        }

        if num_N > 9 {
            for i in 0..hmm::NUM_STATE {
                if i != hmm::R_STATE {
                    alpha[t][i] = f64::INFINITY;
                    path[t][i] = hmm::R_STATE;
                }
            }
        }
    }

    let vpath = backtrack(&alpha, path);
    output(
        &locals[cg],
        head,
        seq,
        aastream,
        metastream,
        dnastream,
        whole_genome,
        formatted,
        vpath,
        gene_len,
        alpha,
    )
}

fn backtrack(alpha: &Vec<[f64; hmm::NUM_STATE]>, path: Vec<[usize; hmm::NUM_STATE]>) -> Vec<usize> {
    // backtrack array to find the optimal path
    let mut vpath: Vec<usize> = vec![0];
    let mut prob = f64::INFINITY;
    for (i, &prob_) in alpha.last().expect("empty seq").iter().enumerate() {
        if prob_ < prob {
            vpath[0] = i;
            prob = prob_;
        }
    }

    // backtrack the optimal path
    for t in (0..=path.len() - 2).rev() {
        vpath.push(path[t + 1][*vpath.last().unwrap()]);
    }
    vpath.reverse();
    vpath
}

fn output<W: Write>(
    local: &hmm::Local,
    head: Vec<u8>,
    seq: Vec<u8>,
    aastream: &mut Option<W>,
    metastream: &mut Option<File>,
    dnastream: &mut Option<File>,
    whole_genome: bool,
    formatted: bool,
    vpath: Vec<usize>,
    gene_len: usize,
    alpha: Vec<[f64; hmm::NUM_STATE]>,
) -> Result<(), Box<dyn Error>> {
    let mut codon_start = 0; // ternaire boolean?
    let mut start_t: isize = -1;
    let mut dna_start_t_withstop: usize = 0;
    let mut dna_start_t: usize = 0;

    let mut dna = vec![];
    let mut _dna1: [u8; 300000] = [0; 300000];
    let mut dna_f = vec![];
    let mut _dna_f1: [u8; 300000] = [0; 300000];
    let mut insert = vec![];
    let mut delete = vec![];

    let mut prev_match = 0;

    let mut start_orf = 0; // initialize?

    for t in 0..seq.len() {
        if codon_start == 0
            && start_t < 0
            && ((vpath[t] >= hmm::M1_STATE && vpath[t] <= hmm::M6_STATE)
                || (vpath[t] >= hmm::M1_STATE_1 && vpath[t] <= hmm::M6_STATE_1)
                || vpath[t] == hmm::S_STATE
                || vpath[t] == hmm::S_STATE_1)
        {
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            start_t = (t + 1) as isize;
            // introduce dna_start_t_withstop YY July 2018
        }

        if codon_start == 0
            && (vpath[t] == hmm::M1_STATE
                || vpath[t] == hmm::M4_STATE
                || vpath[t] == hmm::M1_STATE_1
                || vpath[t] == hmm::M4_STATE_1)
        {
            dna.clear();
            _dna1 = [0; 300000];
            dna_f.clear();
            _dna_f1 = [0; 300000];
            insert.clear();
            delete.clear();

            dna.push(seq[t]);
            dna_f.push(seq[t]);
            dna_start_t_withstop = t + 1;
            dna_start_t = t + 1;
            if vpath[t] == hmm::M1_STATE || vpath[t] == hmm::M4_STATE_1 {
                if t > 2 {
                    dna_start_t_withstop = t - 2;
                }
            }

            start_orf = t + 1;
            prev_match = vpath[t];

            codon_start = if vpath[t] < hmm::M6_STATE { 1 } else { -1 }
        } else if codon_start != 0
            && (vpath[t] == hmm::E_STATE || vpath[t] == hmm::E_STATE_1 || t == seq.len() - 1)
        {
            let mut end_t;
            if vpath[t] == hmm::E_STATE || vpath[t] == hmm::E_STATE_1 {
                end_t = t + 3
            } else {
                let mut temp_t = t;
                while vpath[temp_t] != hmm::M1_STATE
                    && vpath[temp_t] != hmm::M4_STATE
                    && vpath[temp_t] != hmm::M1_STATE_1
                    && vpath[temp_t] != hmm::M4_STATE_1
                {
                    dna.pop();
                    dna_f.pop();
                    temp_t -= 1;
                }
                end_t = temp_t;
            }

            if dna.len() > gene_len {
                let final_score = (alpha[end_t - 4][vpath[end_t - 4]]
                    - alpha[start_t as usize + 2][vpath[start_t as usize + 2]])
                    / (end_t - start_t as usize - 5) as f64;
                let mut frame = start_orf % 3;
                if frame == 0 {
                    frame = 3
                }

                if codon_start == 1 {
                    if start_t == dna_start_t as isize - 3 {
                        dna_start_t -= 3;
                    }

                    if whole_genome {
                        // add refinement of the start codons here
                        let start_old = start_t as usize;
                        let mut codon = &seq[start_old..start_old + 3];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon != b"TAA" || codon != b"TAG" || codon != b"TGA")
                            && start_old >= 1 + s + 35
                        {
                            if codon != b"ATG" || codon != b"GTG" || codon != b"TTG" {
                                let utr = &seq[start_old - 1 - s - 30..start_old - 1 - s - 30 + 63];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    let idx = trinucleotide(utr[j], utr[j + 1], utr[j + 2]);
                                    freq_sum -= local.tr_s[j][idx];
                                }
                                if s == 0 {
                                    e_save = freq_sum;
                                    s_save = 0;
                                } else if freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = -(s as isize); // posivite chain, upstream s_save = -1 * 3
                                }
                            }
                            s += 3;
                            codon = &seq[start_old - 1 - s..start_old - 1 - s + 3];

                            dna_start_t += s_save as usize;
                        }
                    }

                    if let Some(metastream) = metastream {
                        fasta::OwnedRecord {
                            head: head.clone(),
                            seq: format!(
                                "{}\t{}\t+\t{}\t{}\tI:{}\tD:{}",
                                dna_start_t,
                                end_t,
                                frame,
                                final_score,
                                insert
                                    .iter()
                                    .map(|i: &usize| { format!("{},", i) })
                                    .collect::<String>(),
                                delete
                                    .iter()
                                    .map(|i: &usize| { format!("{},", i) })
                                    .collect::<String>()
                            )
                            .into_bytes(),
                        }
                        .write(&mut *metastream)?;
                    }

                    // dna = seq[dna_start_t - 1..end_t].to_vec();
                    if let Some(aastream) = aastream {
                        let mut infohead = head.clone();
                        infohead.append(&mut format!("_{}_{}_+", dna_start_t, end_t).into_bytes());
                        print_protein(&dna, true, whole_genome, infohead, aastream)?;
                    }
                    if let Some(dnastream) = dnastream {
                        let mut infohead = head.clone();
                        infohead.append(&mut format!("_{}_{}_+", dna_start_t, end_t).into_bytes());
                        print_dna(if formatted { &dna_f } else { &dna }, &infohead, dnastream)?;
                    }
                } else if codon_start == -1 {
                    if whole_genome {
                        // add refinement of the start codons here
                        // l 946
                        let end_old = end_t;
                        let mut codon = &seq[end_t - 1 - 2..end_t];
                        let mut s = 0;
                        // find the optimal start codon within 30bp up- and downstream of start codon
                        let mut e_save = 0.0;
                        let mut s_save = 0;
                        while !(codon != b"TAA" || codon != b"TAG" || codon != b"TGA")
                            && end_old - 2 + s + 35 < seq.len()
                        {
                            if codon != b"ATG" || codon != b"GTG" || codon != b"TTG" {
                                let utr = &seq[end_old - 1 - 2 + s - 30..end_old + s + 30];
                                let mut freq_sum = 0.0;
                                for j in 0..utr.len() - 2 {
                                    let idx = trinucleotide(utr[j], utr[j + 1], utr[j + 2]);
                                    freq_sum -= local.tr_e1[j][idx]; // TODO stop1? (their note)
                                }
                                if s == 0 || freq_sum < e_save {
                                    e_save = freq_sum;
                                    s_save = s; // negative chain, s_save = s
                                }
                            }
                            s += 3;
                            codon = &seq[end_old - 1 - 2 + s..end_old];
                        }

                        end_t = end_old + s_save;
                    }

                    if let Some(metastream) = metastream {
                        fasta::OwnedRecord {
                            head: head.clone(),
                            seq: format!(
                                "{}\t{}\t-\t{}\t{}\tI:{}\tD:{}",
                                dna_start_t,
                                end_t,
                                frame,
                                final_score,
                                insert
                                    .iter()
                                    .map(|i: &usize| { format!("{},", i) })
                                    .collect::<String>(),
                                delete
                                    .iter()
                                    .map(|i: &usize| { format!("{},", i) })
                                    .collect::<String>()
                            )
                            .into_bytes(),
                        }
                        .write(&mut *metastream)?;
                    }

                    //dna = seq[dna_start_t_withstop - 1..end_t].to_vec();
                    if let Some(aastream) = aastream {
                        let mut infohead = head.clone();
                        infohead.append(
                            &mut format!("_{}_{}_-", dna_start_t_withstop, end_t).into_bytes(),
                        );
                        print_protein(&dna, false, whole_genome, infohead, aastream)?;
                    }
                    if let Some(dnastream) = dnastream {
                        let mut infohead = head.clone();
                        infohead.append(
                            &mut format!("_{}_{}_-", dna_start_t_withstop, end_t).into_bytes(),
                        );
                        print_dna(
                            &get_rc_dna(if formatted { &dna_f } else { &dna }),
                            &infohead,
                            dnastream,
                        )?;
                    }
                }
            }

            codon_start = 0;
            start_t = -1;
        } else if codon_start != 0
            && ((vpath[t] >= hmm::M1_STATE && vpath[t] <= hmm::M6_STATE)
                || (vpath[t] >= hmm::M1_STATE_1 && vpath[t] <= hmm::M6_STATE_1))
            && vpath[t] < prev_match + 6
        {
            let out_nt = if vpath[t] < prev_match {
                vpath[t] + 6 - prev_match
            } else {
                vpath[t] - prev_match
            };
            for kk in 0..out_nt {
                // for deleted nt in reads
                dna.push(b'N');
                dna_f.push(b'x');
                if kk > 0 {
                    delete.push(t + 1);
                }
            }
            dna.pop();
            dna_f.pop();
            dna.push(seq[t]);
            dna_f.push(seq[t]);
            prev_match = vpath[t];
        } else if codon_start != 0
            && ((vpath[t] >= hmm::I1_STATE && vpath[t] <= hmm::I6_STATE)
                || (vpath[t] >= hmm::I1_STATE_1 && vpath[t] <= hmm::I6_STATE_1))
        {
            dna_f.push(seq[t].to_ascii_lowercase());
            insert.push(t + 1);
        } else if codon_start != 0 && vpath[t] == hmm::R_STATE {
            // for long NNNNNNNN, pretend R state
            codon_start = 0;
            start_t = -1;
        }
    }

    Ok(())
}

fn nt2int(nt: u8) -> Option<usize> {
    match nt {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn trinucleotide(a: u8, b: u8, c: u8) -> usize {
    if let (Some(a_), Some(b_), Some(c_)) = (nt2int(a), nt2int(b), nt2int(c)) {
        16 * a_ + 4 * b_ + c_
    } else {
        0
    }
}

fn trinucleotide_pep(a: u8, b: u8, c: u8) -> usize {
    if let (Some(a_), Some(b_), Some(c_)) = (nt2int(a), nt2int(b), nt2int(c)) {
        16 * a_ + 4 * b_ + c_
    } else {
        64
    }
}

const CODON_CODE: [u8; 65] = [
    b'K', b'N', b'K', b'N', b'T', b'T', b'T', b'T', b'R', b'S', b'R', b'S', b'I', b'I', b'M', b'I',
    b'Q', b'H', b'Q', b'H', b'P', b'P', b'P', b'P', b'R', b'R', b'R', b'R', b'L', b'L', b'L', b'L',
    b'E', b'D', b'E', b'D', b'A', b'A', b'A', b'A', b'G', b'G', b'G', b'G', b'V', b'V', b'V', b'V',
    b'*', b'Y', b'*', b'Y', b'S', b'S', b'S', b'S', b'*', b'C', b'W', b'C', b'L', b'F', b'L', b'F',
    b'X',
];

const ANTI_CODON_CODE: [u8; 65] = [
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'M', b'W', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'F', b'V', b'L', b'I', b'C', b'G', b'R', b'S', b'S', b'A', b'P', b'T', b'Y', b'D', b'H', b'N',
    b'L', b'V', b'L', b'I', b'*', b'G', b'R', b'R', b'S', b'A', b'P', b'T', b'*', b'E', b'Q', b'K',
    b'X',
];

fn print_protein<W: Write>(
    dna: &Vec<u8>,
    forward_strand: bool,
    whole_genome: bool,
    head: Vec<u8>,
    file: &mut W,
) -> Result<(), Box<dyn Error>> {
    let mut protein: Vec<u8> = if forward_strand {
        dna.chunks_exact(3)
            .map(|c| CODON_CODE[trinucleotide_pep(c[0], c[1], c[2])])
            .collect()
    } else {
        dna.rchunks_exact(3)
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
        if forward_strand {
            let s = trinucleotide_pep(dna[0], dna[1], dna[2]);
            if s == trinucleotide_pep(b'G', b'T', b'G') || s == trinucleotide_pep(b'T', b'T', b'G')
            {
                protein[0] = b'M';
            }
        } else {
            let s = trinucleotide_pep(dna[dna.len() - 3], dna[dna.len() - 2], dna[dna.len() - 1]);
            if s == trinucleotide_pep(b'C', b'A', b'C') || s == trinucleotide_pep(b'C', b'A', b'A')
            {
                protein[0] = b'M';
            }
        }
    }

    fasta::OwnedRecord {
        head: head,
        seq: protein,
    }
    .write(file)?;

    Ok(())
}

fn print_dna(dna: &Vec<u8>, head: &Vec<u8>, file: &mut File) -> Result<(), Box<dyn Error>> {
    write!(
        file,
        "> {}\n{}\n",
        std::str::from_utf8(head)?,
        std::str::from_utf8(dna)?
    )?;
    Ok(())
}

fn get_rc_dna(dna: &Vec<u8>) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'N' => b'N',
            b'n' => b'n',
            _ => b'x',
        })
        .collect()
}
