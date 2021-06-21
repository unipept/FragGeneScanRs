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
        matches.is_present("complete"),
    )?;

    Ok(())
}

fn run<R: Read, W: Write>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: &mut fasta::Reader<R, StdPolicy>,
    outputstream: &mut W,
    whole_genome: bool,
) -> Result<(), Box<dyn Error>> {
    while let Some(record) = inputseqs.next() {
        let record = record?;
        viterbi(&global, &locals, record, outputstream, whole_genome)?;
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
    global: &hmm::Global,
    locals: &Vec<hmm::Local>,
    record: fasta::RefRecord,
    outputstream: &mut W,
    whole_genome: bool,
) -> Result<(), Box<dyn Error>> {
    let seq = record.full_seq();
    let cg = count_cg_content(&seq);

    let gene_len = if whole_genome { 120 } else { 60 }; // minimum length to be output

    let mut alpha: Vec<[f64; hmm::NUM_STATE]> = vec![];
    let mut path: Vec<[usize; hmm::NUM_STATE]> = vec![];
    let mut temp_i: [usize; 6] = [0; 6];
    let mut temp_i_1: [usize; 6] = [0; 6];

    for _ in 0..seq.len() {
        alpha.push([0.0; hmm::NUM_STATE]);
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
                if t == 0 { // TODO could t ever be 0 here? (elsewhere, too)
                } else {
                    if i == hmm::M1_STATE {
                        // from M state
                        alpha[t][i] = alpha[t - 1][hmm::M6_STATE]
                            - global.tr[hmm::TR_MM]
                            - locals[cg].e_m[0][from2][to];
                        path[t][i] = hmm::M6_STATE;

                        // from D state
                        if !whole_genome {
                            for j in (hmm::M1_STATE..=hmm::M5_STATE).rev() {
                                let num_d = if j >= i {
                                    (i - j + 6) as i32
                                } else if j + 1 < i {
                                    (i - j) as i32
                                } else {
                                    -10
                                };
                                if num_d > 0 {
                                    let temp_alpha = alpha[t - 1][j] - global.tr[hmm::TR_MD] - locals[cg].e_m[0][from2][to]
                                                   - 0.25_f64.ln() * (num_d - 1) as f64 // TODO num_d - 1 could be negative?
                                                   - global.tr[hmm::TR_DD] * (num_d - 2) as f64 // TODO ^
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
                                    (i - j + 6) as i32
                                } else if j + 1 < i {
                                    (i - j) as i32
                                } else {
                                    -10
                                };
                                if num_d > 0 {
                                    let temp_alpha = alpha[t - 1][j] - global.tr[hmm::TR_MD]
                                                   - locals[cg].e_m[i - hmm::M1_STATE][from2][to]
                                                   - 0.25_f64.ln() * (num_d - 1) as f64 // TODO num_d - 1 could be negative?
                                                   - global.tr[hmm::TR_DD] * (num_d - 2) as f64 // TODO ^
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
                        hmm::I6_STATE + (i - hmm::M1_STATE - 1)
                    };

                    // to avoid stop codon
                    if t < 2 {
                    } else if (i == hmm::M2_STATE || i == hmm::M5_STATE)
                        && (seq[temp_i[j - hmm::I1_STATE]] == b'T')
                        && (
                            (seq[t] == b'A' && seq[t + 1] == b'A') // TODO check t border?
                                   || (seq[t] == b'A' && seq[t + 1] == b'G') // TODO check t border?
                                   || (seq[t] == b'G' && seq[t + 1] == b'A')
                            // TODO check t border?
                        )
                    {
                    } else if (i == hmm::M3_STATE || i == hmm::M6_STATE)
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
        }

        // I state
        for i in hmm::I1_STATE..=hmm::I6_STATE {
            if t == 0 {
            } else {
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
                if t == 0 {
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
                                    (i - j + 6) as i32
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
                                    (i - j + 6) as i32
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
                               && (seq[t + 1] == b'A') // TODO check t border?
                               && (
                                      (seq[t] == b'T' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T')
                                   || (seq[t] == b'T' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'C')
                                   || (seq[t] == b'A' && seq[temp_i_1[j - hmm::I1_STATE_1]] == b'T')
                                  )
                    {
                    } else if (i == hmm::M3_STATE_1 || i == hmm::M6_STATE_1)
                        && (seq[t] == b'A')
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
        }

        // I' state
        for i in hmm::I1_STATE_1..=hmm::I6_STATE_1 {
            if t == 0 {
            } else {
                // from I state
                alpha[t][i] = alpha[t - 1][i] - global.tr[hmm::TR_II] - global.tr_ii[from][to];
                path[t][i] = i;

                // from M state
                if path[t - 3][hmm::S_STATE_1] != hmm::R_STATE
                    && path[t - 4][hmm::S_STATE_1] != hmm::R_STATE
                    && path[t - 5][hmm::S_STATE_1] != hmm::R_STATE
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
        }

        // non_coding state
        if t == 0 {
        } else {
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
        }

        // end state
        if alpha[hmm::E_STATE][t] == 0.0 {
            // TODO could it ever not be?
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

                alpha[hmm::E_STATE][t] = f64::INFINITY;
                alpha[hmm::E_STATE][t + 1] = f64::INFINITY;
                path[hmm::E_STATE][t + 1] = hmm::E_STATE;
                path[hmm::E_STATE][t + 2] = hmm::E_STATE;

                alpha[hmm::M6_STATE][t + 2] = f64::INFINITY;
                alpha[hmm::M5_STATE][t + 1] = f64::INFINITY;
                alpha[hmm::M4_STATE][t] = f64::INFINITY;
                alpha[hmm::M3_STATE][t + 2] = f64::INFINITY;
                alpha[hmm::M2_STATE][t + 1] = f64::INFINITY;
                alpha[hmm::M1_STATE][t] = f64::INFINITY;

                alpha[t + 2][hmm::E_STATE] -= if seq[t + 1] == b'A' && seq[t + 2] == b'A' {
                    0.54_f64.ln()
                } else if seq[t + 1] == b'A' && seq[t + 2] == b'G' {
                    0.16_f64.ln()
                } else {
                    0.30_f64.ln()
                };

                // adjustment based on probability distribution
                let mut start_freq = 0.0;
                let mut freq_id = 0;
                let mut sub_sum = 0.0;
                let mut sub_count = 0.0;

                if t >= 60 {
                    // bug reported by Yu-Wei ------ TODO 60 is the incomplete minimum length? can be merged?
                    for i in (t - 60)..=(t - 3) {
                        if i + 2 < seq.len() {
                            start_freq -=
                                locals[cg].tr_e[i - t + 60][codon(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else {
                    for i in 0..=(t - 3) {
                        if i + 2 < seq.len() {
                            sub_sum +=
                                locals[cg].tr_e[i - t + 60][codon(seq[i], seq[i + 1], seq[i + 2])];
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
                let mut freq_id = 0;

                // TODO needs same 60-edgecase as above?
                for i in 3..=60 {
                    if t + i + 2 < seq.len() {
                        start_freq += locals[cg].tr_s1[i - 3]
                            [codon(seq[t + i], seq[t + i + 1], seq[t + i + 2])];
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
                alpha[t + 2][hmm::S_STATE_1] -=
                    (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln().ln();
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
                let mut freq_id = 0;
                let mut sub_sum = 0.0;
                let mut sub_count = 0.0;

                if t >= 30 {
                    // TODO why 30?
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -=
                                locals[cg].tr_s[i - t + 30][codon(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum +=
                                locals[cg].tr_s[i - t + 30][codon(seq[i], seq[i + 1], seq[i + 2])];
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
                alpha[t + 2][hmm::S_STATE] -= (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln().ln();
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
                let mut freq_id = 0;
                let mut sub_sum = 0.0;
                let mut sub_count = 0.0;

                if t >= 30 {
                    for i in (t - 30)..=(t + 30) {
                        if i + 2 < seq.len() {
                            start_freq -=
                                locals[cg].tr_e1[i - t + 30][codon(seq[i], seq[i + 1], seq[i + 2])];
                        }
                    }
                } else {
                    for i in 0..=(t + 30) {
                        if i + 2 < seq.len() {
                            sub_sum +=
                                locals[cg].tr_e1[i - t + 30][codon(seq[i], seq[i + 1], seq[i + 2])];
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
                alpha[t + 2][hmm::E_STATE_1] -=
                    (h_kd / (h_kd + r_kd)).max(0.01).min(0.99).ln().ln();
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

    // backtrack array to find the optimal path

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

fn codon(a: u8, b: u8, c: u8) -> usize {
    if let (Some(a_), Some(b_), Some(c_)) = (nt2int(a), nt2int(b), nt2int(c)) {
        16 * a_ + 4 * b_ + c_
    } else {
        0
    }
}
