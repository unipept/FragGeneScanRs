#![allow(non_snake_case)] // TODO
#![allow(dead_code)] // TODO

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::str::FromStr;

const NUM_STATE: usize = 29;
const NUM_TRANSITIONS: usize = 14;
const CG: usize = 44;
const ACGT: usize = 4;
const PERIOD: usize = 6;
const CONDITION: usize = 16;

pub struct HmmGlobal {
    pi: [f64; NUM_STATE],
    tr: [f64; NUM_TRANSITIONS],
    trII: [[f64; ACGT]; ACGT],
    trMI: [[f64; ACGT]; ACGT],
}

pub struct HmmCG {
    eM1: [[[f64; 6]; 16]; 4],
    eM: [[[f64; 6]; 16]; 4],

    trRR: [[f64; 4]; 4],

    trS: [[f64; 64]; 61],
    trE: [[f64; 64]; 61],
    trS1: [[f64; 64]; 61],
    trE1: [[f64; 64]; 61],

    distS: [f64; 6],
    distE: [f64; 6],
    distS1: [f64; 6],
    distE1: [f64; 6],
}

pub struct Train {
    trans: [[[[f64; ACGT]; CONDITION]; PERIOD]; CG],
    rtrans: [[[[f64; ACGT]; CONDITION]; PERIOD]; CG],
    noncoding: [[[f64; ACGT]; ACGT]; CG],
    start: [[[f64; 64]; 61]; CG],
    stop: [[[f64; 64]; 61]; CG],
    start1: [[[f64; 64]; 61]; CG],
    stop1: [[[f64; 64]; 61]; CG],

    distS: [[f64; 6]; CG],
    distE: [[f64; 6]; CG],
    distS1: [[f64; 6]; CG],
    distE1: [[f64; 6]; CG],
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(HmmGlobal, Train), Box<dyn Error>> {
    let (tr, trMI, trII, pi) = read_transitions(train_dir.join(filename))?;
    let trans = read_M_transitions(train_dir.join("gene"))?;
    let rtrans = read_M_transitions(train_dir.join("rgene"))?;
    let noncoding = read_noncoding(train_dir.join("noncoding"))?;
    let start = read_start(train_dir.join("start"))?;
    let stop = read_stop(train_dir.join("stop"))?;
    let start1 = read_stop(train_dir.join("start1"))?;
    let stop1 = read_start(train_dir.join("stop1"))?;
    let (distS, distE, distS1, distE1) = read_pwm(train_dir.join("pwm"))?;

    Ok((
        HmmGlobal { pi, tr, trII, trMI },
        Train {
            trans,
            rtrans,
            noncoding,
            start,
            stop,
            start1,
            stop1,
            distS,
            distE,
            distS1,
            distE1,
        },
    ))
}

fn read_transitions(
    filename: PathBuf,
) -> Result<
    (
        [f64; NUM_TRANSITIONS],
        [[f64; 4]; 4],
        [[f64; 4]; 4],
        [f64; NUM_STATE],
    ),
    Box<dyn Error>,
> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut tr = [0.0; NUM_TRANSITIONS];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..NUM_TRANSITIONS {
        line = lines.next().ok_or("Incomplete training file")??;
        let v: Vec<&str> = line.split(char::is_whitespace).collect();
        tr[i] = f64::from_str(v[1])?.ln();
    }

    let mut trMI = [[0.0; ACGT]; ACGT];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            trMI[i][j] = f64::from_str(v[2])?.ln();
        }
    }

    let mut trII = [[0.0; ACGT]; ACGT];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            trII[i][j] = f64::from_str(v[2])?.ln();
        }
    }

    let mut pi = [0.0; NUM_STATE];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..NUM_STATE {
        line = lines.next().ok_or("Incomplete training file")??;
        let v: Vec<&str> = line.split(char::is_whitespace).collect();
        pi[i] = f64::from_str(v[1])?.ln();
    }

    Ok((tr, trMI, trII, pi))
}

fn read_M_transitions(
    filename: PathBuf,
) -> Result<[[[[f64; ACGT]; CONDITION]; PERIOD]; CG], Box<dyn Error>> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut trans = [[[[0.0; ACGT]; CONDITION]; PERIOD]; CG];
    for cg in 0..CG {
        lines.next().ok_or("Incomplete training file")??; // Transition header
        for p in 0..PERIOD {
            for c in 0..CONDITION {
                line = lines.next().ok_or("Incomplete training file")??;
                let v: Vec<&str> = line.split(char::is_whitespace).collect();
                for e in 0..ACGT {
                    trans[cg][p][c][e] = f64::from_str(v[e])?.ln();
                }
            }
        }
    }

    Ok(trans)
}

fn read_noncoding(filename: PathBuf) -> Result<[[[f64; ACGT]; ACGT]; CG], Box<dyn Error>> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut noncoding = [[[0.0; ACGT]; ACGT]; CG];
    for cg in 0..CG {
        lines.next().ok_or("Incomplete training file")??; // Transition header
        for e1 in 0..ACGT {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            for e2 in 0..ACGT {
                noncoding[cg][e1][e2] = f64::from_str(v[e2])?.ln();
            }
        }
    }

    Ok(noncoding)
}

fn read_start(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut start = [[[0.0; 64]; 61]; CG];
    for cg in 0..CG {
        lines.next().ok_or("Incomplete training file")??; // Transition header
        for j in 0..61 {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            for k in 0..64 {
                start[cg][j][k] = f64::from_str(v[k])?.ln();
            }
        }
    }

    Ok(start)
}

fn read_stop(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut stop = [[[0.0; 64]; 61]; CG];
    for cg in 0..CG {
        lines.next().ok_or("Incomplete training file")??; // Transition header
        for j in 0..58 {
            // TODO why only 58?
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            for k in 0..64 {
                stop[cg][j][k] = f64::from_str(v[k])?.ln();
            }
        }
    }

    Ok(stop)
}

fn read_pwm(
    filename: PathBuf,
) -> Result<
    (
        [[f64; 6]; CG],
        [[f64; 6]; CG],
        [[f64; 6]; CG],
        [[f64; 6]; CG],
    ),
    Box<dyn Error>,
> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    let mut distS = [[0.0; 6]; CG];
    let mut distE = [[0.0; 6]; CG];
    let mut distS1 = [[0.0; 6]; CG];
    let mut distE1 = [[0.0; 6]; CG];

    for cg in 0..CG {
        let mut v: Vec<&str>;
        lines.next().ok_or("Incomplete training file")??; // Transition header
        line = lines.next().ok_or("Incomplete training file")??;
        v = line.split(char::is_whitespace).collect();
        for j in 0..6 {
            distS[cg][j] = f64::from_str(v[j])?;
        }
        v = line.split(char::is_whitespace).collect();
        for j in 0..6 {
            distE[cg][j] = f64::from_str(v[j])?;
        }
        v = line.split(char::is_whitespace).collect();
        for j in 0..6 {
            distS1[cg][j] = f64::from_str(v[j])?;
        }
        v = line.split(char::is_whitespace).collect();
        for j in 0..6 {
            distE1[cg][j] = f64::from_str(v[j])?;
        }
    }

    Ok((distS, distE, distS1, distE1))
}
