#![allow(non_snake_case)] // TODO
#![allow(dead_code)] // TODO
#![allow(unused_variables)] // TODO

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::str::FromStr;

const NUM_STATE: usize = 29;
const NUM_TRANSITIONS: usize = 14;
const CG: usize = 44;
const ACGT: usize = 4;

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
    trans: [[[[f64; 4]; 16]; 6]; CG],
    rtrans: [[[[f64; 4]; 16]; 6]; CG],
    noncoding: [[[f64; 4]; 4]; CG],
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
    let rtrans = read_M1_transitions(train_dir.join("rgene"))?;
    let noncoding = read_noncoding(train_dir.join("noncoding"))?;
    let start = read_start(train_dir.join("start"))?;
    let stop = read_stop(train_dir.join("stop"))?;
    let start1 = read_start1(train_dir.join("start1"))?;
    let stop1 = read_stop1(train_dir.join("stop1"))?;
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
    let mut line;

    let mut tr = [0.0; NUM_TRANSITIONS];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..NUM_TRANSITIONS {
        line = lines.next().ok_or("Incomplete training file")??;
        let v: Vec<&str> = line.split(char::is_whitespace).collect();
        tr[i] = f64::from_str(v[1])?;
    }

    let mut trMI = [[0.0; ACGT]; ACGT];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            trMI[i][j] = f64::from_str(v[2])?;
        }
    }

    let mut trII = [[0.0; ACGT]; ACGT];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines.next().ok_or("Incomplete training file")??;
            let v: Vec<&str> = line.split(char::is_whitespace).collect();
            trII[i][j] = f64::from_str(v[2])?;
        }
    }

    let mut pi = [0.0; NUM_STATE];
    lines.next().ok_or("Incomplete training file")??; // Transition header
    for i in 0..NUM_STATE {
        line = lines.next().ok_or("Incomplete training file")??;
        let v: Vec<&str> = line.split(char::is_whitespace).collect();
        pi[i] = f64::from_str(v[1])?;
    }

    Ok((tr, trMI, trII, pi))
}

fn read_M_transitions(filename: PathBuf) -> Result<[[[[f64; 4]; 16]; 6]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_M1_transitions(filename: PathBuf) -> Result<[[[[f64; 4]; 16]; 6]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_noncoding(filename: PathBuf) -> Result<[[[f64; 4]; 4]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_start(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_stop(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_start1(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
}

fn read_stop1(filename: PathBuf) -> Result<[[[f64; 64]; 61]; CG], Box<dyn Error>> {
    Err("Not implemented yet.")?
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
    Err("Not implemented yet.")?
}
