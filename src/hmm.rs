#![allow(non_snake_case)] // TODO

use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::str::FromStr;

extern crate anyhow;
use anyhow::Context;
use anyhow::Result;

extern crate thiserror;
use thiserror::Error;

const NUM_STATE: usize = 29;
const NUM_TRANSITIONS: usize = 14;
const CG: usize = 44;
const ACGT: usize = 4;
const PERIOD: usize = 6;
const CONDITION: usize = 16;

#[derive(Default)]
pub struct HmmGlobal {
    pi: [f64; NUM_STATE],
    tr: [f64; NUM_TRANSITIONS],
    trII: [[f64; ACGT]; ACGT],
    trMI: [[f64; ACGT]; ACGT],
}

pub struct Train {
    eM: [[[f64; ACGT]; CONDITION]; PERIOD],
    eM1: [[[f64; ACGT]; CONDITION]; PERIOD],

    trRR: [[f64; ACGT]; ACGT],

    trS: [[f64; 64]; 61],
    trE: [[f64; 64]; 61],
    trS1: [[f64; 64]; 61],
    trE1: [[f64; 64]; 61],

    distS: [f64; 6],
    distE: [f64; 6],
    distS1: [f64; 6],
    distE1: [f64; 6],
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(Box<HmmGlobal>, Vec<Train>)> {
    let mut hmmGlobal: Box<HmmGlobal> = Box::new(Default::default());
    let mut train: Vec<Train> = Vec::with_capacity(CG);
    for _ in 0..CG {
        train.push(Train {
            eM: [[[0.0; ACGT]; CONDITION]; PERIOD],
            eM1: [[[0.0; ACGT]; CONDITION]; PERIOD],
            trRR: [[0.0; ACGT]; ACGT],
            trS: [[0.0; 64]; 61],
            trE: [[0.0; 64]; 61],
            trS1: [[0.0; 64]; 61],
            trE1: [[0.0; 64]; 61],

            distS: [0.0; 6],
            distE: [0.0; 6],
            distS1: [0.0; 6],
            distE1: [0.0; 6],
        });
    }

    read_transitions(&mut hmmGlobal, train_dir.join(filename))
        .context("Could not read training file")?;
    read_M_transitions(&mut train, train_dir.join("gene"))
        .context("Could not read gene training file")?;
    read_M1_transitions(&mut train, train_dir.join("rgene"))
        .context("Could not read rgene training file")?;
    read_noncoding(&mut train, train_dir.join("noncoding"))
        .context("Could not read noncoding training file")?;
    read_start(&mut train, train_dir.join("start"))
        .context("Could not read start training file")?;
    read_stop(&mut train, train_dir.join("stop")).context("Could not read stop training file")?;
    read_start1(&mut train, train_dir.join("start1"))
        .context("Could not read start1 training file")?;
    read_stop1(&mut train, train_dir.join("stop1"))
        .context("Could not read stop1 training file")?;
    read_pwm(&mut train, train_dir.join("pwm")).context("Could not read pwm training file")?;

    Ok((hmmGlobal, train))
}

fn read_transitions(hmmGlobal: &mut HmmGlobal, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
    for i in 0..NUM_TRANSITIONS {
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        hmmGlobal.tr[i] = f64::from_str(v[1])
            .context("Could not read Transition")?
            .ln();
    }

    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            hmmGlobal.trMI[i][j] = f64::from_str(v[2])
                .context("Could not read TransitionMI")?
                .ln();
        }
    }

    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
    for i in 0..ACGT {
        for j in 0..ACGT {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            hmmGlobal.trII[i][j] = f64::from_str(v[2])
                .context("Could not read TransisionII")?
                .ln();
        }
    }

    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
    for i in 0..NUM_STATE {
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        hmmGlobal.pi[i] = f64::from_str(v[1]).context("Could not read PI")?.ln();
    }

    Ok(())
}

fn read_M_transitions(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for p in 0..PERIOD {
            for c in 0..CONDITION {
                line = lines
                    .next()
                    .ok_or(TrainingDataError::IncompleteTrainingFile)??;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    train[cg].eM[p][c][e] = f64::from_str(v[e])
                        .context("Could not read gene transitions")?
                        .ln();
                }
            }
        }
    }

    Ok(())
}

fn read_M1_transitions(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for p in 0..PERIOD {
            for c in 0..CONDITION {
                line = lines
                    .next()
                    .ok_or(TrainingDataError::IncompleteTrainingFile)??;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    train[cg].eM1[p][c][e] = f64::from_str(v[e])
                        .context("Could not read rgene transitions")?
                        .ln();
                }
            }
        }
    }

    Ok(())
}

fn read_noncoding(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for e1 in 0..ACGT {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            for e2 in 0..ACGT {
                train[cg].trRR[e1][e2] = f64::from_str(v[e2])
                    .context("Could not read nonecoding transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_start(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for j in 0..61 {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..64 {
                train[cg].trS[j][k] = f64::from_str(v[k])
                    .context("Could not read start transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_stop1(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for j in 0..61 {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..64 {
                train[cg].trE1[j][k] = f64::from_str(v[k])
                    .context("Could not read stop1 transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_stop(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for j in 0..61 {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..64 {
                train[cg].trE[j][k] = f64::from_str(v[k])
                    .context("Could not read stop transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_start1(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header
        for j in 0..61 {
            line = lines
                .next()
                .ok_or(TrainingDataError::IncompleteTrainingFile)??;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..64 {
                train[cg].trS1[j][k] = f64::from_str(v[k])
                    .context("Could not read start1 transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_pwm(train: &mut Vec<Train>, filename: PathBuf) -> Result<()> {
    let mut lines = io::BufReader::new(File::open(filename)?).lines();
    let mut line: String;

    for cg in 0..CG {
        lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??; // Transition header

        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            train[cg].distS[j] = f64::from_str(v[j]).context("Could not read distS transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            train[cg].distE[j] = f64::from_str(v[j]).context("Could not read distE transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            train[cg].distS1[j] =
                f64::from_str(v[j]).context("Could not read distS1 transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            train[cg].distE1[j] =
                f64::from_str(v[j]).context("Could not read distE1 transitions")?;
        }
    }

    Ok(())
}

#[derive(Error, Debug)]
pub enum TrainingDataError {
    #[error("incomplete training file")]
    IncompleteTrainingFile,
}
