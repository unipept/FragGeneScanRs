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
pub struct Global {
    pi: [f64; NUM_STATE],
    tr: [f64; NUM_TRANSITIONS],
    tr_ii: [[f64; ACGT]; ACGT],
    tr_mi: [[f64; ACGT]; ACGT],
}

pub struct Local {
    e_m: [[[f64; ACGT]; CONDITION]; PERIOD],
    e_m1: [[[f64; ACGT]; CONDITION]; PERIOD],

    tr_rr: [[f64; ACGT]; ACGT],

    tr_s: [[f64; 64]; 61],
    tr_e: [[f64; 64]; 61],
    tr_s1: [[f64; 64]; 61],
    tr_e1: [[f64; 64]; 61],

    dist_s: [f64; 6],
    dist_e: [f64; 6],
    dist_s1: [f64; 6],
    dist_e1: [f64; 6],
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(Box<Global>, Vec<Local>)> {
    let mut global: Box<Global> = Box::new(Default::default());
    let mut locals: Vec<Local> = Vec::with_capacity(CG);
    for _ in 0..CG {
        locals.push(Local {
            e_m: [[[0.0; ACGT]; CONDITION]; PERIOD],
            e_m1: [[[0.0; ACGT]; CONDITION]; PERIOD],
            tr_rr: [[0.0; ACGT]; ACGT],
            tr_s: [[0.0; 64]; 61],
            tr_e: [[0.0; 64]; 61],
            tr_s1: [[0.0; 64]; 61],
            tr_e1: [[0.0; 64]; 61],

            dist_s: [0.0; 6],
            dist_e: [0.0; 6],
            dist_s1: [0.0; 6],
            dist_e1: [0.0; 6],
        });
    }

    read_transitions(&mut global, train_dir.join(filename))
        .context("Could not read training file")?;
    read_m_transitions(&mut locals, train_dir.join("gene"))
        .context("Could not read gene training file")?;
    read_m1_transitions(&mut locals, train_dir.join("rgene"))
        .context("Could not read rgene training file")?;
    read_noncoding(&mut locals, train_dir.join("noncoding"))
        .context("Could not read noncoding training file")?;
    read_start(&mut locals, train_dir.join("start"))
        .context("Could not read start training file")?;
    read_stop(&mut locals, train_dir.join("stop")).context("Could not read stop training file")?;
    read_start1(&mut locals, train_dir.join("start1"))
        .context("Could not read start1 training file")?;
    read_stop1(&mut locals, train_dir.join("stop1"))
        .context("Could not read stop1 training file")?;
    read_pwm(&mut locals, train_dir.join("pwm")).context("Could not read pwm training file")?;

    Ok((global, locals))
}

fn read_transitions(global: &mut Global, filename: PathBuf) -> Result<()> {
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
        global.tr[i] = f64::from_str(v[1])
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
            global.tr_mi[i][j] = f64::from_str(v[2])
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
            global.tr_ii[i][j] = f64::from_str(v[2])
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
        global.pi[i] = f64::from_str(v[1]).context("Could not read PI")?.ln();
    }

    Ok(())
}

fn read_m_transitions(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                    locals[cg].e_m[p][c][e] = f64::from_str(v[e])
                        .context("Could not read gene transitions")?
                        .ln();
                }
            }
        }
    }

    Ok(())
}

fn read_m1_transitions(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                    locals[cg].e_m1[p][c][e] = f64::from_str(v[e])
                        .context("Could not read rgene transitions")?
                        .ln();
                }
            }
        }
    }

    Ok(())
}

fn read_noncoding(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                locals[cg].tr_rr[e1][e2] = f64::from_str(v[e2])
                    .context("Could not read nonecoding transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_start(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                locals[cg].tr_s[j][k] = f64::from_str(v[k])
                    .context("Could not read start transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_stop1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                locals[cg].tr_e1[j][k] = f64::from_str(v[k])
                    .context("Could not read stop1 transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_stop(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                locals[cg].tr_e[j][k] = f64::from_str(v[k])
                    .context("Could not read stop transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_start1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
                locals[cg].tr_s1[j][k] = f64::from_str(v[k])
                    .context("Could not read start1 transitions")?
                    .ln();
            }
        }
    }

    Ok(())
}

fn read_pwm(locals: &mut Vec<Local>, filename: PathBuf) -> Result<()> {
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
            locals[cg].dist_s[j] =
                f64::from_str(v[j]).context("Could not read distS transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            locals[cg].dist_e[j] =
                f64::from_str(v[j]).context("Could not read distE transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            locals[cg].dist_s1[j] =
                f64::from_str(v[j]).context("Could not read distS1 transitions")?;
        }
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..6 {
            locals[cg].dist_e1[j] =
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
