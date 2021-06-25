use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::str::FromStr;

extern crate anyhow;
use anyhow::Context;
use anyhow::Result;

extern crate thiserror;
use thiserror::Error;

pub const NUM_STATE: usize = 29;

pub const S_STATE: usize = 0;
pub const E_STATE: usize = 1;
pub const R_STATE: usize = 2;
pub const S_STATE_1: usize = 3;
pub const E_STATE_1: usize = 4;
pub const M1_STATE: usize = 5;
pub const M2_STATE: usize = 6;
pub const M3_STATE: usize = 7;
pub const M4_STATE: usize = 8;
pub const M5_STATE: usize = 9;
pub const M6_STATE: usize = 10;
pub const M1_STATE_1: usize = 11;
pub const M2_STATE_1: usize = 12;
pub const M3_STATE_1: usize = 13;
pub const M4_STATE_1: usize = 14;
pub const M5_STATE_1: usize = 15;
pub const M6_STATE_1: usize = 16;
pub const I1_STATE: usize = 17;
pub const I2_STATE: usize = 18;
pub const I3_STATE: usize = 19;
pub const I4_STATE: usize = 20;
pub const I5_STATE: usize = 21;
pub const I6_STATE: usize = 22;
pub const I1_STATE_1: usize = 23;
pub const I2_STATE_1: usize = 24;
pub const I3_STATE_1: usize = 25;
pub const I4_STATE_1: usize = 26;
pub const I5_STATE_1: usize = 27;
pub const I6_STATE_1: usize = 28;
pub const NOSTATE: usize = 30;

const NUM_TRANSITIONS: usize = 14;

pub const TR_MM: usize = 0;
pub const TR_MI: usize = 1;
pub const TR_MD: usize = 2;
pub const TR_II: usize = 3;
pub const TR_IM: usize = 4;
pub const TR_DD: usize = 5;
pub const TR_DM: usize = 6;
pub const TR_GE: usize = 7;
pub const TR_GG: usize = 8;
pub const TR_ER: usize = 9;
pub const TR_RS: usize = 10;
pub const TR_RR: usize = 11;
pub const TR_ES: usize = 12;
pub const TR_ES1: usize = 13;

const CG: usize = 44;
const ACGT: usize = 4;
const PERIOD: usize = 6;
const ACGTACGT: usize = 4 * 4;
#[derive(Default)]
pub struct Global {
    pub pi: [f64; NUM_STATE],
    pub tr: [f64; NUM_TRANSITIONS],
    pub tr_ii: [[f64; ACGT]; ACGT],
    pub tr_mi: [[f64; ACGT]; ACGT],
}

pub struct Local {
    pub e_m: [[[f64; ACGT]; ACGTACGT]; PERIOD],
    pub e_m1: [[[f64; ACGT]; ACGTACGT]; PERIOD],

    pub tr_rr: [[f64; ACGT]; ACGT],

    pub tr_s: [[f64; 64]; 61],
    pub tr_e: [[f64; 64]; 61],
    pub tr_s1: [[f64; 64]; 61],
    pub tr_e1: [[f64; 64]; 61],

    pub dist_s: [f64; 6],
    pub dist_e: [f64; 6],
    pub dist_s1: [f64; 6],
    pub dist_e1: [f64; 6],
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(Box<Global>, Vec<Local>)> {
    let mut global: Box<Global> = Box::new(Default::default());
    let mut locals: Vec<Local> = Vec::with_capacity(CG);
    for _ in 0..CG {
        locals.push(Local {
            e_m: [[[0.0; ACGT]; ACGTACGT]; PERIOD],
            e_m1: [[[0.0; ACGT]; ACGTACGT]; PERIOD],
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
    for _ in 0..NUM_TRANSITIONS {
        line = lines
            .next()
            .ok_or(TrainingDataError::IncompleteTrainingFile)??;
        let v: Vec<&str> = line.split_whitespace().collect();
        let l = match v[0] {
            "MM" => 0,
            "MI" => 1,
            "MD" => 2,
            "II" => 3,
            "IM" => 4,
            "DD" => 5,
            "DM" => 6,
            "GE" => 7,
            "GG" => 8,
            "ER" => 9,
            "RS" => 10,
            "RR" => 11,
            "ES" => 12,
            "ES1" => 13,
            _ => Err(TrainingDataError::UnknownTransitionState)?,
        };
        global.tr[l] = f64::from_str(v[1])
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
            for c in 0..ACGTACGT {
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
            for c in 0..ACGTACGT {
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
    #[error("unknown transition state")]
    UnknownTransitionState,
}
