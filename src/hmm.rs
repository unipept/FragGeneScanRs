use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines};
use std::path::PathBuf;
use std::str::FromStr;

extern crate thiserror;
use thiserror::Error;

use crate::dna::{ACGT, BI_ACGT, CG_MAX, CG_MIN, TRI_ACGT};

pub const NUM_STATE: usize = 29;

// 19 #[derive(Copy, Clone, PartialEq, PartialOrd, Eq, EnumIter, EnumCount)]
// 20 #[allow(non_camel_case_types)]
// 21 pub enum State {
// 22     S = 0,
// 23     E = 1,
// 24     R = 2,

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum State {
    S = 0,
    E = 1,
    R = 2,
    Sr = 3,
    Er = 4,
    M1 = 5,
    M2 = 6,
    M3 = 7,
    M4 = 8,
    M5 = 9,
    M6 = 10,
    M1r = 11,
    M2r = 12,
    M3r = 13,
    M4r = 14,
    M5r = 15,
    M6r = 16,
    I1 = 17,
    I2 = 18,
    I3 = 19,
    I4 = 20,
    I5 = 21,
    I6 = 22,
    I1r = 23,
    I2r = 24,
    I3r = 25,
    I4r = 26,
    I5r = 27,
    I6r = 28,
    N = 30,
}

const NUM_TRANSITIONS: usize = 14;

impl State {
    pub fn next(self) -> Self {
        match self {
            State::S => State::E,
            State::E => State::R,
            State::R => State::Sr,
            State::Sr => State::Er,
            State::Er => State::M1,
            State::M1 => State::M2,
            State::M2 => State::M3,
            State::M3 => State::M4,
            State::M4 => State::M5,
            State::M5 => State::M6,
            State::M6 => State::M1r,
            State::M1r => State::M2r,
            State::M2r => State::M3r,
            State::M3r => State::M4r,
            State::M4r => State::M5r,
            State::M5r => State::M6r,
            State::M6r => State::I1,
            State::I1 => State::I2,
            State::I2 => State::I3,
            State::I3 => State::I4,
            State::I4 => State::I5,
            State::I5 => State::I6,
            State::I6 => State::I1r,
            State::I1r => State::I2r,
            State::I2r => State::I3r,
            State::I3r => State::I4r,
            State::I4r => State::I5r,
            State::I5r => State::I6r,
            State::I6r => State::N,
            State::N => State::S,
        }
    }

    pub fn previous(self) -> Self {
        match self {
            State::S => State::N,
            State::E => State::S,
            State::R => State::E,
            State::Sr => State::R,
            State::Er => State::Sr,
            State::M1 => State::Er,
            State::M2 => State::M1,
            State::M3 => State::M2,
            State::M4 => State::M3,
            State::M5 => State::M4,
            State::M6 => State::M5,
            State::M1r => State::M6,
            State::M2r => State::M1r,
            State::M3r => State::M2r,
            State::M4r => State::M3r,
            State::M5r => State::M4r,
            State::M6r => State::M5r,
            State::I1 => State::M6r,
            State::I2 => State::I1,
            State::I3 => State::I2,
            State::I4 => State::I3,
            State::I5 => State::I4,
            State::I6 => State::I5,
            State::I1r => State::I6,
            State::I2r => State::I1r,
            State::I3r => State::I2r,
            State::I4r => State::I3r,
            State::I5r => State::I4r,
            State::I6r => State::I5r,
            State::N => State::I6r,
        }
    }

    pub fn up_through(self, other: State) -> StateIterator {
        StateIterator {
            start: self,
            finish: other.next(),
        }
    }
}

pub struct StateIterator {
    start: State,
    finish: State,
}

impl Iterator for StateIterator {
    type Item = State;

    fn next(&mut self) -> Option<State> {
        if self.start == self.finish {
            None
        } else {
            self.start = self.start.next();
            Some(self.start.previous())
        }
    }
}

impl DoubleEndedIterator for StateIterator {
    fn next_back(&mut self) -> Option<State> {
        if self.start == self.finish {
            None
        } else {
            self.finish = self.finish.previous();
            Some(self.finish)
        }
    }
}

pub fn states() -> StateIterator {
    State::S.up_through(State::I6r)
}

impl<T> std::ops::Index<State> for [T] {
    type Output = T;

    fn index(&self, idx: State) -> &T {
        &self[idx as usize]
    }
}

impl<T> std::ops::IndexMut<State> for [T] {
    fn index_mut(&mut self, idx: State) -> &mut T {
        &mut self[idx as usize]
    }
}

#[derive(Default)]
pub struct Transition {
    pub mm: f64,
    pub mi: f64,
    pub md: f64,
    pub ii: f64,
    pub im: f64,
    pub dd: f64,
    pub dm: f64,
    pub ge: f64,
    pub gg: f64,
    pub er: f64,
    pub rs: f64,
    pub rr: f64,
    pub es: f64,
    pub es1: f64,
}

const PERIOD: usize = 6;
const WINDOW: usize = 61;

#[derive(Default)]
pub struct Global {
    pub pi: [f64; NUM_STATE],
    pub tr: Transition,
    pub tr_ii: [[f64; ACGT]; ACGT],
    pub tr_mi: [[f64; ACGT]; ACGT],
}

pub struct Local {
    pub e_m: [[[f64; ACGT]; BI_ACGT]; PERIOD],
    pub e_m1: [[[f64; ACGT]; BI_ACGT]; PERIOD],

    pub tr_rr: [[f64; ACGT]; ACGT],

    pub tr_s: [[f64; TRI_ACGT]; WINDOW],
    pub tr_e: [[f64; TRI_ACGT]; WINDOW],
    pub tr_s1: [[f64; TRI_ACGT]; WINDOW],
    pub tr_e1: [[f64; TRI_ACGT]; WINDOW],

    pub dist_s: [f64; PERIOD],
    pub dist_e: [f64; PERIOD],
    pub dist_s1: [f64; PERIOD],
    pub dist_e1: [f64; PERIOD],
}

#[derive(Error, Debug)]
pub enum TrainingDataError {
    #[error("incomplete training file")]
    IncompleteTrainingFile,
    #[error("unknown transition state")]
    UnknownTransitionState,
    #[error("could not read training file {0}")]
    Io(PathBuf, #[source] io::Error),
    #[error("malformed value '{2}' in section {1} of {0}")]
    MalformedValue(PathBuf, String, String, #[source] std::num::ParseFloatError),
}

pub fn get_train_from_file(
    train_dir: PathBuf,
    filename: PathBuf,
) -> Result<(Box<Global>, Vec<Local>), TrainingDataError> {
    let mut global: Box<Global> = Box::new(Default::default());
    let mut locals: Vec<Local> = (0..CG_MAX - CG_MIN)
        .map(|_| Local {
            e_m: [[[0.0; ACGT]; BI_ACGT]; PERIOD],
            e_m1: [[[0.0; ACGT]; BI_ACGT]; PERIOD],
            tr_rr: [[0.0; ACGT]; ACGT],
            tr_s: [[0.0; TRI_ACGT]; WINDOW],
            tr_e: [[0.0; TRI_ACGT]; WINDOW],
            tr_s1: [[0.0; TRI_ACGT]; WINDOW],
            tr_e1: [[0.0; TRI_ACGT]; WINDOW],

            dist_s: [0.0; PERIOD],
            dist_e: [0.0; PERIOD],
            dist_s1: [0.0; PERIOD],
            dist_e1: [0.0; PERIOD],
        })
        .collect();

    read_transitions(&mut global, train_dir.join(filename))?;
    read_m_transitions(&mut locals, train_dir.join("gene"))?;
    read_m1_transitions(&mut locals, train_dir.join("rgene"))?;
    read_noncoding(&mut locals, train_dir.join("noncoding"))?;
    read_start(&mut locals, train_dir.join("start"))?;
    read_stop(&mut locals, train_dir.join("stop"))?;
    read_start1(&mut locals, train_dir.join("start1"))?;
    read_stop1(&mut locals, train_dir.join("stop1"))?;
    read_pwm(&mut locals, train_dir.join("pwm"))?;

    Ok((global, locals))
}

fn lines_from_file(filename: &PathBuf) -> Result<Lines<BufReader<File>>, TrainingDataError> {
    Ok(io::BufReader::new(
        File::open(filename).map_err(|e| TrainingDataError::Io(filename.to_owned(), e))?,
    )
    .lines())
}

fn next_line(
    filename: &PathBuf,
    lines: &mut Lines<BufReader<File>>,
) -> Result<String, TrainingDataError> {
    lines
        .next()
        .ok_or(TrainingDataError::IncompleteTrainingFile)?
        .map_err(|e| TrainingDataError::Io(filename.to_path_buf(), e))
}

fn parse_float(
    filename: &PathBuf,
    part: &String,
    line: &String,
    string: &str,
) -> Result<f64, TrainingDataError> {
    f64::from_str(string).map_err(|e| {
        TrainingDataError::MalformedValue(filename.to_owned(), part.to_owned(), line.to_owned(), e)
    })
}

fn parse_float_col(
    filename: &PathBuf,
    part: &String,
    line: String,
    column: usize,
) -> Result<f64, TrainingDataError> {
    let v: Vec<&str> = line.split_whitespace().collect();
    parse_float(filename, part, &line, v[column])
}

fn read_transitions(global: &mut Global, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    let mut header;

    header = next_line(&filename, &mut lines)?;
    for _ in 0..NUM_TRANSITIONS {
        let line = next_line(&filename, &mut lines)?;
        let v = line.split_whitespace().collect::<Vec<&str>>();
        let value = parse_float(&filename, &header, &line, v[1])?.ln();
        match v[0] {
            "MM" => global.tr.mm = value,
            "MI" => global.tr.mi = value,
            "MD" => global.tr.md = value,
            "II" => global.tr.ii = value,
            "IM" => global.tr.im = value,
            "DD" => global.tr.dd = value,
            "DM" => global.tr.dm = value,
            "GE" => global.tr.ge = value,
            "GG" => global.tr.gg = value,
            "ER" => global.tr.er = value,
            "RS" => global.tr.rs = value,
            "RR" => global.tr.rr = value,
            "ES" => global.tr.es = value,
            "ES1" => global.tr.es1 = value,
            _ => Err(TrainingDataError::UnknownTransitionState)?,
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..ACGT {
        for j in 0..ACGT {
            global.tr_mi[i][j] =
                parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 2)?.ln();
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..ACGT {
        for j in 0..ACGT {
            global.tr_ii[i][j] =
                parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 2)?.ln();
        }
    }

    header = next_line(&filename, &mut lines)?;
    for i in 0..NUM_STATE {
        global.pi[i] =
            parse_float_col(&filename, &header, next_line(&filename, &mut lines)?, 1)?.ln();
    }

    Ok(())
}

fn read_m_transitions(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for p in 0..PERIOD {
            for c in 0..BI_ACGT {
                let line = next_line(&filename, &mut lines)?;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    locals[cg].e_m[p][c][e] = parse_float(&filename, &header, &line, v[e])?.ln();
                }
            }
        }
    }

    Ok(())
}

fn read_m1_transitions(
    locals: &mut Vec<Local>,
    filename: PathBuf,
) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for p in 0..PERIOD {
            for c in 0..BI_ACGT {
                let line = next_line(&filename, &mut lines)?;
                let v: Vec<&str> = line.split_whitespace().collect();
                for e in 0..ACGT {
                    locals[cg].e_m1[p][c][e] = parse_float(&filename, &header, &line, v[e])?.ln();
                }
            }
        }
    }

    Ok(())
}

fn read_noncoding(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for e1 in 0..ACGT {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for e2 in 0..ACGT {
                locals[cg].tr_rr[e1][e2] = parse_float(&filename, &header, &line, v[e2])?.ln();
            }
        }
    }

    Ok(())
}

fn read_start(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_s[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_stop1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_e1[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_stop(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_e[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_start1(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        for j in 0..WINDOW {
            let line = next_line(&filename, &mut lines)?;
            let v: Vec<&str> = line.split_whitespace().collect();
            for k in 0..TRI_ACGT {
                locals[cg].tr_s1[j][k] = parse_float(&filename, &header, &line, v[k])?.ln();
            }
        }
    }

    Ok(())
}

fn read_pwm(locals: &mut Vec<Local>, filename: PathBuf) -> Result<(), TrainingDataError> {
    let mut lines = lines_from_file(&filename)?;
    let mut line: String;
    for cg in 0..(CG_MAX - CG_MIN) {
        let header = next_line(&filename, &mut lines)?;
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_s[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_e[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_s1[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
        line = next_line(&filename, &mut lines)?;
        let v: Vec<&str> = line.split_whitespace().collect();
        for j in 0..PERIOD {
            locals[cg].dist_e1[j] = parse_float(&filename, &header, &line, v[j])?;
            // no ln
        }
    }

    Ok(())
}
