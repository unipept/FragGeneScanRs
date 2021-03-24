use std::path::PathBuf;
use std::error::Error;

const NUM_STATE: usize = 29;
const NUM_TRANSITIONS: usize = 14;

pub struct HmmGlobal {
	pi: [f64; NUM_STATE],
	tr: [f64; NUM_TRANSITIONS],
	trII: [[f64; 4]; 4],
	trMI: [[f64; 4]; 4],
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
	trans: [[[[f64; 4]; 16]; 6]; 44],
	rtrans: [[[[f64; 4]; 16]; 6]; 44],
	noncoding: [[[f64; 4]; 4]; 44],
	start: [[[f64; 64]; 61]; 44],
	stop: [[[f64; 64]; 61]; 44],
	start1: [[[f64; 64]; 61]; 44],
	stop1: [[[f64; 64]; 61]; 44],

	distS: [[f64; 6]; 44],
	distE: [[f64; 6]; 44],
	distS1: [[f64; 6]; 44],
	distE1: [[f64; 6]; 44],
}

pub fn get_train_from_file(train_dir: PathBuf, filename: PathBuf) -> Result<(HmmGlobal, Train), Box<dyn Error>> {
	let (tr, trMI, trII, pi) = read_transitions(train_dir.join(filename))?;
	let trans = read_M_transitions(train_dir.join("gene"))?;
	let rtrans = read_M1_transitions(train_dir.join("rgene"))?;
	let noncoding = read_noncoding(train_dir.join("noncoding"))?;
	let start = read_start(train_dir.join("start"))?;
	let stop = read_stop(train_dir.join("stop"))?;
	let start1 = read_start1(train_dir.join("start1"))?;
	let stop1 = read_stop1(train_dir.join("stop1"))?;
	let (distS, distE, distS1, distE1) = read_pwm(train_dir.join("pwm"))?;

	Ok((HmmGlobal { pi, tr, trII, trMI },
	    Train { trans, rtrans, noncoding, start, stop, start1, stop1, distS, distE, distS1, distE1 }
	))
}

fn read_transitions(filename: PathBuf) -> Result<([f64; NUM_TRANSITIONS], [[f64; 4]; 4], [[f64; 4]; 4], [f64; NUM_STATE]), Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_M_transitions(filename: PathBuf) -> Result<[[[[f64; 4]; 16]; 6]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_M1_transitions(filename: PathBuf) -> Result<[[[[f64; 4]; 16]; 6]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_noncoding(filename: PathBuf) -> Result<[[[f64; 4]; 4]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_start(filename: PathBuf) -> Result<[[[f64; 64]; 61]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_stop(filename: PathBuf) -> Result<[[[f64; 64]; 61]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_start1(filename: PathBuf) -> Result<[[[f64; 64]; 61]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_stop1(filename: PathBuf) -> Result<[[[f64; 64]; 61]; 44], Box<dyn Error>> {
	Err("Not implemented yet.")?
}

fn read_pwm(filename: PathBuf) -> Result<([[f64; 6]; 44], [[f64; 6]; 44], [[f64; 6]; 44], [[f64; 6]; 44]), Box<dyn Error>> {
	Err("Not implemented yet.")?
}
