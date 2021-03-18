/// FragGeneScanRs executable

extern crate clap;
use clap::{Arg, App};

fn main() {
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
}
