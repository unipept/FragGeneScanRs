# FragGeneScanRs

## Installation

### From release

Download the build of the [latest release][release] for your platform
and extract it somewhere in your path.

[release]: https://github.com/unipept/FragGeneScanRs/releases/latest

### From source

FragGeneScanRs is written in Rust, so first head over to their
[installation instructions][Rust]. After, clone this repository or
download the source code of the [latest release][release]. In this
directory, run `cargo install --path .` to install. The installation
progress may prompt you to add a directory to your path so you can
easily execute it.

[Rust]: https://www.rust-lang.org/tools/install

## Usage

You can use FragGeneScanRs with the short options of FragGeneScan but
it also provides some additional options and long-form options. It
reads from and writes to standard input and standard output by default,
allowing shorter calls in case you only need the predicted proteins.

```sh
# get predictions for 454 pyrosequencing reads with about 1% error rate
FragGeneScanRs -t 454_10 < example/NC_000913-454.fna > example/NC_000913-454.faa

# get predictions for complete reads
FragGeneScanRs -t complete -w 1 < example/NC_000913.fna > example/NC_000913.faa
```

### Backwards compatible mode

```sh
FragGeneScanRs -s seq_file_name -o output_file_name -w [0 or 1] -t train_file_name -p num_threads
```

where:

* `seq_file_name` is the (FASTA) sequence file name including the full
  path

* `output_file_name` is the base name for the 3 outputfiles, including
  the full path. A `.out`, `.faa` and `.ffn` file will be created
  containing the gene prediction metadata, the predicted proteins, and
  the corresponding DNA reads.

* `0 or 1` for short sequence reads or complete genomic sequences.

* `train_file_name` is used to select the training file for one of the following types:

  - `complete` for complete genomic sequences or short sequence reads without sequencing error
  - `sanger_5` for Sanger sequencing reads with about 0.5% error rate
  - `sanger_10` for Sanger sequencing reads with about 1% error rate
  - `454_5` for 454 pyrosequencing reads with about 0.5% error rate
  - `454_10` for 454 pyrosequencing reads with about 1% error rate
  - `454_30` for 454 pyrosequencing reads with about 3% error rate
  - `illumina_5` for Illumina sequencing reads with about 0.5% error rate
  - `illumina_10` for Illumina sequencing reads with about 1% error rate

  The corresponding file should be in the `train` directory below the
  working directory. Other files can be added and selected here.

* `num_threads` is the number of threads to be used. Defaults to 1.

### Additional options

* `-m meta_file`, `-d dna_file` and `-e aa_file` can be used to write
  output to specific files, instead of having the program create filenames
  with predetermined extentions. These take precedence over the `-o`
  option.

* Leaving out the `-o` option or using the name `stdout` causes
  FragGeneScanRs to only write the predicted proteins to standard output.
  The other files can still be requested with the specific options above.

* Leaving out the `-s` options causes FragGeneScanRs to read the
  sequences from standard input.

* `-r train_file_dir` can change the directory containing the training
  files, so you can put it anywhere on your system.

The complete list of options will be printed when running
`FragGeneScanRs --help`.

## Execution time

Benchmarks were done using the `meta/benchmark.sh` script on a 4-core
Intel(R) Core(TM) i5-7600K CPU @ 3.80GHz with 16GB RAM. The datasets
used are the example datasets provided by FragGeneScan. The table
below shows the average execution time of 10 runs. Detailed results
may be found in `meta/benchmark.csv`. For the short reads (80bp),
FragGeneScanRs is about 12 times faster than FragGeneScan and 1.2
times as fast as FragGeneScanPlus. For the long reads (1328bp) and the
complete genome (Escherichia coli str. K-12 substr. MG1655, 4639675bp),
FragGeneScanRs is 3.5 and 2.4 times faster than FragGeneScan and 1.2 and
many times faster than FGS+.

| Short reads      |  1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|----------:|----------:|----------:|----------:|
| FragGeneScan     |  11.2361s |   6.0421s |   4.5947s |   3.5992s |
| FragGeneScanPlus |   1.1971s |   0.6623s |   0.4772s |   0.4035s |
| FragGeneScanRs   |   0.9974s |   0.5423s |   0.3887s |   0.3058s |

| Long reads       |  1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|----------:|----------:|----------:|----------:|
| FragGeneScan     |  31.7019s |  17.0807s |   12.453s |   9.7152s |
| FragGeneScanPlus |  10.5704s |   5.6104s |     3.95s |    2.956s |
| FragGeneScanRs   |    8.643s |   4.4444s |   3.0728s |   2.3115s |

| complete genome  |  1 thread |
|:-----------------|----------:|
| FragGeneScan     |   3.5526s |
| FragGeneScanPlus | 386.1725s |
| FragGeneScanRs   |   1.3605s |

The commands and arguments used for this benchmarks were:

```sh
./FragGeneScan -t 454_10 -s example/NC_000913-454.fna -o example/NC_000913-454 -w 0
./FGS+ -t 454_10 -s example/NC_000913-454.fna -o example/NC_000913-454 -w 0
./FragGeneScanRs -t 454_10 -s example/NC_000913-454.fna -o example/NC_000913-454 -w 0

./FragGeneScan -t complete -s example/contigs.fna -o example/contigs -w 1
./FragGeneScanRs -t complete -s example/contigs.fna -o example/contigs -w 1

./FragGeneScan -t complete -s example/NC_000913.fna -o example/NC_000913 -w 1
./FragGeneScanRs -t complete -s example/NC_000913.fna -o example/NC_000913 -w 1
```

By default, FragGeneScanPlus outputs only the predicted genes, not the
metadata and DNA files. Below are measurements taken when those files
aren't output by FragGeneScanRs either.

| Short reads      |  1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|----------:|----------:|----------:|----------:|
| FragGeneScanPlus |    1.194s |   0.6606s |   0.4764s |   0.4018s |
| FragGeneScanRs   |   0.9321s |   0.4994s |   0.3589s |   0.2817s |

The commands used here are:

```sh
./FGS+ -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 > /dev/null
./FragGeneScanRs -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 > /dev/null
```
