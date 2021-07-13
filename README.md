# FragGeneScanRs

## Installation

### From source

Run `cargo install --path .` in this directory after cloning it.

### From release

To do.

## Usage

You can use FragGeneScanRs with the short options of FragGeneScan, as
listed below. It also provides some additional options and long-form
options, listed below that.

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
Intel(R) Core(TM) i7-6500U CPU @ 2.50GHz with 8GB RAM. The datasets
used are the example datasets provided by FragGeneScan. The table
below shows the average execution time of 10 runs. Detailed results
may be found in `meta/benchmark.csv`. For the short reads (80bp),
FragGeneScanRs is about 10 times faster than FragGeneScan and twice
as fast as FragGeneScanPlus. For the long reads (1328bp) and the
complete genome (Escherichia coli str. K-12 substr. MG1655, 4639675bp),
FragGeneScanRs is 3.5 and 2.4 times faster than FragGeneScan (FGS+
crashes on complete reads).

| Short reads      | 1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|---------:|----------:|----------:|----------:|
| FragGeneScan     | 17.4531s |   9.2082s |   9.3818s |   7.3247s |
| FragGeneScanPlus |  1.7709s |   0.9410s |   0.8632s |   0.8170s |
| FragGeneScanRs   |  1.4469s |   0.7865s |   0.7256s |   0.6668s |

| Long reads       | 1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|---------:|----------:|----------:|----------:|
| FragGeneScan     | 46.4709s |  25.4329s |  27.2826s |  22.0585s |
| FragGeneScan     | 13.1824s |   6.7994s |   6.2426s |   6.0322s |

| Complete genome  | 1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|---------:|----------:|----------:|----------:|
| FragGeneScan     |  4.9669s |   5.0641s |   5.0420s |   5.0722s |
| FragGeneScanRs   |  2.0493s |   2.0635s |   2.0761s |   2.1108s |

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

By default, FragGeneScanPlus outputs only the predicted genes, not the metadata and DNA files. Below are measurements taken when those files aren't output by FragGeneScanRs either.

| Short reads      | 1 thread | 2 threads | 3 threads | 4 threads |
|:-----------------|---------:|----------:|----------:|----------:|
| FragGeneScanPlus |  1.7911s |   0.9610s |   0.8821s |   0.8283s |
| FragGeneScanRs   |  1.3573s |   0.7333s |   0.6949s |   0.6268s |

The commands used here are:

```sh
./FGS+ -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 > /dev/null
./FragGeneScanRs -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 > /dev/null
```
