# Evaluation of FGS, FGSrs, FGS+, Prodigal on whole genomes

Source assembly: https://www.ebi.ac.uk/ena/browser/view/GCA_001628815?show=chromosomes

The 'FASTA' download `ena_data_20210917-1328.fasta` is the complete assembly.

The 'TEXT' download `ena_data_20210917-1328.txt` also contains annotated genes.

## Simulating reads 

Uses https://github.com/seqan/seqan/tree/master/apps/mason2 (cmake from repository root, takes 2 hours)

```yaml
authors:
- family: Holtgrewe
  given: M.
title: "Mason - A Read Simulator for Second Generation Sequencing Data"
url: "http://publications.imp.fu-berlin.de/962/"
year: 2010
container-title: "Technical Report FU Berlin"
```

```sh
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 100 --fragment-mean-size 2000 --fragment-min-size 1000 -o 100.fasta -oa 100.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 200 --fragment-mean-size 2000 --fragment-min-size 1000 -o 200.fasta -oa 200.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 300 --fragment-mean-size 2000 --fragment-min-size 1000 -o 300.fasta -oa 300.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 400 --fragment-mean-size 2000 --fragment-min-size 1000 -o 400.fasta -oa 400.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 500 --fragment-mean-size 2000 --fragment-min-size 1000 -o 500.fasta -oa 500.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 600 --fragment-mean-size 2000 --fragment-min-size 1000 -o 600.fasta -oa 600.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 700 --fragment-mean-size 2000 --fragment-min-size 1000 -o 700.fasta -oa 700.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 800 --fragment-mean-size 2000 --fragment-min-size 1000 -o 800.fasta -oa 800.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 900 --fragment-mean-size 2000 --fragment-min-size 1000 -o 900.fasta -oa 900.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 1000 --fragment-mean-size 2000 --fragment-min-size 1000 -o 1000.fasta -oa 1000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 2000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 2000.fasta -oa 2000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 3000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 3000.fasta -oa 3000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 4000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 4000.fasta -oa 4000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 5000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 5000.fasta -oa 5000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 6000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 6000.fasta -oa 6000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 7000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 7000.fasta -oa 7000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 8000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 8000.fasta -oa 8000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 9000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 9000.fasta -oa 9000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 10000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 10000.fasta -oa 10000.sam
```

# Create the FGS* files

(swap directories to execute these)

```sh
FragGeneScanRs -s 2000.fasta -a FGSrs2000.faa -t illumina_10 -w 0
cd path/to/FGS
./FragGeneScan -s ~-/2000.fasta -o ~-/FGS2000 -t illumina_10 -w 0
./FGS+ -s ~-/2000.fasta -o ~-/FGS+2000 -t illumina_10 -w 0
cd -
rm FGS2000.out FGS2000.ffn
```

## Print comparison table

You'll need to `pip install pysam`. Execute `evaluate.py`.

```
tool            TP      FP      TN      FN    prec    sens    spec     NPV     MCC
FGS         63.59%  28.80%   4.63%   2.98%  68.83%  95.52%  13.86%  60.85%    0.17 2.91%
FGS+        63.92%  29.12%   4.31%   2.64%  68.70%  96.03%  12.89%  61.99%    0.17 2.91%
FGSrs       63.59%  28.79%   4.64%   2.98%  68.84%  95.52%  13.88%  60.89%    0.17 2.91%
```

## Timings for these predictions using [hyperfine](https://github.com/sharkdp/hyperfine)

```
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 100 --fragment-mean-size 2000 --fragment-min-size 1000 -o 100.fasta -oa 100.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 200 --fragment-mean-size 2000 --fragment-min-size 1000 -o 200.fasta -oa 200.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 300 --fragment-mean-size 2000 --fragment-min-size 1000 -o 300.fasta -oa 300.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 400 --fragment-mean-size 2000 --fragment-min-size 1000 -o 400.fasta -oa 400.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 500 --fragment-mean-size 2000 --fragment-min-size 1000 -o 500.fasta -oa 500.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 600 --fragment-mean-size 2000 --fragment-min-size 1000 -o 600.fasta -oa 600.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 700 --fragment-mean-size 2000 --fragment-min-size 1000 -o 700.fasta -oa 700.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 800 --fragment-mean-size 2000 --fragment-min-size 1000 -o 800.fasta -oa 800.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 900 --fragment-mean-size 2000 --fragment-min-size 1000 -o 900.fasta -oa 900.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 1000 --fragment-mean-size 2000 --fragment-min-size 1000 -o 1000.fasta -oa 1000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 2000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 2000.fasta -oa 2000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 3000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 3000.fasta -oa 3000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 4000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 4000.fasta -oa 4000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 5000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 5000.fasta -oa 5000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 6000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 6000.fasta -oa 6000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 7000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 7000.fasta -oa 7000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 8000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 8000.fasta -oa 8000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 9000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 9000.fasta -oa 9000.sam
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 10000 --fragment-mean-size 20000 --fragment-min-size 1000 -o 10000.fasta -oa 10000.sam
```

Run in the FGS or FGS+ directory (for the training files).

```
hyperfine -L length 100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000 \
    './FragGeneScan -s meta/simulated-reads/{length}.fasta -o meta/simulated-reads/FGS{length} -t illumina_10 -w 0' \
    './FGS+ -s meta/simulated-reads/{length}.fasta -o meta/simulated-reads/FGS+{length} -t illumina_10 -w 0' \
    'FragGeneScanRs -s meta/simulated-reads/{length}.fasta -a meta/simulated-reads/FGSrs{length} -t illumina_10 -w 0' \
    --export-csv hyperfine-output.csv --export-json hyperfine-output.json --export-markdown hyperfine-output.md
```
