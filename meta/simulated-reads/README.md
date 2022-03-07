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
# mason_variator -ir ena_data_20210917-1328.fasta -ov variants.vcf
mason_simulator -ir ena_data_20210917-1328.fasta -n 10000 --illumina-read-length 1000 --fragment-mean-size 2000 --fragment-min-size 1000 -o reads.fasta -oa alignments.sam
```

# Create the FGS* files

(swap directories to execute these)

```sh
FragGeneScanRs -s reads.fasta -a FGSrs.faa -t illumina_10 -w 0
cd path/to/FGS
./FragGeneScan -s ~-/reads.fasta -o ~-/FGS -t illumina_10 -w 0
./FGS+ -s ~-/reads.fasta -o ~-/FGS+ -t illumina_10 -w 0
cd -
rm FGS.out FGS.ffn
```

## Print comparison table

You'll need to `pip install pysam`.

Execute `evaluate.py`.

```
tool            TP      FP      TN      FN    prec    sens    spec     NPV     MCC
FGS         63.59%  28.80%   4.63%   2.98%  68.83%  95.52%  13.86%  60.85%    0.17 2.91%
FGS+        63.92%  29.12%   4.31%   2.64%  68.70%  96.03%  12.89%  61.99%    0.17 2.91%
FGSrs       63.59%  28.79%   4.64%   2.98%  68.84%  95.52%  13.88%  60.89%    0.17 2.91%
```

## Timings for these predictions using [hyperfine](https://github.com/sharkdp/hyperfine)

Run in the FGS or FGS+ directory (for the training files).

```sh
hyperfine 'FragGeneScan -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGS -t complete -w 1' \
          'FGS+ -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGS+ -t complete -w 1' \
          'FragGeneScanRs -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGSrs -t complete -w 1' \
          'prodigal -i meta/evaluation/ena_data_20210917-1328.fasta -f gff -o meta/evaluation/prodigal.gff'
```

```
Benchmark #1: ./FragGeneScan -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGS -t complete -w 1
  Time (mean ± σ):      3.797 s ±  0.006 s    [User: 3.413 s, System: 0.348 s]
  Range (min … max):    3.792 s …  3.807 s    5 runs

Benchmark #2: ./FGS+ -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGS+ -t complete -w 1
  Time (mean ± σ):     369.979 s ± 25.774 s    [User: 367.679 s, System: 0.517 s]
  Range (min … max):   353.713 s … 415.649 s    5 runs

Benchmark #1: FragGeneScanRs -s meta/evaluation/ena_data_20210917-1328.fasta -o meta/evaluation/FGSrs -t complete -w 1
  Time (mean ± σ):      1.703 s ±  0.014 s    [User: 1.395 s, System: 0.275 s]
  Range (min … max):    1.684 s …  1.719 s    5 runs

Benchmark #4: prodigal -i meta/evaluation/ena_data_20210917-1328.fasta -f gff -o meta/evaluation/prodigal.gff
  Time (mean ± σ):      8.533 s ±  0.038 s    [User: 8.453 s, System: 0.047 s]
  Range (min … max):    8.493 s …  8.573 s    5 runs
```
