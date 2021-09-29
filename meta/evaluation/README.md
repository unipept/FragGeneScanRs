# Evaluation of FGS, FGSrs, FGS+, Prodigal on whole genomes

Source assembly: https://www.ebi.ac.uk/ena/browser/view/GCA_001628815?show=chromosomes

The 'FASTA' download `ena_data_20210917-1328.fasta` is the complete assembly.

The 'TEXT' download `ena_data_20210917-1328.txt` also contains annotated genes.

## Create the annotations and lengths files

Execute `annotations.py`.

## Create the FGS/FGS+ files (from .aa)

(swap directories to execute these)

```sh
cd path/to/FGS
./FragGeneScan -s ~-/ena_data_20210917-1328.fasta -o ~-/FGS -t complete -w 1
./FGS+ -s ~-/ena_data_20210917-1328.fasta -o ~-/FGS+ -t complete -w 1
cd -
rm FGS.out FGS.ffn
sed -n 's/^>ENA|\([^|]*\)|.*_\([0-9]*\)_\([0-9]*\)_\([+-]\)$/\1,\2,\3,\4/p' FGS.faa > FGS.csv
sed -n 's/^>ENA|\([^|]*\)|.*_\([0-9]*\)_\([0-9]*\)_\([+-]\)$/\1,\2,\3,\4/p' FGS+.faa > FGS+.csv
```

## Create the FGSrs/Prodigal files (from .gff)

```sh
FragGeneScanRs -s ena_data_20210917-1328.fasta -g FGSrs.gff -t complete -w 1
prodigal -i ena_data_20210917-1328.fasta -f gff -o prodigal.gff
grep -v '^#' FGSrs.gff | tr '\t' ',' | cut -d, -f1,4,5,7 | sed 's/ENA|//;s/|[^,]*,/,/' > FGSrs.csv
grep -v '^#' prodigal.gff | tr '\t' ',' | cut -d, -f1,4,5,7 | sed 's/ENA|//;s/|[^,]*,/,/' > prodigal.csv
```

## Print comparison table

Execute `rates.py`.

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
