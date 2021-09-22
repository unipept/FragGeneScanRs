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
prodigal -i ena_data_20210917-1328.fasta -p meta -f gff -o prodigal.gff
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
          'prodigal -i meta/evaluation/ena_data_20210917-1328.fasta -p meta -f gff -o meta/evaluation/prodigal.gff'
```
