#!/bin/sh

d="$(mktemp -d)"

for flags in "-t 454_10 -s example/NC_000913-454.fna -o $d/NC_000913-454 -w 0" \
             "-t complete -s example/contigs.fna -o $d/contigs -w 1" \
             "-t complete -s example/NC_000913.fna -o $d/NC_000913 -w 1" \
do
	for threadnum in $(seq 1 4); do
		for program in ./FragGeneScan ./target/release/FragGeneScanRs ./FGS+; do
			printf "%s,%s,%d" "$flags" "$program" "$threadnum"
			for run in $(seq 1 10); do
				before="$(date +"%s%3N")"
				$program $flags -p "$threadnum" > /dev/null 2>&1
				after="$(date +"%s%3N")"
				printf ",%d" "$((after - before))"
			done
			printf "\n"
		done
	done
done

rm -r "$d"

for threadnum in $(seq 1 4); do
	for program in ./target/release/FragGeneScanRs ./FGS+; do
		printf "NC_000913-454 faa,%s,%d" "$program" "$threadnum"
		for run in $(seq 1 10); do
			before="$(date +"%s%3N")"
			$program -w 0 -t 454_10 -s example/NC_000913-454.fna -o stdout -p "$threadnum" > /dev/null 2>&1
			after="$(date +"%s%3N")"
			printf ",%d" "$((after - before))"
		done
		printf "\n"
	done
done
