#!/bin/sh

RUNS=5

printf "tool,type,reads,threads,time,all memory,resmem\n"
printf "tool,type,reads,threads,run,time\n" >&2

go() {
	name="$1"
	command="$2"
	sum=0
	for run in $(seq 1 "$RUNS"); do
		before="$(date +"%s%3N")"
		$command > /dev/null 2>&1
		after="$(date +"%s%3N")"
		printf "%s,%d,%d\n" "$name" "$run" "$((after - before))" >&2
		sum="$((sum + after - before))"
	done
	f="$(mktemp)"
	valgrind --tool=massif --pages-as-heap=yes --massif-out-file="$f" $command > /dev/null 2>&1
	mem="$(sed '/heap_tree=peak/,$d' < "$f" | sed -n 's/mem_heap_B=//p' | sed -n '$p')"
	valgrind --tool=massif --stacks=yes --massif-out-file="$f" $command > /dev/null 2>&1
	heap="$(sed '/heap_tree=peak/,$d' < "$f" | sed -n 's/mem_heap_B=//p' | sed -n '$p')"
	extra="$(sed '/heap_tree=peak/,$d' < "$f" | sed -n 's/mem_heap_extra_B=//p' | sed -n '$p')"
	stack="$(sed '/heap_tree=peak/,$d' < "$f" | sed -n 's/mem_stacks_B=//p' | sed -n '$p')"
	rm "$f"
	printf "%s,%d,%d,%d\n" "$name" "$((sum / RUNS))" "$mem" "$((heap + extra + stack))"
}

d="$(mktemp -d)"

go "FGS,complete,1,1" "./FragGeneScan -t complete -s example/NC_000913.fna -o $d/NC_000913 -w 1"
go "FGS+,complete,1,1" "./FGS+ -t complete -s example/NC_000913.fna -o $d/NC_000913 -w 1"
go "FGSrs,complete,1,1" "./target/release/FragGeneScanRs -t complete -s example/NC_000913.fna -o $d/NC_000913 -w 1"
go "FGSrsu,complete,1,1" "./target/release/FragGeneScanRs -u -t complete -s example/NC_000913.fna -o $d/NC_000913 -w 1"

for threadnum in $(seq 1 16); do
	go "FGS,long,19349,$threadnum" "./FragGeneScan -t complete -s example/contigs.fna -o $d/contigs -w 1 -p $threadnum"
	if [ "$threadnum" -lt 11 ]; then
		go "FGS+,long,19349,$threadnum" "./FGS+ -t complete -s example/contigs.fna -o $d/contigs -w 1 -p $threadnum"
	done
	go "FGSrs,long,19349,$threadnum" "./target/release/FragGeneScanRs -t complete -s example/contigs.fna -o $d/contigs -w 1 -p $threadnum"
	go "FGSrsu,long,19349,$threadnum" "./target/release/FragGeneScanRs -u -t complete -s example/contigs.fna -o $d/contigs -w 1 -p $threadnum"
done

for threadnum in $(seq 1 16); do
	go "FGS,short,23373,$threadnum" "./FragGeneScan -t 454_10 -s example/NC_000913-454.fna -o $d/NC_000913-454 -w 0 -p $threadnum"
	if [ "$threadnum" -lt 11 ]; then
		go "FGS+,short,23373,$threadnum" "./FGS+ -t 454_10 -s example/NC_000913-454.fna -o $d/NC_000913-454 -w 0 -p $threadnum"
	fi
	go "FGSrs,short,23373,$threadnum" "./target/release/FragGeneScanRs -t 454_10 -s example/NC_000913-454.fna -o $d/NC_000913-454 -w 0 -p $threadnum"
	go "FGSrsu,short,23373,$threadnum" "./target/release/FragGeneScanRs -u -t 454_10 -s example/NC_000913-454.fna -o $d/NC_000913-454 -w 0 -p $threadnum"
done

for threadnum in $(seq 1 10); do
	go "stdout FGS+,short,23373,$threadnum" "./FGS+ -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 -p $threadnum"
	go "stdout FGSrsu,short,23373,$threadnum" "./target/release/FragGeneScanRs -u -t 454_10 -s example/NC_000913-454.fna -o stdout -w 0 -p $threadnum"
done

rm -r "$d"
