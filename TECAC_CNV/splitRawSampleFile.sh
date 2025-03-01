#!/bin/bash
#This script can be used to split the large raw signal intensity file.
#"n" can be adjusted to the actual number of samples in a file.
#Depends on the computation strength can split runs by adjusting n.
for ((n=1; n<=14159; n++))
do
SPL=$n
let start=3*$n+1
let end=3*$n+3
cut -d "," -f 1-3,$start-$end /project/knathanslab/TECAC/CNV/tecac_all_cnv.csv > /project/knathanslab/TECAC/CNV/ByBID/sample${SPL}.txt &
if (( (n % 1000) == 0 )); then
wait
fi
done
