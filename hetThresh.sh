#!/bin/bash
# pluta 12/13/16

# bash script to select subjects with excessive heterozygosity
# any subjects with het rate > mean +/- 3 * SD are flagged


# get the heterozygosity rate and attach subject ids
# store in temp file X
tail -n +2 hetrate.het | awk '{if( $5 + 0 != 0) print $2, ($5 - $3)/$5}' > X

# hideous awk code to get the high and low thresholds from
# mean and standard deviation
STAT=($(awk '{for(i=1;i<=NF;i++) {sum += $i; sumsq += ($i)^2}}
	END { printf "%f %f \n", (sum/NR) - 3 * sqrt(sumsq/NR - (sum/NR)^2), (su
m/NR) + 3 * sqrt(sumsq/NR - (sum/NR)^2)}' X))

lo=${STAT[0]}
hi=${STAT[1]}

# find subjects that exceed these thresholds and write to het-rm.txt
awk '$2 < "'"$lo"'" || $2 > "'"$hi"'" {print 0,$1}' X > het-rm.txt

# cleanup
rm X
