

#!/bin/bash
# pluta 8/13/2024

# script to reduce 2 plink files to the common set of snps and merge them
module load plink/1.9-20210416

NPARAM=$#

if [ $NPARAM -lt 3 ]
then
	echo ""
	echo "USAGE :: "
	echo "FILE1 = the root name of plink file 1 (no .bed/bim/fam)"
	echo "FILE2 = the root name of plink file 2"
	echo "OUTNAME = the root name of the output file"
	echo ""
	exit
fi

FILE1=$1
FILE2=$2
OUTNAME=$3

more ${FILE1}.bim | awk '{print $2}' | sort > snps1
more ${FILE2}.bim | awk '{print $2}' | sort > snps2

comm -12 snps1 snps2 > common-snps

plink --bfile ${FILE1} --extract common-snps --make-bed --out ${FILE1}_comm
plink --bfile ${FILE2} --extract common-snps --make-bed --out ${FILE2}_comm

plink --bfile ${FILE1}_comm --bmerge ${FILE2}_comm --make-bed --out merged

plink --bfile ${FILE1}_comm --exclude merged-merge.missnp --make-bed --out ${FILE1}_comm_qc
plink --bfile ${FILE2}_comm --exclude merged-merge.missnp --make-bed --out ${FILE2}_comm_qc

plink --bfile ${FILE1}_comm_qc --bmerge ${FILE2}_comm_qc --make-bed --out ${OUTNAME}
