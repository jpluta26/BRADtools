#!/bin/bash
module load bcftools/1.20
module load htslib/1.20


NPARAM=$#

if [ $NPARAM -lt 3 ]
then
	echo ""
	echo "USAGE :: "
	echo "./extractDosage.sh INFILE SNPLST OUTFILE ... "
	echo "INFILE - the imputed file with dosage, usually .dose.vcf.gz"
	echo "SNPLST - the snps of interest with tab separated columns for Chr and position"
	echo "OUTFILE - the name of the output file"
	echo ""
	exit
fi

INFILE=$1
SNPLST=$2
OUTFILE=$3


if [ ! -e $INFILE ]
then
	echo " $INFILE :: file not found."
	exit 1
fi

if [ ! -e $SNPLST ]
then
	echo " $SNPLST :: file not found."
	exit 1
fi

bcftools view ${INFILE} -R $SNPLST | bcftools query -f '%CHROM\t%POS\t%ID\t[%DS\t]\n' > $OUTFILE
