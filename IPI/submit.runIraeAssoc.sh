#!/bin/bash
# pluta 10/4/20

# wrapper script for runAssocIrae.R

NPARAM=$#

BEDFILE=$1
COVARS=$2
COVFILE=$3
CHR=$4
MDL=$5
MAF=$6
OUTPREFIX=$7

if [ $NPARAM -lt 7 ]
then
	echo ""
	echo "USAGE :: "
	echo "./submitRunIPI.AssociationTests.sh BEDFILE COVARS COVFILE CHR MDL OUTPREFIX"
	echo "BEDFILE - the bedfile, corresponding .fam and .bim must exist"
	echo "COVARS - covariates in comma seperated format, eg ECOG,Stage,NDose"
	echo "COVFILE - covariate file, eg bms2.pheno.txt"
	echo "CHR - chromosome (integer valued)"
	echo "MDL - the association model"
	echo "1: irae, no prior; 2: irae, prior int; 3: cox ph, no prior; 4: cox ph, prior int"
	echo "MAF - MAF file corresponding to the BED file"
	echo "OUTPREFIX - prefix to attach to output files"
	echo ""
	exit
fi

Rscript ~/IPI/runIraeAssoc.R $1 $2 $3 $4 $5 $6 $7
