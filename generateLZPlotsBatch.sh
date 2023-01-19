#!/bin/bash

# get the number of input parameters
NPARAM=$#

# if user enters less than the required number of arguments, print the usage
if [ $NPARAM -lt 2 ]
then
	echo ""
	echo "USAGE :: "
	echo "./generateLZPlotsBatch.sh SNPLST WINDOWSIZE"
	echo "SNPLST is a text file containing the list of SNPs to plot."
	echo "SNPLST should be two columns: CHR BP"
	echo ""
	exit
fi


WINDOWSIZE=$2
WINDOWSIZE=${WINDOWSIZE}kb
#if [ {$WINDOWSIZE: -2} -ne "kb" ]
#then
#	WINDOWSIZE=$(echo ${WINDOWSIZE}kb)
#fi

# input is the list of snps to plot, in two columns- CHR and BP
while IFS='' read -r line || [[ -n "$line" ]]
do
	CHR=$(echo $line | awk '{print $1}')
	if [ $CHR -eq 23 ]
	then
		CHR="X"
	fi
	
	echo "LINE = $line"
	BP=$(echo $line | awk '{print $2}')
	REFSNP=$(echo $line | awk '{print $1 ":" $2}')
#	METAFILE=../Meta/meta-chr${CHR}.txt
	METAFILE=/project/knathanslab/TECAC/Meta/newVal_8site_chr${CHR}_ggman.txt	
	
	echo "BP = $BP"
	echo "REFSNP = $REFSNP"
	echo "METAFILE = $METAFILE"
	echo "CHR = $CHR"
	echo "WINDOWSIZE = $WINDOWSIZE"

	if [ ! -e $METAFILE ]
	then
		echo "$METAFILE not found! aborting"
		exit 1
	fi

	echo "plotting $REFSNP from $METAFILE..."
#	bsub -o ${CHR}.${BP}_LZ.o -e ${CHR}.${BP}_LZ.e sh ../scripts/plotLZ-LD.sh $METAFILE $REFSNP $WINDOWSIZE
	bsub -o ${CHR}.${BP}_LZ.o -e ${CHR}.${BP}_LZ.e sh /project/knathanslab/TECAC/scripts/plotLZ.sh $METAFILE $REFSNP $WINDOWSIZE
done < "$1"
