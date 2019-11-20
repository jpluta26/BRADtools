#!/bin/bash
# pluta 11/20/19

# input is a file of CRVs in chr:pos format in hg19
# script converts them into the format liftOver likes, and transforms to hg38

NPARAMS=$#

if [ $NPARAM -lt 1 2 ]
then
  echo ""
  echo "USAGE :: 
  echo "./convertCRVstoGRCH38.sh SNPLST OUTNAME"
  echo "SNPLST is a list of snps in chr:pos format, in hg19"
  echo "OUTNAME is the root output name"
  echo ""
  exit 1
fi

SNPLST=$1
OUTNAME=$2
HG19=${OUTNAME}_hg19.bed
HG38=${OUTNAME}_hg38.bed
more $SNPLST | tr ":" " " | awk '{print "chr" $1, $2 - 1, $2}' > ${HG19}
~/liftOver ${HG19} ~/hg19Tohg38.over.chain ${HG38} unmapped
