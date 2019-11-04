#!/bin/bash

# pluta 10/14/19
# run customAnnotation on hpc

NPARAM=$#

if [ $NPARAM -lt 1 ]
then
        echo "USAGE :: "
        echo "garfieldCustomAnnotation.R BEDFILE"
        echo "BEDFILE is the file to draw annotation from"
        echo ""
        exit 1
fi


# dont run on the head node
HEADNODE="consign.hpc.local"
if [ "$HOSTNAME" == "$HEADNODE" ]
then
        echo "dont run scripts on the head node!"
        exit 1
fi

INFILE=$1

if [ ! -e $INFILE ]
then
        echo "$INFILE :: file not found"
        exit 1
fi

Rscript garfieldCustomAnnotation.R $INFILE
