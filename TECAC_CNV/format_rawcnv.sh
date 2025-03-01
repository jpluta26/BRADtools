#!/bin/bash

# script to prepare .rawcnv file for processing
# dos2unix; spaces instead of tabs; the prefix before
# each subject id needs to be uniform (cannot be cases and controls)
# HandyCNV expects an 8th column, but PennCNV 1.05 only produces 7 columns,
# add a dummy 8th column

# pluta 7/8/21


RAWCNV=$1

dos2unix $RAWCNV
more $RAWCNV | tr "\t" " " | awk '{print $0, "conf=1"}' > tmp

if grep  -q "cases" tmp
then
	more tmp | sed 's/cases/controls/g' > tmp2
	mv tmp2 tmp
fi


FNAME=$(basename $RAWCNV .rawcnv)

mv tmp ${FNAME}_fix.rawcnv
