#!/bin/bash
#  pluta  7/13/21

#  parse a .rawcnv file into del and dup files

RAWCNV=$1
NAME=`basename $RAWCNV .rawcnv`

more $RAWCNV | awk -F " " '$4 == "state1,cn=0" || $4  == "state2,cn=1"' > ${NAME}.del.rawcnv
more $RAWCNV | awk -F " " '$4 == "state5,cn=3" || $4  == "state6,cn=4"' > ${NAME}.dup.rawcnv