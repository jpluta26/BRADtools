#!/bin/bash

# Run penncnv commands by specifying .pfb file, .list file and the output name. Please upload this script to the penncnv-1.0.5 folder, and run it from there.
#The list file should containing paths to the sample signal intensity files that will be analyzed. Each line should ontain only one path.
#After finishing the penncnv detection command, this script will combine segments in the detection results that are likely split by the detection algorithm.
#Both the joined by snp (default) option and the joined by bp option will be used.
#It will also run a command to generate the .qcsum file that will be used later.
#Example: ./runPennCNV.sh Infinium24HumanCore_KLN.pfb ListOfSamples.txt penncnv_run3
PFB=$1
LIST=$2
OUTNAME=$3


perl detect_cnv.pl -test -hmm lib/exome.hmm -pfb $PFB --list $LIST -log ${OUTNAME}.log -out ${OUTNAME}.rawcnv

perl clean_cnv.pl combineseg --signalfile $PFB ${OUTNAME}.rawcnv > ${OUTNAME}._joinedBySnp.rawcnv

perl clean_cnv.pl combineseg --signalfile $PFB ${OUTNAME}.rawcnv --bp > ${OUTNAME}._joinedBySnp.rawcnv

perl filter_cnv.pl ${OUTNAME}.rawcnv -qclogfile ${OUTNAME}.log 
