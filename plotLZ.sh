
#!/bin/sh

# get the number of input parameters

module load python/2.7.9
NPARAM=$#

# if user enters less than the required number of arguments, print the usage
if [ $NPARAM -lt 3 ]
then
	echo ""
	echo "USAGE :: "
	echo "./plotLZ.sh PLOTFILE REFSNP WINDOWSIZE "
	echo "PLOTFILE is the file with snps and p-values"
	echo "REFSNP is the reference snp"
	echo "WINDOWSIZE: flanking region around ref snp (500kb)"
	echo ""
	exit
fi

# assign inputs to variable names
# $1 is the first variable entered, $2 is the second, etc
PLOTFILE=$1
REFSNP=$2
WINDOWSIZE=$3

# STOP if any of the input files are not found
if [ ! -e $PLOTFILE ]
then
	echo " $PLOTFILE :: file not found."
	exit 1
fi

echo "Running locuszoom with the following parameters"
echo "Data File: $PLOTFILE"
echo "Reference SNP: $REFSNP"
echo "Window Size: $WINDOWSIZE"
echo ""

# setup a command
CMD="/home/jpluta/locuszoom/bin/locuszoom --metal $PLOTFILE 
					  --cache None
					  --build hg19 
					  --pop EUR 
                                          --source 1000G_March2012 
                                          --refsnp $REFSNP
                                          --delim ' ' 
					  --flank $WINDOWSIZE"
                                         # --plotonly"

eval $CMD

if [ $? -ne 0 ]
then
	echo "$CMD failed. Aborting."
	exit 1
fi

echo "Done!"
