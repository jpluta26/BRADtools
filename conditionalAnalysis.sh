#!/bin/sh

# pluta 3/21/18
# script to run conditional analysis using GCTA

NPARAM=$#

if [ $NPARAM -lt 2 ]
then
	echo ""
	echo "USAGE :: "
	echo "conditionalAnalysis.sh REFSNP WINDOW"
	echo "REFSNP: snp to condition on in chr:bp format"
	echo "WINDOW: flanking length in kbp"
	echo ""
	exit 
fi


# in chr:pos format
REFSNP=$1
WINDOW=$2

### CONSTANTS ###
TECACDIR="/project/knathanslab/TECAC"
METATBL="${TECACDIR}/Meta/TECAC_META_newVal_8site1.tbl"

# summary stats, this comes from the METATBL
COJOFILE="${TECACDIR}/Meta/GCTA-input3.ma"
#################

echo "Conditional analysis on snp $REFSNP"
echo "BEGIN"

echo "1. Set up variables..."

# GCTA needs text file input
# write reference snp to a temporary file
echo "$REFSNP" > cond.snplist

# avoid colons in filenames
CHR=$(echo $REFSNP | tr ":" " " | awk '{print $1}')
BP=$(echo $REFSNP | tr ":" " " | awk '{print $2}')
OUTNAME=${CHR}_${BP}


# imputed gwas data- this is the replication set only
# TODO: check if snp is actually in the replication data
# necessary to calculate LD
BFILE=${TECACDIR}/chr${CHR}/HRC/Imputed/chr${CHR} 

# get the p-value of the reference snp in the meta analysis
# but first have to make sure it exists in the data
tmp=$(LC_ALL=C grep -w "$REFSNP" $METATBL)

if [ $? != 0 ]
then
	echo "$REFSNP not found in $METATBL."
	echo "exiting"
	exit 1
fi

PVAL=$(echo $tmp | awk '{print $8}')

echo "done!"
echo ""



echo "2. Extracting snp list..."
# extract the list of snps with WINDOW flanking the reference snp
awk -v bp="$BP" -v win="$WINDOW" '$4 > bp - win * 1000 && $4 < bp + win * 1000 {print $2}' ${BFILE}.bim > SNPLST

echo "done!"
echo ""


echo "3. Conditional analysis..."

if [[ ${CHR} -eq "X" ]]
then
	CHR=23
fi

# run independence test with GCTA
CMD="/home/jpluta/gcta64 --bfile $BFILE
			 --chr $CHR 
			 --cojo-file $COJOFILE 
			 --cojo-cond cond.snplist 
			 --extract SNPLST
			 --cojo-collinear 0.8	
			 --out $OUTNAME"

eval $CMD > ${OUTNAME}_cojo.log

if [ $? != 0 ]
then
	echo "$CMD failed!"
	echo "exiting"
	exit 1
fi

echo "done!"
echo ""




# get ld
echo "4. Get LD with reference snp..."
CMD="/home/jpluta/gcta64 --bfile $BFILE 
			 --chr $CHR 
			 --ld cond.snplist 
			 --ld-sig 0.05 
			 --ld-wind 10000 
			 --extract SNPLST 
			 --out $OUTNAME"

eval $CMD > ${OUTNAME}_ld.log
if [ $? != 0 ]
then
	echo "$CMD failed!"
	echo "exiting"
	exit 1
fi

echo "done!"
echo ""



echo "5. Plotting data..."
CMD="Rscript ${TECACDIR}/Meta/plotCond.R $REFSNP $PVAL"

eval $CMD > ${OUTNAME}_Rplot.log
if [ $? != 0 ]
then
	echo "$CMD failed!"
	echo "exiting"
	exit 1
fi
 
echo "done!"

echo ""
echo "conditional analysis - successfully complete."

# clean up
rm cond.snplist
rm SNPLST
