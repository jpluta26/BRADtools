#!/bin/bash

# pluta 1/24/2024

# script to split plink file by chromosome, convert into vcf format, and check some parameters
# against HRC data (essentially pre phasing)

# for use on the sanger server, which required updating several plink commands
module load htslib/1.9
module load plink/1.9-20210416
module load vcftools/0.1.16
module load perl/5.20.2 
# phasing is handled by the HRC

# get user defined input
NPARAM=$#

INFILE=$1
CHR=$2


if [ $NPARAM -lt 2 ]
then
	echo ""
	echo "USAGE :: "
	echo "./preImputeCheck.sh INFILE CHR"
	echo "this program uses a lot of memory. Set the -M >= 25000. This is the minimum, for chr22. Larger chromosomes might need much more memory."
	echo "run this from inside the Phasing directory"
	echo "TO DO: update the python scripts to support dynamically defined directories."
	echo ""
	exit
fi



# ---------------- constants -------------------------------------- #
# location of the perl and python scripts called in this program
REF=/project/knathans_tecac/REF/HRC
HRCDIR=chr${CHR}/HRC
# ----------------------------------------------------------------- #

BED=${INFILE}.bed
BIM=${INFILE}.bim
FAM=${INFILE}.fam

# do the necessary plink files exist?


if [ ! -d $HRCDIR ]
then
	mkdir $HRCDIR
fi

for i in $BED $BIM $FAM
do
	if [ ! -f ${i} ]
	then
		echo "$i not found!"
		exit 1	
	fi
done

echo "PLINK files found"
echo "extracting chr${CHR}..."

# extract chromosome

CMD="plink --bfile $INFILE 
                   --chr $CHR 
                   --recode 
                   --make-bed 
		   
                   --out ${HRCDIR}/chr${CHR}"
eval $CMD

if [ $? -ne 0 ]
then
	echo "$CMD failed! Aborting."
	exit 1
else
	echo "Done!"
fi


cd $HRCDIR

# get MAF
plink --bfile chr${CHR} --freq --out maf 

echo "align to HRC ..."

# compare phased data to HRC and apply corrections
# this part of the script is time and memory intensive
# the smallest chromosome, chr22, used 22244 MB of memory
# chr1 is the largest, and is about 6 times larger than chr22
CMD="perl ${REF}/HRC-1000G-check-bim-NoReadKey.pl -b chr${CHR}.bim 
                                              -f maf.frq 
                                              -r ${REF}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"
eval $CMD

if [ $? -ne 0 ]
then
	echo "alignment failed!"
	echo "Aborting."
	exit 1
else
	echo "Done!"
fi

if [ ! -f Exclude-chr${CHR}-HRC.txt ]
then
	echo "HRC-1000G-check-bim.pl output missing."
	echo "Aborting"
	exit 1
fi


echo "Correcting PLINK files..."

# these files are all output from the above script
plink --bfile chr${CHR} --exclude Exclude-chr${CHR}-HRC.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-chr${CHR}-HRC.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-chr${CHR}-HRC.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-chr${CHR}-HRC.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-chr${CHR}-HRC.txt --make-bed --out chr${CHR}-fix --output-chr M

echo "done!"
rm TEMP*
rm temp*


# use this call for imputing on the Sanger server
# for michigan i have been using the same line without the --a2-allele flag
echo "converting to VCF...."
CMD="plink --bfile chr${CHR}-fix --recode vcf-iid --a2-allele Force-Allele1-chr${CHR}-HRC.txt --output-chr M --out chr${CHR}"
eval $CMD

if [ $? -ne 0 ]
then 
	echo "$CMD failed!"
	echo "Aborting"
	exit 1
fi

vcf-sort chr${CHR}.vcf | bgzip -c > chr${CHR}.vcf.gz

