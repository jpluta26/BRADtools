#!/bin/bash

# pluta 4/3/2017

# script to split plink file by chromosome, convert into vcf format, and check some parameters
# against HRC data (essentially pre phasing)

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
        echo ""
        exit
fi



# ---------------- constants -------------------------------------- #
# location of the perl and python scripts called in this program
SCRIPTDIR=/project/knathanslab/TECAC/scripts
REF=/project/knathanslab/REF
HRCDIR=/project/knathanslab/TECAC/chr${CHR}/HRC
# ----------------------------------------------------------------- #

BED=/project/knathanslab/TECAC/${INFILE}.bed
BIM=/project/knathanslab/TECAC/${INFILE}.bim
FAM=/project/knathanslab/TECAC/${INFILE}.fam

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
CMD="plink1.09 --bfile /project/knathanslab/TECAC/$INFILE
                   --chr $CHR
                   --recode
                   --make-bed
                   --out ${HRCDIR}/chr${CHR} --noweb"
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
plink1.09 --bfile chr${CHR} --freq --out maf --noweb
echo "align to HRC ..."

# compare phased data to HRC and apply corrections
# this part of the script is time and memory intensive
# the smallest chromosome, chr22, used 22244 MB of memory
# chr1 is the largest, and is about 6 times larger than chr22
CMD="perl ${SCRIPTDIR}/HRC-1000G-check-bim.pl -b chr${CHR}.bim
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
plink1.09 --bfile chr${CHR} --exclude Exclude-chr${CHR}-HRC.txt --make-bed --out TEMP1 --noweb
plink1.09 --bfile TEMP1 --update-map Chromosome-chr${CHR}-HRC.txt --update-chr --make-bed --out TEMP2
plink1.09 --bfile TEMP2 --update-map Position-chr${CHR}-HRC.txt --make-bed --out TEMP3 --noweb
plink1.09 --bfile TEMP3 --flip Strand-Flip-chr${CHR}-HRC.txt --make-bed --out TEMP4 --noweb
plink1.09 --bfile TEMP4 --reference-allele Force-Allele1-chr${CHR}-HRC.txt --make-bed --out chr${CHR}-fix --noweb

echo "done!"
rm TEMP*


echo "converting to VCF...."
CMD="plink1.09 --bfile chr${CHR}-fix --recode vcf-iid --out chr${CHR}"
eval $CMD

if [ $? -ne 0 ]
then
        echo "$CMD failed!"
        echo "Aborting"
        exit 1
fi

/home/jpluta/vcftools_0.1.13/bin/vcf-sort chr${CHR}.vcf | bgzip -c > chr${CHR}.vcf.gz


