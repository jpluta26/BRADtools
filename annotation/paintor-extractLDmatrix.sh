#!/bin/bash
# pluta 9/5/19

# step 1 of PAINTOR annotation pipeline

# given the top hit/reference snp, this script will extract the LD matrix of
# that snp's predefined credible set. it will also account for redundant snps
# (pairs of snps with LD=1)- this will make the annotation matrix singular
# and cause a crash

# the output is the ld matrix, and a .credSet file, which is used in the next steps
# of annotation
NPARAM=$#

if [ $NPARAM -lt 1 ]
then
        echo ""
        echo "USAGE :: "
        echo "./extraLDmatrix.sh SNPNAME"
        echo "where SNPNAME is the reference snp (top hit) of some credible set"
        echo "its assumed this snp is in the annotated CRV file"
        exit 1
fi

HEADNODE="consign.hpc.local"

if [ "$HOSTNAME" == "$HEADNODE" ]
then
  echo "dont run scripts on the head node"
  exit 1
fi

SNPNAME=$1

# use either chr:pos or chr_pos
SNPNAME=$(echo $SNPNAME | tr ":" "_")
CHR=$(echo $SNPNAME | tr "_" " " | awk '{print $1}')

# credset files include a header, but dont include the reference snp
# so the number of snps in the credset including reference should be
# the same size as the credSet file
NCREDSET=$(wc -l ${SNPNAME}_credSet.txt | awk '{print $1}')

# these files are already predefined to include snps with LD >= 0.8
# with reference snp; this comes from prior gwas steps
tail -n +2 ${SNPNAME}_credSet.txt | awk '{print $2}' > snplst1

# make sure the reference snp is added to the credible set
echo $SNPNAME >> snplst1
sed -i 's/_/:/g' snplst1

NSNPLST1=$(wc -l snplst1 | awk '{print $1}')

if [[ $NCREDSET -ne $NSNPLST1 ]]
then
                echo "number of snps in credset does not equal number of snps found"
                echo "there were $NCREDSET snps in the credset"
                echo "there were $NSNPLST1 snps found by in novelHits_CRVs_annot"
                echo "something went wrong"
                exit 1
fi

# extract the ld matrix from the LD data
echo "extracting ld matrix..."
CHRDIR=/project/knathanslab/TECAC/chr${CHR}/HRC/Imputed

# remove any duplicate snps
CMD="~/bcftools-1.3.1/bcftools norm -d snps --threads=12 chr${CHR}.qc.vcf.gz -O z -o tmp"
eval $CMD

# really dont want this command to file and overwrite the imputed data
if [ $? -ne 0 ]
then
  echo "$CMD failed. aborting"
  exit 1
fi

mv tmp chr${CHR}.qc.vcf.gz

CMD="plink1.09 --vcf ${CHRDIR}/chr${CHR}.qc.vcf.gz --extract snplst1 --double-id --r2 square --out $SNPNA
ME"
eval $CMD

if [ $? -ne 0 ]
then
  echo "$CMD failed. aborting"
  exit 1
fi

echo "done!"

echo "checking for redundant snps..."
Rscript removeRedundantSnps.R ${SNPNAME}.ld
echo "done!"

# compare the two credible sets before and after correcting for redundancy
NSNPLST2=$(wc -l snplst2 | awk '{print $1}')

# if some snps were removed due to redundancy,
if [[ $NSNPLST1 -gt $NSNPLST2 ]]
then
        SNPLST=snplst2
        NRM=$(expr $NSNPLST1 - $NSNPLST2)
        echo "$NRM snps removed"
else
        SNPLST=snplst1
fi

# novelHits_CRVs_annot.txt is the summary statistics of all top snps and corresponding
# credible sets
echo "writing credSet file for annotation..."
tail -n +2 novelHits_CRVs_annot.txt | LC_ALL=C grep -f $SNPLST  > tmp
cat novelHits-header tmp > ${SNPNAME}.credSet
rm tmp
echo "done"
