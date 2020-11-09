#!/bin/bash

# pluta 2/7/16

# run PCA on genotype data. Data should be truncated to independent snps.
# output is vectors of eigen values
NPARAM=$#

if [ $NPARAM -lt 2 ]
then
    	echo ""
        echo "USAGE :: "
        echo "./runPCA plink.ROOT OUTNAME"
        echo "this should be the plink file containing only SNPs in ld
        echo "e.g. the output from runLD.sh"
        echo ""
        exit 1
fi



BFILE=$1
OUTNAME=$2
EIGPATH=/home/jpluta/EIG-6.1.4/bin
OUTDIR=/project/knathanslab/TECAC/PCA


BED=${1}.bed
BIM=${1}.bim
FAM=${1}.fam

# setup files for EIG
echo "converting to eigenstrat format..."

# create a parameter file for convertf
echo "genotypename:     $BED"   > ${OUTDIR}/conv.par
echo "snpname:          $BIM"   >> ${OUTDIR}/conv.par
echo "indivname:        $FAM"   >> ${OUTDIR}/conv.par
echo "outputformat:     EIGENSTRAT"             >> ${OUTDIR}/conv.par
echo "genotypeoutname:  ${OUTDIR}/${OUTNAME}.eigenstratgeno" >> ${OUTDIR}/conv.par
echo "snpoutname:	${OUTDIR}/${OUTNAME}.snp"     >> ${OUTDIR}/conv.par
echo "indivoutname:     ${OUTDIR}/${OUTNAME}.ind"     >> ${OUTDIR}/conv.par
echo "familynames:	NO"                     >> ${OUTDIR}/conv.par
echo "pordercheck:	NO"                     >> ${OUTDIR}/conv.par

CMD='${EIGPATH}/convertf -p ${OUTDIR}/conv.par'
eval $CMD > ${OUTDIR}/pca.log

if [ $? -ne 0 ]
then
    	echo "$CMD failed! Aborting."
        exit
fi

echo "Done!"
echo ""

# make sure the conversion actually worked
if [ ! -e ${OUTDIR}/${OUTNAME}.eigenstratgeno ];
then
    	echo "file '${OUTNAME}.eigenstratgeno' not found!"
        exit 1
fi

echo "Calculating eigenvectors..."


# 
NEIGVEC=3

echo "genotypename: ${OUTDIR}/${OUTNAME}.eigenstratgeno"  > ${OUTDIR}/pca.par
echo "snpname:      ${OUTDIR}/${OUTNAME}.snp"            >> ${OUTDIR}/pca.par
echo "indivname:    ${OUTDIR}/${OUTNAME}.ind"            >> ${OUTDIR}/pca.par
echo "evecoutname:  ${OUTDIR}/${OUTNAME}.evec"           >> ${OUTDIR}/pca.par
echo "evaloutname:  ${OUTDIR}/${OUTNAME}.eval"           >> ${OUTDIR}/pca.par
echo "altnormstyle: NO"                            >> ${OUTDIR}/pca.par
echo "numoutevec:   $NEIGVEC"                      >> ${OUTDIR}/pca.par
echo "numoutlieriter: 0"                           >> ${OUTDIR}/pca.par
echo "numoutlierevec: 10"                          >> ${OUTDIR}/pca.par
echo "outliersigmathresh: 6.0"                     >> ${OUTDIR}/pca.par
echo "qtmode: 0"                                   >> ${OUTDIR}/pca.par
echo "nsnpldregress: 2"                            >> ${OUTDIR}/pca.par
echo "outlieroutname: pca-outlier"                 >> ${OUTDIR}/pca.par
echo "familynames:  NO"                            >> ${OUTDIR}/pca.par
echo "grmoutname:   ${OUTDIR}/${OUTNAME}-grmout"          >> ${OUTDIR}/pca.par


CMD='${EIGPATH}/smartpca -p ${OUTDIR}/pca.par' > ${OUTDIR}/pca.log

# arguments:
# -i :the genotype file from convertf
# -a :SNP names
# -b :subject names
# -o :output eigenvectors
# -p :plots the output, requires gnuplot
# -e :output eigenvalues
# -l :log, including individuals defined as outliers
# -m :number of outlier removal iterations
# -t :number of components from which outliers should be removed
# -k :number of components to compute
# -s :minimum number of standard deviations from mean to be considered an outli$

eval $CMD >> ${OUTDIR}/pca_output.log

if [ $? -ne 0 ]
then
    	echo "$CMD failed! Aborting"
        exit 1
fi

echo "Done!"
echo ""



