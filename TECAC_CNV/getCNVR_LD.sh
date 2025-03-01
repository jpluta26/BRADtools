
#!/bin/bash

NPARAM=$#

if [[ "$NPARAM" -lt 3 ]]
then
	echo ""
	echo "USAGE :: "
	echo "./getCNVR_LD.sh CNVRFILE SNP GENOVCF CNVRID "
	echo "CNVRFILE: the matrix of subject x CNVR, created in runCNVRassoc.R"
	echo "SNP: the snp of interest in chr:pos format"
	echo "GENOVCF: vcf file containing genotypes"
	echo "CNVRID: CNVR of interest in the format CNVR_###"
	echo ""
	exit
fi

# make sure you're not on the head node
HEADNODE="consign.hpc.local"

if [ "$HOSTNAME" == "$HEADNODE" ]
then
  echo "dont run scripts on the head node"
  exit 1
fi

CNVRFILE=$1
SNP=$2
GENOVCF=$3
CNVRID=$4

# CHR=$(echo $SNP | awk -F ":" '{print $1}')
echo $SNP > snplst

more $CNVRFILE | awk '{print $NF "_" $NF}' | awk '{print $1,  $1}' > cnv_subs

plink1.09 --vcf $GENOVCF --keep cnv_subs --recodeAD --double-id --extract snplst --out tmp.geno

more tmp.geno.raw | awk '{print $1, $7}' > geno


Rscript getCNVR_LD.R $CNVRFILE $CNVRID geno
