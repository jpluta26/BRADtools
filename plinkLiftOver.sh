#!/bin/bash

# pluta 1/12/26
# run liftOver on plink files
module load plink/1.9-20210416 
module load liftOver

# currently only goes from b37 to b38
IN=$1

# remove duplicated variants
more ${IN}.bim | awk '{print $2}' | sort | uniq -d > duplicated-vars
plink --bfile ${IN} --exclude duplicated-vars --make-bed --out ${IN}

# attach "chr" to match b38 designations in chain file
awk 'BEGIN{OFS="\t"}{chr=$1; pos=$4; print "chr"chr, pos-1, pos, "chr"$2}' ${IN}.bim > ${IN}.b37.bed
liftOver ${IN}.b37.bed /project/knathans_tecac/REF/hg19ToHg38.over.chain ${IN}.b38.bed ${IN}.unmapped.b
ed
awk 'BEGIN{OFS="\t"}{
  chr=$1; gsub(/^chr/,"",chr);
  snp=$4;
  bp=$3;
  print snp, chr, bp
}' ${IN}.b38.bed > ${IN}.lift.map

awk '{print $4}' ${IN}.unmapped.bed > ${IN}.unmapped.snps

plink --bfile ${IN} \
  --exclude ${IN}.unmapped.snps \
  --update-chr ${IN}.lift.map 2 \
  --update-map ${IN}.lift.map 3 \
  --make-bed \
  --out ${IN}.b38
