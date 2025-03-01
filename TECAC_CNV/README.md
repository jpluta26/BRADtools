# TECAC_CNV
This project aims at identifying germline Copy Number Variations (CNVs) asscociated with TCGT.

## Project Flow Chart ##

## Processing Steps ##

# CNVR Association
* **1. PennCNV** to generate .rawcnv files <br />
* **2. format_rawcnv.sh** to format .rawcnv file for further processing <br />
* **3. subsetRawcnv.sh** to split .rawcnv file into duplication and deletion files <br />
* **4. runCNVRassoc.R** map CNVs to CNVRs; run association testing and gene/cnvr mapping <br />
* **5. plotGeneInCNVR.R** to visualize CNVR/gene overlap <br />
<br />

# Correlation of CNV/Gene Expression
* **1. getTCGAMatrices.R** create matrices for gene-level CNV, cpg-methylation, and gene expression for all genes of interest, from TCGA data <br />
* **2. plotExpressionByCNV.R** create plots of log2(FPKM) adjusted for cpg methylation against CNV

## Tools ##
* Shell <br />
* **PennCNV** <br /> 
* **HandyCNV** <br />


## Supporting Files ##
* **TECAC_CNV_SUBLIST.txt**: The file containing BID and SSID for all samples that will be retained for CNV Calling <br />
<br />

* **TECAC_CNV_PHENO.txt**: Contains BID, SSID, and phenotype for all TECAC samples <br /> 
<br />

* **tecac.evec**: principle components for TECAC data, used as covariates <br />
<br />

* **header_all.txt**: The file containing all column names of the large, unsplit raw signal intensity file. <br />
<br />
