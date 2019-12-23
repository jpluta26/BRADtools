# BRADtools
Biomedical Research and Development (BRAD) tools for genetics analysis. Code used in GWAS analysis. Descriptions:

MAIN/

*bash_script_template.sh*: template for bash scripts to be run on the HPC. handles obtaining user input and basic command error checking.

*bayes_tutorial.R*: runs some simple examples and plots explaining the fundamentals of Bayesian statistics.

*calcImputationValidation.R*: Calculate the concordance between select imputed and genotyped SNPs, for the TECAC replication paper.

*callerHeatmap.R*: Read in a file of caller summary statistics. For each caller, each variant is either PASS, REJECT, or NO CALL. Visualize these in a heatmap, and compare each caller to each other via chisq test. Currently assumes that the 4 callers used are MUTECT2, STRELKA2, VARDICTJAVA, VARSCAN2. Can be further generalized, if Brad needs me to.

*compare2ROC.R*: function to plot and statistically compare two ROC curves. the inputs come from the function "getROCstats"

*conditionalAnalysis.sh*: Run conditional analysis on the meta-analysis summary statistics using GCTA. The .tbl output comes from METAL. Invokes plotCond.R.

*convertGenoToPed.R*: Convert a raw genotype file from the CAG to plink PED format.

*getGeneCoordinates.R*: get gene coordinates for a given gene

*getHeritability.R*: calculate heritability for a given list of snps. The input is OR (or logOR) and corresponding MAF for a list of SNPs. lamba values are defined for TGCT. This method is used in Wang et al. 2017.

*getROCstats.R*: generate ROC curve w/ LOOCV

*hrd_wrapper.R*: either update this or delete- moved to HRDex

*manPlot_Fast.R*: subsample SNP data below a certain threshold to speed-up creation of manhattan plots

*plotCond.R*: Create conditional analysis plots. Invoked by conditionalAnalysis.sh

*plotPRS_OR.R*: Make a plot of Odds ratios with error bars over quantiles, from PRS scores

*preImputeCheck.sh*: verify that a file is valid before submitting to HHRC server

*reverseOddsRatio.R*: reverse odds ratio and confidence interval so that it is aligned to the alternate allele

*reverseOddsRatio.vba*: same as above, but works in Excel

*runPCA.sh*: PCA of genetic data to determine ancestry

*wakefield_bayes.R*: implement Wakefiled Baye's Factor
