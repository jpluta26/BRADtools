#This script is used to sort the previously filtered 10249 samples
#based on phenotype. 
#It also randomly pick 1000 case samples and generate a file path .txt file.
#Later this path file will be used to generate .pfb file by PennCNV.
#----Load Libraries
library(dplyr)
library(tidyverse)
library(filesstrings)
library(pbapply)
library(vroom)

#---------Reference Data------------
#Read in the PHENO files containing phenotype information
pheno_all <- read.csv(file = "/project/knathanslab/TECAC/CNV/TECAC_CNV_PHENO.csv")
pheno_all <- as_tibble(pheno_all)
#Separate cases from controls.
cases_all <- pheno_all %>% dplyr::filter(PHENO == 2)
controls_all <- pheno_all %>% dplyr::filter(PHENO == 1)

#----Function------------------------
#Read in sample BID and generate a path for the final destination of the sample file.
read_BID_and_move <- function(BID){
  old_path <- file.path("/project/knathanslab/TECAC/CNV/GenoCN_Samples", 
                        paste(BID, ".txt", sep = ""))
  if (BID %in% cases_all$BID){
    file.move(old_path, "/project/knathanslab/TECAC/CNV/GenoCN_Samples/cases")
    new_path <- file.path("/project/knathanslab/TECAC/CNV/GenoCN_Samples/cases", 
                          paste(BID, ".txt", sep = ""))
  }else{
    file.move(old_path, "/project/knathanslab/TECAC/CNV/GenoCN_Samples/controls")
    new_path <- file.path("/project/knathanslab/TECAC/CNV/GenoCN_Samples/controls", 
                          paste(BID, ".txt", sep = ""))
  }
  return(new_path)
}

#A function to read in comma deliminated files, and write out as tab deliminated.
#Argument: path <- the path of the comma-separated file to be transformed
#Output: A new path pointing to the transformed tab-deliminated files.
comma_to_tab <- function(path){
  buffer <- vroom::vroom(file = path, delim = ",")
  if (nrow(buffer) > 1){
    new_path <- file.path("/project/knathanslab/TECAC/CNV/GenoCN_Samples/For_PFB",
                          basename(path))
    vroom_write(buffer, new_path,
                delim = "\t", col_names = T)
  }
  return(new_path)
}

#-----------------------------------------------------------------------
#Sort all signal intensity files based on sample phenotype.
#Also produce a file storing the new path of every sample.
new_paths_list <- unlist(pblapply(pheno_all$BID, FUN = read_BID_and_move))
write.table(new_paths_list, file = "/project/knathanslab/TECAC/CNV/GenoCN/all_paths.txt",
            row.names = F, sep = "\n")

#Generate a list of 1000 controls for running th compile_pfb.pl program in PennCNV
#Write the list into a path .txt file
controls_path_all <- paste("/project/knathanslab/TECAC/CNV/GenoCN_Samples/controls/",
                           controls_all$BID, ".txt", sep = "")
controls_path_all <- data.frame(seq(1,4967), controls_path_all)
colnames(controls_path_all) <- c("ControlID", "Paths")
control_id_sub <- sample.int(4967, size = 1000)
controls_paths_sub <- controls_path_all[controls_path_all$ControlID %in% control_id_sub]

#Because compile_pfb.pl only reads tab deliminated files, change all selected control files 
#To tab-deliminated format, and copy to a new directory For_PFB.
control_paths_sub <- unlist(pblapply(control_paths_sub, FUN = comma_to_tab))
write.table(controls_paths_sub, "/project/knathanslab/TECAC/CNV/GenoCN_Samples/paths_for_pfb.txt",
            row.names = F, col.names = F, sep = "\n")







