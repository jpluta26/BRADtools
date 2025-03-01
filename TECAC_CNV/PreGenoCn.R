##Create this script to filter the 14519 samples that failed QC.
#Rename the filtered samples with sample BID
#Set working directory
setwd("~/project/knathanslab/TECAC/CNV")
#Load Library----
library(vroom)
library(tidyverse)
library(dplyr)
library(stringr)
library(stringi)
library(pbapply)
library(filesstrings)

#------------------Functions-----------------------------------------
#A function for parsing over the headers and extract sampleBID in original order
extract_BID <- function(string){
  if (length(grep("^(.+)(\\.[a-zA-Z]{5})$", string, value = T) > 0)){
    return (unlist(str_split(string, "\\."))[1])
  }
}

#A function for matching BIDs of the total samples to the list of retained BIDs,
#And return a list of splitIDs
match_and_find_splitID <- function(bid){
  if(bid %in% retained_list$V1){
    return((total_sample %>% dplyr::filter(BID == bid))$SplitID)
  }
}

#Write splitIDs into filepaths pointing to the corresponding split results on HPC server
splitID_to_path <- function(id){
  filename <- paste(id, ".txt", sep = "")
  if(id %in% paste("sample", 1:2900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample1_2900", 
                     filename))
  }else if(id %in% paste("sample", 2901:3900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample2901_3900", filename))
  }else if(id %in% paste("sample", 3901:4900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample3901_4900", filename))
  }else if(id %in% paste("sample", 4901:5900,  sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample4901_5900", filename))
  }else if(id %in% paste("sample", 5901:6900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample5901_6900", filename))
  }else if(id %in% paste("sample", 6901:7900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample6901_7900", filename))
  }else if(id %in% paste("sample", 7901:8900,  sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample7901_8900", filename))
  }else if(id %in% paste("sample", 8901:9900, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample8901_9900", filename))
  }else if(id %in% paste("sample", 9901:10000, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample9901_10000", filename))
  }else if(id %in% paste("sample", 10001:11000, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample10001_11000", filename))
  }else if(id %in% paste("sample", 11001:12000,  sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample11001_12000", filename))
  }else if(id %in% paste("sample", 12001:13000, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample12001_13000", filename))
  }else if(id %in% paste("sample", 13001:14000, sep = "")){
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample13001_14000", filename))
  }else{
    return(file.path("/project/knathanslab/TECAC/CNV/ByBID_Sample14001_14519", filename))
  }
}

#A function for renameing retained sample data files with BIDs
rename_and_move <- function(filepath){
  splitid <- unlist(strsplit(basename(filepath), split = "\\."))[1]
  new_path <- file.path(dirname(filepath),
                                paste((total_sample %>% dplyr::filter(SplitID == splitid))$BID,
                   ".txt", sep = ""))
  file.rename(from = filepath, to = new_path)
  file.copy(new_path, "/project/knathanslab/TECAC/CNV/GenoCN_Samples")
  return(file.path("/project/knathanslab/TECAC/CNV/GenoCN_Samples", basename(new_path)))
}


#Read in the list of sample BIDs that will be retained 
retained_list <- as_tibble(read.table("TECAC_CNV_SUBLIST.txt", header = F, sep = "\t"))

#Read in the header of the tecac_all_cnv.csv file that contains all 14592 samples

headers <- scan(file = "header_all.txt", what = character(), sep = ",")
headers <- headers[c(4:length(headers))]
headers <- as.list(headers)
#Extract sample BIDs from the header_all.txt, and match splitID to BID
total_sample_BID <- unlist(sapply(headers, FUN = extract_BID))
total_sample <- as_tibble(data.frame(paste("sample",1:14519, sep = ''), total_sample_BID))
names(total_sample) <- c("SplitID", "BID")
#Acquire the list of retained splitID, and rewrite into filepaths pointing to files on HPC server.
retained_SplitID <- unlist(pblapply(total_sample$BID, 
                             FUN = match_and_find_splitID))
retained_paths <- unlist(pblapply(retained_SplitID, FUN = splitID_to_path))

#Rename selected sample data files and move them to the GenoCN sample folder
genocn_paths <- unlist(pblapply(retained_paths, FUN = rename_and_move))












