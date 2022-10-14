# script dividing SCE into tissues
# Luise A. Seeker
# 20210118


###
# 



# load libraries

library(scater)
library(scran)
library(Seurat)
library(devtools)
library(dplyr)


setwd("/Users/lseeker/Documents/Work/HumanCellAtlas")

###### PART 4
###### SPLIT DATASET IN TISSUES AND SCRAN NORMALISATION

SCE_comb <- readRDS(file.path(getwd(), 
                              "splice_control_out",
                              "datasets",
                              "03_combined_matrices",
                              "combined_QC_filtered.RDS"))


# I noticed that there are 5 levels in the column "Tissue" instead of three, 
# because CB and BA4 may have an empty space behind their last letter or not,
# I am correcting this here

SCE_comb$Tissue <- ifelse(SCE_comb$Tissue == "BA4  ", "BA4",
                          ifelse(SCE_comb$Tissue == "CB ", "CB",
                                 SCE_comb$Tissue))



normaliseSCE <- function(sce_dataset, split_id_column){
  for(i in 1: length(levels(as.factor(sce_dataset[[split_id_column]])))){
    tissue <- levels(as.factor(sce_dataset[[split_id_column]]))[i]
    split_boul <- sce_dataset[[split_id_column]] == tissue
    split_sce <- sce_dataset[, split_boul]
    
    # scran normalisation
    ## calculation of deconvolution factors
    set.seed(100)
    clust_sce <- quickCluster(split_sce) 
    #deconv_sf <- calculateSumFactors(split_sce, cluster=clust_sce)
    split_sce <- computeSumFactors(split_sce, 
                                   cluster=clust_sce, 
                                   min.mean=0.1)
    
    split_sce <- logNormCounts(split_sce)
    
    dir.create(file.path(getwd(), 
                         "splice_control_out",
                         "datasets",
                         "04_scran_normalised",
                         paste(tissue)))
               
               
               saveRDS(split_sce, file.path(getwd(), 
                                          "splice_control_out",
                                          "datasets",
                                          "04_scran_normalised",
                                          paste(tissue),
                                          paste(tissue,
                                                "sce_fil_norm.RDS",
                                                sep = "_")))
  }
  
}


normaliseSCE(SCE_comb, "Tissue")



print("done")




