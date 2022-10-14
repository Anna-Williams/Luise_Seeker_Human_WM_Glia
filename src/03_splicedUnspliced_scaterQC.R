# script for the quality control of 
# spliced and unspliced matrices seperately
# combining of both matrices
# normalisation with scran
# Luise A. Seeker
# 20200418


###
# use scripts in following order
# 1) saveSCE_splicedUnspliced.R for converting loom to SCE
# 2) splicedUnsplicedQC.R (optional) (generated QC plots)
# 3) splicedUnspliced_scaterQC.R for QC


# load libraries

library(scater)
library(scran)
library(Seurat)
library(devtools)
library(dplyr)


# set working directory




setwd("/Users/lseeker/Documents/Work/HumanCellAtlas")


# Define functions

# originally I looked for which genes to keep intividually 
# for each sample, but I think this is not neccessary. 
# simply using the combined splided and unspliced matrices
# containing all samples should be fine. The filtering threshold 
# may be increased for that a bit. 

generateKeepGeneList<-function(splicedSCE, 
                               unsplicedSCE, 
                               combExprThres = 200){
  
  # combined expression in spliced and unspliced matrices is considered 
  rowSumsCombined <- nexprs(splicedSCE, 
                            byrow= TRUE) + 
    nexprs(unsplicedSCE, byrow = TRUE)
  
  keep_feature <- rowSumsCombined > combExprThres
  
  # The next step is done with splicedSCE only, because rownames of both
  # matrices are the same (as long as no filter was applied upstream
  # which would have caused to throw and error when rowSumsCombined
  # were calculated.
  
  # generate list of genes to keep downstream:
  keep_genes <- rownames(splicedSCE)[keep_feature]
  
  return(keep_genes)
}



#set paths to input data. Those input directories should have
#been created by previous script (saveSCE_splicedUnspliced.R)

pathToSplicedSCE <- file.path(getwd(), 
                              "splice_control_out",
                              "datasets",
                              "raw_matrices",
                              "spliced")



pathToUnsplicedSCE <-  file.path(getwd(), 
                                 "splice_control_out",
                                 "datasets",
                                 "raw_matrices",
                                 "unspliced") 



# Create output folders


dir.create(file.path(getwd(), 
                     "splice_control_out", 
                     "qc_plots", 
                     "scater_plots"), 
           recursive = TRUE)

dir.create(file.path(getwd(), 
                     "splice_control_out",
                     "datasets",
                     "01_gene_filtered"))

dir.create(file.path(getwd(), 
                     "splice_control_out",
                     "datasets",
                     "02_scater_cell_filtered"))

dir.create(file.path(getwd(), 
                     "splice_control_out",
                     "datasets",
                     "03_combined_matrices"))

dir.create(file.path(getwd(), 
                     "splice_control_out",
                     "datasets",
                     "04_scran_normalised"))

dir.create(file.path(getwd(),
                     "splice_control_out",
                     "qc_data",
                     "scater"))

dir.create(file.path(getwd(),
                     "splice_control_out",
                     "processing_plots"))




#####################
#####################
#####################
## PART 1#############
# perform gene QC and save filtered matrices to be used 
# in future downstream analysis

# read in matrices that contain all spliced and all unspliced 
# data respectively and create a list of genes that are to be kept. 

path_list <- c(pathToSplicedSCE, pathToUnsplicedSCE)


spliced_filename <- list.files(path = path_list[1], pattern = ".RDS")
unspliced_filename <- list.files(path = path_list[2], pattern = ".RDS")


spliced_sce <- readRDS(file.path(path_list[1], spliced_filename))
unspliced_sce <- readRDS(file.path(path_list[2], unspliced_filename))


keep_gene_list <- generateKeepGeneList(spliced_sce, 
                                       unspliced_sce, 
                                       combExprThres = 200)


dim(spliced_sce)# number of genes before filtering
length(keep_gene_list) # number of genes after filtering

## use the keep_gene_list to subset spliced and unspliced 
# datasets for genes that are expressed more than the combined
# threshold 

filter_boulean <- rownames(spliced_sce) %in% keep_gene_list

spliced_sce<-spliced_sce[filter_boulean,]
unspliced_sce<-unspliced_sce[filter_boulean,]




# save SCE objects that were filtered for not/ very lowly
# expressed genes

spliced_path <- file.path(getwd(),
                          "splice_control_out",
                          "datasets",
                          "01_gene_filtered",
                          paste("spliced_SCE_gene_filtered.RDS"))
unspliced_path <- file.path(getwd(),
                            "splice_control_out",
                            "datasets",
                            "01_gene_filtered",
                            paste("unspliced_SCE_gene_filtered.RDS"))


saveRDS(spliced_sce, spliced_path)


saveRDS(unspliced_sce, unspliced_path)




######## PART 2
######## Scater Quality control to remove low quality nuclei

file_path_list <- c(spliced_path, unspliced_path)

for(i in 1: length(file_path_list)){
  # saving and then loading the file again is slow. As an alternative
  # I tried to loop through a list of SCE objects, however this does not
  # seem to work. The result is alwacs a SCE with only 1 gene.
  
  SCE <- readRDS(file_path_list[i])
  if(i == 1){
    splicedUnspliced <- "spliced"
    print("looking at spliced data")
  }else if (i == 2){
    splicedUnspliced <- "unspliced"
    print("looking at unspliced data")
  }
  
  
  
  location <- rowRanges(SCE)
  is_mito <- any(seqnames(location) == "MT")
  
  #df <- perCellQCMetrics(SCE, subsets=list(Mito = is_mito))
  
  
  #generate QC data frame
  df <- perCellQCMetrics(SCE)
  df$mitoPerc<-SCE$percent.mt
  #df
  
  #add QC information to metadata
  SCE <- addPerCellQC(SCE)
  #colnames(colData(SCE))
  
  #calculate thresholds. CAVE! Parts here are hard coded. 
  #Some thresholds are too permissive for nuclei, so I set 
  #them manually. Would be good to pull this part out of 
  #the script and let the user specify cutoffs in the future. 
  
  #calculate if cells are outliers for low UMI count
  qc_lib2 <- isOutlier(df$sum, log=TRUE, type="lower")
  
  lowUMI <- as.numeric(attr(qc_lib2, "thresholds")[2])
  
  if(lowUMI == Inf){
    lowUMI <-0
  }
  
  if(min(df$sum) >= lowUMI){
    print("Filter for low UMI counts are too permissive. 
          Apply manual set threshold. Default is 100")
    qc_lib2 <- df$sum < 100 #changed from 20
    lowUMI <- 100 #changed from 20
  }else{
    print("Calculated thresholds seem to filter for 
          low UMI counts")
  }
  
  #calculate if cells are outliers for high UMI count 
  #(could be doublets or something similar)
  
  qc_lib2_higher<-isOutlier(df$sum, 
                            log=TRUE, 
                            type="higher") 
  #for the human cell atlas, the above is not filtering out much
  
  higherUMI <- as.numeric(attr(qc_lib2_higher, "thresholds")[2])
  
  if(max(df$sum)<= higherUMI){
    print("Filter for high UMI counts too permissive. 
          Apply manual set threshold.")
    qc_lib2_higher <- df$sum > 30000
    higherUMI <- 30000
  }else{
    print("Calculated thresholds seem to filter for high UMI counts")
  }
  
  
  
  qc_nexprs2_higher <- isOutlier(df$detected, log=TRUE, type="higher")
  
  higherFeatures <- as.numeric(attr(qc_nexprs2_higher, "thresholds")[2])
  
  if(max(df$detected)<= higherFeatures){
    print("Filter for high gene counts too permissive. 
          Apply manual set threshold.")
    qc_lib2_higher<-(df$detected)> 8000
    higherFeatures<-8000
  }else{
    print("Calculated thresholds seem to filter for high gene counts")
  }
  
  qc_nexprs2 <- isOutlier(df$detected, log=TRUE, type="lower")
  
  
  lowerFeatures <- as.numeric(attr(qc_nexprs2, "thresholds")[1])
  
  
  if(min(df$detected) >= lowerFeatures | lowerFeatures < 10){
    print("Filter for high gene counts too permissive. 
          Apply manual set threshold.")
    qc_nexprs2<-(df$detected)<10
    lowerFeatures<-10
  }else{
    print("Calculated thresholds seem to filter for high gene counts")
  }
  
  
  
  qc_mito2 <- isOutlier(df$mitoPerc, type="higher")
  attr(qc_mito2, "thresholds")
  
  percMitoThres<- as.numeric(attr(qc_mito2, "thresholds")[2])
  #no mitochondrial genes in unspliced matrices...
  
  
  discard2 <- qc_lib2 | qc_nexprs2 | qc_mito2 | qc_lib2_higher | qc_nexprs2_higher
  
  
  filterDF<-DataFrame(LibSize_low=sum(qc_lib2),
                      LibSize_low_threshold = lowUMI,  
                      LibSize_higher=sum(qc_lib2_higher), 
                      LibSize_higher_threshold = higherUMI, 
                      NExprs_low=sum(qc_nexprs2), 
                      NExprs_low_threshold = lowerFeatures,
                      NExprs_higher=sum(qc_nexprs2_higher), 
                      NExprs_higher_threshold = higherFeatures,
                      MitoPerc=sum(qc_mito2), 
                      MitoPeec_threshold = percMitoThres,
                      Total=sum(discard2))
  
  write.csv(filterDF, file.path(getwd(),
                                "splice_control_out",
                                "qc_data",
                                "scater",
                                paste(
                                  splicedUnspliced,
                                  "cell_qc_threshold_table.csv",
                                  sep = "_"
                                )))
  
  
  
  
  datafr<-data.frame(low_lib_size = qc_lib2, 
                     large_lib_size = qc_lib2_higher, 
                     low_n_features = qc_nexprs2, 
                     high_n_features = qc_nexprs2_higher,
                     high_subsets_mito_percent = qc_mito2, 
                     discard = discard2)
  
  
  
  colData(SCE) <- cbind(colData(SCE), datafr)
  SCE$ProcessNumber <- factor(SCE$process_number)
  SCE$X10XBatch <- factor(SCE$X10XBatch)
  
  
  SCE$ScaterQC_failed <- datafr$discard
  
  tag <- paste("10X batch ", splicedUnspliced, sep= "")
  
  
  UMIplot <- plotColData(SCE, 
                         x="X10XBatch", 
                         y="sum", 
                         colour_by="ScaterQC_failed") +  
    scale_y_log10() + 
    ggtitle("Total UMI count per nucleus") + 
    xlab(tag)+
    guides(fill=guide_legend(title= "Failed Scater QC"))
  
  
  
  FeaturePlotQC <- plotColData(SCE, 
                               x="X10XBatch", 
                               y="detected", 
                               colour_by="ScaterQC_failed") +
    scale_y_log10() + 
    ggtitle("Detected total features per nucleus") + 
    xlab(tag)+
    guides(fill=guide_legend(title= "Failed Scater QC"))
  
  
  
  MitoPlot<-plotColData(SCE, 
                        x="X10XBatch", 
                        y="percent.mt", 
                        colour_by="ScaterQC_failed") + 
    ggtitle("Percent mitochondrial genes per cell") + 
    xlab(tag) +
    guides(fill=guide_legend(title= "Failed Scater QC"))
  
  
  myQCPlotList<- list()
  myQCPlotList[[1]]<-UMIplot
  myQCPlotList[[2]]<-FeaturePlotQC
  myQCPlotList[[3]]<-MitoPlot
  
  
  
  mp <-  multiplot(plotlist = myQCPlotList, cols = 1)
  
  
  pdf(file.path(getwd(),
                "splice_control_out",
                "qc_plots",
                "scater_plots",
                paste(splicedUnspliced,
                      "_Scater_FilterPlot.pdf",
                      sep="")), width = 7, height = 6)
  print(mp)
  dev.off()
  
  
  scatter1 <- plotColData(SCE, 
                          x = "sum", 
                          y = "percent.mt", 
                          colour_by = "ScaterQC_failed")
  
  pdf(file.path(getwd(),
                "splice_control_out",
                "qc_plots",
                "scater_plots",
                paste(splicedUnspliced,
                      "_Scater_scatterP.pdf",
                      sep="")), width = 5, height = 5)
  print(scatter1)
  dev.off()
  
  
  
  
  #discard low quality nuclei
  SCE <- SCE[,!datafr$discard]
  
  dir.create(file.path(getwd(),
                       "splice_control_out",
                       "datasets",
                       "02_scater_cell_filtered",
                       splicedUnspliced))
  
  saveRDS(SCE, file.path(getwd(),
                         "splice_control_out",
                         "datasets",
                         "02_scater_cell_filtered",
                         splicedUnspliced,
                         paste("SCE", 
                               splicedUnspliced,
                               "gene_cell_filtered.RDS",
                               sep = "_")))
  
  
}


######## PART 3
######## Add spliced and unspliced counts for nuclei that passed
######## QC for both spliced and unspliced counts


spliced_sce <- readRDS(file.path(getwd(),
                                 "splice_control_out",
                                 "datasets",
                                 "02_scater_cell_filtered",
                                 "spliced",
                                 "SCE_spliced_gene_cell_filtered.RDS"))

unspliced_sce <- readRDS(file.path(getwd(),
                                   "splice_control_out",
                                   "datasets",
                                   "02_scater_cell_filtered",
                                   "unspliced",
                                   "SCE_unspliced_gene_cell_filtered.RDS"))


#Scater was used independently on spliced and unspliced matrices. 
#This means different cells may have been removed based on each input
#matrix. Nuclei that failed the Scater QC for one mode (spliced or 
#unspliced), will now also be removed for the other mode. 

summary(rownames(spliced_sce) == rownames(unspliced_sce))
#head(colnames(spliced))
#head(colnames(unspliced))
splicedInUnspliced <- colnames(spliced_sce) %in% 
  colnames(unspliced_sce)

summary(splicedInUnspliced)

spliced_sce<-spliced_sce[,splicedInUnspliced]

unsplicedInSpliced<-colnames(unspliced_sce) %in% 
  colnames(spliced_sce)

summary(unsplicedInSpliced)

unspliced_sce<-unspliced_sce[, unsplicedInSpliced]

#check if colnames are the same now

summary(colnames(spliced_sce) == colnames(unspliced_sce))


metadata <- colData(spliced_sce)

splicedM <- counts(spliced_sce)
unsplicedM <- counts(unspliced_sce)


combMatrix <- splicedM + unsplicedM

SCE_comb <- SingleCellExperiment(assays = list(counts = combMatrix), 
                                 colData = metadata)


saveRDS(SCE_comb, file.path(getwd(), 
                            "splice_control_out",
                            "datasets",
                            "03_combined_matrices",
                            "combined_QC_filtered.RDS"))

print("Done")