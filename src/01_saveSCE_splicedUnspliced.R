# Script for converting loom files to single cell experiments
# for all spliced and all unspliced matrices respectively
# LASeeker
# 20200418

# LOAD LIBRARIES
library(lattice)
library(matrix)
library(KernSmooth)
library(MASS)
library(nlme)
library(cluster)
library(survival)
library(devtools)
library(Seurat)
library(sctransform)
library(hdf5r)
library(pcaMethods)
library(loomR)
library(velocyto.R)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)



e <- 40000 * 1024^2
options(future.globals.maxSize = e)


####### MANUAL INPUT REQUIRED HERE #####################
####### Set working directory############################

setwd("/exports/eddie/scratch/lseeker")


########### Find input file locations####################
# loom files should be saved in a directory called "loom" and
# the metadata in a directory called "metadata" within
# the working directory

loom_names <- list.files(path = "loom", pattern = ".loom")
met_data <- list.files(path = "metadata", pattern = ".csv")

############ Set variables for creating Seurat object########
min_cells <- 0 # you can set QC thresholds here if you like
# (for example enter a 3 here and remove all genes that were
# expressed in less cells)

min_features <- 0 # if you want to filter some more, you can
# remove here genes that are not frequently expressed
# (for example all genes below 200 copies detected across all
# cells)

# I decided to so all filtering later.


############# create output folders#########################

dir.create(file.path(getwd(),
  "splice_control_out",
  "datasets",
  "raw_matrices",
  "unspliced",
  sep = ""
), recursive = TRUE)

dir.create(file.path(getwd(),
  "splice_control_out",
  "datasets",
  "raw_matrices",
  "spliced",
  sep = ""
))

dir.create(file.path(
  "splice_control_out",
  "dupl_genes"
), recursive = TRUE)



##############################################################
#################### START ###################################


## read in raw metadata
raw_met_data <- read.csv(paste(getwd(), "/metadata/",
  met_data,
  sep = ""
))


for (i in 1:length(loom_names)) {
  process_number <- sapply(strsplit(loom_names[i], "__"), `[`, 2)
  process_number <- sapply(strsplit(process_number, "_"), `[`, 2)

  # generate a sample name that will be used for saving data
  # and naming Seurat/SCE objects
  sample_name <- sapply(strsplit(loom_names[i], "__"), `[`, 2)
  sample_name <- sub(
    pattern = "(.*)\\..*$",
    replacement = "\\1", basename(sample_name)
  )


  # read in .loom file
  loom_matrices <- read.loom.matrices(file.path(
    getwd(),
    "loom",
    loom_names[i]
  ))

  # extract a matrix with counts for spliced mRNAs
  loom_spliced <- loom_matrices$spliced
  # extract a matrix with counts for unspliced mRNAs
  loom_unspliced <- loom_matrices$unspliced



  ############## Save duplicated gene names######################

  # depending on the used reference genome, it is common that some
  # gene names are duplicated. It is not completely clear to us why
  # this happens. We tested several ways of dealing with this
  # issue and decided to remember duplicated genes so that they
  # can be checked further downstream. If they happen to be genes
  # of interest, differen solutions need to be explored. The combined
  # matrix is only created here for the sake of saving gene counts
  # for the duplicated genes. make.unique maked names unique by
  # appending an ".1" to the name of the second item.

  comb_matrix <- loom_spliced + loom_unspliced

  rownames(comb_matrix) <- make.unique(rownames(comb_matrix))

  d <- rownames(loom_spliced)[duplicated(rownames(loom_spliced))] 
  # gene names in all matrices within a .loom file are the same


  d <- c(d, paste(d, ".1", sep = ""))
  s <- as.data.frame(rowSums(comb_matrix)[which(rownames(comb_matrix)
  %in% d)])
  
  s$genes <- rownames(s)
  rownames(s) <- NULL
  names(s)[1] <- paste("count",
                     sample_name,
                     sep = "_")

  
  if(i == 1){
    t <- s
  } else {
    t <- merge(t, s, by = "genes", all = TRUE)
  }




  ###################################################################
  #### here we continue with the spliced and unspliced matrices #######

  rownames(loom_unspliced) <- make.unique(rownames(loom_unspliced))
  rownames(loom_spliced) <- make.unique(rownames(loom_spliced))


  # prepare to link count matrix to metadata:

  df <- data.frame(
    Barcode = colnames(loom_spliced),
    process_number = as.integer(paste(process_number))
  )

  # Link barcodes to metadata
  metadata <- merge(df, raw_met_data, by = "process_number", all.x = TRUE)

  # rownames must be the barcodes to allow merging

  rownames(metadata) <- metadata$Barcode


  # Convert the combined matrix to a Seurat object
  y_unspliced <- CreateSeuratObject(
    counts = loom_unspliced,
    min.cells = 0,
    min.features = 0,
    project = sample_name,
    assay = "RNA",
    names.field = 1,
    names.delim = ".",
    meta.data = metadata
  )

  y_unspliced[["percent.mt"]] <-
    PercentageFeatureSet(y_unspliced, pattern = "^MT-")


  # same for spliced
  y_spliced <- CreateSeuratObject(
    counts = loom_spliced,
    min.cells = 0,
    min.features = 0,
    project = sample_name,
    assay = "RNA", 
    names.field = 1,
    names.delim = ".",
    meta.data = metadata
  )

  y_spliced[["percent.mt"]] <-
    PercentageFeatureSet(y_spliced, pattern = "^MT-")

  # Save the first one and then apend the others to it
  if (i == 1) {
    x_unspliced <- y_unspliced
    x_spliced <- y_spliced
  } else {
    x_unspliced <- merge(
      x = x_unspliced,
      y = y_unspliced, project = "HCA_raw_unspliced"
    )
    x_spliced <- merge(
      x = x_spliced,
      y = y_spliced, project = "HCA_raw_spliced"
    )
  }
  
  write.csv(t, file.path(
    getwd(),
    "splice_control_out",
    "dupl_genes",
    paste("dupl_gene_counts_",
          ".csv",
          sep = ""
    )
  ))
}



sample_count <- i




sing_cell_unspliced <- as.SingleCellExperiment(x_unspliced)
sing_cell_spliced <- as.SingleCellExperiment(x_spliced)

saveRDS(sing_cell_unspliced, file.path(
  getwd(),
  "splice_control_out",
  "datasets",
  "raw_matrices",
  "unspliced",
  paste("Combined_sce_",
    sample_count,
    "_unspliced_samples.RDS",
    sep = ""
  )
))

saveRDS(sing_cell_spliced, file.path(
  getwd(),
  "splice_control_out",
  "datasets",
  "raw_matrices",
  "spliced",
  paste("Combined_sce_",
    sample_count,
    "_spliced_samples.RDS",
    sep = ""
  )
))




# save session info to file
sink(file.path(
  getwd(),
  "splice_control_out",
  paste(Sys.Date(),
    "_session_info.txt",
    sep = ""
  )
))
sessionInfo()
sink()


# End of Script
