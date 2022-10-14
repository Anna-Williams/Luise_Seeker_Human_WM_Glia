# Script for generating quality control plots on spliced
# and unspliced matrices separately and for saving a list
# of gene names that were duplicated (in the regerence
# genome?).
# 20200413
# Luise A. Seeker
# Human Cell Atlas project: Human oligodendrocyte
# heterogeneity with CNS region, age and sex


# Introduction
# This script is only generating QC plots 

# data has previously been aligned using cellranger and
# run through velocyto to retrieve spliced and unspliced
# feature count matrices. Then loom matrices were saved as 
# single cell experiment object using the script 
# "saveSCE_splicedUnspliced.R".
# The present script is not performing any 
# Quality control steps but different filters can be set to 
# observe how this affects plots and summary stats. 
# To actually filter data, use the script 
# splicedUnspliced_scaterQC.R.


####### Load Libraries:

library(Seurat)
library(ggplot2)
library(scater)
library(ggsci)
library(gridExtra)



#### set working directory###

setwd("/exports/eddie/scratch/lseeker")



##### Set paths to input files

###################################################################
# Set Quality control thresholds or leave them as they have
# been used for the HCA

# Set QC paramters that will be used for filtering the data:
min_umi <- 1 # 1 less stict than scater becasue scater filters
# were done on combined spliced and unspliced counts

max_umi <- 15000 # 15000 previously determined using Scater
min_feature <- 1 # 1 less stict than scater becasue scater
# filters were done on combined spliced and unspliced counts

max_feature <- 8000 # 8000 previously determined using scater
max_mt_pc <- 11

######### Set variables for creating Seurat object##############
min_cells <- 0 
# you can set QC thresholds here if you like
# for example enter a 3 here and remove
# all genes that were expressed in less cells

min_features <- 0 # if you want to filter some more,
# you can remove here genes that are not
# frequently expressed (for example all genes below
# 200 copies detected across all cells)


########## create QC plot directories ######################

dir.create(file.path(
  "splice_control_out", "qc_plots"))

dir.create(file.path("splice_control_out", "qc_data"))

###### Pick colour pallet ####################################


colour_palette <- c(
  pal_npg("nrc", alpha = 0.7)(10),
  pal_tron("legacy", alpha = 0.7)(7),
  pal_lancet("lanonc", alpha = 0.7)(9),
  pal_simpsons(
    palette = c("springfield"),
    alpha = 0.7
  )(16)
)

################################################################
################################################################
### specify functions

determine_plot_tag <- function(filter) {
  if (filter == "YES") {
    plot_tag <<- "filtered"
  } else if (filter == "NO") {
    plot_tag <<- "unfiltered"
  } else {
    stop('invalid value for filter')
  }
}

generate_individual_cq_plots <- function(spli_unspli, filter) {
  
  if (spli_unspli == "spliced") {
    # determine which data is processed and set plot colour and 
    # x-axis max.
    
    # split seurat object into list of single samples
    seur_list <<- SplitObject(seur_intron,
                             split.by = "process_number")
    plot_colour <- colour_palette[5]
    max_x <- 4000 # this is a hard-coded limit for x axes of plots
    #would be nice to find a more variable measure
    print("spliced data will be investigated")
    
    # Define the variables to be for the unspliced
  } else if (spli_unspli == "unspliced") {
    seur_list <<- SplitObject(seur_exon,
                             split.by = "process_number"
    )
    plot_colour <- colour_palette[6]
    max_x <- 6000# this is a hard-coded limit for x axes of plots
    #would be nice to find a more variable measure
    print("unspliced data will be investigated")
  }
  
  
  
  for (i in 1:length(seur_list)) {
    # look at individual sample within a list of seurat objects,
    # generate a sample_id that can be used to name files.
    sample_seurat <<- seur_list[[i]]
    sample_id <<- paste(sample_seurat$uniq_id[1],
                       "_",
                       sample_seurat$caseNO[1],
                       "_",
                       sample_seurat$Tissue[1],
                       sep = ""
    )
    sample_id <<- gsub("/", "_", sample_id)
    
    if (spli_unspli == "spliced") {
      if (filter == "YES") {
        sample_seurat <<- subset(
          x = sample_seurat,
          subset = nCount_RNA > min_umi &
            nCount_RNA < max_umi &
            nFeature_RNA > min_feature &
            nFeature_RNA < max_feature &
            percent.mt > -Inf &
            percent.mt < max_mt_pc
        )
      } else if (filter == "NO") {
      }
      }else if (spli_unspli == "unspliced") {
        if (filter == "YES") {
          sample_seurat <<- subset(
            x = sample_seurat,
            subset = nCount_RNA > min_umi &
              nCount_RNA < max_umi &
              nFeature_RNA > min_feature &
              nFeature_RNA < max_feature
          )
        } 
      }
        
        
        
        p1 <- ggplot(sample_seurat@meta.data, aes(x = caseNO, 
                                                 y = nCount_RNA)) + 
          geom_violin(fill= colour_palette[6])+ 
          geom_jitter(shape=16, position=position_jitter(0.2)) + 
          theme_minimal()+
          labs(title= paste(sample_seurat@meta.data$process_number[1],
                            sample_seurat@meta.data$caseNO[1],
                            sample_seurat@meta.data$Tissue[1],
                            sep = "_"),
               x="Individual ID", y = "RNA molecule count")
          
        
        p2 <- ggplot(sample_seurat@meta.data, aes(x = caseNO, 
                                                    y = nFeature_RNA)) + 
          geom_violin(fill= colour_palette[7])+ 
          geom_jitter(shape=16, position=position_jitter(0.2)) + 
          theme_minimal()+
          labs(title= paste(sample_seurat@meta.data$process_number[1],
                            sample_seurat@meta.data$caseNO[1],
                            sample_seurat@meta.data$Tissue[1],
                            sep = "_"),
               x="Individual ID", y = "Unique gene count")
        
        p3 <- ggplot(sample_seurat@meta.data, aes(x = caseNO, 
                                                  y = percent.mt)) + 
          geom_violin(fill= colour_palette[8])+ 
          geom_jitter(shape=16, position=position_jitter(0.2)) + 
          theme_minimal()+
          labs(title= paste(sample_seurat@meta.data$process_number[1],
                            sample_seurat@meta.data$caseNO[1],
                            sample_seurat@meta.data$Tissue[1],
                            sep = "_"),
               x="Individual ID", 
               y = "Mitochondrial genes (%)")
        
        
        v1 <- grid.arrange(p1, p2, p3, ncol = 3)
        
        # apparently there are no unspliced mitochondrial genes in the
        # human cell nucleus
        
        seur_plot_l[[i]] <<- v1
        
        
        v2 <- FeatureScatter(sample_seurat,
                             feature1 = "nCount_RNA",
                             feature2 = "nFeature_RNA"
        )
        
        seur_cor_plot_l[[i]] <<- v2
        
        
        # calculate some QC stats and save them to a data frame.
        # Append list of dataframes in consecutive iterations and bind
        # them together after loop. Same as for spliced genes but without
        # mitochondrial info.
        
        plot_data <- data.frame(
          sample_id = sample_id,
          processNumber =
            sample_seurat$process_number[1],
          RinValue = sample_seurat$RINvalue[1],
          matrix_mod = spli_unspli,
          filtered = plot_tag,
          average_umi =
            mean(sample_seurat$nCount_RNA),
          median_umi =
            median(sample_seurat$nCount_RNA),
          umi_98_percentile =
            quantile(sample_seurat$nCount_RNA, 0.98),
          averageGenes =
            mean(sample_seurat$nFeature_RNA),
          medianGenes =
            median(sample_seurat$nFeature_RNA), 
          average_pc_mt = 
            mean(sample_seurat$percent.mt),
          median_pc_mt = 
            median(sample_seurat$percent.mt)
        )
        
        
        rownames(plot_data) <- ""
        
        datalist[[i]] <<- plot_data
        
        
        # Prepare histograms
        
        umi_98_percentile <- quantile(sample_seurat$nCount_RNA, 0.98)
        
        cell_count <- ncol(sample_seurat)
        
        hist_data <- data.frame(numberOfUMIs = sample_seurat$nCount_RNA)
        
        intercept_mean <- mean(sample_seurat$nCount_RNA)
        intercept_median <- median(sample_seurat$nCount_RNA)
        intercept98 <- plot_data$umi_98_percentile
        
        
        hist_plot <- ggplot(hist_data, aes(x = numberOfUMIs)) +
          geom_histogram(
            col = "black", fill = plot_colour,
            binwidth = 100
          ) +
          theme_classic() +
          ggtitle(paste(sample_id, spli_unspli, "RIN",
                        sample_seurat$RINvalue[1], "Cells:",
                        cell_count,
                        sep = " "
          )) +
          xlab("UMI count") +
          ylab("Number of cells") +
          xlim(0, max_x) +
          geom_vline(aes_(
            xintercept = intercept_mean,
            colour = "royalblue", show.legend = T
          )) +
          geom_vline(aes_(
            xintercept = intercept_median,
            colour = "red", show.legend = T
          )) +
          geom_vline(aes_(
            xintercept = intercept98,
            colour = "cyan", show.legend = T
          )) +
          scale_color_identity(
            name = "Vertical lines inticate",
            breaks = c("cyan", "red", "royalblue"),
            labels = c(
              "98th UMI percentile",
              "Median UMI", "mean UMI"
            ),
            guide = "legend"
          )
        
        myplots[[i]] <<- hist_plot
    }
  
  #setting number of olumns for plots
  col_count <- floor(sqrt(length(myplots)))
  
  
  #save histrograms to file
  #I decided to always have four plots per row, but the
  #number of rows clearly varies with the number of plots
  #in the list. Therefore the plot height will be adjusted
  #accordingly.
  
  pdf(file.path(
    "splice_control_out",
    "qc_plots",
    paste(spli_unspli,
          "_",
          plot_tag,
          "_Histograms.pdf",
          sep = ""
      )
    ),
    width = 18, height = (length(myplots)/4)* 2.7
  )
  
  
  #mp <- multiplot(plotlist = myplots[1:4], cols = 4)
  mp <- do.call("grid.arrange", c(myplots, ncol=4))  
  
  print(mp)
  
  dev.off()
  
  
  
  #save seurat plots to file
  
  pdf(file.path(
    "splice_control_out",
    "qc_plots",
    paste(spli_unspli,
          "_",
          plot_tag,
          "_SeuratQC1.pdf",
          sep = ""
      )
    ),
    width = 18, height = (length(myplots)/3)* 3.5
  )
  
  
  seur_mp1 <- multiplot(plotlist = seur_plot_l, cols = 3)
  
  
  print(seur_mp1)
  dev.off()
  
  
  pdf(file.path(
    "splice_control_out",
    "qc_plots",
    paste(spli_unspli,
          "_",
          plot_tag,
          "_SeuratQC2.pdf",
          sep = ""
      )
    ),
    width = 18, height = (length(myplots)/4)* 2.7
  )
  
  seur_mp2 <- multiplot(plotlist = seur_cor_plot_l, cols = 4)
  print(seur_mp2)
  dev.off()
  
  
  #####################
}

generate_overview_cq_plots <- function(spli_unspli) {
  big_plot_data <- do.call(rbind, datalist)
  
  
  head(big_plot_data)
  
  big_plot_data$num_rin <-
    as.numeric(paste(big_plot_data$RinValue))
  
  big_plot_data$rin_group <-
    ifelse(big_plot_data$num_rin >= 6.5, "6.5 - 7.5",
           ifelse(big_plot_data$num_rin < 6.5 &
                    big_plot_data$num_rin >= 5, "5.0 - 6.5",
                  ifelse(big_plot_data$num_rin < 5 &
                           big_plot_data$num_rin >= 4, "4.0 - 5.0",
                         ifelse(big_plot_data$num_rin == "NA",
                                "Not yet measured", "2.1 - 2.8"
                         )
                  )
           )
    )
  
  
  
  big_plot_data$rin_group <- as.factor(big_plot_data$rin_group)
  
  big_plot_data$qc_check <-
           ifelse(big_plot_data$average_umi < 580,
                  "Low average UMI",
                  "Pass"
           )
  
  
  
  
  
  p1 <- ggplot(
    big_plot_data,
    aes(x = average_umi, y = median_umi, col = RIN_GGroup)
  ) +
    geom_point(aes(colour = rin_group), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("Median UMI") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  
  p2 <- ggplot(big_plot_data, aes(x = average_umi, y = umi_98_percentile)) +
    geom_point(aes(colour = rin_group), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("UMI 98th percentile") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p3 <- ggplot(big_plot_data, aes(x = median_umi, y = umi_98_percentile,
                                  col = rin_group)) +
    geom_point(aes(colour = rin_group), size = 3) +
    theme_classic() +
    xlab("Median UMI") +
    ylab("UMI 98th percentile") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  
  p4 <- ggplot(big_plot_data, aes(x = average_umi, y = median_umi)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("Median UMI") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  
  p5 <- ggplot(big_plot_data, aes(x = average_umi, y = umi_98_percentile,
                                  col = qc_check)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("UMI 98th percentile") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p6 <- ggplot(big_plot_data, aes(x = median_umi, y = umi_98_percentile,
                                  col = qc_check)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Median UMI") +
    ylab("UMI 98th percentile") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  
  pdf(file.path(
    "splice_control_out",
    "qc_plots",
    paste(spli_unspli, "_",
          plot_tag,
          "_overviewQCplots.pdf",
          sep = ""
    )
  ))
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  
  dev.off()
  
  return(big_plot_data)
}

generate_mito_plots <- function(spli_unspli) {
  keep_data <- big_plot_data
  
  p7 <- ggplot(keep_data, aes(x = average_umi, y = average_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("Average % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p8 <- ggplot(keep_data, aes(x = median_umi, y = average_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Median UMI") +
    ylab("Average % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p9 <- ggplot(keep_data, aes(x = umi_98_percentile, y = average_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("UMI 98th percentile") +
    ylab("Average % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p10 <- ggplot(keep_data, aes(x = average_umi, y = median_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Average UMI") +
    ylab("Median % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p11 <- ggplot(keep_data, aes(x = median_umi, y = median_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("Median UMI") +
    ylab("Median % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  p12 <- ggplot(keep_data, aes(x = umi_98_percentile, y = median_pc_mt)) +
    geom_point(aes(colour = qc_check), size = 3) +
    theme_classic() +
    xlab("UMI 98th percentile") +
    ylab("Median % mt") +
    ggtitle(spli_unspli) +
    scale_colour_manual(values = colour_palette, na.value = "black")
  
  pdf(file.path(
    "splice_control_out", "qc_plots",
    paste(spli_unspli, "_",
          plot_tag,
          "overviewQCplots_mitochondria.pdf",
          sep = ""
    )
  ))
  
  
  print(p7)
  print(p8)
  print(p9)
  print(p10)
  print(p11)
  print(p12)
  
  dev.off()
}


#### if single cell experiment objects generated in the previous
# step are still in environment, use. If not load
# separately for spliced and unspliced SCE


if (exists("sing_cell_intron") == FALSE) {
  intron_input <-
    list.files(
      path =
        file.path(
          getwd(),
          "splice_control_out",
          "datasets",
          "raw_matrices",
          "unspliced"
        ),
      pattern = ".RDS"
    )
  sing_cell_intron <-
    readRDS(file.path(
      getwd(),
      "splice_control_out",
      "datasets",
      "raw_matrices",
      "unspliced",
      intron_input
    ))
} else {
  print("Taking input from environment")
}


if (exists("sing_cell_exon") == FALSE) {
  exon_input <-
    list.files(
      path =
        file.path(
          getwd(),
          "splice_control_out",
          "datasets",
          "raw_matrices",
          "spliced"
        ),
      pattern = ".RDS"
    )
  sing_cell_exon <-
    readRDS(file.path(
      getwd(),
      "splice_control_out",
      "datasets",
      "raw_matrices",
      "spliced",
      exon_input
    ))
} else {
  print("Taking input from environment")
}


#convert to Seurat object
seur_intron <- as.Seurat(sing_cell_intron)
seur_exon <- as.Seurat(sing_cell_exon)



# initialise list for collecting sample QC metrics and
# qc plots

datalist <- list()
myplots <- list()
seur_plot_l <- list()
seur_cor_plot_l <- list()




# both spliced and unspliced data will be plotted once without
# QC filter being applied, once with

matrix_mod <- c("spliced", "unspliced")

qual_contr <- c("YES", "NO")

plot_tag <- "NA"



# For spliced and unspliced
for (k in 1 : length(matrix_mod)){
  # For apply QC and not apply QC
  for (n in 1:length(qual_contr)) {
    determine_plot_tag(qual_contr[n])
    generate_individual_cq_plots(matrix_mod[k], qual_contr[n])
  }

  big_plot_data <- generate_overview_cq_plots(matrix_mod[k])
    
  if (matrix_mod[k] == "spliced") {
    generate_mito_plots(matrix_mod[k])
  }

  # save data to csv
  write.csv(big_plot_data, file.path(
    getwd(),
    "splice_control_out",
    "qc_data",
    paste(matrix_mod[k],
      "_",
      plot_tag,
      "_summaryData.csv",
      sep = "")))
}
  



print("DONE")




#To do: fix that at the momebt big_plot_data is not
# saved for filtered data. Why??