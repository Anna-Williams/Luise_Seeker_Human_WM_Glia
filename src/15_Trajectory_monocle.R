library(dplyr)
library(monocle)
library(Seurat)

#Extract data, phenotype data, and feature data from the SeuratObject

split_sr_by_t <- c("CSC", "CB", "BA4")


run_monocle_split <- function(seur_obj, 
                        split_by = "NULL",
                        split_sr_by,
                        mod_form_string, 
                        save_dir = getwd(),
                        cores_use = 8){
   
  if(split_by != "NULL"){
    
    for(i in 1: length(split_sr_by)){
     
    Idents(seur_obj) <- split_by
    seur_spli_obj <- subset(seur_obj, ident = split_sr_by[i])
    mon_data <- as(as.matrix(seur_spli_obj@assays$RNA@data), 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = seur_spli_obj@meta.data)
    f_data <- data.frame(gene_short_name = row.names(mon_data), 
                         row.names = row.names(mon_data))
    fd <- new('AnnotatedDataFrame', data = f_data)
    
    #Construct monocle cds
    ol_cds <- newCellDataSet(mon_data,
                             phenoData = pd,
                             featureData = fd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())
    ol_cds <- estimateSizeFactors(ol_cds)
    ol_cds <- estimateDispersions(ol_cds)
    clustering_oligos <- differentialGeneTest(ol_cds,
                                              fullModelFormulaStr = mod_form_string,
                                              cores = cores_use, verbose = T)
    
    # select genes that are significant at an FDR <10%
    genes_ol_cds <- subset(clustering_oligos, qval < 0.1) 
    
    
    # my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
    ol_cds <- setOrderingFilter(ol_cds, ordering_genes = genes_ol_cds)
    ol_cds <- reduceDimension(ol_cds, method = 'DDRTree', cores=cores_use, verbose = TRUE)
    ol_cds <- orderCells(ol_cds)
    my_pseudotime_de <- differentialGeneTest(ol_cds,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                             cores = cores_use)
    
    saveRDS(ol_cds, paste(save_dir, "/monocle_obj_", split_sr_by[i], ".RDS", sep = ""))
    saveRDS(my_pseudotime_de, paste(save_dir, "/pseudotime_de_", split_sr_by[i], ".RDS", sep = ""))
    }
  }
  print("Done")
}

dir.create("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle")
run_monocle_split(nad_ol, split_by = "Tissue",
                  split_sr_by = split_sr_by_t,
                  mod_form_string = '~ol_clusters_named', 
                  save_dir = "/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle",
                  cores_use = 8)

# For me this created pseudotimes with the wrong root

ol_cds_csc <-readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/wrong_start/monocle_obj_CSC.RDS")

ol_cds_csc <- orderCells(ol_cds_csc, root_state = 2)

saveRDS(ol_cds_csc,"/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/monocle_obj_CSC.RDS" )

ol_cds_cb <-readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/wrong_start/monocle_obj_CB.RDS")

ol_cds_cb <- orderCells(ol_cds_cb, root_state = 4)

saveRDS(ol_cds_cb,"/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/monocle_obj_CB.RDS" )

ol_cds_ba4 <-readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/wrong_start/monocle_obj_BA4.RDS")

ol_cds_ba4 <- orderCells(ol_cds_ba4, root_state = 3)
saveRDS(ol_cds_ba4,"/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/monocle_obj_BA4.RDS" )


ol_cds<- ol_cds_csc

ol_cds<- ol_cds_ba4



#whole data

nad_ol <- readRDS("/exports/eddie/scratch/lseeker/srt_oligos_and_opcs_LS.RDS")

ol_data <- as(as.matrix(nad_ol@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = nad_ol@meta.data)

f_data <- data.frame(gene_short_name = row.names(ol_data), row.names = row.names(ol_data))
fd <- new('AnnotatedDataFrame', data = f_data)

#Construct monocle cds
ol_cds <- newCellDataSet(ol_data,
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())

pData(ol_cds)


ol_cds <- estimateSizeFactors(ol_cds)
ol_cds <- estimateDispersions(ol_cds)

clustering_oligos <- differentialGeneTest(ol_cds,
                                          fullModelFormulaStr = '~ol_clusters_named',
                                          cores = 1, verbose = T)

# select genes that are significant at an FDR <10%
genes_ol_cds <- subset(clustering_oligos, qval < 0.1) 


# my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
ol_cds <- setOrderingFilter(ol_cds, ordering_genes = genes_ol_cds)
ol_cds <- reduceDimension(ol_cds, method = 'DDRTree', cores=1, verbose = TRUE)
ol_cds <- orderCells(ol_cds)

saveRDS(ol_cds, "/exports/eddie/scratch/lseeker/ol_cds.RDS")

ol_cds <- readRDS("/exports/eddie/scratch/lseeker/ol_cds.RDS")

my_pseudotime_de <- differentialGeneTest(ol_cds,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)")

saveRDS(ol_cds, "/exports/eddie/scratch/lseeker/monocle_all_cds.RDS")
saveRDS(my_pseudotime_de, "/exports/eddie/scratch/lseeker/monocle_all_pseudotime.RDS")




plot_cell_trajectory(ol_cds, color_by = "State") + facet_wrap(~State) + 
  scale_colour_manual(values=mycoloursP[20:40])

plot_cell_trajectory(ol_cds, color_by = "State")  + 
  scale_colour_manual(values=mycoloursP[20:40])


plot_cell_trajectory(ol_cds, color_by = "ol_clusters_named") + 
  facet_wrap(~ol_clusters_named) +
  scale_colour_manual(values=mycoloursP[6:40])

#plot_cell_trajectory(ol_cds, color_by = "Tissue") + facet_wrap(~Tissue)
plot_cell_trajectory(ol_cds, color_by = "gender", cell_size = 0.2)+
  scale_colour_manual(values=mycoloursP[10:40])

plot_cell_trajectory(ol_cds, color_by = "Pseudotime")
plot_cell_trajectory(ol_cds, color_by = "ol_clusters_named")+
  scale_colour_manual(values=mycoloursP[6:40])


head(pData(ol_cds))
#the pseudotime and its state is now in the pData

my_pseudotime_de <- differentialGeneTest(ol_cds,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 18)
top_genes <- my_pseudotime_de %>% arrange(qval) %>% head(10)

top_genes
# select the 10 highest genes and use them for plotting
plot_genes_in_pseudotime(ol_cds[rownames(top_genes)])


plot_genes_in_pseudotime(ol_cds[c("DSCAM","TNR","PTPRZ1","LHFPL3", "PCDH15")])
plot_genes_in_pseudotime(ol_cds[c("SOX6","VCAN","MMP16","FGF14", "CA10")])
plot_genes_in_pseudotime(ol_cds[c("DPYD","COL11A1","TNR","BRINP3", "BCAN")])

plot_genes_in_pseudotime(ol_cds[c("PDGFRA","PCDH15", "PLP1","MBP")], 
                         color_by = "ol_clusters_named") +  
  scale_colour_manual(values=mycoloursP[6:40])


# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$gene_short_name
gene_to_cluster <- c(gene_to_cluster, "MAG", "MOG", "PLP1", "SPARC", "RBFOX1", 
                     "OPALIN", "KLK6", "SPARCL1", "SLC22A3", "PAX3", "NELL1", 
                     "AFF3", "LGALS1", "FMN1", "HHIP", "TUBB2A")
my_pseudotime_cluster <- plot_pseudotime_heatmap(ol_cds[gene_to_cluster,],
                                                 num_clusters = 3,
                                                 cores = 8,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_cluster <- cutree(my_pseudotime_cluster$tree_row, 5)
my_cluster

my_pseudotime_de[names(my_cluster[my_cluster == 1]),"gene_short_name"]



#  I am surprised about how the pseudotime trajectories look like. I expected 
# one cluster on one side of a fork and the other onthe other side. So if this
# trajectory analysis is not picking up on cluster differences, what is driving
# the fork then? 

# A differential gene expression of the two end points may be helpful

ol_cds$subs_pt <- ifelse(ol_cds$State == 1 & ol_cds$Pseudotime > 45, "state1_end_point", "FALSE")
ol_cds$subs_pt <- ifelse(ol_cds$State == 3 & ol_cds$Pseudotime > 28, "state3_end_point", ol_cds$subs_pt)
plot_cell_trajectory(ol_cds, color_by = "subs_pt")+
  scale_colour_manual(values=mycoloursP[36:40])

# extract barcode , pseudotime, state and and subs_pt and merge with seur metadata

mon_df <- data.frame(Barcode = ol_cds$Barcode, 
                     Pseudotime = ol_cds$Pseudotime, 
                     monocle_state = ol_cds$State, 
                     monocle_end_points = ol_cds$subs_pt)

# Do the above for all the different tissues, create master mon_df, bring barcodes
# in same order as in seurat obj, add columns to seurat object metadata. 
# Plot Dimplots that show states for all nuclei and separated by tissue. 
# Do differential gene expression analysis of end states. 


ol_cds_cb<- readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/monocle_obj_CB.RDS")
ol_cds_cb$subs_pt <- ifelse(ol_cds_cb$State == 1 & ol_cds_cb$Pseudotime > 40, "state1_end_point", "FALSE")
ol_cds_cb$subs_pt <- ifelse(ol_cds_cb$State == 5 & ol_cds_cb$Pseudotime > 32, "state5_end_point", ol_cds_cb$subs_pt)
plot_cell_trajectory(ol_cds_cb, color_by = "subs_pt")+
  scale_colour_manual(values=mycoloursP[36:40])

# extract barcode , pseudotime, state and and subs_pt and merge with seur metadata

mon_df_cb <- data.frame(Barcode = ol_cds_cb$Barcode, 
                     Pseudotime = ol_cds_cb$Pseudotime, 
                     monocle_state = ol_cds_cb$State, 
                     monocle_end_points = ol_cds_cb$subs_pt)

#And for BA4

ol_cds_ba4<- readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/2021_oligos_out/monocle/monocle_obj_BA4.RDS")
ol_cds_ba4$subs_pt <- ifelse(ol_cds_ba4$State == 1 & ol_cds_ba4$Pseudotime > 46, "state1_end_point", "FALSE")
plot_cell_trajectory(ol_cds_ba4, color_by = "subs_pt")+
  scale_colour_manual(values=mycoloursP[36:40])

# extract barcode , pseudotime, state and and subs_pt and merge with seur metadata

mon_df_ba4 <- data.frame(Barcode = ol_cds_ba4$Barcode, 
                        Pseudotime = ol_cds_ba4$Pseudotime, 
                        monocle_state = ol_cds_ba4$State, 
                        monocle_end_points = ol_cds_ba4$subs_pt)


merged_mon_df <- rbind(mon_df, mon_df_cb, mon_df_ba4)


# sonrt barcodes so that they are the same in this dataframe and in the seurat
# object, so that data can be added to metadata of seurat object

nad_ol <- readRDS("/Users/lseeker/Documents/Work/HumanCellAtlas/srt_oligos_Nadine/srt_oligos_and_opcs_LS.RDS")

summary(nad_ol@meta.data$Barcode == merged_mon_df$Barcode)

nrow(nad_ol@meta.data) == nrow(merged_mon_df)

merged_mon_df <- merged_mon_df[match(nad_ol@meta.data$Barcode,
                                     merged_mon_df$Barcode),]

summary(nad_ol@meta.data$Barcode == merged_mon_df$Barcode)

nad_ol$mon_pseudotime <- merged_mon_df$Pseudotime
nad_ol$monocle_state <- merged_mon_df$monocle_state
nad_ol$monocle_end_points <- merged_mon_df$monocle_end_points

DimPlot(nad_ol, group.by = "monocle_end_points", split.by = "Tissue")


# looking at differential gene expression of end points seems to make more 
# sense with in the tissue than across tissues. 

Idents(nad_ol) <- "Tissue"
csc <- subset(nad_ol, ident = "CSC")

Idents(csc) <- "monocle_end_points"

end1_markers <- FindMarkers(object = csc, ident.1 = "state1_end_point", 
                            ident.2 = "state3_end_point", min.pct = 0.25)

print(x = head(x = end1_markers, n = 5))

end3_markers <- FindMarkers(object = csc, ident.1 = "state3_end_point", 
                            ident.2 = "state1_end_point", min.pct = 0.25, 
                            only.pos = TRUE)

print(x = head(x = end3_markers, n = 5))


fil_end1_mark <- subset(end1_markers, end1_markers$avg_log2FC > 1.5)
fil_end3_mark <- subset(end3_markers, end3_markers$avg_log2FC > 1.5)
