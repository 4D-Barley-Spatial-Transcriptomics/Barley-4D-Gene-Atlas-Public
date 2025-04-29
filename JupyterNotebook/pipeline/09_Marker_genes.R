library(cowplot)
library(rainbow)
library(patchwork)
library(dplyr)
library(Seurat)
library(gridExtra)

find_marker_genes <- function(data, output_folder, name, res, timepoints) {
  print("Marker genes")
  
  output_folder <- paste0(output_folder, "09_Marker_genes/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  title <- paste0(name, " res=", res)
  all_markers <- FindAllMarkers(data, min.pct = 0.05, logfc.threshold = 0.25, 
                                only.pos = TRUE, test.use = "wilcox", res = res)
  
  filtered_markers <- all_markers[all_markers$p_val_adj < 0.05, ] 
  
  write.csv(filtered_markers, file = paste0(output_folder, title, " all_marker_genes.csv"))
  
  
}