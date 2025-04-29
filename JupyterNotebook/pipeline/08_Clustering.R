library(rainbow)
library(patchwork)
library(dplyr)
library(Seurat)

perform_clustering_and_marker_genes <- function(seu_obj, output_folder, resolution_list, umap_clusters, dims, 
                                                custom_palette_org, name, alpha, timepoints, levels_org) {
  
  image_clustering_folder_path <- paste0(output_folder, "08_Clustering/")
  if(!file.exists(image_clustering_folder_path))
    dir.create(image_clustering_folder_path, recursive = TRUE)
  
  image_folder_path_spatial_plots <- paste0(image_clustering_folder_path, "Spatial_plots/")
  if(!file.exists(image_folder_path_spatial_plots))
    dir.create(image_folder_path_spatial_plots, recursive = TRUE)
  
  image_folder_path_UMAPs <- paste0(image_clustering_folder_path, "UMAPs/")
  if(!file.exists(image_folder_path_UMAPs))
    dir.create(image_folder_path_UMAPs, recursive = TRUE)
  
  for(res in resolution_list){
    print(paste0("RESSSSSS: ", res))
    
    seu_obj <- FindClusters(seu_obj, verbose = TRUE, assay = "SCT", resolution = res)
    if(!umap_clusters) {
      seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims=1:dims)
    }
    
    levels <- levels(seu_obj@meta.data[["seurat_clusters"]])
    no_clusters <- max(as.numeric(levels)) + 1
    
    custom_palette <- custom_palette_org[1:no_clusters]
    mapped_colors <- setNames(custom_palette, as.character(levels))
    
    source("./08_01_Draw_clusters.R")
    
    draw_spatial_plot(seu_obj, image_folder_path_spatial_plots, name, FALSE, 
                      res, no_clusters, mapped_colors, alpha)
    draw_spatial_plot(seu_obj, image_folder_path_spatial_plots, name, TRUE, 
                      res, no_clusters, mapped_colors, alpha)
    
    draw_UMAP_clusters(seu_obj, image_folder_path_UMAPs, name, FALSE, res, mapped_colors)
    draw_UMAP_clusters(seu_obj, image_folder_path_UMAPs, name, TRUE, res, mapped_colors)
    
    # Marker genes
    source("./09_Marker_genes.R")
    find_marker_genes(seu_obj, output_folder, name, res, timepoints)
    
  }
  
  draw_UMAP_timepoints(seu_obj, image_folder_path_UMAPs, name, timepoints)
  
  custom_palette <- custom_palette_org[1:length(levels_org)]
  mapped_colors <- setNames(custom_palette, as.character(levels_org))
  draw_UMAP_timepoints_sections(seu_obj, image_folder_path_UMAPs, name, levels_org, mapped_colors, TRUE)
  draw_UMAP_timepoints_sections(seu_obj, image_folder_path_UMAPs, name, levels_org, mapped_colors, FALSE)
  
  return(seu_obj)
}
