library(cowplot)
library(rainbow)
library(patchwork)
library(dplyr)
library(Seurat)
library(gridExtra)
library(ggplot2)


draw_spatial_plot <- function(seu_obj, output_folder, name, label, res, no_clusters, mapped_colors, alpha) {
  title_integrated <- paste0(name," Merged Bright Field Image and Spatial Clusters - res=", res, if(label==TRUE) " labeled")
  print(title_integrated)
  
  plot_integrated <- SpatialDimPlot(seu_obj, alpha = c(alpha, alpha), label = label, label.box = TRUE,  
                                    label.size = 2, repel = TRUE, stroke = 0.0, cols = mapped_colors, crop = TRUE)
  
  spatial_plot_list <- list()
  
  for (i in 1:length(plot_integrated)) {
    title <- plot_integrated[[i]][["labels"]][["title"]]
    coord <- GetTissueCoordinates(object = seu_obj@images[[title]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    spatial_plot_list[[i]] <- plot_integrated[[i]] + 
      theme(plot.background = element_rect(fill="white"), legend.position = "None", aspect.ratio = myratio)
    
    spatial_plot_list[[i]][["labels"]][["title"]] <- substr(title, start = 2, stop = nchar(title))
  }
  
  ncol_grid <- floor(sqrt(length(spatial_plot_list)))
  
  grid <- grid.arrange(grobs = spatial_plot_list, ncol = ncol_grid)
  
  title <- ggdraw() + draw_label(title_integrated, size = 12)
  
  plot_umap <- DimPlot(seu_obj, reduction = "umap", cols = mapped_colors) + 
    theme(legend.position = "top", legend.key.width = unit(0.4, "inches"), legend.justification = "center") +
    guides(color=guide_legend(override.aes = list(size=5),  nrow=ceiling(no_clusters/8), byrow = TRUE))
  print(plot_umap)
  
  # legend <- get_legend(plot_umap)
  legend <- get_plot_component(plot_umap, 'guide-box-top', return_all = TRUE)
  print(legend)
  
  grid_plot <- plot_grid(title, legend, grid, nrow = 3, rel_heights = c(0.10, 0.25, 1.5)) +
    theme(panel.background = element_rect(fill="white")) 
  
  ggsave(filename=paste0(output_folder, title_integrated, ".tiff"), plot = grid_plot, 
         dpi = 300, width = 8.3, height = 11.7)
}


draw_UMAP_clusters <- function(seu_obj, output_folder, name, label, res, mapped_colors) {
  title_umap <- paste0(name, " UMAP res=", res, " - clusters", if(label==TRUE) " labeled")
  
  plot_umap <- DimPlot(seu_obj, reduction = "umap", label = label, label.box = FALSE,  
                       label.size = 4, repel = TRUE, cols = mapped_colors) + 
    ggtitle(title_umap) + 
    theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position = "top", legend.justification = "center") +
    guides(color=guide_legend(override.aes = list(size = 5),  nrow = ceiling(length(mapped_colors)/12), byrow = TRUE))
  
  ggsave(filename = paste0(output_folder, title_umap, ".tiff"), plot = plot_umap, 
         dpi = 300, width = 8, height = 10)
}


draw_UMAP_timepoints <- function(seu_obj, output_folder, name, timepoints, label = FALSE) {
  seu_obj@meta.data$timepoint <- factor(seu_obj@meta.data$timepoint, levels = timepoints)
  
  title_umap <- paste0(name, " UMAP - timepoints", if(label==TRUE) " labeled")
  
  plot_umap <- DimPlot(seu_obj, reduction = "umap", group.by = "timepoint",
                       label = label, label.box = FALSE,  label.size = 4, repel = TRUE) + 
    ggtitle(title_umap) + 
    theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = "top", legend.justification = "center", 
          legend.text = element_text(size = 7), legend.key.size = unit(0.5, "lines"), legend.spacing.y = unit(0, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 5), nrow = ceiling(length(timepoints)/8), byrow = TRUE))
  
  ggsave(filename=paste0(output_folder, title_umap, ".tiff"), plot = plot_umap, 
         dpi = 300, width = 8, height = 8)
}


draw_UMAP_timepoints_sections <- function(seu_obj, output_folder, name, levels, mapped_colors, label = FALSE) {
  seu_obj@meta.data$orig.ident <- factor(seu_obj@meta.data$orig.ident, levels = levels)
  title_umap <- paste0(name, " UMAP - per timepoint per section", if(label==TRUE) " labeled")
  
  plot_umap <- DimPlot(seu_obj, reduction = "umap", group.by = "orig.ident",
                       label = label, label.box = FALSE,  label.size = 4, repel = TRUE, 
                       cols = mapped_colors) + 
    ggtitle(title_umap) + 
    theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = "top", legend.justification = "center", 
          legend.text = element_text(size = 7), legend.key.size = unit(0.5, "lines"), legend.spacing.y = unit(0, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 3), nrow = ceiling(length(levels)/8), byrow = TRUE))
  
  ggsave(filename = paste0(output_folder, title_umap, ".tiff"), plot = plot_umap, 
         dpi = 300, width = 8, height = 8)
  
  print(plot_umap)
  
}
