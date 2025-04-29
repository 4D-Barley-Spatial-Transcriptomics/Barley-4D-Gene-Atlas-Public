library(ggplot2)
library(patchwork)
library(cowplot)
library(gridExtra)

# Draw Violin Plots, Spatial Plots for nCount_Spatial (nUMIs) and nFeature_Spatial (nGenes) and 
# Histograms for Library size

inspect_data <- function(seu_obj, output_folder, name, alpha) {
  
  output_folder <- paste0(output_folder, "01_Inspect_data/")
  
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  coord <- GetTissueCoordinates(object = seu_obj@images[[paste0("X", name)]])
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  
  color <- "#0D98BA"
  
  # Plot 1 - VlnPlot nCount_Spatial
  # nCount_Spatial = nUMIs 
  title_nCount <- paste0(name, " nCount_Spatial")
  vln_plot_nCount <- VlnPlot(seu_obj, features = "nCount_Spatial", pt.size = 0.1) + 
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    guides(fill = FALSE) +
    ggtitle(title_nCount) + 
    NoLegend() 
  
  spatial_plot_nCount <- SpatialFeaturePlot(seu_obj, features = "nCount_Spatial", crop = TRUE, 
                                            alpha = c(alpha, alpha)) + 
    ggtitle(title_nCount) +
    theme(legend.position = "top", legend.key.width = unit(0.5, "inches"),
          plot.title = element_text(size = 12, hjust=0.5), aspect.ratio = myratio)
  
  plot_nCount <- wrap_plots(vln_plot_nCount, spatial_plot_nCount)
  ggsave(paste0(output_folder, title_nCount, ".tiff"), plot_nCount, dpi = 300, 
         width = 11.7, height = 8.3)
  
  # Plot 2 - VlnPlot nFeature_Spatial
  # nFeature_Spatial = nGenes
  title_nFeature <- paste0(name, " nFeature_Spatial")
  vln_plot_nFeature <- VlnPlot(seu_obj, features = "nFeature_Spatial", pt.size = 0.1) + 
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    labs(x = NULL) +
    guides(fill = FALSE) +
    ggtitle(title_nFeature) + 
    NoLegend()
  
  spatial_plot_nFeature <- SpatialFeaturePlot(seu_obj, features = "nFeature_Spatial", crop = TRUE, 
                                              alpha = c(alpha, alpha)) + 
    ggtitle(title_nFeature) + 
    theme(legend.position = "top", legend.key.width = unit(0.5, "inches"),
          plot.title = element_text(size = 12, hjust = 0.5), aspect.ratio = myratio)
  
  plot_nFeature <- wrap_plots(vln_plot_nFeature, spatial_plot_nFeature)
  ggsave(paste0(output_folder, title_nFeature, ".tiff"), plot_nFeature, dpi = 300, 
         width = 11.7, height = 8.3)
  
  # Plot 3 - FeatureScatter
  # library size is nCount_Spatial = nUMIs
  df_data <- data.frame(nCount_Spatial = seu_obj$nCount_Spatial, nFeature_Spatial = seu_obj$nFeature_Spatial)
  
  title_plot <- paste0(name, " FeatureScatter - nUMIs (nCount_Spatial) vs nGenes (nFeature_Spatial)")
  plot <- ggplot(df_data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
    geom_point(color = color) +
    labs(title = title_plot, x = "nUMIs (nCount_Spatial)", y = "nGenes (nFeature_Spatial)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_plot, ".tiff"),
         plot = plot, dpi = 300, width = 11.7, height = 8.3)
  
  # library size is nCount_Spatial = nUMIs
  # Plot 3 - Histogram 500 library size
  df_data <- data.frame(library_size = seu_obj$nCount_Spatial)
  
  title_hist_500 <- paste0(name, " histogram library size (binwidth = 500)")
  hist_500 <- ggplot(df_data, aes(x = library_size)) +
    geom_histogram(binwidth = 500, fill = color, color = "black") +
    labs(title = title_hist_500, x = "Library size", y = "Counts") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_hist_500, ".tiff"),
         plot = hist_500, dpi = 300, width = 11.7, height = 8.3)
  
  # Plot 3 - Histogram 1200 library size
  title_hist_1200 <- paste0(name, " histogram library size (binwidth = 1200)")
  hist_1200 <- ggplot(df_data, aes(x = library_size)) +
    geom_histogram(binwidth = 1200, fill = color, color = "black") +
    labs(title = title_hist_1200, x = "Library size", y = "Counts") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_hist_1200, ".tiff"),
         plot = hist_1200, dpi = 300, width = 11.7, height = 8.3)
  
}


inspect_data_merged <- function(seu_obj, output_folder, name, alpha, levels) {
  
  output_folder <- paste0(output_folder, "01_Inspect_data/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  color <- "#0D98BA"
  
  seu_obj$orig.ident <- factor(x = seu_obj$orig.ident, levels = levels)
  
  # Plot 1 - VlnPlot nCount_Spatial (nUMIs)
  # nCount_Spatial = nUMIs 
  title_vln_plot_nCount <- paste0(name, " Vln plot nCount_Spatial (nUMIs)")
  vln_plot_nCount <- VlnPlot(seu_obj, features = "nCount_Spatial", pt.size = 0.1, group.by = "orig.ident") + 
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    labs(x = NULL) +
    guides(fill = FALSE) +
    ggtitle(title_vln_plot_nCount) + 
    NoLegend()
  
  ggsave(paste0(output_folder, title_vln_plot_nCount, ".tiff"), vln_plot_nCount, 
         dpi = 300, width = 11.7, height = 8.3)
  
  # Plot 2 - VlnPlot nCount_Spatial (nUMIs)
  title_spatial_plot_nCount <- paste0(name, " SpatialFeature plot nCount_Spatial (nUMIs)")
  spatial_plots_nCount <- SpatialFeaturePlot(seu_obj, features = "nCount_Spatial", crop = TRUE,
                                             alpha = c(alpha, alpha))
  
  spatial_plot_list <- list()
  
  for (i in 1:length(spatial_plots_nCount)) {
    title <- spatial_plots_nCount[[i]][["labels"]][["title"]]
    coord <- GetTissueCoordinates(object = seu_obj@images[[title]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    spatial_plot_list[[i]] <- spatial_plots_nCount[[i]] + 
      theme(plot.background = element_rect(fill="white"), legend.position = "None", aspect.ratio = myratio)
    
    spatial_plot_list[[i]][["labels"]][["title"]] <- substr(title, start=2, stop=nchar(title))
  }
  
  ncol_grid <- floor(sqrt(length(spatial_plot_list)))
  grid <- grid.arrange(grobs = spatial_plot_list, ncol = ncol_grid)  
  
  title <- ggdraw() + draw_label(title_spatial_plot_nCount)
  
  temp_plot <- SpatialFeaturePlot(seu_obj, features = "nCount_Spatial") + 
    theme(legend.position = "top", legend.key.width = unit(10, "inches"), legend.justification = "center") +
    guides(color=guide_legend(nrow=1, byrow = TRUE))
  
  legend <- get_legend(temp_plot)
  
  grid_plot <- plot_grid(title, legend, grid, nrow=3, rel_heights = c(0.1, 0.15, 1)) +
    theme(panel.background = element_rect(fill="white"))
  
  ggsave(paste0(output_folder, title_spatial_plot_nCount, ".tiff"), grid_plot, 
         dpi = 300, width = 8.3, height = 11.7)
  
  # Plot 3 - VlnPlot nFeature_Spatial (nGenes)
  # nFeature_Spatial
  title_vln_plot_nFeature <- paste0(name, " Vln plot nFeature_Spatial (nGenes)")
  vln_plot_nFeature <- VlnPlot(seu_obj, features = "nFeature_Spatial", pt.size = 0.1, group.by = "orig.ident") + 
    theme(plot.title = element_text(size = 12, hjust = 0.5)) +
    labs(x = NULL) +
    guides(fill = FALSE) +
    ggtitle(title_vln_plot_nFeature) + 
    NoLegend()
  
  ggsave(paste0(output_folder, title_vln_plot_nFeature, ".tiff"), vln_plot_nFeature, 
         dpi = 300, width = 11.7, height = 8.3)
  
  # Plot 4 - SpatialFeature nFeature_Spatial (nGenes)
  title_spatial_plot_nFeature <- paste0(name, " SpatialFeature plot nFeature_Spatial (nGenes)")
  spatial_plots_nFeature <- SpatialFeaturePlot(seu_obj, features = "nFeature_Spatial", crop = TRUE,
                                               alpha = c(alpha, alpha))
  
  spatial_plot_list <- list()
  
  for (i in 1:length(spatial_plots_nFeature)) {
    title <- spatial_plots_nFeature[[i]][["labels"]][["title"]]
    coord <- GetTissueCoordinates(object = seu_obj@images[[title]])
    myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
    spatial_plot_list[[i]] <- spatial_plots_nFeature[[i]] + 
      theme(plot.background = element_rect(fill="white"), legend.position = "None", aspect.ratio = myratio)
    
    spatial_plot_list[[i]][["labels"]][["title"]] <- substr(title, start = 2, stop=nchar(title))
  }
  
  ncol_grid <- floor(sqrt(length(spatial_plot_list)))
  grid <- grid.arrange(grobs = spatial_plot_list, ncol = ncol_grid)  
  
  title <- ggdraw() + draw_label(title_spatial_plot_nFeature)
  
  temp_plot <- SpatialFeaturePlot(seu_obj, features = "nFeature_Spatial") + 
    theme(legend.position = "top", legend.key.width = unit(10, "inches"), legend.justification = "center") +
    guides(color=guide_legend(nrow=1, byrow = TRUE))
  
  legend <- get_legend(temp_plot)
  
  grid_plot <- plot_grid(title, legend, grid, nrow=3, rel_heights = c(0.1, 0.15, 1)) +
    theme(panel.background = element_rect(fill="white"))
  
  ggsave(paste0(output_folder, title_spatial_plot_nFeature, ".tiff"), grid_plot, 
         dpi = 300, width = 8.3, height = 11.7)
 
  # Plot 5 - FeatureScatter nFeature_Spatial (nGenes)
  df_data <- data.frame(nCount_Spatial = seu_obj$nCount_Spatial, nFeature_Spatial = seu_obj$nFeature_Spatial)
  
  title_plot <- paste0(name, " FeatureScatter - nUMIs (nCount_Spatial) vs nGenes (nFeature_Spatial)")
  plot <- ggplot(df_data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
    geom_point(color = color) +
    labs(title = title_plot, x = "nUMIs (nCount_Spatial)", y = "nGenes (nFeature_Spatial)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_plot, ".tiff"),
         plot = plot, dpi = 300, width = 11.7, height = 8.3)

  # Plot 6 - Histogram 500 library size
  df_data <- data.frame(library_size = seu_obj$nCount_Spatial)
  
  title_hist_500 <- paste0(name, " histogram library size (binwidth = 500)")
  hist_500 <- ggplot(df_data, aes(x = library_size)) +
    geom_histogram(binwidth = 500, fill = color, color = "black") +
    labs(title = title_hist_500, x = "Library size", y = "Counts") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_hist_500, ".tiff"),
         plot = hist_500, dpi = 300, width = 11.7, height = 8.3)

  # Plot 7 - Histogram 1200 library size
  title_hist_1200 <- paste0(name, " histogram library size (binwidth = 1200)")
  hist_1200 <- ggplot(df_data, aes(x = library_size)) +
    geom_histogram(binwidth = 1200, fill = color, color = "black") +
    labs(title = title_hist_1200, x = "Library size", y = "Counts") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(output_folder, title_hist_1200, ".tiff"),
         plot = hist_1200, dpi = 300, width = 11.7, height = 8.3)
  
}

