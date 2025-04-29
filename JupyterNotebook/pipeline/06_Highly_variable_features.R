library(cowplot)

find_highly_variable_features <- function(seu_obj, output_folder, name) {
  output_folder <- paste0(output_folder, "06_Highly_variable_features/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  output_folder <- paste0(output_folder, "FindVariableFeatures/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  title <- paste0(name, " Top 10 Highly Variable Features before norm")
  data_highly_var <- FindVariableFeatures(seu_obj, assay = "Spatial", 
                                          selection.method = "vst", nfeatures = 2000)
  
  top10 <- head(VariableFeatures(data_highly_var), 10)
  
  plot1 <- VariableFeaturePlot(data_highly_var, assay = "Spatial")
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  ggsave(paste0(output_folder, title, ".tiff"), plot2, dpi = 300, width = 8.3, height = 11.7)
  
  for (top in top10) {
    plot_integrated <- SpatialFeaturePlot(data_highly_var, features = c(top), ncol = 1, 
                                          pt.size.factor = 1.7, crop = TRUE)
    title_ <- paste0(name, " Spatial Feature Plot ", top)
    
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
    
    plot <- plot_integrated[[1]] +
      theme(legend.position = "top", legend.key.width = unit(0.4, "inches"))
    
    print(plot)
    # legend <- get_legend(plot)
    legend <- get_plot_component(plot, 'guide-box-top', return_all = TRUE)
    
    title <- ggdraw() + draw_label(title_)
    grid_plot <- plot_grid(title, legend, grid, nrow = 3, rel_heights = c(0.1, 0.1, 1)) +
      theme(panel.background = element_rect(fill = "white"))
    
    ggsave(paste0(output_folder, title_, ".tiff"), grid_plot, dpi = 300, width = 8.3, height = 11.7)
  }
}

find_spatially_variable_feature_after_norm <- function(seu_obj, output_folder, name) {
  output_folder <- paste0(output_folder, "06_Highly_variable_features/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  output_folder <- paste0(output_folder, "FindSpatiallyVariableFeatures_after_norm/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  title <- paste0(name, " Top 10 Spatially Highly Variable Features after norm")
  data_highly_var <- FindSpatiallyVariableFeatures(seu_obj, assay = "SCT", 
                                                   selection.method = "moransi", nfeatures = 2000)
  
  top10 <- head(VariableFeatures(data_highly_var), 10)
  
  plot1 <- VariableFeaturePlot(data_highly_var, assay = "SCT")
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  ggsave(paste0(output_folder, title, ".tiff"), plot2, dpi = 300, width = 8.3, height = 11.7)
  
  for (top in top10) {
    plot_integrated <- SpatialFeaturePlot(data_highly_var, features = c(top), ncol = 1, 
                                          pt.size.factor = 1.7, crop = TRUE)
    title_ <- paste0(name, " Spatial Feature Plot ", top)
    
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
    
    plot <- plot_integrated[[1]] +
      theme(legend.position = "top", legend.key.width = unit(0.4, "inches"))
    # legend <- get_legend(plot)
    legend <- get_plot_component(plot, 'guide-box-top', return_all = TRUE)
    
    title <- ggdraw() + draw_label(title_)
    grid_plot <- plot_grid(title, legend, grid, nrow = 3, rel_heights = c(0.1, 0.1, 1)) +
      theme(panel.background = element_rect(fill = "white"))
    
    ggsave(paste0(output_folder, title_, ".tiff"), grid_plot, dpi = 300, width = 8.3, height = 11.7)
    
  }
}