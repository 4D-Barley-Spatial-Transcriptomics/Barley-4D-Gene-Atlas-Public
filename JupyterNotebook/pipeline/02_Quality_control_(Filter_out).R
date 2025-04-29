filter_data_graphs <- function(seu_obj, output_folder, cut_off_nUMIs_list, cut_off_nGenes_list, 
                               alpha, color_red, color_green, dataset_name) {
  
  output_folder <- paste0(output_folder, "02_Quality_control_(Filter_out)/")
  
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  for(i in seq(length(cut_off_nUMIs_list))) {
    spatial_coordinates <- seu_obj@images[[paste0("X", dataset_name)]]@coordinates
    
    nUMIs_min <- cut_off_nUMIs_list[i]
    nGenes_min <- cut_off_nGenes_list[i]
    
    qc_nUMIs <-  seu_obj@meta.data[["nCount_Spatial"]] >= nUMIs_min
    qc_nGenes <-  seu_obj@meta.data[["nFeature_Spatial"]] >= nGenes_min
    qc_nUMIs_nGenes <- qc_nUMIs & qc_nGenes
    
    percent_discarded <- round(sum(!qc_nUMIs_nGenes) / length(qc_nUMIs_nGenes) * 100, 2)
    
    print(paste0(dataset_name, " nUMIs cut of ", nUMIs_min, " & nGenes cut off ", nGenes_min, 
                                "\n", paste(capture.output(table(qc_nUMIs_nGenes)), collapse = "\n"), 
                                "\nPercentage of discarded spots: ", percent_discarded, "%"))
    seu_obj$qc_nUMIs_nGenes <- as.numeric(qc_nUMIs_nGenes)
    
    title <- paste0(dataset_name, " nUMIs above ", nUMIs_min, ", nGenes above ", nGenes_min)
    
    plot_qc_nUMIs_nGenes <- SpatialFeaturePlot(seu_obj, features = "qc_nUMIs_nGenes", pt.size.factor = 1,
                                               alpha = c(alpha, alpha), crop = FALSE) + 
      ggtitle(title) +
      scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 1.0), 
                           low = color_red, high = color_green, midpoint = 0.5) +
      theme(legend.position = "right", plot.title = element_text(size = 12, hjust = 0.5))
    
    ggsave(paste0(output_folder, title, ".tiff"), plot_qc_nUMIs_nGenes, dpi = 300,
           width = 8.3, height = 11.7)
  }
}
