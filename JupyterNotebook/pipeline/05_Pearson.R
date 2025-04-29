library(cowplot)
library(patchwork)
library(dplyr)
library(Seurat)
library(gridExtra)
library(ggplot2)

calculate_Pearson_correlation <- function(seu_obj, output_folder, timepoints, sections) {
  output_folder <- paste0(output_folder, "05_Pearson_correlation/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  # discared genes that are expressed in less than 5 spots
  DefaultAssay(seu_obj) <- "SCT"
  more_than_5_spots_expressed <- rowSums(GetAssayData(object = seu_obj, slot = "counts") > 0) >= 5
  num_false <- sum(!more_than_5_spots_expressed)
  seu_obj.filtered <- seu_obj[more_than_5_spots_expressed, ]
  
  section_pairs <- combn(sections, 2)
  
  for(timepoint in timepoints) {
    
    col_list <- NULL
    plot_list <- list()
    plot_num <- 0
    
    for(section in sections) {
      subset <- subset(seu_obj.filtered, subset = orig.ident == paste0(timepoint, "_", section))
      col_list[[section]] <- rowMeans(subset)
    }
    
    data <- as.data.frame(col_list)
    
    for (i in seq(ncol(section_pairs))) {
      pair <- section_pairs[, i]
      print(pair)
      if(pair[1] != pair[2]) {
        
        correlation <- cor(data[[pair[1]]], data[[pair[2]]])
        title <- paste0("R = ", round(correlation, 5))
        
        plot <- ggplot(data, aes(x = .data[[pair[1]]], y = .data[[pair[2]]])) +
          geom_point(size = 1) +
          geom_smooth(method = "lm", se = FALSE, color = "black") +
          ggtitle(title) +
          labs(x = paste0("Section ", pair[1]), y = paste0("Section ", pair[2])) +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
                plot.background = element_blank(),
                panel.background = element_rect(fill="white"), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(color= "black"), 
                plot.margin = margin(15, 15, 15, 15))
        
        plot_num = plot_num + 1
        plot_list[[plot_num]] <- plot
      }
    }
    
    ncol_grid <- ceiling(sqrt(length(plot_list)))
    grid <- grid.arrange(grobs = plot_list, ncol = ncol_grid)
    
    title_integrated <- paste0(timepoint, " Pearson correlation")
    title <- ggdraw() + draw_label(title_integrated, size = 15)
    
    grid_plot <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.10, 1.5)) +
      theme(panel.background = element_rect(fill="white"))
    
    ggsave(filename = paste0(output_folder, title_integrated, ".tiff"), plot = grid_plot, 
           dpi = 300, width = 8, height = 7)
  }
}