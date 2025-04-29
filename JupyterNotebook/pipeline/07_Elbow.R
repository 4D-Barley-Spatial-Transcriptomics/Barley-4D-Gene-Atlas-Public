
draw_elbow_plot <- function(seu_obj, name, output_folder, dims) {
  
  output_folder <- paste0(output_folder, "07_Elbow/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  title <- paste0(name, " Elbow Plot")
  elbow_plot <- ElbowPlot(seu_obj, ndims = 2*dims) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5), plot.background = element_rect(fill = "white"))
  
  ggsave(paste0(output_folder, title, ".tiff"), elbow_plot, dpi = 300, width = 11.7, height = 8.3)
  
}
