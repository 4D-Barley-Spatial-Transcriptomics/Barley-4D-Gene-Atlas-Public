library(clustree)

draw_cluster_tree <- function(data, output_folder, name) {
  
  output_folder <- paste0(output_folder, "08_Clustering/")
  if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)
  
  title <- paste0(name, " Cluster tree")
  cluster_tree <- clustree(data, prefix = "SCT_snn_res.")
  ggsave(filename = paste0(output_folder, title, ".pdf"), plot = cluster_tree)
  
}