library(Seurat)
library(data.table)

SEED <- 314
set.seed(SEED)

run_name <- "BARLEY GER 0HAI example 30th June 2024"
input_folder <- "/media/user/DATA_SSD_2/Barley_germination/spaceranger_count_2_1_1_ref_IBSC_v2_51_json_loupe/"
output_folder <- paste0("/home/user/Documents/Barley/GER/", run_name, "/")

name <- "0HAI"
sections <- c("A", "B", "C", "D")
timepoints <- c("0HAI")
levels <- c("0HAI_A", "0HAI_B", "0HAI_C", "0HAI_D")

cut_off_nUMIs_list <- c(20, 20, 30, 30, 30, 100, 150, 200, 200, 200, 300, 300)
cut_off_nGenes_list <- c(10, 20, 10, 20, 30, 100, 150, 150, 200, 250, 250, 300)
cut_off_nUMIs <- list("0HAI" = 30)
cut_off_nGenes <- list("0HAI" = 10)

regress_out <- TRUE
umap_clusters <- TRUE
dims <- 30
resolution_list <- seq(0.1, 2, 0.1)

custom_palette_org <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2",
                    "#7E6148B2", "#B09C85B2", "#FF410DCC", "#6EE2FFCC", "#F7C530CC", "#95CC5ECC", "#D0DFE6CC", "#F79D1ECC",
                    "#748AA6CC", "#00468BB2", "#ED0000B2", "#42B540B2", "#0099B4B2", "#925E9FB2", "#FDAF91B2", "#AD002AB2",
                    "#ADB6B6B2", "#1B1919B2", "#FED439CC", "#709AE1CC", "#8A9197CC", "#D2AF81CC", "#FD7446CC", "#D5E4A2CC",
                    "#197EC0CC", "#F05C3BCC", "#46732ECC", "#71D0F5CC", "#370335CC", "#075149CC", "#C80813CC", "#91331FCC",
                    "#1A9993CC", "#FD8CC1CC", "#FAFD7CCC", "#82491ECC", "#24325FCC", "#B7E4F9CC", "#FB6467CC", "#526E2DCC",
                    "#FF6F00B2", "#C71000B2", "#008EA0B2", "#8A4198B2", "#5A9599B2", "#CC0C00B2", "#5C88DAB2", "#84BD00B2",
                    "#FFCD00B2", "#7C878EB2", "#171717B2", "#7D0226B2", "#300049B2", "#165459B2", "#3F2327B2", "#0B1948B2",
                    "#E71012B2", "#555555B2", "#193006B2", "#A8450CB2", "#D51317B2", "#F39200B2", "#EFD500B2", "#95C11FB2",
                    "#007B3DB2", "#31B7BCB2", "#0094CDB2", "#164194B2", "#6F286AB2", "#706F6FB2", "#F9CA24CC", "#F0932BCC",
                    "#EB4D4BCC", "#6AB04CCC", "#C7ECEECC", "#22A6B3CC", "#BE2EDDCC", "#4834D4CC", "#130F40CC", "#535C68CC")


alpha <- 0.75
color_red <- "#DF2E38"
color_green <- "#A7DF97"

seu_obj_unfiltered_list <- NULL
seu_obj_filtered_list <- NULL

for(timepoint in timepoints) {
  for(section in sections) {
    dataset_name <- paste0(timepoint, "_", section)
    data_folder <- paste0(input_folder, dataset_name, "/outs")
    
    # STEP 0 - Load the data
    filtered_feature_bc_matrix_path <- "./filtered_feature_bc_matrix.h5"
    
    seu_obj <- Load10X_Spatial(
      data.dir = data_folder, 
      filename = filtered_feature_bc_matrix_path,
      assay = "Spatial",
      slice = dataset_name, 
      filter.matrix = TRUE, 
    )
    
    seu_obj <- RenameCells(object = seu_obj, add.cell.id = dataset_name)
    seu_obj@meta.data$section <- section
    seu_obj@meta.data$timepoint <- timepoint
    seu_obj@meta.data$orig.ident <- dataset_name
    seu_obj$section <- dataset_name
    Idents(seu_obj) <- dataset_name
    
    feature_tsv_file_path <- paste0(data_folder, "/filtered_feature_bc_matrix/features.tsv.gz")
    feature_df_data <- fread(feature_tsv_file_path, header=FALSE, sep='\t')
    seu_obj@assays[["Spatial"]]@counts@Dimnames[[1]] <- feature_df_data$V1
    seu_obj@assays[["Spatial"]]@data@Dimnames[[1]] <- feature_df_data$V1
    
    # STEP 1 - INSPECT DATA
    source("./01_Inspect_data.R")
    inspect_data(seu_obj, output_folder, dataset_name, alpha)
    
    seu_obj_unfiltered_list <- append(seu_obj_unfiltered_list, seu_obj)
    
    # STEP 2 - FILTER DATA
    source("./02_Quality_control_(Filter_out).R")
    filter_data_graphs(seu_obj, output_folder, cut_off_nUMIs_list, cut_off_nGenes_list, 
                       alpha, color_red, color_green, dataset_name)
    seu_obj <- subset(seu_obj, subset = nFeature_Spatial >= cut_off_nGenes & nCount_Spatial >= cut_off_nUMIs)  
    seu_obj_filtered_list <- append(seu_obj_filtered_list, seu_obj)
  }
  
  # STEP 3 - INSPECT MERGED UNFILTERED DATA
  seu_obj_unfiltered_merged <- Reduce(merge, seu_obj_unfiltered_list)
  inspect_data_merged(seu_obj_unfiltered_merged, output_folder, name, alpha)
}

# MERGE FILTERED DATA 
seu_obj_merged <- Reduce(merge, seu_obj_filtered_list)

# STEP 4 - NORMALISATION
if(regress_out) {
  seu_obj_merged <- SCTransform(seu_obj_merged, assay = "Spatial", vars.to.regress =  c('section'),
                      variable.features.n = NULL, variable.features.rv.th = 1, return.only.var.genes = FALSE)
} else { 
  seu_obj_merged <- SCTransform(seu_obj_merged, assay = "Spatial", variable.features.n = NULL, return.only.var.genes = FALSE)
}

# STEP 5 - Pearson correlation coefficient
source("./05_Pearson.R")
calculate_Pearson_correlation(seu_obj_merged, output_folder, timepoints, sections)

# STEP 6 - Highly variable features
source("./06_Highly_variable_features.R")
find_highly_variable_features(seu_obj_merged, output_folder, name)
find_spatially_variable_feature_after_norm(seu_obj_merged, output_folder, name)

seu_obj_merged <- RunPCA(seu_obj_merged, assay = "SCT", npcs = 2*dims)

source("./07_Elbow.R")
draw_elbow_plot(seu_obj_merged, name, output_folder, dims)

if(umap_clusters)
{ 
  seu_obj_merged <- RunUMAP(seu_obj_merged, reduction = "pca", dims = 1:dims)
}

seu_obj_merged <- FindNeighbors(seu_obj_merged, reduction = "pca", dims = 1:dims)

source("./08_Clustering.R")
seu_obj_merged <- perform_clustering_and_marker_genes(seu_obj_merged, output_folder, resolution_list, umap_clusters, dims, 
                                                      custom_palette_org, name, alpha, timepoints, levels)

source("./08_02_Cluster_tree.R")
draw_cluster_tree(seu_obj_merged, output_folder, name)


