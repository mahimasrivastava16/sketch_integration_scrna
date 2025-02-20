remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
update.packages(oldPkgs = c("withr", "rlang"))
if(!require(dplyr)) install.packages("dplyr")
if(!require(Seurat)) install.packages("Seurat")
if(!require(HGNChelper)) install.packages("HGNChelper")
install.packages("openxlsx")
install.packages("hdf5r")
BiocManager::install("SingleR", force = TRUE)
BiocManager::install("scRNAseq")  # For built-in references

library(SingleR)
library(SummarizedExperiment)
library(hdf5r)
install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", "https://cloud.r-project.org"))
library(Seurat)
#library(Azimuth)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
plan("sequential")
library(SeuratData)
library(SeuratWrappers)
options(future.globals.maxSize = 60 * 1024^3)  # Set to 20 GB
library(dplyr)
gc() # Garbage collection


metadata <- read.csv('/Users/mahimasrivastava/Downloads/Single Cell Metadata Annotation - Sheet7.csv')
# Function to add stage-wise metadata to Seurat objects 
add_stage_metadata <- function(seurat_obj, metadata) {   # Extract the sample name from Seurat object's orig.ident   
  sample_name <- seurat_obj@meta.data$orig.ident[1]      # Filter the metadata for the corresponding sample   
  sample_info <- dplyr::filter(metadata, GSM.ID == sample_name)      # Add the Tumor.Stage metadata if the sample is found   
  if (nrow(sample_info) > 0) {     
    seurat_obj@meta.data$Tumor_Stage <- sample_info$Tumor.Grade   } 
  else {     
    warning(paste("No metadata found for sample:", sample_name))   }      
  return(seurat_obj)    }  

# Load 10X data and create Seurat objects 
main_dir <- '/Users/mahimasrivastava/Downloads/GSE1615492'
sub_dirs <- list.dirs(main_dir, recursive = FALSE) 
seurat_list <- list()  
for (dir in sub_dirs) {   
  subfolder_name <- basename(dir)  
  data <- Read10X(data.dir = dir)   
  seurat_list[[subfolder_name]] <- CreateSeuratObject(counts = data, project = subfolder_name)   
  cat("Created Seurat object for:", subfolder_name, "\n") } 

# Add stage-wise metadata to all Seurat objects 
seurat_list <- lapply(seurat_list, function(obj) add_stage_metadata(obj, metadata))  
rm(metadata)
rm(data)
# Add mitochondrial gene percentage 
seurat_list <- lapply(seurat_list, function(seurat_obj) {   
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")   
  return(seurat_obj) }) 



# Directory to save pre-filter QC plots 
#pre_filter_dir <- "/Users/mahimasrivastava/scrna_seurat/QC_plots/new_small_dataset" 
#dir.create(pre_filter_dir, showWarnings = FALSE)  
# Generate and save pre-filter QC plots 
#lapply(names(seurat_list), function(sample_name) { 
# Access the Seurat object for the current sample 
# seurat_obj <- seurat_list[[sample_name]]  
# Generate the violin plot for nFeature_RNA, nCount_RNA, and percent.mt 
#vln_plot_pre_filter <- VlnPlot(seurat_obj, 
#                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#                             ncol = 3) +  #ggtitle(paste("Pre-filter QC Metrics for", sample_name))  
# Generate individual feature scatter plot 
#feature_scatter <- FeatureScatter(seurat_obj,  
#         feature1 = "nCount_RNA",  #        
#feature2 = "nFeature_RNA") + 
#    ggtitle(paste("Feature Scatter for", sample_name))  
# Save the individual scatter plot 
#  ggsave(filename = paste0(pre_filter_dir, "/", sample_name, "_feature_scatter.png"), 
#        plot = feature_scatter, width = 10, height = 5) 
# Save the violin plot #  ggsave(filename = paste0(pre_filter_dir, "/", sample_name, "_pre_filter_qc_plot.png"), 
#        plot = vln_plot_pre_filter, width = 10, height = 5) #})  



# Apply filtering criteria to Seurat objects (example filters applied) 
seurat_list[["GSM4909253"]] <- subset(seurat_list[["GSM4909253"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6500 &  percent.mt < 20) 
seurat_list[["GSM4909254"]] <- subset(seurat_list[["GSM4909254"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &  percent.mt < 20) 
seurat_list[["GSM4909257"]] <- subset(seurat_list[["GSM4909257"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 30) 
seurat_list[["GSM4909261"]] <- subset(seurat_list[["GSM4909261"]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 &  percent.mt < 20) 
seurat_list[["GSM4909263"]] <- subset(seurat_list[["GSM4909263"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 20) 
seurat_list[["GSM4909265"]] <- subset(seurat_list[["GSM4909265"]], subset = nFeature_RNA > 500 & nFeature_RNA < 9000 &  percent.mt < 20) 
seurat_list[["GSM4909266"]] <- subset(seurat_list[["GSM4909266"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
#seurat_list[["GSM4909267"]] <- subset(seurat_list[["GSM4909267"]], subset = nFeature_RNA > 500 & nFeature_RNA < 4500 &  percent.mt < 20)  
seurat_list[["GSM4909268"]] <- subset(seurat_list[["GSM4909268"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6500 &  percent.mt < 20) 
#seurat_list[["GSM4909269"]] <- subset(seurat_list[["GSM4909269"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 20) 
seurat_list[["GSM4909270"]] <- subset(seurat_list[["GSM4909270"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
seurat_list[["GSM4909271"]] <- subset(seurat_list[["GSM4909271"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  percent.mt < 20) 
seurat_list[["GSM4909272"]] <- subset(seurat_list[["GSM4909272"]], subset = nFeature_RNA > 500 & nFeature_RNA < 2700 &  percent.mt < 20) 
#seurat_list[["GSM4909273"]] <- subset(seurat_list[["GSM4909273"]], subset = nFeature_RNA > 500 & nFeature_RNA < 4500 &  percent.mt < 20) 
seurat_list[["GSM4909274"]] <- subset(seurat_list[["GSM4909274"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6500 &  percent.mt < 20) 
#seurat_list[["GSM4909275"]] <- subset(seurat_list[["GSM4909275"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 20) 
seurat_list[["GSM4909276"]] <- subset(seurat_list[["GSM4909276"]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 &  percent.mt < 20) 
seurat_list[["GSM4909277"]] <- subset(seurat_list[["GSM4909277"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  percent.mt < 20) 
seurat_list[["GSM4909278"]] <- subset(seurat_list[["GSM4909278"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3800 &  percent.mt < 10) 
seurat_list[["GSM4909279"]] <- subset(seurat_list[["GSM4909279"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 20) 
seurat_list[["GSM4909280"]] <- subset(seurat_list[["GSM4909280"]], subset = nFeature_RNA > 500 & nFeature_RNA < 4500 &  percent.mt < 20) 
seurat_list[["GSM4909281"]] <- subset(seurat_list[["GSM4909281"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
seurat_list[["GSM4909282"]] <- subset(seurat_list[["GSM4909282"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 20) 
seurat_list[["GSM4909283"]] <- subset(seurat_list[["GSM4909283"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
seurat_list[["GSM4909284"]] <- subset(seurat_list[["GSM4909284"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 30) 
seurat_list[["GSM4909285"]] <- subset(seurat_list[["GSM4909285"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 20) 
seurat_list[["GSM4909286"]] <- subset(seurat_list[["GSM4909286"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 30) 
seurat_list[["GSM4909287"]] <- subset(seurat_list[["GSM4909287"]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 &  percent.mt < 40) 
seurat_list[["GSM4909288"]] <- subset(seurat_list[["GSM4909288"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909289"]] <- subset(seurat_list[["GSM4909289"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 25) 
seurat_list[["GSM4909290"]] <- subset(seurat_list[["GSM4909290"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  percent.mt < 25) 
seurat_list[["GSM4909291"]] <- subset(seurat_list[["GSM4909291"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3300 &  percent.mt < 20) 
seurat_list[["GSM4909292"]] <- subset(seurat_list[["GSM4909292"]], subset = nFeature_RNA > 400 & nFeature_RNA < 4000 &  percent.mt < 25) 
seurat_list[["GSM4909293"]] <- subset(seurat_list[["GSM4909293"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 30) 
seurat_list[["GSM4909294"]] <- subset(seurat_list[["GSM4909294"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5500 &  percent.mt < 30) 
seurat_list[["GSM4909295"]] <- subset(seurat_list[["GSM4909295"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 30) 
seurat_list[["GSM4909296"]] <- subset(seurat_list[["GSM4909296"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
seurat_list[["GSM4909297"]] <- subset(seurat_list[["GSM4909297"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &  percent.mt < 20) 
seurat_list[["GSM4909298"]] <- subset(seurat_list[["GSM4909298"]], subset = nFeature_RNA > 300 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909299"]] <- subset(seurat_list[["GSM4909299"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 20) 
seurat_list[["GSM4909300"]] <- subset(seurat_list[["GSM4909300"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  percent.mt < 20) 
seurat_list[["GSM4909301"]] <- subset(seurat_list[["GSM4909301"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &  percent.mt < 20) 
seurat_list[["GSM4909302"]] <- subset(seurat_list[["GSM4909302"]], subset = nFeature_RNA > 300 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909303"]] <- subset(seurat_list[["GSM4909303"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 30) 
seurat_list[["GSM4909304"]] <- subset(seurat_list[["GSM4909304"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &  percent.mt < 30) 
seurat_list[["GSM4909305"]] <- subset(seurat_list[["GSM4909305"]], subset = nFeature_RNA > 300 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909306"]] <- subset(seurat_list[["GSM4909306"]], subset = nFeature_RNA > 300 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909307"]] <- subset(seurat_list[["GSM4909307"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3800 &  percent.mt < 20) 
#seurat_list[["GSM4909308"]] <- subset(seurat_list[["GSM4909308"]], subset = nFeature_RNA > 500 & nFeature_RNA < 4700 &  percent.mt < 20) 
seurat_list[["GSM4909309"]] <- subset(seurat_list[["GSM4909309"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  percent.mt < 40) 
#seurat_list[["GSM4909310"]] <- subset(seurat_list[["GSM4909310"]], subset = nFeature_RNA > 300 & nFeature_RNA < 5500 &  percent.mt < 40) 
seurat_list[["GSM4909311"]] <- subset(seurat_list[["GSM4909311"]], subset = nFeature_RNA > 300 & nFeature_RNA < 1200 &  percent.mt < 20) 
#seurat_list[["GSM4909312"]] <- subset(seurat_list[["GSM4909312"]], subset = nFeature_RNA > 500 & nFeature_RNA < 2500 &  percent.mt < 20) 
seurat_list[["GSM4909313"]] <- subset(seurat_list[["GSM4909313"]], subset = nFeature_RNA > 500 & nFeature_RNA < 3000 &  percent.mt < 20) 
#seurat_list[["GSM4909314"]] <- subset(seurat_list[["GSM4909314"]], subset = nFeature_RNA > 300 & nFeature_RNA < 2000 &  percent.mt < 20) 
seurat_list[["GSM4909315"]] <- subset(seurat_list[["GSM4909315"]], subset = nFeature_RNA > 500 & nFeature_RNA < 6500 &  percent.mt < 30) 
#seurat_list[["GSM4909316"]] <- subset(seurat_list[["GSM4909316"]], subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &  percent.mt < 40) 
seurat_list[["GSM4909317"]] <- subset(seurat_list[["GSM4909317"]], subset = nFeature_RNA > 500 & nFeature_RNA < 5800 &  percent.mt < 30) 
#seurat_list[["GSM4909318"]] <- subset(seurat_list[["GSM4909318"]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 &  percent.mt < 40) 
seurat_list[["GSM4909319"]] <- subset(seurat_list[["GSM4909319"]], subset = nFeature_RNA > 300 & nFeature_RNA < 1300 &  percent.mt < 30) 
seurat_list[["GSM4909320"]] <- subset(seurat_list[["GSM4909320"]], subset = nFeature_RNA > 300 & nFeature_RNA < 1300 &  percent.mt < 20) 
#seurat_list[["GSM4909321"]] <- subset(seurat_list[["GSM4909321"]], subset = nFeature_RNA > 300 & nFeature_RNA < 1700 &  percent.mt < 20)    

# Modify cell names in each Seurat object to make them unique
seurat_list_unique <- lapply(names(seurat_list), function(sample_name) {
  seurat_obj <- seurat_list[[sample_name]]
  # Append the sample name to the cell names to ensure uniqueness
  seurat_obj$cell_names <- paste0(sample_name, "_", colnames(seurat_obj))
  colnames(seurat_obj) <- seurat_obj$cell_names
  return(seurat_obj)
})


rm(seurat_list)


merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
saveRDS(merged_seurat,"merged_seurat.rds") #split rna

rm(seurat_list_unique)
library(future)
plan("multisession", workers = 8) # Adjust workers to your CPU cores

merged_seurat <- readRDS('/Users/mahimasrivastava/merged_seurat.rds')
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$Tumor_Stage)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- SketchData(merged_seurat, ncells = 15000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(merged_seurat) <- "sketch"
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
# integrate the datasets
merged_seurat <- IntegrateLayers(merged_seurat, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                          dims = 1:30, k.anchor = 20, 
                          reference = which(Layers(merged_seurat, search = "data") == "data.None"))
#cluster the integrated data
merged_seurat <- FindNeighbors(merged_seurat, reduction = "integrated.rpca", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.7)
merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.rpca", dims = 1:30, return.model = T)

merged_seurat[["sketch"]] <- JoinLayers(merged_seurat[["sketch"]])
bc_markers <- FindMarkers(merged_seurat, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
head(bc_markers)
DimPlot(merged_seurat, group.by = "Tumor_Stage", reduction = "umap")

DimPlot(merged_seurat, group.by = "seurat_clusters", reduction = "umap")


all_markers = FindAllMarkers(merged_seurat, 
                             only.pos = TRUE, # genes more expressed in the cluster compared
                             min.pct = 0.25, # % of cell expressing the marker
                             logfc.threshold = 0.25) 

# Number of markers identified by cluster
table(all_markers$cluster)

# Save in a table 3 genes the most differentially expressed in one cluster VS all the other clusters
top3_markers = as.data.frame(all_markers %>% 
                               group_by(cluster) %>% 
                               top_n(n = 3, wt = avg_log2FC))

# Create a dotplot the vidualise the expression of genes by cluster
Seurat::DotPlot(merged_seurat, features = unique(top3_markers$gene)) +
  # this second part of the code is just for esthetics :
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     vjust = 1,
                                                     size = 8, 
                                                     hjust = 1)) +
  Seurat::NoLegend()


# Top 10 markers for each cluster
top3 <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# Visualize marker gene expression for each cluster
DotPlot(merged_seurat, features = unique(top3$gene)) + RotatedAxis()


hpca <- HumanPrimaryCellAtlasData()














# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
merged_seurat[["sketch"]] <- split(merged_seurat[["sketch"]], f = merged_seurat$Tumor_Stage)

merged_seurat <- ProjectIntegration(merged_seurat, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap",
                  reduction.key = "UMAP")

p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters", alpha = 0.1)
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "Tumor_Stage", alpha = 0.1)
p1
p2

merged_seurat[["sketch"]] <- JoinLayers(merged_seurat[["sketch"]])




