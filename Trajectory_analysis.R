# install monocle3
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils',
#                        'HDF5Array', 'terra', 'ggrastr'))
# install.packages("devtools")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
# devtools::install_github('cole-trapnell-lab/monocle3')

# install SeuratWrappers
# install.packages("remotes")
# remotes::install_github("satijalab/seurat-wrappers")

# data source: http://221.239.103.203:3838/abc/
# Instruction source: https://www.youtube.com/watch?v=iq4T_uzMFcY

setwd("C:/Users/zhang/Desktop/R")

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(magrittr)

# For trajectory analysis, in addition to expression profile, the known annotation for each cell is also required.
# We have to make sure that trajectory is exist in our data.

# 1. Import data
markers = read.delim('ABC_Marker.txt', header = T)
cell = read.delim('ABC_Meta.txt', header = T)
expression = read.csv('ABC_umi_matrix_7551_cells.csv', header = T)

# 2. Create seurat object and filter the useful information
# We require gene as rows and cell as columns
expression_t = t(expression)
remove(expression)
seurat_obj = CreateSeuratObject(counts = expression_t)
remove(expression_t)
# Merge cell information to the metadata
seurat_obj@meta.data = merge(seurat_obj@meta.data, cell, by.x = 'row.names', by.y = 'cell_id')
# Set row names
seurat_obj@meta.data = seurat_obj@meta.data %>% column_to_rownames(var = 'Row.names')
#Calculate percentage of mitochondrial DNA and ribosome
seurat_obj$mito = PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# Filter data with low quality
seurat_obj = subset(seurat_obj, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mito < 10)


unique(seurat_obj@meta.data$population)
# There are 7 types of cells in the data set
#  "sp" "t"  "mo" "nk" "e"  "b"  "n" 
# I'm interested in monocyte, "mo"
Idents(seurat_obj) = seurat_obj@meta.data$population
mo_seurat = subset(seurat_obj, idents = 'mo')
unique(mo_seurat@meta.data$redefined_cluster)
# There are 6 sub types of monocyte
# Classical monocyte, Intermediate monocyte, hMDP/cMoP, Pre-monocyte, hMDP, Non-classical monocyte

# 3. Pre-processing by using Seurat
mo_seurat = NormalizeData(mo_seurat)
mo_seurat = FindVariableFeatures(mo_seurat)
mo_seurat = ScaleData(mo_seurat)
mo_seurat = RunPCA(mo_seurat)
mo_seurat = FindNeighbors(mo_seurat, dims = 1:30)
mo_seurat = FindClusters(mo_seurat, resolution = 0.9)
# This step will add a column named seurat_clusters in metedata
mo_seurat = RunUMAP(mo_seurat, dims = 1:30, n.neighbors = 50)

# Create plot by seurat
p1 = DimPlot(mo_seurat, reduction = 'umap', group.by = 'redefined_cluster', label = T)
p2 = DimPlot(mo_seurat, reduction = 'umap', group.by = 'seurat_clusters', label = T)
p1 | p2

# 4. Prepare data for trajectory analysis
# Convert seurat object to cell_data_set object
monocyte = as.cell_data_set(mo_seurat)

colData(monocyte)
# It contains the metadata of seurat object

fData(monocyte)
# It's just a list of feature(gene) names as row names
# Add gene name as content
fData(monocyte)$gene_short_name = rownames(fData(monocyte))

# Assign cell name to partition
names(monocyte@clusters$UMAP$partitions) = monocyte@colData@rownames


# Create plot by monocle
# Since the object is derived from seurat, the plots are the same as before.
cluster_before_trajectory = plot_cells(monocyte,
                                       color_cells_by = 'cluster',
                                       label_groups_by_cluster = F,
                                       group_label_size = 5) + 
  theme(legend.position = 'right')
cluster.names = plot_cells(monocyte,
                           color_cells_by = 'redefined_cluster',
                           label_groups_by_cluster = F,
                           group_label_size = 5) + scale_color_manual(values = c('red', 'blue', 'green','maroon', 'yellow', 'green', 'cyan')) +
  theme(legend.position = 'right')

cluster_before_trajectory | cluster.names

# 5. Learn trajectory graph
monocyte = learn_graph(monocyte, use_partition = F)
plot_cells(monocyte,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)

# 6. Order the cells in pseudotime
# Based on cluster result, we can see that hMDP belong to cluster 7
monocyte = order_cells(monocyte, reduction_method = 'UMAP', root_cells = colnames(monocyte[, clusters(monocyte) == 7]))
plot_cells(monocyte,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)

# Create a data frame for ggplot2
monocyte$pseudotime = pseudotime(monocyte)
data_pseudotime = as.data.frame(colData(monocyte))
ggplot(data_pseudotime, aes(pseudotime, reorder(redefined_cluster, pseudotime, median), fill = redefined_cluster)) +
  geom_boxplot()

# 7. Find genes change as a function of pseudotime
deg = graph_test(monocyte, neighbor_graph = 'principal_graph', cores = 4)
deg %>% arrange(q_value) %>% filter(status == 'OK') %>% head()
# 'TNFRSF1B', 'STMN1', 'SH3BGRL3', 'RPS8' are top4 markers
# Create feature plot for markers
FeaturePlot(mo_seurat, features = c('TNFRSF1B', 'STMN1', 'SH3BGRL3', 'RPS8'))

# 8. Visualize pseudotime in seurat
mo_seurat$pseudotime = pseudotime(monocyte)
Idents(mo_seurat) = mo_seurat$redefined_cluster
FeaturePlot(mo_seurat, features = 'pseudotime', label = T)
