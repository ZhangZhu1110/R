# BiocManager::install('Seurat')
# BiocManager::install('tidyverse')
# BiocManager::install('hdf5r')
library(Seurat)
library(tidyverse)
library(hdf5r)
setwd("C:/Users/zhang/Desktop/R")


# Data source: https://www.10xgenomics.com/resources/datasets/40-k-mixture-of-nsclc-dt-cs-from-7-donors-3-ht-v-3-1-3-1-high-6-1-0

# Load the NSCLC dataset
nsclc.sparse.m = Read10X_h5(filename = '40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
cts = nsclc.sparse.m$`Gene Expression`

# Initialize the Seurat object with the raw (non-normalized data)
# Reserve features expressed in at least 3 cells
# Reserve cells expressing at least 200 features
nsclc.seurat.obj = CreateSeuratObject(counts = cts, project = 'NSCLC', min.cells = 3, min.features = 200)

# 1.QC

# Define mitochondrial gene and ribosome protein
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
nsclc.seurat.obj[["percent.rb"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^RP[SL]")
remove(nsclc.sparse.m)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering
nsclc.seurat.obj = subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & percent.rb < 10)

# 3. Normalization
nsclc.seurat.obj = NormalizeData(nsclc.seurat.obj)

# 4. Identify highly variable features
nsclc.seurat.obj = FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 variable features
top10 = head(VariableFeatures(nsclc.seurat.obj), 10)

plot1 = VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# 5. Z-score normalization
all.genes = rownames(nsclc.seurat.obj)
nsclc.seurat.obj = ScaleData(nsclc.seurat.obj, features = all.genes)

# 6. Perform Linear dimension reduction
nsclc.seurat.obj = RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# Visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(nsclc.seurat.obj, reduction = "pca")
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500,balanced = TRUE)

# Find how many PCs can be used without much information loss
ElbowPlot(nsclc.seurat.obj)

# 7. Clustering
nsclc.seurat.obj = FindNeighbors(nsclc.seurat.obj, dims = 1:10)

# Understanding resolution
# Higher resolution leads to more clusters
nsclc.seurat.obj = FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# 8. Non-linear dimension reduction - UMAP and t-SNE
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:10)
DimPlot(nsclc.seurat.obj, reduction = "umap", label = T)

nsclc.seurat.obj <- RunTSNE(nsclc.seurat.obj, dims = 1:10)
DimPlot(nsclc.seurat.obj, reduction = "tsne", label = T)

# 9.Feature plot
FeaturePlot(nsclc.seurat.obj, features = c("IGHG2", "IGLV2-14", "HBA2", "IGKV4-1"))

FeaturePlot(nsclc.seurat.obj, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))
FeaturePlot(nsclc.seurat.obj, features = "nCount_RNA") & theme(plot.title = element_text(size=10))

VlnPlot(nsclc.seurat.obj,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))

# 10.Find markers
DefaultAssay(nsclc.seurat.obj) <- "RNA"
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
all.markers <- FindAllMarkers(nsclc.seurat.obj, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
head(all.markers)
table(all.markers$cluster)

# Fine top3 markers for each cluster
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers
