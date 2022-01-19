library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reticulate)
set.seed(42)

# c1.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/SRR9123032/outs/raw_feature_bc_matrix")
# c2.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C2/outs/raw_feature_bc_matrix")
# c3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C3/outs/raw_feature_bc_matrix")
# c4.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C4/outs/raw_feature_bc_matrix")
# c5.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C5/outs/raw_feature_bc_matrix")
# c6.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C6/outs/raw_feature_bc_matrix")
# c7.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C7/outs/raw_feature_bc_matrix")
# c8.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C8/outs/raw_feature_bc_matrix")
# c9.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/C9/outs/raw_feature_bc_matrix")
# 
# p1.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS1/outs/raw_feature_bc_matrix")
# p2.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS2/outs/raw_feature_bc_matrix")
# p3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS3/outs/raw_feature_bc_matrix")
# p4.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS4/outs/raw_feature_bc_matrix")
# p5.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS5/outs/raw_feature_bc_matrix")
# p6.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS6/outs/raw_feature_bc_matrix")
# p7.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS7/outs/raw_feature_bc_matrix")
# p8.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS8/outs/raw_feature_bc_matrix")
# #p9.data <- Read10X(data.dir ="~/Documents/Thesis/IMPx/IMPx5/Matrices/MS/Lesion/P9/")
# p10.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS10/outs/raw_feature_bc_matrix")
# p11.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS11/outs/raw_feature_bc_matrix")
# p12.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/MS/MS12/outs/raw_feature_bc_matrix")
# 
# #Set up first object
# c1 <- CreateSeuratObject(counts = c1.data, project = "CONTROL", min.cells = 3)
# c1[["percent.mt"]] <- PercentageFeatureSet(c1, pattern = "^MT-")
# c1 <- subset(c1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c1$stim <- "Control"
# 
# #Set up first object
# c2 <- CreateSeuratObject(counts = c2.data, project = "CONTROL", min.cells = 3)
# c2[["percent.mt"]] <- PercentageFeatureSet(c2, pattern = "^MT-")
# c2 <- subset(c2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c2$stim <- "Control"
# 
# #Set up first object
# c3 <- CreateSeuratObject(counts = c3.data, project = "CONTROL", min.cells = 3)
# c3[["percent.mt"]] <- PercentageFeatureSet(c3, pattern = "^MT-")
# c3 <- subset(c3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c3$stim <- "Control"
# 
# #Set up first object
# c4 <- CreateSeuratObject(counts = c4.data, project = "CONTROL", min.cells = 3)
# c4[["percent.mt"]] <- PercentageFeatureSet(c4, pattern = "^MT-")
# c4 <- subset(c4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c4$stim <- "Control"
# 
# #Set up first object
# c5 <- CreateSeuratObject(counts = c5.data, project = "CONTROL", min.cells = 3)
# c5[["percent.mt"]] <- PercentageFeatureSet(c5, pattern = "^MT-")
# c5 <- subset(c5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c5$stim <- "Control"
# 
# #Set up first object
# c6 <- CreateSeuratObject(counts = c6.data, project = "CONTROL", min.cells = 3)
# c6[["percent.mt"]] <- PercentageFeatureSet(c6, pattern = "^MT-")
# c6 <- subset(c6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c6$stim <- "Control"
# 
# #Set up first object
# c7 <- CreateSeuratObject(counts = c7.data, project = "CONTROL", min.cells = 3)
# c7[["percent.mt"]] <- PercentageFeatureSet(c7, pattern = "^MT-")
# c7 <- subset(c7, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c7$stim <- "Control"
# 
# #Set up first object
# c8 <- CreateSeuratObject(counts = c8.data, project = "CONTROL", min.cells = 3)
# c8[["percent.mt"]] <- PercentageFeatureSet(c8, pattern = "^MT-")
# c8 <- subset(c8, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c8$stim <- "Control"
# 
# #Set up first object
# c9 <- CreateSeuratObject(counts = c9.data, project = "CONTROL", min.cells = 3)
# c9[["percent.mt"]] <- PercentageFeatureSet(c9, pattern = "^MT-")
# c9 <- subset(c9, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# c9$stim <- "Control"
# 
# #Set up first object
# p1 <- CreateSeuratObject(counts = p1.data, project = "Patient", min.cells = 3)
# p1[["percent.mt"]] <- PercentageFeatureSet(p1, pattern = "^MT-")
# p1 <- subset(p1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p1$stim <- "Patient"
# 
# #Set up first object
# p2 <- CreateSeuratObject(counts = p2.data, project = "Patient", min.cells = 3)
# p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")
# p2 <- subset(p2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p2$stim <- "Patient"
# 
# #Set up first object
# p3 <- CreateSeuratObject(counts = p3.data, project = "Patient", min.cells = 3)
# p3[["percent.mt"]] <- PercentageFeatureSet(p3, pattern = "^MT-")
# p3 <- subset(p3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p3$stim <- "Patient"
# 
# #Set up first object
# p4 <- CreateSeuratObject(counts = p4.data, project = "Patient", min.cells = 3)
# p4[["percent.mt"]] <- PercentageFeatureSet(p4, pattern = "^MT-")
# p4 <- subset(p4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p4$stim <- "Patient"
# 
# #Set up first object
# p5 <- CreateSeuratObject(counts = p5.data, project = "Patient", min.cells = 3)
# p5[["percent.mt"]] <- PercentageFeatureSet(p5, pattern = "^MT-")
# p5 <- subset(p5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p5$stim <- "Patient"
# 
# #Set up first object
# p6 <- CreateSeuratObject(counts = p6.data, project = "Patient", min.cells = 3)
# p6[["percent.mt"]] <- PercentageFeatureSet(p6, pattern = "^MT-")
# p6 <- subset(p6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p6$stim <- "Patient"
# 
# #Set up first object
# p7 <- CreateSeuratObject(counts = p7.data, project = "Patient", min.cells = 3)
# p7[["percent.mt"]] <- PercentageFeatureSet(p7, pattern = "^MT-")
# p7 <- subset(p7, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p7$stim <- "Patient"
# 
# #Set up first object
# p8 <- CreateSeuratObject(counts = p8.data, project = "Patient", min.cells = 3)
# p8[["percent.mt"]] <- PercentageFeatureSet(p8, pattern = "^MT-")
# p8 <- subset(p8, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p8$stim <- "Patient"
# 
# # #Set up first object
# # p9 <- CreateSeuratObject(counts = p9.data, project = "Patient", min.cells = 3)
# # p9[["percent.mt"]] <- PercentageFeatureSet(p9, pattern = "^MT-")
# # p9 <- subset(p9, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# # p9$stim <- "Patient"
# 
# #Set up first object
# p10 <- CreateSeuratObject(counts = p10.data, project = "Patient", min.cells = 3)
# p10[["percent.mt"]] <- PercentageFeatureSet(p10, pattern = "^MT-")
# p10 <- subset(p10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p10$stim <- "Patient"
# 
# #Set up first object
# p11 <- CreateSeuratObject(counts = p11.data, project = "Patient", min.cells = 3)
# p11[["percent.mt"]] <- PercentageFeatureSet(p11, pattern = "^MT-")
# p11 <- subset(p11, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p11$stim <- "Patient"
# 
# #Set up first object
# p12 <- CreateSeuratObject(counts = p12.data, project = "Patient", min.cells = 3)
# p12[["percent.mt"]] <- PercentageFeatureSet(p12, pattern = "^MT-")
# p12 <- subset(p12, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
# p12$stim <- "Patient"
# 
# anchor.set <- FindIntegrationAnchors(object.list = list(c1, c2, c3, c4, c5, c6, c7, c8, c9,
#                                                         p1, p2, p3, p4, p5, p6, p7, p8, p10,
#                                                         p11, p12), dims = 1:20)
# #p9
# ms.combined <- IntegrateData(anchorset= anchor.set, dims = 1:20)
# 
# 
# DefaultAssay(ms.combined) <- "integrated"
# 
# VlnPlot(ms.combined, c("PLP1"))
# 
# ms.combined <- SCTransform(ms.combined, vars.to.regress = "percent.mt", verbose = TRUE)
# 
# ms.combined <- RunPCA(ms.combined, npcs = 30, verbose = TRUE)
# 
# #Determine appropriate number of PCA
# print(ms.combined[["pca"]], dims = 1:9, nfeatures = 5)
# ElbowPlot(ms.combined)
# 
# # t-SNE and Clustering
# ms.combined <- RunUMAP(ms.combined, reduction = "pca", dims = 1:14)
# ms.combined <- FindNeighbors(ms.combined, reduction = "pca", dims = 1:14)
# 
# #Change resolution to modify number or clusters
# ms.combined <- FindClusters(ms.combined, resolution = 0.6)
# 
# #Find markers for every cluster compared to all remaining cells, report only the positive ones
# ms.combined.marker <- FindAllMarkers(ms.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ms.combined.marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# 
# # Visualization
# pan1 <- DimPlot(ms.combined, reduction = "umap")
# pan2 <- DimPlot(ms.combined, reduction = "umap", label = TRUE)
# plot_grid(pan1, pan2)
# DimPlot(ms.combined, reduction = "umap", split.by = "stim")
load("~/Ms.RData")
#save.image("~/Ms.RData")
# my.data=fetch.data(ms.combined,c("ident","PC1","nGene","orig.ident","PLP1"))
# write.csv(my.data, file="./plp_levels.csv")
Idents(ms.combined, WhichCells(object = ms.combined, expression = PLP1 > 1.5, slot = 'data')) <- 'plp.pos'
Idents(ms.combined, WhichCells(object = ms.combined, expression = PLP1 <= 1.5, slot = 'data')) <- 'plp.neg'
Idents(ms.combined) <- paste0(Idents(ms.combined), ms.combined$stim)
plp.markers <- FindMarkers(ms.combined, ident.1 = 'plp.posPatient', ident.2 = 'plp.posControl')
write.csv(plp.markers, file="./plp_DEG.csv")
#ms.combined[["plp"]] <- PercentageFeatureSet(ms.combined, pattern = "^PLP-")
VlnPlot(ms.combined, c("PLP1"))
ms.combined <- subset(ms.combined, subset = PLP1 > 1.5)
scd.counts <- ms.combined[["integrated"]]["SCD"]
write.csv(scd.counts, file="./scd_counts.csv")

scd5.counts <- ms.combined[["integrated"]]["SCD5"]
write.csv(scd5.counts, file="./scd5_counts.csv")

# Visualization
pan1 <- DimPlot(ms.combined, reduction = "umap")
pan2 <- DimPlot(ms.combined, reduction = "umap", label = TRUE)
plot_grid(pan1, pan2)
DimPlot(ms.combined, reduction = "umap", split.by = "stim")



# orig.levels <- levels(ms.combined)
# Idents(ms.combined) <- gsub(pattern = " ", replacement = "_", x = Idents(ms.combined))
# orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
# levels(ms.combined) <- orig.levels
# cluster.averages <- AverageExpression(ms.combined, return.seurat = TRUE)
# VlnPlot(cluster.averages, features = c("SCD", "SCD5", "QKI"))
# RidgePlot(cluster.averages, features = c("SCD", "SCD5", "QKI"))


DefaultAssay(ms.combined) <- "RNA"
FeaturePlot(ms.combined, features = c("SCD", "SCD5", "QKI"), split.by = "stim")
VlnPlot(ms.combined, features = c("SCD", "SCD5", "QKI"), pt.size = 0)
RidgePlot(ms.combined, features = c("SCD", "SCD5", "QKI"), split.by = "stim")



#save png
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="/rsrch3/home/genomic_med/dzamler/")










