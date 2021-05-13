library(Seurat)
library(dplyr)
library(cowplot)
library(reticulate)
library(ggplot2)
set.seed(42)

gbm.data <- Read10X(data.dir ="~/Documents/Thesis/IMPx/IMPx5/Matrices/HuFre_Frz/")
pbmc.data <- Read10X(data.dir ="~/Documents/Rezvani/filtered_gene_bc_matrices/hg19/")
immune.gbm <- agg
#Set up GBM object
immune.gbm <- CreateSeuratObject(counts = gbm.data, project = "NKCompare", min.cells = 3)
immune.gbm[["percent.mt"]] <- PercentageFeatureSet(immune.gbm, pattern = "^MT-")
immune.gbm <- subset(immune.gbm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
immune.gbm$stim <- "GBM"

#Set up PBMC object
immune.pbmc <- CreateSeuratObject(counts = pbmc.data, project = "NKCompare", min.cells = 3)
immune.pbmc[["percent.mt"]] <- PercentageFeatureSet(immune.pbmc, pattern = "^MT-")
immune.pbmc <- subset(immune.pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
immune.pbmc$stim <- "PBMC"

immune.gbm <- SCTransform(immune.gbm, vars.to.regress = "percent.mt", verbose = TRUE)

immune.gbm <- RunPCA(immune.gbm, npcs = 30, verbose = TRUE)

#Determine appropriate number of PCA
print(immune.gbm[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(immune.gbm)

# t-SNE and Clustering
immune.gbm <- RunUMAP(immune.gbm, reduction = "pca", dims = 1:15)
immune.gbm <- FindNeighbors(immune.gbm, reduction = "pca", dims = 1:15)

#Change resolution to modify number or clusters
immune.gbm <- FindClusters(immune.gbm, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
gimmune.markers <- FindAllMarkers(immune.gbm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gimmune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



# Visualization
p1 <- DimPlot(immune.gbm, reduction = "umap", split.by = "stim")
p2 <- DimPlot(immune.gbm, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(immune.gbm, reduction = "umap", split.by = "stim")

#High resolution figs
FeaturePlot(immune.gbm, features = c("KLRD1", "NKG7", "NKTR"))
VlnPlot(immune.gbm, features = c("KLRD1", "NKG7", "NKTR"))
nk.gbm <- subset(immune.gbm, subset =  seurat_clusters == 2)



# What proportion of cells are in each cluster?
counts.cluster <- prop.table(table(Idents(immune.gbm)))
# How does cluster membership vary by replicate?
counts.cluster.sample <- table(Idents(immune.gbm), immune.gbm$stim)


immune.pbmc <- SCTransform(immune.pbmc, vars.to.regress = "percent.mt", verbose = TRUE)

immune.pbmc <- RunPCA(immune.pbmc, npcs = 30, verbose = TRUE)

#Determine appropriate number of PCA
print(immune.pbmc[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(immune.pbmc)

# t-SNE and Clustering
immune.pbmc <- RunUMAP(immune.pbmc, reduction = "pca", dims = 1:15)
immune.pbmc <- FindNeighbors(immune.pbmc, reduction = "pca", dims = 1:15)

#Change resolution to modify number or clusters
immune.pbmc <- FindClusters(immune.pbmc, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
immune.markers <- FindAllMarkers(immune.pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/FreshPatients0.6.csv")


# Visualization
p1 <- DimPlot(immune.pbmc, reduction = "umap", split.by = "stim")
p2 <- DimPlot(immune.pbmc, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(immune.pbmc, reduction = "umap", split.by = "stim")

#High resolution figs
FeaturePlot(immune.pbmc, features = c("KLRD1", "NKG7", "NKTR"))
VlnPlot(immune.pbmc, features = c("KLRD1", "NKG7", "NKTR"))
nk.pbmc <- subset(immune.pbmc, subset =  seurat_clusters == 2)



# What proportion of cells are in each cluster?
counts.cluster <- prop.table(table(Idents(immune.pbmc)))
# How does cluster membership vary by replicate?
counts.cluster.sample <- table(Idents(immune.pbmc), immune.pbmc$stim)

nk.anchors <- FindIntegrationAnchors(object.list = list(nk.gbm, nk.pbmc), dims = 1:30)
nk.combined <- IntegrateData(anchorset = nk.anchors, dims = 1:30)

DefaultAssay(nk.combined) <- "SCT"
nk.combined <- SCTransform(nk.combined, vars.to.regress = "percent.mt", verbose = TRUE)
nk.combined <- RunPCA(nk.combined, npcs = 30, verbose = TRUE)

#Determine appropriate number of PCA
print(nk.combined[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(nk.combined)

# t-SNE and Clustering
nk.combined <- RunUMAP(nk.combined, reduction = "pca", dims = 1:15)
nk.combined <- FindNeighbors(nk.combined, reduction = "pca", dims = 1:15)

#Change resolution to modify number or clusters
nk.combined <- FindClusters(nk.combined, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
nimmune.markers <- FindAllMarkers(nk.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nimmune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/FreshPatients0.6.csv")


# Visualization
p1 <- DimPlot(nk.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(nk.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(nk.combined, reduction = "umap", split.by = "stim")

DefaultAssay(nk.combined) <- "RNA"
#Activation
FeaturePlot(nk.combined, features = c("CD226", "KLRK1", "TNFSF10"), split.by = "stim")
VlnPlot(nk.combined, features = c("CD226", "KLRK1", "TNFSF10"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CD69", "NCR3", "NCR2"))
VlnPlot(nk.combined, features = c("CD69", "NCR3", "NCR2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("NCR1", "CD244", "FAS"))
VlnPlot(nk.combined, features = c("NCR1", "CD244", "FAS"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("FASLG", "GZMB", "GZMA"))
VlnPlot(nk.combined, features = c("FASLG", "GZMB", "GZMA"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GZMH", "GZMK", "TBX21"))
VlnPlot(nk.combined, features = c("GZMH", "GZMK", "TBX21"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("PRF1", "MKI67", "KIT"))
VlnPlot(nk.combined, features = c("PRF1", "MKI67", "KIT"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KLRC2", "SELL", "IL2RA"))
VlnPlot(nk.combined, features = c("KLRC2", "SELL", "IL2RA"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("FCGR3A", "EOMES", "ZAP70"))
VlnPlot(nk.combined, features = c("FCGR3A", "EOMES", "ZAP70"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SYK", "KLRB1", "CD247"))
VlnPlot(nk.combined, features = c("SYK", "KLRB1", "CD247"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CD27", "LAMP1", "IFNG"))
VlnPlot(nk.combined, features = c("CD27", "LAMP1", "IFNG"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("TNF", "LTA", "TGFB1"))
VlnPlot(nk.combined, features = c("TNF", "LTA", "TGFB1"), group.by = "stim", pt.size = 0)


#Inhibition
FeaturePlot(nk.combined, features = c("KLRG1", "KLRC1", "LAG3"))
VlnPlot(nk.combined, features = c("KLRG1", "KLRC1", "LAG3"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("HAVCR2", "PDCD1", "TIGIT"))
VlnPlot(nk.combined, features = c("HAVCR2", "PDCD1", "TIGIT"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CD96", "B3GAT1", "CD8A"))
VlnPlot(nk.combined, features = c("CD96", "B3GAT1", "CD8A"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SIGLEC7", "KLRD1", "KIR2DL4"))
VlnPlot(nk.combined, features = c("SIGLEC7", "KLRD1", "KIR2DL4"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KIR2DL3", "KIR3DL2", "KIR2DL1"))
VlnPlot(nk.combined, features = c("KIR2DL3", "KIR3DL2", "KIR2DL1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KIR2DL5A", "KIR2DL2", "LILRB1"))
VlnPlot(nk.combined, features = c("KIR2DL5A", "KIR2DL2", "LILRB1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("IL1RAPL1", "CD9", "ITGAE"))
VlnPlot(nk.combined, features = c("IL1RAPL1", "CD9", "ITGAE"), group.by = "stim", pt.size = 0)

#TGF-B pathway
FeaturePlot(nk.combined, features = c("ACVR1", "ACVR1B", "ACVR1C"))
VlnPlot(nk.combined, features = c("ACVR1", "ACVR1B", "ACVR1C"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("ACVR2A", "ACVR2B", "ACVRL1"))
VlnPlot(nk.combined, features = c("ACVR2A", "ACVR2B", "ACVRL1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("AMHR2", "ATF2", "BAMBI"))
VlnPlot(nk.combined, features = c("AMHR2", "ATF2", "BAMBI"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("BMP1", "BMP10", "BMP15"))
VlnPlot(nk.combined, features = c("BMP1", "BMP10", "BMP15"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("BMP2", "BMP3", "BMP4"))
VlnPlot(nk.combined, features = c("BMP2", "BMP3", "BMP4"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("BMP5", "BMP6", "BMP7"))
VlnPlot(nk.combined, features = c("BMP5", "BMP6", "BMP7"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("BMP8A", "BMP8B", "BMPR1A"))
VlnPlot(nk.combined, features = c("BMP8A", "BMP8B", "BMPR1A"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("BMPR1B", "BMPR2", "CITED1"))
VlnPlot(nk.combined, features = c("BMPR1B", "BMPR2", "CITED1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CITED2", "CREBBP", "DCP1B"))
VlnPlot(nk.combined, features = c("CITED2", "CREBBP", "DCP1B"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("EP300", "FKBP1A", "FOSL1"))
VlnPlot(nk.combined, features = c("EP300", "FKBP1A", "FOSL1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("FOXH1", "GDF1", "GDF10"))
VlnPlot(nk.combined, features = c("FOXH1", "GDF1", "GDF10"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GDF11", "GDF15", "GDF2"))
VlnPlot(nk.combined, features = c("GDF11", "GDF15", "GDF2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GDF3", "GDF5", "GDF6"))
VlnPlot(nk.combined, features = c("GDF3", "GDF5", "GDF6"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GDF7", "GDF9", "GDNF"))
VlnPlot(nk.combined, features = c("GDF7", "GDF9", "GDNF"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("HRAS", "INHBA", "INHBB"))
VlnPlot(nk.combined, features = c("HRAS", "INHBA", "INHBB"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("INHBC", "INHBE", "JUN"))
VlnPlot(nk.combined, features = c("INHBC", "INHBE", "JUN"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("JUNB", "JUND", "LEFTY1"))
VlnPlot(nk.combined, features = c("JUNB", "JUND", "LEFTY1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("LEFTY2", "MAP3K7", "MAP3K7CL"))
VlnPlot(nk.combined, features = c("LEFTY2", "MAP3K7", "MAP3K7CL"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("MAPK1", "MAPK10", "MAPK11"))
VlnPlot(nk.combined, features = c("MAPK1", "MAPK10", "MAPK11"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("MAPK12", "MAPK13", "MAPK14"))
VlnPlot(nk.combined, features = c("MAPK12", "MAPK13", "MAPK14"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("MAPK3", "MAPK8", "MAPK9"))
VlnPlot(nk.combined, features = c("MAPK3", "MAPK8", "MAPK9"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("MSTN", "NODAL", "NRAS"))
VlnPlot(nk.combined, features = c("MSTN", "NODAL", "NRAS"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("RRAS", "SKI", "SKIL"))
VlnPlot(nk.combined, features = c("RRAS", "SKI", "SKIL"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD1", "SMAD2", "SMAD3"))
VlnPlot(nk.combined, features = c("SMAD1", "SMAD2", "SMAD3"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD4", "SMAD6", "SMAD7"))
VlnPlot(nk.combined, features = c("SMAD4", "SMAD6", "SMAD7"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD9", "SMURF1", "SMURF2"))
VlnPlot(nk.combined, features = c("SMAD9", "SMURF1", "SMURF2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SNIP1", "TAB1", "TGFB2"))
VlnPlot(nk.combined, features = c("SNIP1", "TAB1", "TGFB2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("TGFB3", "TGFBR1", "TGFBR2"))
VlnPlot(nk.combined, features = c("TGFB3", "TGFBR1", "TGFBR2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("TLL", "TLL2", "ZFYVE9"))
VlnPlot(nk.combined, features = c("TLL", "TLL2", "ZFYVE9"), group.by = "stim", pt.size = 0)


FeaturePlot(nk.combined, features = c("NT5E"), split.by = "stim")


#save png
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="~/Documents/Rezvani/TGF_Pathway")

nkg2d.counts <- nk.combined[["SCT"]]["NKG2D"]
write.csv(immune.markers, file="./Top_pbmc_genes.csv")

jund.counts <- nk.combined[["SCT"]]["JUND"]
p.dat <- rbind(ncr3.counts, jund.counts)
p.data <- rbind(nk.combined[["SCT"]]["CD226"], nk.combined[["SCT"]]["CD69"], nk.combined[["SCT"]]["NCR3"], nk.combined[["SCT"]]["GZMB"], nk.combined[["SCT"]]["GZMA"], nk.combined[["SCT"]]["GZMH"],  nk.combined[["SCT"]]["GZMK"], nk.combined[["SCT"]]["PRF1"],  nk.combined[["SCT"]]["SELL"], nk.combined[["SCT"]]["KLRC2"],  nk.combined[["SCT"]]["FCGR3A"], nk.combined[["SCT"]]["CD247"],  nk.combined[["SCT"]]["KLRD1"], nk.combined[["SCT"]]["SIGLEC7"],  nk.combined[["SCT"]]["KIR2DL4"], nk.combined[["SCT"]]["KIR2DL1"], nk.combined[["SCT"]]["JUND"],  nk.combined[["SCT"]]["SMAD3"], nk.combined[["SCT"]]["SMAD4"], nk.combined[["SCT"]]["SMAD7"], nk.combined[["SCT"]]["SMURF2"] )
write.csv(p.data, file="./p_vals.csv")

 HeatGenes <- c("CD226", "KLRK1", "TNFSF10","CD69", "NCR3", "NCR2", "NCR1", 
               "CD244", "FAS", "FASLG", "GZMB", "GZMA", "GZMH", "GZMK", "TBX21",
               "PRF1", "MKI67", "KIT", "KLRC2", "SELL", "IL2RA", "FCGR3A", "EOMES",
               "ZAP70", "SYK", "KLRB1", "CD247", "CD27", "LAMP1", "IFNG",
               "TNF", "LTA", "TGFB1")
DoHeatmap(nk.combined, features = HeatGenes, group.by = "stim")
HeatGenes[HeatGenes %in% rownames(nk.combined@scale.data)]


DimHeatmap(nk.combined, dims = 1:3)

FeaturePlot(nk.combined, features = c("CD226"))
VlnPlot(nk.combined, features = c("CD226"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CD69"))
VlnPlot(nk.combined, features = c("CD69"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("NCR3"))
VlnPlot(nk.combined, features = c("NCR3"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GZMB"))
VlnPlot(nk.combined, features = c("GZMB"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GZMA"))
VlnPlot(nk.combined, features = c("GZMA"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GZMH"))
VlnPlot(nk.combined, features = c("GZMH"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("GZMK"))
VlnPlot(nk.combined, features = c("GZMK"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("PRF1"))
VlnPlot(nk.combined, features = c("PRF1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SELL"))
VlnPlot(nk.combined, features = c("SELL"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KLRC2"))
VlnPlot(nk.combined, features = c("KLRC2"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("FCGR3A"))
VlnPlot(nk.combined, features = c("FCGR3A"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("CD247"))
VlnPlot(nk.combined, features = c("CD247"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KLRD1"))
VlnPlot(nk.combined, features = c("KLRD1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SIGLEC7"))
VlnPlot(nk.combined, features = c("SIGLEC7"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KIR2DL4"))
VlnPlot(nk.combined, features = c("KIR2DL4"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("KIR2DL1"))
VlnPlot(nk.combined, features = c("KIR2DL1"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("JUND"))
VlnPlot(nk.combined, features = c("JUND"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD3"))
VlnPlot(nk.combined, features = c("SMAD3"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD4"))
VlnPlot(nk.combined, features = c("SMAD4"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMAD7"))
VlnPlot(nk.combined, features = c("SMAD7"), group.by = "stim", pt.size = 0)
FeaturePlot(nk.combined, features = c("SMURF2"))
VlnPlot(nk.combined, features = c("SMURF2"), group.by = "stim", pt.size = 0)
