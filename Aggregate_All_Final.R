install.packages('dplyr')
install.packages('Seurat')
library(Seurat)
library(dplyr)
library(cowplot)
library(reticulate)
library(ggplot2)
set.seed(42)
#######################################################
####REVIEWER TOKEN TO ACCESS DATA: opejqqwmtvwxhwb ####
#######################################################
impagg.data <- Read10X(data.dir ="~/Documents/Thesis/IMPx/IMPx5/Matrices/ImpAgg/")
spoagg.data <- Read10X(data.dir ="~/Documents/Thesis/IMPx/IMPx5/Matrices/SpoAgg/")
agg.data <- Read10X(data.dir ="~/Documents/Thesis/IMPx/IMPx5/Matrices/HuFFAgg/")

#Set up first object
#Set up first object
agg <- CreateSeuratObject(counts = agg.data, project = "AGGREGATE", min.cells = 3)
agg[["percent.mt"]] <- PercentageFeatureSet(agg, pattern = "^MT-")
agg <- subset(agg, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
gemgroup <- sapply(strsplit(rownames(agg@meta.data), split="-"), "[[", 2) 
map = setNames(c("10", "9", "8", "7", "2", "5", "3", "1", "6", "5", "4", "3", "2" ), c("13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"))
gemgroup[] <- map[unlist(gemgroup)]
agg <- AddMetaData(object=agg, metadata=data.frame(gemgroup=gemgroup, row.names=rownames(agg@meta.data)))
agg <- SCTransform(agg, vars.to.regress = "percent.mt", verbose = TRUE)
VlnPlot(agg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(agg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(agg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Set up second object
imp <- CreateSeuratObject(counts = impagg.data, project = "IMPLANT_CELL", min.cells = 3)
imp[["percent.mt"]] <- PercentageFeatureSet(imp, pattern = "^mt-")
imp <- subset(imp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
igemgroup <- sapply(strsplit(rownames(imp@meta.data), split="-"), "[[", 2) 
imp <- AddMetaData(object=imp, metadata=data.frame(gemgroup=igemgroup, row.names=rownames(imp@meta.data)))
imp$stim <- "IMP"
imp <- SCTransform(imp, vars.to.regress = "percent.mt", verbose = TRUE)
VlnPlot(imp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(imp, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(imp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Set up third objimp
spo <- CreateSeuratObject(counts = spoagg.data, project = "SPONT_CELL", min.cells = 3)
spo[["percent.mt"]] <- PercentageFeatureSet(spo, pattern = "^mt-")
spo <- subset(spo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
sgemgroup <- sapply(strsplit(rownames(spo@meta.data), split="-"), "[[", 2) 
spo <- AddMetaData(object=spo, metadata=data.frame(gemgroup=sgemgroup, row.names=rownames(spo@meta.data)))
spo$stim <- "SPO"
spo <- SCTransform(spo, vars.to.regress = "percent.mt", verbose = TRUE)
VlnPlot(spo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(spo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(spo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Low resolution Figures

#Spontaneous Mice Alone
spo <- RunPCA(spo, npcs = 30, verbose = FALSE)

#Determine appropriate number of PCA
print(spo[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(spo)
pca.mouse.spo <- spo[["pca"]]
# t-SNE and Clustering
spo <- RunUMAP(spo, reduction = "pca", dims = 1:9)
spo <- FindNeighbors(spo, reduction = "pca", dims = 1:9)

#Change resolution to modify number or clusters
spo <- FindClusters(spo, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
spomune.markers <- FindAllMarkers(spo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
spomune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(spomune.markers, file="../Thesis/IMPx/IMPx5/MoSpo_0.1Deg.csv")


# Visualization
p1 <- DimPlot(spo, reduction = "umap", group.by = "stim")
p2 <- DimPlot(spo, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(spo, reduction = "umap", split.by = "gemgroup")

#Implanted Mice Alone
imp <- RunPCA(imp, npcs = 30, verbose = FALSE)

#Determine appropriate number of PCA
print(imp[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(imp)
pca.mouse.imp <- imp[["pca"]]
# t-SNE and Clustering
imp <- RunUMAP(imp, reduction = "pca", dims = 1:9)
imp <- FindNeighbors(imp, reduction = "pca", dims = 1:9)

#Change resolution to modify number or clusters
imp <- FindClusters(imp, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
impmune.markers <- FindAllMarkers(imp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
impmune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(impmune.markers, file="../Thesis/IMPx/IMPx5/MoImp_0.1Deg.csv")


# Visualization
p1 <- DimPlot(imp, reduction = "umap", group.by = "stim")
p2 <- DimPlot(imp, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(imp, reduction = "umap", split.by = "gemgroup")

#Combined Mice
immune.anchors <- FindIntegrationAnchors(object.list = list(imp, spo), dims = 1:20)

immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)

immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)


#Determine appropriate number of PCA
print(immune.combined[["pca"]], dims = 1:9, nfeatures = 5)
ElbowPlot(immune.combined)
pca.mouse.combined <- immune.combined[["pca"]]
write.csv(pca.mouse.combined, file="../Thesis/IMPx/IMPx5/MoAgg_PCA.csv")

# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:9)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:9)

#Change resolution to modify number or clusters
immune.combined <- FindClusters(immune.combined, resolution = 0.1)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
immune.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/MoAgg_0.1Deg.csv")


# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(immune.combined, reduction = "umap", split.by = "stim")

#Human Alone
##Run this Portion if you want to assign grade subsets for separate analysis
# lgg <- subset(agg, subset = gemgroup == 1 | gemgroup == 7 | gemgroup == 9)
# lgg$grade <- "LGG"
# gbm <- subset(agg, subset = gemgroup == 2 | gemgroup == 3 | gemgroup == 4 | gemgroup == 5 | gemgroup == 6 | gemgroup == 8 | gemgroup == 10)
# gbm$grade <- "GBM"
# 
# immune.anchors <- FindIntegrationAnchors(object.list = list(lgg, gbm), dims = 1:20)
# agg <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
# 
# DefaultAssay(agg) <- "integrated"
agg <- RunPCA(agg, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
agg <- RunUMAP(agg, reduction = "pca", dims = 1:15)
agg <- FindNeighbors(agg, reduction = "pca", dims = 1:15)

#Change resolution to modify number or clusters
agg <- FindClusters(agg, resolution = 0.1)

# Visualization
p1 <- DimPlot(agg, reduction = "umap", group.by = "gemgroup")
p2 <- DimPlot(agg, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(agg, reduction = "umap", label = TRUE )

#Find markers for every cluster compared to all remaining cells, report only the positive ones
agg.markers <- FindAllMarkers(agg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
agg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/HuAgg_0.1Deg.csv")



#Plots
DefaultAssay(imp) <- "RNA"
DefaultAssay(spo) <- "RNA"
DefaultAssay(agg) <- "RNA"
DefaultAssay(immune.combined) <- "RNA"

FeaturePlot(agg, features = c("CD3D", "CD3E", "CD3G"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD3D", "CD3E", "CD3G"))
FeaturePlot(agg, features = c("HCST", "KLRD1", "NKG7"), min.cutoff = "q9")
VlnPlot(agg, features = c("HCST", "KLRD1", "NKG7"))
FeaturePlot(agg, features = c("C1QA", "C1QB", "C1QC"), min.cutoff = "q9")
VlnPlot(agg, features = c("C1QA", "C1QB", "C1QC"))
FeaturePlot(agg, features = c("CD14", "CD68", "GBP2"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "GBP2"))
FeaturePlot(agg, features = c("CD14", "CD68", "GRN"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "GRN"))
FeaturePlot(agg, features = c("CD14", "CD68", "MRC1"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "MRC1"))
FeaturePlot(agg, features = c("CD14", "CD68", "ITGAM"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "ITGAM"))
FeaturePlot(agg, features = c("ITGAX", "CD83", "CD14"), min.cutoff = "q9")
VlnPlot(agg, features = c("ITGAX", "CD83", "CD14"))
FeaturePlot(agg, features = c("CD74", "HLA-DQA1", "HLA-DRB5"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD74", "HLA-DQA1", "HLA-DRB5"))
FeaturePlot(agg, features = c("CD68", "MARCKS", "IL1B"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD68", "MARCKS", "IL1B"))
FeaturePlot(agg, features = c("ITGAX", "S100A8", "S100A9"), min.cutoff = "q9")
VlnPlot(agg, features = c("ITGAX", "S100A8", "S100A9"))
FeaturePlot(agg, features = c("NCF2", "S100A8", "S100A9"), min.cutoff = "q9")
VlnPlot(agg, features = c("NCF2", "S100A8", "S100A9"))
FeaturePlot(agg, features = c("CD14", "LILRB1", "LILRB2"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "LILRB1", "LILRB2"))
FeaturePlot(agg, features = c("CD68", "CX3CR1", "TMEM119"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD68", "CX3CR1", "TMEM119"))
FeaturePlot(agg, features = c("CD79A", "CD79B", "MS4A3"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD79A", "CD79B", "MS4A7"))


FeaturePlot(imp, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(imp, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(imp, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(imp, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(imp, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(imp, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(imp, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)
FeaturePlot(imp, features = c("Ly6a", "Ly6c2", "Ly6d"), min.cutoff = "q9")
VlnPlot(imp, features = c("Ly6a", "Ly6c2", "Ly6d"), pt.size = 0)


FeaturePlot(spo, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(spo, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(spo, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(spo, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(spo, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(spo, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(spo, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)




FeaturePlot(immune.combined, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Ly6a", "Ly6c2", "Ly6d"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Ly6a", "Ly6c2", "Ly6d"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Ms4a3"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("Cd1", "Cd1b", "Cd1c", "Cd1d"))


#Percentage contribution
agg.meta.data <- agg@meta.data

agg.counts <- group_by(agg.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(agg.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Percentage contribution
mo.meta.data <- immune.combined@meta.data

mo.counts <- group_by(mo.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(mo.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Dot plots
markers.to.plot <- c("CD19","CD79A","JCHAIN", "KLRD1", "NKG7", "NKTR",
                     "CD3D", "CD3E", "THY1", "CD4", "CD8A", "FOXP3", "IL2RA", "GATA3", 
                     "RORA", "NCR1",
                     "ITGAM","CD14", "GRN", "C1QA", "LILRB1", "IL1B", "MRC1", "ITGAX", 
                     "ITGAE", "CCXCR1", "SIGLEC9", "LY6D", "S100A8", "S100A9", "LY6H", 
                     "CSF1R", "ADGRE1", "CX3CR1", "TGFB1", "CD274", "IDO1", "PTGES2", "FCGR3A",
                     "FUT4", "CD33", "SELL", "ARG1", "CEACAM1", "HLA-DRB5", "LY6G6D", "FCGR1A", "NCAM1")
DotPlot(agg, features = rev(markers.to.plot), dot.scale = 10) + RotatedAxis()


mouskers.to.plot <- c("Cd19","Cd79a","Jchain","Klrd1", "Nkg7", "Nktr",
                     "Cd3d", "Cd3e", "Thy1", "Cd4", "Cd8a", "Foxp3", "Il2ra", "Gata3", 
                     "Rora", "Ncr1",
                     "Itgam", "Cd14", "Grn", "C1qa", "Pirb", "Il1b", "Mrc1", "Itgax", 
                     "Itgae", "Xcr1", "Siglech", "Ly6d", "S100a8", "S100a9", "Ly6c1", 
                     "Csf1r", "Adgre1", "Cx3cr1", "Tgfb1", "Cd274", "Ido1", "Ptges2", "Cd16",
                     "Fut4", "Cd33", "Sell", "Arg1", "Ceacam1", "H2-Eb1", "Ly6g", "Fcgr1", "Ncam1")
           
DotPlot(immune.combined, features = rev(mouskers.to.plot), dot.scale = 10) + RotatedAxis()

#High Resolution Figures

#Change resolution to modify number or clusters
spo <- FindClusters(spo, resolution = 0.65)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
spomune.markers <- FindAllMarkers(spo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(spomune.markers, file="../Thesis/IMPx/IMPx5/MoSpo_0.65Deg.csv")


# Visualization
p1 <- DimPlot(spo, reduction = "umap", group.by = "stim")
p2 <- DimPlot(spo, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(spo, reduction = "umap", split.by = "gemgroup")

#Change resolution to modify number or clusters
imp <- FindClusters(imp, resolution = 0.65)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
impmune.markers <- FindAllMarkers(imp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(impmune.markers, file="../Thesis/IMPx/IMPx5/MoImp_0.65Deg.csv")


# Visualization
p1 <- DimPlot(imp, reduction = "umap", group.by = "stim")
p2 <- DimPlot(imp, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(imp, reduction = "umap", split.by = "gemgroup")

#Change resolution to modify number or clusters
immune.combined <- FindClusters(immune.combined, resolution = 0.65)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
immune.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/MoAgg_0.65Deg.csv")

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(immune.combined, reduction = "umap", split.by = "stim")

#Change resolution to modify number or clusters
agg <- FindClusters(agg, resolution = 0.65)

# Visualization
p1 <- DimPlot(agg, reduction = "umap", group.by = "gemgroup")
p2 <- DimPlot(agg, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(agg, reduction = "umap", label = TRUE )

#Find markers for every cluster compared to all remaining cells, report only the positive ones
agg.markers <- FindAllMarkers(agg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
agg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/HuAgg_0.65Deg.csv")

#Plots
DefaultAssay(imp) <- "RNA"
DefaultAssay(spo) <- "RNA"
DefaultAssay(agg) <- "RNA"
DefaultAssay(immune.combined) <- "RNA"

FeaturePlot(agg, features = c("CD3D", "CD3E", "CD3G"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD3D", "CD3E", "CD3G"))
FeaturePlot(agg, features = c("HCST", "KLRD1", "NKG7"), min.cutoff = "q9")
VlnPlot(agg, features = c("HCST", "KLRD1", "NKG7"))
FeaturePlot(agg, features = c("C1QA", "C1QB", "C1QC"), min.cutoff = "q9")
VlnPlot(agg, features = c("C1QA", "C1QB", "C1QC"))
FeaturePlot(agg, features = c("CD14", "CD68", "GBP2"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "GBP2"))
FeaturePlot(agg, features = c("CD14", "CD68", "GRN"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "GRN"))
FeaturePlot(agg, features = c("CD14", "CD68", "MRC1"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "MRC1"))
FeaturePlot(agg, features = c("CD14", "CD68", "ITGAM"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "CD68", "ITGAM"))
FeaturePlot(agg, features = c("ITGAX", "CD83", "CD14"), min.cutoff = "q9")
VlnPlot(agg, features = c("ITGAX", "CD83", "CD14"))
FeaturePlot(agg, features = c("CD74", "HLA-DQA1", "HLA-DRB5"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD74", "HLA-DQA1", "HLA-DRB5"))
FeaturePlot(agg, features = c("CD68", "MARCKS", "IL1B"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD68", "MARCKS", "IL1B"))
FeaturePlot(agg, features = c("ITGAX", "S100A8", "S100A9"), min.cutoff = "q9")
VlnPlot(agg, features = c("ITGAX", "S100A8", "S100A9"))
FeaturePlot(agg, features = c("NCF2", "S100A8", "S100A9"), min.cutoff = "q9")
VlnPlot(agg, features = c("NCF2", "S100A8", "S100A9"))
FeaturePlot(agg, features = c("CD14", "LILRB1", "LILRB2"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD14", "LILRB1", "LILRB2"))
FeaturePlot(agg, features = c("CD68", "CX3CR1", "TMEM119"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD68", "CX3CR1", "TMEM119"))
FeaturePlot(agg, features = c("CD79A", "CD79B", "MS4A3"), min.cutoff = "q9")
VlnPlot(agg, features = c("CD79A", "CD79B", "MS4A7"))


FeaturePlot(imp, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(imp, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(imp, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(imp, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(imp, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(imp, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(imp, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(imp, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)
FeaturePlot(imp, features = c("Ly6a", "Ly6c2", "Ly6d"), min.cutoff = "q9")
VlnPlot(imp, features = c("Ly6a", "Ly6c2", "Ly6d"), pt.size = 0)


FeaturePlot(spo, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(spo, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(spo, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(spo, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(spo, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(spo, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(spo, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(spo, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)




FeaturePlot(immune.combined, features = c("Cd24a", "S100a8", "S100a9"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd24a", "S100a8", "S100a9"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd3d", "Cd3e", "Cd3g"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd74", "H2-Eb1", "H2-Aa"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd74", "H2-Eb1", "H2-Aa"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd68", "Cx3cr1", "Tmem119"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd68", "Cx3cr1", "Tmem119"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Cd79a", "Cd79b", "Ms4a7"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Cd79a", "Cd79b", "Ms4a7"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Ly6a", "Ly6c2", "Ly6d"), min.cutoff = "q9")
VlnPlot(immune.combined, features = c("Ly6a", "Ly6c2", "Ly6d"), pt.size = 0)
FeaturePlot(immune.combined, features = c("Ms4a3"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("Cd1", "Cd1b", "Cd1c", "Cd1d"))


#Percentage contribution
agg.meta.data <- agg@meta.data

agg.counts <- group_by(agg.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(agg.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Percentage contribution
mo.meta.data <- immune.combined@meta.data

mo.counts <- group_by(mo.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(mo.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Dot plots
markers.to.plot <- c("CD19","CD79A","JCHAIN", "KLRD1", "NKG7", "NKTR",
                     "CD3D", "CD3E", "THY1", "CD4", "CD8A", "FOXP3", "IL2RA", "GATA3", 
                     "RORA", "NCR1",
                     "ITGAM","CD14", "GRN", "C1QA", "LILRB1", "IL1B", "MRC1", "ITGAX", 
                     "ITGAE", "CCXCR1", "SIGLEC9", "LY6D", "S100A8", "S100A9", "LY6H", 
                     "CSF1R", "ADGRE1", "CX3CR1", "TGFB1", "CD274", "IDO1", "PTGES2", "FCGR3A",
                     "FUT4", "CD33", "SELL", "ARG1", "CEACAM1", "HLA-DRB5", "LY6G6D", "FCGR1A", "NCAM1")
DotPlot(agg, features = rev(markers.to.plot), dot.scale = 10) + RotatedAxis()


mouskers.to.plot <- c("Cd19","Cd79a","Jchain","Klrd1", "Nkg7", "Nktr",
                      "Cd3d", "Cd3e", "Thy1", "Cd4", "Cd8a", "Foxp3", "Il2ra", "Gata3", 
                      "Rora", "Ncr1",
                      "Itgam", "Cd14", "Grn", "C1qa", "Pirb", "Il1b", "Mrc1", "Itgax", 
                      "Itgae", "Xcr1", "Siglech", "Ly6d", "S100a8", "S100a9", "Ly6c1", 
                      "Csf1r", "Adgre1", "Cx3cr1", "Tgfb1", "Cd274", "Ido1", "Ptges2", "Cd16",
                      "Fut4", "Cd33", "Sell", "Arg1", "Ceacam1", "H2-Eb1", "Ly6g", "Fcgr1", "Ncam1")

DotPlot(immune.combined, features = rev(mouskers.to.plot), dot.scale = 10) + RotatedAxis()


#save png
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="~/Documents/Thesis/IMPx/IMPx5/")

