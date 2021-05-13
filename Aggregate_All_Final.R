install.packages('dplyr')
install.packages('Seurat')
install.packages('batchbench')
install.packages('devtools')
devtools::install_github('xzhoulab/iDEA')
install.packages("lattice")
install.packages("clusterProfiler")
require(devtools)
devtools::install_version('flexmix', '2.3-13')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
library(Seurat)
library(dplyr)
library(cowplot)
library(reticulate)
library(ggplot2)
library(tictoc)
library(future)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)
library(iDEA)
library(vegan)
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

plot1 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


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
DimPlot(immune.combined, reduction = "umap", label = FALSE) + NoLegend()

immune.combined[["percent.hsp"]] <- PercentageFeatureSet(immune.combined, pattern = "^Hsp")
imp[["percent.hsp"]] <- PercentageFeatureSet(imp, pattern = "^Hsp")
spo[["percent.hsp"]] <- PercentageFeatureSet(spo, pattern = "^Hsp")


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
ElbowPlot(agg)

# t-SNE and Clustering
agg <- RunUMAP(agg, reduction = "pca", dims = 1:15)
agg <- FindNeighbors(agg, reduction = "pca", dims = 1:15)

#Change resolution to modify number or clusters
agg <- FindClusters(agg, resolution = 0.65)

# Visualization
p1 <- DimPlot(agg, reduction = "umap", group.by = "gemgroup")
p2 <- DimPlot(agg, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(agg, reduction = "umap", label = TRUE )

agg[["percent.hsp"]] <- PercentageFeatureSet(agg, pattern = "^HSP")


#Find markers for every cluster compared to all remaining cells, report only the positive ones
agg.markers <- FindAllMarkers(agg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
agg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(immune.markers, file="../Thesis/IMPx/IMPx5/HuAgg_0.1Deg.csv")



#Plots
DefaultAssay(imp) <- "RNA"
DefaultAssay(spo) <- "RNA"
DefaultAssay(agg) <- "RNA"
DefaultAssay(immune.combined) <- "RNA"

FeaturePlot(agg, features = c("IL18"))
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
FeaturePlot(agg, features = c("QKI"), min.cutoff = "q9")


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
FeaturePlot(imp, features = c("percent.mt"), min.cutoff = "q9")
FeaturePlot(imp, features = c("percent.hsp"), min.cutoff = "q9")
FeaturePlot(imp, features = c("nFeature_RNA"), min.cutoff = "q9")
FeaturePlot(imp, features = c("nCount_RNA"), min.cutoff = "q9")


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
FeaturePlot(spo, features = c("percent.mt"), min.cutoff = "q9")
FeaturePlot(spo, features = c("percent.hsp"), min.cutoff = "q9")
FeaturePlot(spo, features = c("nFeature_RNA"), min.cutoff = "q9")
FeaturePlot(spo, features = c("nCount_RNA"), min.cutoff = "q9")





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
FeaturePlot(immune.combined, features = c("Arg1"))
FeaturePlot(immune.combined, features = c("Cd1", "Cd1b", "Cd1c", "Cd1d"))
FeaturePlot(immune.combined, features = c("percent.mt"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("percent.hsp"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("nFeature_RNA"), min.cutoff = "q9")
FeaturePlot(immune.combined, features = c("nCount_RNA"), min.cutoff = "q9")


#Percentage contribution
agg.meta.data <- agg@meta.data

agg.counts <- group_by(agg.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(agg.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Percentage contribution
mo.meta.data <- immune.combined@meta.data

mo.counts <- group_by(mo.meta.data, stim, seurat_clusters) %>% summarise(count = n())

ggplot(mo.counts, aes(stim, count, fill = seurat_clusters)) +
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
DimPlot(agg, reduction = "umap", label = FALSE)

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
FeaturePlot(agg, features = c("VTCN1"), min.cutoff = "q9")



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

#Percentage contribution
spo.meta.data <- spo@meta.data

spo.counts <- group_by(spo.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(spo.counts, aes(gemgroup, count, fill = seurat_clusters)) +
  geom_bar(stat = 'identity')

#Percentage contribution
imp.meta.data <- imp@meta.data

imp.counts <- group_by(imp.meta.data, gemgroup, seurat_clusters) %>% summarise(count = n())

ggplot(imp.counts, aes(gemgroup, count, fill = seurat_clusters)) +
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
file.copy(from=plots.png.paths, to="~/Documents/Watowich/GO_Violins_Features")


geneset.list<-list(
"Microglial_function"=c("CASP1","CX3CR1","NOD2","TLR6","TLR7"),
"Complement" = c("C1QA","C1QB","C1QC")

)

agg<-AddModuleScore(agg, 
                    features =geneset.list, 
                    name=c("Microglial_function", "Complement"))

agg.matrix <- as.matrix(GetAssayData(agg, slot = "counts"))
agg.matrix <- rbind(agg.matrix, agg$seurat_clusters)

markers <- himmune.markers
# Filter
markers <- markers[ markers$p_val_adj<0.01, ]
# Sort
markers <- markers[ order(markers$cluster, -markers$avg_logFC), ]

#Set your number of clusters
Groups <- c("0":"19")

# Set number of upregulated genes
NTop <-50

#Entrez4Enrichment- convert genes to EntrezID for GO conversion
for (group in Groups)
{
  GeneSymbol4Enrichment<-markers[markers$cluster==group,] %>% group_by(cluster) %>% 
    top_n(n = NTop, wt = avg_logFC)
  GeneSymbol4Enrichment.group<-GeneSymbol4Enrichment$gene
  converted <- bitr(GeneSymbol4Enrichment.group, fromType="SYMBOL", toType="ENTREZID", 
                    OrgDb="org.Hs.eg.db")
  Entrez4Enrichment <- converted$ENTREZID
  assign(paste0("Entrez4Enrichment", group), converted$ENTREZID)
}


# background gene selection
NTopBkgGenes<-2000
ExpressionVector<-rowMeans(as.matrix(GetAssayData(agg, slot = "data")))
ExpressionVector<-sort(ExpressionVector, decreasing = T)
GeneSymbol4Bkg<-names(ExpressionVector[1:NTopBkgGenes])
converted <- bitr(GeneSymbol4Bkg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Entrez4Bkg <- converted$ENTREZID

#choose the gene ontologies you want to check
Onts <- c("CC","MF","BP")

for (group in Groups)
{
  Entrez_counter <- paste("Entrez4Enrichment", group, sep = "")
  Entrez_counter <- eval(parse(text = Entrez_counter))
  for (ont in Onts)
  {
    print("Working on")
    print(group)
    print(ont)
    enrich.res<-enrichGO(gene= Entrez_counter, universe = Entrez4Bkg,
                         ont=ont,OrgDb='org.Hs.eg.db' ,readable=T,
                         minGSSize = 20, maxGSSize = 500, qvalueCutoff = 0.05)
    if(is.null(enrich.res) == T){
      
    }
    else{
      print("computed")
      assign(paste0("enrich.res", group, ont), enrich.res)
    }
  }
}

options(repr.plot.width=10, repr.plot.height=5) # control figure size
# dot plot of enrichment analysis
dotplot(enrich.res11CC, showCategory=20) +ggtitle("Significantly Enriched Pathways")

# bar plot of enrichment analysis
barplot(enrich.res11BP, showCategory=20)+ggtitle("Significantly Enriched Pathways")

# table of enrichment analysis results
enrich.df <- as.data.frame(enrich.res)
enrich.df$cluster <-Group
(head(enrich.df))

# After the first step of over-representation analysis, we observed some meaningful gene set terms:
# Glucosamine metabolic process, cartilage development. 
# To compute their scores, we get their gene lists from MSigDB and calculate their scores with Seurat.
# We could also evaluate these scores with some other R packages like GSVA or AUCell.

mo.geneset.list<-list(
  "Adhesion"=c("Ada","Alcam","Amica1","Angpt2","Cd164","Cd22","Cd33","Cd6","Cd63","Cd84","Cd9","Cd96","Cd97","Cdh5","Col3a1","Csf3r","Cyfip2","Fn1","Glycam1","Icam2","Icam4","Irf2","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgax","Itgb1","Itgb4","Jam3","Kdr","Klra4","Klra5","Klra6","Klra7","Ly9","Lyve1","Map2k1","Mcam","Mertk","Mfge8","Mmp9","Msln","Ncam1","Nrp1","Pecam1","Plau","Pnma1","Pvrl2","Saa1","Sell","Selplg","Siglec1","Spink5","Spn","Spp1","Tek","Thy1","Tnfrsf12a","Vcam1","Vwf"), 
  "Adaptive_immunity"=c("C1qbp","C3ar1","Camp","Ccl1","Ccl11","Ccl12","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccrl2","Cd28","Cd4","Cd40","Cd40lg","Cd80","Cd86","Cd8a","Cd97","Cklf","Cma1","Creb5","Crp","Csf2","Cxcl1","Cxcl10","Cxcl13","Cxcl14","Cxcl2","Cxcl3","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Fasl","Fcer1a","Fcer2a","Foxp3","Fpr2","Gata3","H2-Q10","Hmgb1","Icam1","Ido1","Ifna2","Ifnar1","Ifnb1","Ifng","Ifngr1","Il10","Il13","Il17a","Il18","Il1a","Il1b","Il1r1","Il2","Il22","Il23a","Il24","Il25","Il2ra","Il4","Il5","Il6","Il6st","Il9","Irf3","Irf7","Itgam","Itk","Jak2","Jam3","Klrb1","Lta","Mapk8","Mbl2","Ms4a2","Mx1","Nfatc4","Nfkb1","Nod2","Nos2","Pnma1","Rag1","Rorc","S100a8","Sele","Slc11a1","Stat1","Stat3","Stat4","Stat6","Tbx21","Thbs1","Tlr4","Tlr6","Tnf","Tnfrsf11a","Txk","Xcr1"),
  "Antigen_processing"=c("Ccr7","Cd1d1","Cd1d2","Cd74","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Ea-ps","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Icam1","Mr1","Nod1","Nod2","Psmb8","Psmb9","Relb","Slc11a1","Tap1","Tap2","Tapbp"),
  "Apoptosis"=c("Abl1","Adora2a","Angpt1","Apoe","App","Atg5","Atg7","Atm","Bax","Bcl2l1","Bid","Birc5","Btk","C6","C9","Casp3","Casp8","Cd38","Cd3g","Cd5","Cd59b","Cdk1","Clec5a","Clu","Ctsh","Cyfip2","Cyld","Dusp6","Egfr","Egr3","Ep300","Ets1","Fadd","Fn1","Gpi1","Gzma","Gzmb","Hif1a","Hmgb1","Ifih1","Igf2r","Ikbke","Il19","Il24","Il3","Inpp5d","Itga1","Itga6","Jun","Kdr","Lck","Lcn2","Litaf","Lrp1","Ltbr","Ltk","Map2k4","Map3k1","Map3k5","Map3k7","Mapk1","Mapk3","Mapk8","Mef2c","Mertk","Mfge8","Mmp9","Muc1","Myc","Nefl","Nlrp3","Nos2","Osm","Pdcd1","Pik3cg","Plaur","Pml","Prkcd","Psen1","Ptgs2","Pycard","Rps6","Runx3","S100b","Sell","Spn","Spp1","Tcf7","Tek","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tmem173","Tnfaip3","Tnfrsf10b","Tnfrsf11b","Tnfrsf12a","Tnfrsf18","Tnfrsf8","Tnfsf10","Tnfsf12","Tnfsf14","Tnfsf15","Traf2","Traf3","Trp53","Twist1","Txnip","Vegfa","Xaf1"),
  "Autophagy"=c("Atg10","Atg12","Atg16l1","Atg5","Atg7","Lamp1"),
  "B_activation"=c("Cd28","Cd4","Cd69","Cd70","Cd83","Cd86","Cr2","Cxcr5","Dpp4","Fcer2a","Icosl","Il2ra","Il6","Ms4a1","Tgfb1"),
  "B_differentiation"=c("Ada","Cd79a","Flt3","Il10","Il11","Ptprc","Rag1"),
  "B_function"=c("Ada","Atm","Bcl10","Bcl2","Bcl6","Blk","Blnk","Bmi1","Btla","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd27","Cd28","Cd37","Cd38","Cd4","Cd40","Cd40lg","Cd69","Cd70","Cd74","Cd79a","Cd79b","Cd81","Cd83","Cd86","Cdkn1a","Cr2","Ctla4","Cxcl13","Cxcr5","Dpp4","Fas","Fcer2a","Fcgr2b","Flt3","Foxp3","Gpr183","H2-Q10","Icosl","Ikbkb","Ikbkg","Ikzf1","Il10","Il11","Il13","Il13ra1","Il1r2","Il21","Il2ra","Il2rg","Il4","Il5","Il6","Il7","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Ms4a1","Nt5e","Prdm1","Prkcd","Ptprc","Rag1","Sh2b2","Stat5b","Syk","Tgfb1","Ticam1","Tirap","Tnfaip3","Tnfrsf13b","Tnfrsf13c","Tnfrsf4","Tnfsf13b","Tnfsf4"),
  "B_proliferation"=c("Bcl2","Cd38","Cd81","Cdkn1a","Ctla4","Il7","Prkcd","Tnfrsf13b","Tnfrsf13c","Tnfsf13b"),
  "B_regulation"=c("Btla","Cd27","Cd40","Cd40lg","Foxp3","Il4"),
  "Bacterial_response"=c("Ccl2","Cd14","Fos","Hmgb1","Il10","Il1b","Irak1","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Ripk2","Tlr4","Tlr6","Tnfrsf1a"),
  "Basic_cell_function"=c("Arg1","Arg2","Chil3","Chit1","Cmpk2","Ctsg","Ctsl","Ddx60","Dock9","Dusp4","Emr1","Epsti1","Ewsr1","F12","F13a1","Fpr2","Gbp2b","Hamp","Hcst","Herc6","Hsd11b1","Isg15","Isg20","Map2k2","Mapk11","Mpo","Mpped1","Ncf4","Notch1","Oas2","Oas3","Oasl1","Pdgfc","Pla2g1b","Pla2g6","Pmch","Pou2af1","Prg2","Psmb7","Raet1c","Reps1","Rrad","Rsad2","S100a8","St6gal1","Stat2","Tab1","Tank","Timd4","Trem2","Ubc","Usp18","Usp9y","Ythdf2","Zfp13"),
  "Cancer_progression"=c("Akt3","Angpt1","Angpt2","Apoe","C3","C3ar1","Camp","Casp8","Ccl11","Ccl5","Ccl7","Ccl8","Ccr2","Ccr3","Cd163","Cd34","Cd36","Cd44","Cd46","Cdh1","Cdkn1a","Ceacam1","Cfp","Clu","Cma1","Col1a1","Col3a1","Col4a1","Crebbp","Csf2rb","Cspg4","Dll4","Egfr","Erbb2","Fap","Hif1a","Hspb2","Kdr","Mmp9","Msln","Psma2","Sele","Smad2","Smad3","Smad4","Smn1","Snai1","Tdo2","Tek","Tgfbr1","Tgfbr2","Tie1","Twist1","Vegfc","Vhl","Vim","Vwf"),
  "CD4_T_function"=c("Cd27","Cd4","Crebbp","Ctla4","Il15","Il7","Jak2","Mapk8","Ptprc","Socs3","Tgfb3","Tnfrsf4","Tnfrsf8","Tnfsf4","Tyk2","Yy1"),
  "CD_molecules"=c("Abcb1a","Alcam","Bst1","Bst2","Btla","C1qbp","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr9","Cd14","Cd160","Cd163","Cd164","Cd180","Cd19","Cd1d1","Cd1d2","Cd2","Cd200","Cd207","Cd209e","Cd22","Cd244","Cd247","Cd27","Cd274","Cd276","Cd28","Cd33","Cd34","Cd36","Cd37","Cd38","Cd3d","Cd3e","Cd3eap","Cd3g","Cd4","Cd40","Cd40lg","Cd46","Cd47","Cd48","Cd5","Cd53","Cd55","Cd6","Cd63","Cd68","Cd7","Cd70","Cd74","Cd79a","Cd79b","Cd80","Cd81","Cd83","Cd84","Cd86","Cd8a","Cd8b1","Cd9","Cd96","Cd97","Cd99","Cdh1","Cdh5","Ceacam1","Cr2","Csf1r","Csf2rb","Csf3r","Ctla4","Ctsw","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Dpp4","Eng","Entpd1","Epcam","Fas","Fcer2a","Fcgr1","Fcgr2b","Fcgr3","Fcgr4","Flt3","Icam1","Icam2","Icam4","Icos","Icosl","Ifitm1","Ifngr1","Igf1r","Igf2r","Igll1","Il10ra","Il12rb1","Il13ra1","Il13ra2","Il15ra","Il17ra","Il18r1","Il18rap","Il1r1","Il1r2","Il21r","Il2ra","Il2rb","Il2rg","Il3ra","Il4ra","Il5ra","Il6ra","Il6st","Il7r","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgam","Itgax","Itgb1","Itgb2","Itgb3","Itgb4","Kit","Klrb1","Klrc1","Klrc2","Klrk1","Lag3","Lamp1","Lamp2","Lamp3","Lilra5","Lrp1","Lrrn3","Ly9","Mcam","Mme","Mrc1","Ms4a1","Msr1","Mst1r","Muc1","Ncam1","Ncr1","Nrp1","Nt5e","Pdcd1","Pdcd1lg2","Pdgfrb","Pecam1","Plaur","Psmd7","Ptgdr2","Ptprc","Pvr","Sele","Sell","Selplg","Siglec1","Slamf1","Slamf6","Slamf7","Spn","Tfrc","Thbd","Thy1","Tlr1","Tlr2","Tlr3","Tlr4","Tlr6","Tlr8","Tlr9","Tnfrsf10b","Tnfrsf11a","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf4","Tnfsf8","Trem1","Vcam1"),
  "Cell_cycle"=c("Abcb1a","Abl1","Anp32b","Anxa1","App","Atm","Bid","Birc5","Casp3","Ccnd3","Cdk1","Cdkn1a","Cxcl15","Cyld","Ets1","Il12a","Il12b","Itgb1","Map2k1","Mapk3","Muc1","Myc","Nfatc1","Pin1","Pml","Prkce","Rps6","Runx3","Smpd3","Stat5b","Tal1","Tgfb1","Tgfb2","Thbs1","Trp53","Txnip"),
  "Chemokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cx3cr1","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Xcr1"),
  "Chemokines"=c("Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Cklf","Csf1r","Cx3cl1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Elane","Hc","Il18","Il1rl1","Il22ra1","Il4ra","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Sigirr","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnfsf4","Xcl1"),
  "Chemokines_chemokine_receptors"=c("C5ar1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cklf","Cmklr1","Csf1r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Elane","Hc","Ifng","Il16","Il18","Il1b","Il1rl1","Il22ra1","Il4","Il4ra","Il6","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Ppbp","Sigirr","Tgfb1","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnf","Tnfsf4","Xcl1","Xcr1"),
  "Chemotaxis"=c("C5ar1","Cmklr1","Ifng","Il16","Il1b","Il4","Il6","Ppbp","Tgfb1","Tnf"),
  "Complement_pathway"=c("A2m","C1qa","C1qb","C1qbp","C1ra","C1s1","C2","C3","C3ar1","C4b","C5ar1","C6","C7","C8a","C8b","C8g","C9","Cd46","Cd55","Cd59b","Cfb","Cfd","Cfh","Cfi","Cfp","Cr2","Crp","Hc","Masp1","Masp2","Mbl2","Serping1"),
  "Cytokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Csf1r","Cxcr1","Cxcr4","Erbb2","Flt3","Lyn"),
  "Cytokines"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr2","Cxcr3","Ebi3","Elane","F2rl1","Fasl","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytokines_cytokine_receptors"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf1r","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Ebi3","Elane","Erbb2","F2rl1","Fasl","Flt3","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Lyn","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytotoxicity"=c("Cd1d2","Cd8a","Ctsh","Fcgr1","Fcgr3","Gzmb","Gzmk","Gzmm","H2-D1","H2-K1","H2-M3","H2-T23","H60a","Il12a","Il21","Il23a","Il7r","Klrb1c","Klrk1","Lag3","Prf1","Ptprc","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Stat5b","Tap1","Ulbp1","Xcl1"),
  "DC_function"=c("Ccl19","Ccl5","Ccr1","Ccr2","Ccr5","Cd40","Cd40lg","Cd83","Cd86","Cr2","Cxcr1","Cxcr4","Il10","Lyn","Rag1","Relb","Tgfb1"),
  "Humoral_immune_rersponse"=c("Aire","Blnk","Bmi1","Bst1","Bst2","C4b","C5ar1","Ccl12","Ccl2","Ccl22","Ccl3","Ccl7","Ccr2","Ccr6","Ccr7","Cd28","Cd37","Cd40","Cd83","Crp","Cxcl13","Ebi3","Fcer2a","Foxj1","Gpi1","Gpr183","Ifnb1","Ifng","Il10","Il1b","Il6","Il7","Itgb2","Lta","Ltf","Ly86","Ly96","Mbl2","Mef2c","Mnx1","Ms4a1","Ms4a2","Nfkb1","Nod2","Pax5","Pdcd1","Pou2af1","Pou2f2","Psmb10","Sh2d1a","St6gal1","Tfe3","Tfeb","Tnf","Trem1","Trem2","Ythdf2"),
  "Inflammation"=c("Adora2a","Anxa1","Axl","Bcl6","C1qbp","C3","C3ar1","C4b","Camp","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr3","Ccr4","Ccr7","Ccrl2","Cd14","Cd163","Cd180","Cd276","Cd28","Cd40","Cd40lg","Cd47","Cd97","Cebpb","Cklf","Clec7a","Cma1","Creb5","Crp","Csf1","Cspg4","Cxcl1","Cxcl10","Cxcl11","Cxcl13","Cxcl14","Cxcl15","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Elane","F2rl1","Fas","Fasl","Fcer1a","Fcer2a","Fcgr2b","Fos","Foxp3","Fpr2","Hc","Hck","Hmgb1","Ido1","Il10","Il13","Il17a","Il17b","Il17f","Il18","Il1a","Il1b","Il1r1","Il1rap","Il1rl1","Il1rn","Il22","Il23r","Il24","Il25","Il27","Il2ra","Il34","Il4","Il4ra","Il5ra","Il6","Il6st","Il9","Itgb2","Jam3","Klrb1","Lbp","Lta","Ly86","Ly96","Lyn","Mapkapk2","Mefv","Mif","Ms4a2","Myd88","Nfatc4","Nfkb1","Nlrp3","Nos2","Nt5e","Pik3cd","Pik3cg","Pnma1","Pparg","Ptgs2","Pycard","Ripk2","S100a8","Sbno2","Sele","Stat3","Tgfb1","Thbs1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnf","Tnfaip3","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfsf4","Tollip","Xcl1","Xcr1"),
  "Innate_immune_response"=c("A2m","Abca1","Abcg1","Abl1","Aire","Angpt1","App","Atf1","Atf2","Atg12","Atg5","Axl","Bcl10","Bcl2","Bcl2l1","Bid","Bst2","Btk","C1qa","C1qb","C1ra","C1s1","C2","C3","C4b","C6","C7","C8a","C8b","C8g","C9","Camp","Card9","Casp1","Casp8","Ccl17","Ccl2","Ccl5","Ccr1","Ccr3","Ccr6","Ccr9","Cd14","Cd180","Cd1d1","Cd1d2","Cd36","Cd4","Cd46","Cd55","Cd74","Cd8a","Cd97","Cdk1","Cebpb","Cfb","Cfd","Cfh","Cfi","Cfp","Chuk","Clec4a2","Clec4n","Clec5a","Clec7a","Clu","Colec12","Cr2","Creb1","Crebbp","Crp","Csf1","Csf1r","Csf2","Ctss","Cxcl10","Cxcl11","Cxcl16","Cxcl2","Cxcl9","Cxcr2","Cxcr3","Cxcr6","Cybb","Cyfip2","Cyld","Ddx58","Defb1","Dmbt1","Dusp4","Dusp6","Ecsit","Elk1","Ep300","F12","F2rl1","Fadd","Fcgr1","Fos","Gzmk","Gzmm","Hamp","Hc","Hck","Hmgb1","Ifih1","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifnb1","Ifngr1","Ikbkb","Ikbke","Ikbkg","Il18r1","Il18rap","Il1a","Il1b","Il1r1","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il23a","Il23r","Il27","Il34","Il4","Irak1","Irak2","Irak3","Irak4","Irf3","Irf7","Irgm2","Isg15","Isg20","Itch","Itga5","Itgam","Itgax","Jak1","Jak2","Jak3","Klrg1","Lbp","Lcn2","Lgals3","Lilra5","Ly86","Ly96","Lyn","Map2k1","Map2k2","Map2k4","Map3k1","Map3k5","Map3k7","Map4k2","Mapk1","Mapk11","Mapk14","Mapk3","Mapk8","Mapkapk2","Marco","Masp1","Masp2","Mavs","Mbl2","Mefv","Mif","Mst1r","Mx1","Mx2","Myd88","Ncf4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pik3cd","Pik3cg","Pin1","Pparg","Pycard","Rela","Ripk2","S100b","Saa1","Serping1","Sigirr","Slamf7","Stat1","Syk","Tab1","Tank","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tmem173","Tnfaip3","Tollip","Traf2","Traf3","Traf6","Txnip","Tyk2","Ubc","Xcl1","Zbp1"),
  "Interferon_response"=c("Ccr7","Cd3e","Ciita","Cxcl16","Ddx58","Eomes","Fadd","Gbp5","H2-Aa","H2-Ab1","H60a","Ifi27","Ifi35","Ifi44","Ifi44l","Ifih1","Ifit1","Ifit2","Ifit3","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifna4","Ifnar1","Ifnar2","Ifngr1","Ifnl2","Il12rb2","Irf7","Irf8","Irgm2","Mavs","Nlrc5","Nos2","Runx3","Sh2d1b1","Tbk1","Tmem173","Ulbp1"),
  "Interleukins"=c("A2m","Bcl10","Card11","Card9","Casp1","Ccl2","Ccr2","Ccr7","Cd1d1","Cd1d2","Cd276","Cd28","Cd34","Cd3e","Cd40lg","Cd83","Cma1","Cmklr1","Csf2","Cxcl15","Cxcr1","Cxcr2","Ebi3","Egr1","Elane","Fadd","Fcer1a","Fcer1g","Fcgr2b","Foxj1","Foxp3","Gfi1","H2-Eb1","H2-Q1","Icosl","Ido1","Ifng","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il13","Il13ra1","Il15","Il16","Il17a","Il17f","Il17ra","Il18","Il18rap","Il19","Il1a","Il1b","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il3","Il4","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il9","Inpp5d","Irak2","Irak3","Irak4","Irf1","Irf4","Irf8","Itk","Jak2","Jak3","Lag3","Ltb","Mapk3","Mapkapk2","Mavs","Myd88","Nfkb1","Nlrp3","Nod1","Nod2","Pycard","Rel","Rela","Sele","Stat5b","Syk","Tcf7","Ticam1","Ticam2","Tigit","Tirap","Tlr1","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnfrsf11a","Tnfrsf1a","Tnfsf4","Tollip","Traf2","Traf6","Txk"),
  "Leukocyte_function"=c("Ccl19","Ccl25","Ccl4","Ccr1","Ccr7","Cd34","Cklf","Clec7a","Cxcl10","Cxcl12","Cxcl2","Cxcl5","Elane","F2rl1","Foxj1","Fut7","Hc","Icam1","Il16","Il23r","Il3","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgal","Itgam","Itgb1","Itgb2","Itgb3","Jam3","Lbp","Pecam1","Psen1","Psen2","Sele","Selplg","Syk","Tlr2","Vcam1"),
  "Macrophage_function"=c("C3ar1","Casp8","Ccl2","Ccl5","Ccr2","Ccr5","Ccr7","Cd1d1","Cd69","Cklf","Cmklr1","Crp","Csf1","Csf1r","Csf2","Cx3cl1","Cx3cr1","Eng","Fcer1a","Fcer2a","H60a","Hc","Il13","Il17f","Il18","Il1b","Il1rl1","Il23a","Il34","Il4","Il4ra","Itgb3","Lbp","Msr1","Nfkbia","Pparg","Prkce","Rora","Saa1","Sbno2","Syk","Thbs1","Tlr1","Ulbp1","Vegfa"),
  "Mast_cell_function"=c("C5ar1","Cd48","Fcer1a","Tpsab1"),
  "Mature_B_function"=c("Atm","Bcl10","Bcl6","Blk","Blnk","Bmi1","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd37","Cd74","Cd79b","Cxcl13","Fas","Fcgr2b","Gpr183","H2-Q10","Ikbkb","Ikbkg","Ikzf1","Il13","Il13ra1","Il1r2","Il21","Il2rg","Il5","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Nt5e","Prdm1","Sh2b2","Stat5b","Syk","Ticam1","Tirap","Tnfaip3","Tnfrsf4","Tnfsf4"),
  "Mature_T_function"=c("Anxa1","Bcl10","Card11","Ccl19","Ccl2","Ccl3","Ccr2","Ccr7","Cd1d2","Cd247","Cd274","Cd276","Cd48","Cd5","Cd59b","Cd80","Cd83","Cd86","Ctsh","Cxcl13","Eomes","Fasl","Fcgr4","Foxj1","Fut7","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Icam1","Ido1","Ifnar1","Ifnb1","Ikzf1","Il12a","Il12rb1","Il23a","Il2rg","Il6st","Il7r","Itgal","Itgam","Itgax","Itgb2","Itk","Lcp1","Mapk1","Mill2","Psen1","Psen2","Psmb10","Pvr","Pvrl2","Rag1","Rorc","Rps6","Spn","Stat5b","Syk","Tap1","Tcf7","Tgfb1","Tgfb2","Tigit","Tnfsf18","Tnfsf8","Traf2","Txk","Vcam1","Xcl1","Zap70"),
  "MHC1_MHC2"=c("Cd160","Cd1d1","Cd40lg","Cd74","Ciita","Ctsh","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Il10","Klrk1","Lag3","Mr1","Nlrc5","Pml","Tap1","Tap2","Tapbp"),
  "Microglial_function"=c("Casp1","Cx3cr1","Nod2","Tlr6","Tlr7"),
  "NK_function"=c("Axl","Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Cd2","Cd244","Cd247","Cd7","Cd96","Flt3l","H2-M3","H60a","Ifi27","Ikzf1","Il11ra1","Il12a","Il12b","Il12rb1","Il15","Il15ra","Il21","Il23a","Il2rb","Itgb2","Klra1","Klra15","Klra17","Klra2","Klra20","Klra21","Klra27","Klra3","Klra4","Klra5","Klra6","Klra7","Klrb1c","Klrc1","Klrd1","Klrk1","Lag3","Mertk","Mill2","Ncam1","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Slamf7","Stat5b","Ulbp1"),
  "Pathogen_response"=c("Ccl2","Cd14","Fos","Hmgb1","Ifnb1","Ifng","Il10","Il12a","Il1b","Il6","Irak1","Irf3","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Rela","Ripk2","Tbk1","Ticam1","Tlr3","Tlr4","Tlr6","Tlr7","Tlr8","Tnf","Tnfrsf1a"),
  "Phagocytosis"=c("Anxa1","C3","Cd14","Cd36","Cd44","Cd47","Clec7a","Colec12","Crp","Csf1","Csf2","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Ifng","Il1rl1","Itgam","Itgb2","Marco","Mbl2","Mfge8","Mif","Myd88","Nod1","Pecam1","Siglec1","Ticam1","Tlr3","Tlr9","Tnf","Tnfsf11"),
  "Regulation_inflammatory_response"=c("Bcl6","C3","C3ar1","C4b","Ccl1","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl3","Ccl4","Ccl7","Ccl8","Ccr2","Ccr3","Ccr4","Ccr7","Cd14","Cd40","Cd40lg","Cebpb","Crp","Csf1","Cxcl1","Cxcl10","Cxcl11","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Fos","Il10","Il1r1","Il1rap","Il22","Il9","Itgb2","Ly96","Myd88","Nfkb1","Nos2","Ripk2","Sele","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr9","Tollip"),
  "Treg_function"=c("Ccr6","Foxp3","Ikzf2","Il9","Irf4","Irf8","Pou2f2","Rel","Tnfsf11"),
  "Senescence"=c("Abl1","Atm","Bmi1","Cdkn1a","Egr1","Ets1","Hras","Ifng","Igf1r","Irf3","Irf5","Irf7","Map2k1","Mapk14","Myc","Nfkb1","Plau","Prkcd","Serpinb2","Tgfb1","Trp53","Twist1"),
  "T_cell_function"=c("Ada","Anxa1","Bcl10","Bcl2","Bcl6","Btla","Card11","Casp3","Ccl11","Ccl19","Ccl2","Ccl3","Ccl5","Ccl7","Ccnd3","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Cd1d1","Cd1d2","Cd2","Cd247","Cd27","Cd274","Cd276","Cd28","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd40lg","Cd47","Cd48","Cd5","Cd59b","Cd74","Cd80","Cd83","Cd86","Cd8a","Cd8b1","Cebpb","Cma1","Crebbp","Csf2","Ctla4","Ctsh","Cxcl12","Cxcl13","Cxcr3","Cxcr4","Dpp4","Egr1","Eomes","Fas","Fasl","Fcgr4","Flt3","Foxj1","Foxp3","Fut7","Gata3","Gfi1","Gpr44","Gzmb","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Havcr2","Icam1","Icos","Icosl","Ido1","Ifnar1","Ifnb1","Ifng","Ikzf1","Ikzf2","Il10","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il15","Il17a","Il18","Il18r1","Il18rap","Il1b","Il1r1","Il2","Il21","Il23a","Il25","Il27","Il2ra","Il2rg","Il4","Il4ra","Il5","Il6","Il6st","Il7","Il7r","Il9","Irf1","Irf4","Irf8","Itch","Itga1","Itgal","Itgam","Itgax","Itgb2","Itk","Jak1","Jak2","Jak3","Lag3","Lck","Lcp1","Lgals3","Maf","Map3k7","Mapk1","Mapk8","Mill2","Nfatc1","Nfatc2","Nfkb1","Nos2","Pdcd1","Pdcd1lg2","Pou2f2","Psen1","Psen2","Psmb10","Ptprc","Pvr","Pvrl2","Rag1","Rel","Relb","Ripk2","Rora","Rorc","Rps6","Sell","Socs1","Socs3","Spn","Spp1","Stat1","Stat4","Stat5b","Stat6","Syk","Tap1","Tbx21","Tcf7","Tgfb1","Tgfb2","Tgfb3","Thy1","Tigit","Tlr4","Tlr6","Tmed1","Tnf","Tnfrsf13c","Tnfrsf14","Tnfrsf4","Tnfrsf8","Tnfsf11","Tnfsf13b","Tnfsf14","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf6","Trp53","Txk","Tyk2","Vcam1","Vegfa","Xcl1","Yy1","Zap70"),
  "T_cell_anergy"=c("Cma1","Ctla4","Gzmb","Icos","Itch","Itga1","Jak1","Lgals3","Pdcd1","Sell"),
  "T_cell_differentiation"=c("Ada","Bcl2","Cd1d1","Cd3d","Cd4","Cd74","Egr1","Flt3","Il27","Il7","Irf4","Nos2","Socs1"),
  "T_cell_proliferation"=c("Casp3","Ccnd3","Cd3e","Cxcl12","Cxcr4","Icosl","Il15","Il1b","Pdcd1lg2","Ptprc","Ripk2","Spp1","Tnfrsf13c","Tnfsf13b","Tnfsf14","Traf6","Trp53"),
  "T_cell_regulation"=c("Btla","Cd27","Cd3g","Cd47","Cd8a","Cd8b1","Dpp4","Fas","Foxp3","Icam1","Il2","Lag3","Lck","Map3k7","Tgfb1","Thy1","Tnfrsf14"),
  "Th1_Th2_differentiation"=c("Cd2","Cd28","Cd40","Cd40lg","Ifng","Il12b","Il18","Il2ra","Il4","Il6","Relb"),
  "Th1_function"=c("Ccr5","Csf2","Cxcr3","Havcr2","Il12rb2","Il18r1","Il18rap","Il2","Il27","Irf1","Nfkb1","Socs1","Stat1","Stat4","Tbx21","Tlr4","Tlr6","Tnf","Vegfa"),
  "Th2_function"=c("Bcl6","Ccl11","Ccl5","Ccl7","Ccr3","Ccr4","Cebpb","Gata3","Gfi1","Gpr44","Icos","Il10","Il13","Il13ra1","Il1r1","Il25","Il4ra","Il5","Jak1","Jak3","Maf","Nfatc1","Nfatc2","Stat6","Tmed1"),
  "Th17_function"=c("Il17a","Il21","Rora","Rorc"),
  "TLR"=c("Cd86","Gfi1","Irak1","Irak2","Irf3","Irf4","Map3k7","Mapkapk2","Myd88","Nfkbia","Prkce","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Traf3","Traf6"),
  "TNF_superfamily"=c("Cd27","Cd40","Cd40lg","Cd70","Fas","Fasl","Lta","Ltb","Ltbr","Tnf","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_members"=c("Cd40lg","Cd70","Fasl","Lta","Ltb","Tnf","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_receptors"=c("Cd27","Cd40","Fas","Ltbr","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9"),
  "Transcription_factors"=c("Atf1","Atf2","Batf","Bcl10","Bcl2","Bcl6","Card11","Card9","Cd34","Cd36","Cd40","Cdh1","Cebpb","Ciita","Clu","Cmklr1","Creb1","Creb5","Crebbp","Cyld","Ddx58","Ecsit","Egr1","Egr2","Elk1","Eomes","Ep300","Ets1","Fos","Foxj1","Foxp3","Gata3","Gfi1","Gtf3c1","Hck","Hmgb1","Icam1","Ikbkb","Ikbkg","Ikzf1","Il10","Il4","Il5","Irak1","Irak2","Irak3","Irf1","Irf2","Irf3","Irf4","Irf5","Irf7","Irf8","Itch","Itgb2","Jun","Maf","Mapk1","Mapk3","Mavs","Mef2c","Mnx1","Myc","Nfatc2","Nfatc3","Nfatc4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pax5","Pou2f2","Pparg","Pycard","Rel","Rela","Relb","Rora","Rorc","Runx1","Runx3","Sigirr","Smad4","Snai1","Stat1","Stat3","Stat4","Stat6","Tal1","Tbx21","Tcf7","Tfe3","Tfeb","Ticam1","Tlr3","Tlr9","Tmem173","Tnfrsf11a","Tnfrsf4","Tnfrsf8","Tnfsf18","Traf2","Traf3","Traf6","Twist1","Xbp1","Yy1"),
  "Transmembrane_transporter"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Apoe","App","Atg10","Atg16l1","Atg7","Axl","Bax","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd36","Cd3g","Cmah","Col1a1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fez1","Fyn","Icam1","Ifit1","Igf2r","Il13","Il1b","Il4","Itgb3","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Mertk","Msr1","Nefl","Nfatc1","Nup107","Pparg","Prkcd","Prkce","Psen1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Tmed1"),
  "Transporter_function"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Anxa1","Apoe","App","Atg10","Atg12","Atg16l1","Atg5","Atg7","Axl","Bax","C3","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd14","Cd36","Cd3g","Cd44","Cd47","Clec7a","Cmah","Col1a1","Colec12","Crp","Csf1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Fez1","Fyn","Icam1","Ifit1","Ifng","Igf2r","Il13","Il1b","Il1rl1","Il4","Itgam","Itgb2","Itgb3","Lamp1","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Marco","Mbl2","Mertk","Mfge8","Mif","Msr1","Myd88","Nefl","Nfatc1","Nod1","Nup107","Pecam1","Pparg","Prkcd","Prkce","Psen1","Siglec1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Ticam1","Tlr3","Tlr9","Tmed1","Tnf","Tnfsf11")
)

hu.geneset.list<-list(
  "Adhesion"=c("Ada","Alcam","Amica1","Angpt2","Cd164","Cd22","Cd33","Cd6","Cd63","Cd84","Cd9","Cd96","Cd97","Cdh5","Col3a1","Csf3r","Cyfip2","Fn1","Glycam1","Icam2","Icam4","Irf2","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgax","Itgb1","Itgb4","Jam3","Kdr","Klra4","Klra5","Klra6","Klra7","Ly9","Lyve1","Map2k1","Mcam","Mertk","Mfge8","Mmp9","Msln","Ncam1","Nrp1","Pecam1","Plau","Pnma1","Pvrl2","Saa1","Sell","Selplg","Siglec1","Spink5","Spn","Spp1","Tek","Thy1","Tnfrsf12a","Vcam1","Vwf"), 
  "Adaptive_immunity"=c("C1qbp","C3ar1","Camp","Ccl1","Ccl11","Ccl12","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccrl2","Cd28","Cd4","Cd40","Cd40lg","Cd80","Cd86","Cd8a","Cd97","Cklf","Cma1","Creb5","Crp","Csf2","Cxcl1","Cxcl10","Cxcl13","Cxcl14","Cxcl2","Cxcl3","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Fasl","Fcer1a","Fcer2a","Foxp3","Fpr2","Gata3","H2-Q10","Hmgb1","Icam1","Ido1","Ifna2","Ifnar1","Ifnb1","Ifng","Ifngr1","Il10","Il13","Il17a","Il18","Il1a","Il1b","Il1r1","Il2","Il22","Il23a","Il24","Il25","Il2ra","Il4","Il5","Il6","Il6st","Il9","Irf3","Irf7","Itgam","Itk","Jak2","Jam3","Klrb1","Lta","Mapk8","Mbl2","Ms4a2","Mx1","Nfatc4","Nfkb1","Nod2","Nos2","Pnma1","Rag1","Rorc","S100a8","Sele","Slc11a1","Stat1","Stat3","Stat4","Stat6","Tbx21","Thbs1","Tlr4","Tlr6","Tnf","Tnfrsf11a","Txk","Xcr1"),
  "Antigen_processing"=c("Ccr7","Cd1d1","Cd1d2","Cd74","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Ea-ps","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Icam1","Mr1","Nod1","Nod2","Psmb8","Psmb9","Relb","Slc11a1","Tap1","Tap2","Tapbp"),
  "Apoptosis"=c("Abl1","Adora2a","Angpt1","Apoe","App","Atg5","Atg7","Atm","Bax","Bcl2l1","Bid","Birc5","Btk","C6","C9","Casp3","Casp8","Cd38","Cd3g","Cd5","Cd59b","Cdk1","Clec5a","Clu","Ctsh","Cyfip2","Cyld","Dusp6","Egfr","Egr3","Ep300","Ets1","Fadd","Fn1","Gpi1","Gzma","Gzmb","Hif1a","Hmgb1","Ifih1","Igf2r","Ikbke","Il19","Il24","Il3","Inpp5d","Itga1","Itga6","Jun","Kdr","Lck","Lcn2","Litaf","Lrp1","Ltbr","Ltk","Map2k4","Map3k1","Map3k5","Map3k7","Mapk1","Mapk3","Mapk8","Mef2c","Mertk","Mfge8","Mmp9","Muc1","Myc","Nefl","Nlrp3","Nos2","Osm","Pdcd1","Pik3cg","Plaur","Pml","Prkcd","Psen1","Ptgs2","Pycard","Rps6","Runx3","S100b","Sell","Spn","Spp1","Tcf7","Tek","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tmem173","Tnfaip3","Tnfrsf10b","Tnfrsf11b","Tnfrsf12a","Tnfrsf18","Tnfrsf8","Tnfsf10","Tnfsf12","Tnfsf14","Tnfsf15","Traf2","Traf3","Trp53","Twist1","Txnip","Vegfa","Xaf1"),
  "Autophagy"=c("ATG10","ATG12","ATG16L1","ATG5","ATG7","LAMP1"),
  "B_activation"=c("Cd28","Cd4","Cd69","Cd70","Cd83","Cd86","Cr2","Cxcr5","Dpp4","Fcer2a","Icosl","Il2ra","Il6","Ms4a1","Tgfb1"),
  "B_differentiation"=c("Ada","Cd79a","Flt3","Il10","Il11","Ptprc","Rag1"),
  "B_function"=c("Ada","Atm","Bcl10","Bcl2","Bcl6","Blk","Blnk","Bmi1","Btla","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd27","Cd28","Cd37","Cd38","Cd4","Cd40","Cd40lg","Cd69","Cd70","Cd74","Cd79a","Cd79b","Cd81","Cd83","Cd86","Cdkn1a","Cr2","Ctla4","Cxcl13","Cxcr5","Dpp4","Fas","Fcer2a","Fcgr2b","Flt3","Foxp3","Gpr183","H2-Q10","Icosl","Ikbkb","Ikbkg","Ikzf1","Il10","Il11","Il13","Il13ra1","Il1r2","Il21","Il2ra","Il2rg","Il4","Il5","Il6","Il7","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Ms4a1","Nt5e","Prdm1","Prkcd","Ptprc","Rag1","Sh2b2","Stat5b","Syk","Tgfb1","Ticam1","Tirap","Tnfaip3","Tnfrsf13b","Tnfrsf13c","Tnfrsf4","Tnfsf13b","Tnfsf4"),
  "B_proliferation"=c("Bcl2","Cd38","Cd81","Cdkn1a","Ctla4","Il7","Prkcd","Tnfrsf13b","Tnfrsf13c","Tnfsf13b"),
  "B_regulation"=c("Btla","Cd27","Cd40","Cd40lg","Foxp3","Il4"),
  "Bacterial_response"=c("Ccl2","Cd14","Fos","Hmgb1","Il10","Il1b","Irak1","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Ripk2","Tlr4","Tlr6","Tnfrsf1a"),
  "Basic_cell_function"=c("Arg1","Arg2","Chil3","Chit1","Cmpk2","Ctsg","Ctsl","Ddx60","Dock9","Dusp4","Emr1","Epsti1","Ewsr1","F12","F13a1","Fpr2","Gbp2b","Hamp","Hcst","Herc6","Hsd11b1","Isg15","Isg20","Map2k2","Mapk11","Mpo","Mpped1","Ncf4","Notch1","Oas2","Oas3","Oasl1","Pdgfc","Pla2g1b","Pla2g6","Pmch","Pou2af1","Prg2","Psmb7","Raet1c","Reps1","Rrad","Rsad2","S100a8","St6gal1","Stat2","Tab1","Tank","Timd4","Trem2","Ubc","Usp18","Usp9y","Ythdf2","Zfp13"),
  "Cancer_progression"=c("Akt3","Angpt1","Angpt2","Apoe","C3","C3ar1","Camp","Casp8","Ccl11","Ccl5","Ccl7","Ccl8","Ccr2","Ccr3","Cd163","Cd34","Cd36","Cd44","Cd46","Cdh1","Cdkn1a","Ceacam1","Cfp","Clu","Cma1","Col1a1","Col3a1","Col4a1","Crebbp","Csf2rb","Cspg4","Dll4","Egfr","Erbb2","Fap","Hif1a","Hspb2","Kdr","Mmp9","Msln","Psma2","Sele","Smad2","Smad3","Smad4","Smn1","Snai1","Tdo2","Tek","Tgfbr1","Tgfbr2","Tie1","Twist1","Vegfc","Vhl","Vim","Vwf"),
  "CD4_T_function"=c("Cd27","Cd4","Crebbp","Ctla4","Il15","Il7","Jak2","Mapk8","Ptprc","Socs3","Tgfb3","Tnfrsf4","Tnfrsf8","Tnfsf4","Tyk2","Yy1"),
  "CD_molecules"=c("Abcb1a","Alcam","Bst1","Bst2","Btla","C1qbp","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr9","Cd14","Cd160","Cd163","Cd164","Cd180","Cd19","Cd1d1","Cd1d2","Cd2","Cd200","Cd207","Cd209e","Cd22","Cd244","Cd247","Cd27","Cd274","Cd276","Cd28","Cd33","Cd34","Cd36","Cd37","Cd38","Cd3d","Cd3e","Cd3eap","Cd3g","Cd4","Cd40","Cd40lg","Cd46","Cd47","Cd48","Cd5","Cd53","Cd55","Cd6","Cd63","Cd68","Cd7","Cd70","Cd74","Cd79a","Cd79b","Cd80","Cd81","Cd83","Cd84","Cd86","Cd8a","Cd8b1","Cd9","Cd96","Cd97","Cd99","Cdh1","Cdh5","Ceacam1","Cr2","Csf1r","Csf2rb","Csf3r","Ctla4","Ctsw","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Dpp4","Eng","Entpd1","Epcam","Fas","Fcer2a","Fcgr1","Fcgr2b","Fcgr3","Fcgr4","Flt3","Icam1","Icam2","Icam4","Icos","Icosl","Ifitm1","Ifngr1","Igf1r","Igf2r","Igll1","Il10ra","Il12rb1","Il13ra1","Il13ra2","Il15ra","Il17ra","Il18r1","Il18rap","Il1r1","Il1r2","Il21r","Il2ra","Il2rb","Il2rg","Il3ra","Il4ra","Il5ra","Il6ra","Il6st","Il7r","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgam","Itgax","Itgb1","Itgb2","Itgb3","Itgb4","Kit","Klrb1","Klrc1","Klrc2","Klrk1","Lag3","Lamp1","Lamp2","Lamp3","Lilra5","Lrp1","Lrrn3","Ly9","Mcam","Mme","Mrc1","Ms4a1","Msr1","Mst1r","Muc1","Ncam1","Ncr1","Nrp1","Nt5e","Pdcd1","Pdcd1lg2","Pdgfrb","Pecam1","Plaur","Psmd7","Ptgdr2","Ptprc","Pvr","Sele","Sell","Selplg","Siglec1","Slamf1","Slamf6","Slamf7","Spn","Tfrc","Thbd","Thy1","Tlr1","Tlr2","Tlr3","Tlr4","Tlr6","Tlr8","Tlr9","Tnfrsf10b","Tnfrsf11a","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf4","Tnfsf8","Trem1","Vcam1"),
  "Cell_cycle"=c("Abcb1a","Abl1","Anp32b","Anxa1","App","Atm","Bid","Birc5","Casp3","Ccnd3","Cdk1","Cdkn1a","Cxcl15","Cyld","Ets1","Il12a","Il12b","Itgb1","Map2k1","Mapk3","Muc1","Myc","Nfatc1","Pin1","Pml","Prkce","Rps6","Runx3","Smpd3","Stat5b","Tal1","Tgfb1","Tgfb2","Thbs1","Trp53","Txnip"),
  "Chemokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cx3cr1","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Xcr1"),
  "Chemokines"=c("Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Cklf","Csf1r","Cx3cl1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Elane","Hc","Il18","Il1rl1","Il22ra1","Il4ra","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Sigirr","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnfsf4","Xcl1"),
  "Chemokines_chemokine_receptors"=c("C5ar1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cklf","Cmklr1","Csf1r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Elane","Hc","Ifng","Il16","Il18","Il1b","Il1rl1","Il22ra1","Il4","Il4ra","Il6","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Ppbp","Sigirr","Tgfb1","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnf","Tnfsf4","Xcl1","Xcr1"),
  "Chemotaxis"=c("C5ar1","Cmklr1","Ifng","Il16","Il1b","Il4","Il6","Ppbp","Tgfb1","Tnf"),
  "Complement_pathway"=c("A2m","C1qa","C1qb","C1qbp","C1ra","C1s1","C2","C3","C3ar1","C4b","C5ar1","C6","C7","C8a","C8b","C8g","C9","Cd46","Cd55","Cd59b","Cfb","Cfd","Cfh","Cfi","Cfp","Cr2","Crp","Hc","Masp1","Masp2","Mbl2","Serping1"),
  "Cytokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Csf1r","Cxcr1","Cxcr4","Erbb2","Flt3","Lyn"),
  "Cytokines"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr2","Cxcr3","Ebi3","Elane","F2rl1","Fasl","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytokines_cytokine_receptors"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf1r","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Ebi3","Elane","Erbb2","F2rl1","Fasl","Flt3","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Lyn","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytotoxicity"=c("Cd1d2","Cd8a","Ctsh","Fcgr1","Fcgr3","Gzmb","Gzmk","Gzmm","H2-D1","H2-K1","H2-M3","H2-T23","H60a","Il12a","Il21","Il23a","Il7r","Klrb1c","Klrk1","Lag3","Prf1","Ptprc","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Stat5b","Tap1","Ulbp1","Xcl1"),
  "DC_function"=c("Ccl19","Ccl5","Ccr1","Ccr2","Ccr5","Cd40","Cd40lg","Cd83","Cd86","Cr2","Cxcr1","Cxcr4","Il10","Lyn","Rag1","Relb","Tgfb1"),
  "Humoral_immune_rersponse"=c("Aire","Blnk","Bmi1","Bst1","Bst2","C4b","C5ar1","Ccl12","Ccl2","Ccl22","Ccl3","Ccl7","Ccr2","Ccr6","Ccr7","Cd28","Cd37","Cd40","Cd83","Crp","Cxcl13","Ebi3","Fcer2a","Foxj1","Gpi1","Gpr183","Ifnb1","Ifng","Il10","Il1b","Il6","Il7","Itgb2","Lta","Ltf","Ly86","Ly96","Mbl2","Mef2c","Mnx1","Ms4a1","Ms4a2","Nfkb1","Nod2","Pax5","Pdcd1","Pou2af1","Pou2f2","Psmb10","Sh2d1a","St6gal1","Tfe3","Tfeb","Tnf","Trem1","Trem2","Ythdf2"),
  "Inflammation"=c("Adora2a","Anxa1","Axl","Bcl6","C1qbp","C3","C3ar1","C4b","Camp","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr3","Ccr4","Ccr7","Ccrl2","Cd14","Cd163","Cd180","Cd276","Cd28","Cd40","Cd40lg","Cd47","Cd97","Cebpb","Cklf","Clec7a","Cma1","Creb5","Crp","Csf1","Cspg4","Cxcl1","Cxcl10","Cxcl11","Cxcl13","Cxcl14","Cxcl15","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Elane","F2rl1","Fas","Fasl","Fcer1a","Fcer2a","Fcgr2b","Fos","Foxp3","Fpr2","Hc","Hck","Hmgb1","Ido1","Il10","Il13","Il17a","Il17b","Il17f","Il18","Il1a","Il1b","Il1r1","Il1rap","Il1rl1","Il1rn","Il22","Il23r","Il24","Il25","Il27","Il2ra","Il34","Il4","Il4ra","Il5ra","Il6","Il6st","Il9","Itgb2","Jam3","Klrb1","Lbp","Lta","Ly86","Ly96","Lyn","Mapkapk2","Mefv","Mif","Ms4a2","Myd88","Nfatc4","Nfkb1","Nlrp3","Nos2","Nt5e","Pik3cd","Pik3cg","Pnma1","Pparg","Ptgs2","Pycard","Ripk2","S100a8","Sbno2","Sele","Stat3","Tgfb1","Thbs1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnf","Tnfaip3","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfsf4","Tollip","Xcl1","Xcr1"),
  "Innate_immune_response"=c("A2m","Abca1","Abcg1","Abl1","Aire","Angpt1","App","Atf1","Atf2","Atg12","Atg5","Axl","Bcl10","Bcl2","Bcl2l1","Bid","Bst2","Btk","C1qa","C1qb","C1ra","C1s1","C2","C3","C4b","C6","C7","C8a","C8b","C8g","C9","Camp","Card9","Casp1","Casp8","Ccl17","Ccl2","Ccl5","Ccr1","Ccr3","Ccr6","Ccr9","Cd14","Cd180","Cd1d1","Cd1d2","Cd36","Cd4","Cd46","Cd55","Cd74","Cd8a","Cd97","Cdk1","Cebpb","Cfb","Cfd","Cfh","Cfi","Cfp","Chuk","Clec4a2","Clec4n","Clec5a","Clec7a","Clu","Colec12","Cr2","Creb1","Crebbp","Crp","Csf1","Csf1r","Csf2","Ctss","Cxcl10","Cxcl11","Cxcl16","Cxcl2","Cxcl9","Cxcr2","Cxcr3","Cxcr6","Cybb","Cyfip2","Cyld","Ddx58","Defb1","Dmbt1","Dusp4","Dusp6","Ecsit","Elk1","Ep300","F12","F2rl1","Fadd","Fcgr1","Fos","Gzmk","Gzmm","Hamp","Hc","Hck","Hmgb1","Ifih1","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifnb1","Ifngr1","Ikbkb","Ikbke","Ikbkg","Il18r1","Il18rap","Il1a","Il1b","Il1r1","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il23a","Il23r","Il27","Il34","Il4","Irak1","Irak2","Irak3","Irak4","Irf3","Irf7","Irgm2","Isg15","Isg20","Itch","Itga5","Itgam","Itgax","Jak1","Jak2","Jak3","Klrg1","Lbp","Lcn2","Lgals3","Lilra5","Ly86","Ly96","Lyn","Map2k1","Map2k2","Map2k4","Map3k1","Map3k5","Map3k7","Map4k2","Mapk1","Mapk11","Mapk14","Mapk3","Mapk8","Mapkapk2","Marco","Masp1","Masp2","Mavs","Mbl2","Mefv","Mif","Mst1r","Mx1","Mx2","Myd88","Ncf4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pik3cd","Pik3cg","Pin1","Pparg","Pycard","Rela","Ripk2","S100b","Saa1","Serping1","Sigirr","Slamf7","Stat1","Syk","Tab1","Tank","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tmem173","Tnfaip3","Tollip","Traf2","Traf3","Traf6","Txnip","Tyk2","Ubc","Xcl1","Zbp1"),
  "Interferon_response"=c("Ccr7","Cd3e","Ciita","Cxcl16","Ddx58","Eomes","Fadd","Gbp5","H2-Aa","H2-Ab1","H60a","Ifi27","Ifi35","Ifi44","Ifi44l","Ifih1","Ifit1","Ifit2","Ifit3","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifna4","Ifnar1","Ifnar2","Ifngr1","Ifnl2","Il12rb2","Irf7","Irf8","Irgm2","Mavs","Nlrc5","Nos2","Runx3","Sh2d1b1","Tbk1","Tmem173","Ulbp1"),
  "Interleukins"=c("A2m","Bcl10","Card11","Card9","Casp1","Ccl2","Ccr2","Ccr7","Cd1d1","Cd1d2","Cd276","Cd28","Cd34","Cd3e","Cd40lg","Cd83","Cma1","Cmklr1","Csf2","Cxcl15","Cxcr1","Cxcr2","Ebi3","Egr1","Elane","Fadd","Fcer1a","Fcer1g","Fcgr2b","Foxj1","Foxp3","Gfi1","H2-Eb1","H2-Q1","Icosl","Ido1","Ifng","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il13","Il13ra1","Il15","Il16","Il17a","Il17f","Il17ra","Il18","Il18rap","Il19","Il1a","Il1b","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il3","Il4","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il9","Inpp5d","Irak2","Irak3","Irak4","Irf1","Irf4","Irf8","Itk","Jak2","Jak3","Lag3","Ltb","Mapk3","Mapkapk2","Mavs","Myd88","Nfkb1","Nlrp3","Nod1","Nod2","Pycard","Rel","Rela","Sele","Stat5b","Syk","Tcf7","Ticam1","Ticam2","Tigit","Tirap","Tlr1","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnfrsf11a","Tnfrsf1a","Tnfsf4","Tollip","Traf2","Traf6","Txk"),
  "Leukocyte_function"=c("Ccl19","Ccl25","Ccl4","Ccr1","Ccr7","Cd34","Cklf","Clec7a","Cxcl10","Cxcl12","Cxcl2","Cxcl5","Elane","F2rl1","Foxj1","Fut7","Hc","Icam1","Il16","Il23r","Il3","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgal","Itgam","Itgb1","Itgb2","Itgb3","Jam3","Lbp","Pecam1","Psen1","Psen2","Sele","Selplg","Syk","Tlr2","Vcam1"),
  "Macrophage_function"=c("C3ar1","Casp8","Ccl2","Ccl5","Ccr2","Ccr5","Ccr7","Cd1d1","Cd69","Cklf","Cmklr1","Crp","Csf1","Csf1r","Csf2","Cx3cl1","Cx3cr1","Eng","Fcer1a","Fcer2a","H60a","Hc","Il13","Il17f","Il18","Il1b","Il1rl1","Il23a","Il34","Il4","Il4ra","Itgb3","Lbp","Msr1","Nfkbia","Pparg","Prkce","Rora","Saa1","Sbno2","Syk","Thbs1","Tlr1","Ulbp1","Vegfa"),
  "Mast_cell_function"=c("C5ar1","Cd48","Fcer1a","Tpsab1"),
  "Mature_B_function"=c("Atm","Bcl10","Bcl6","Blk","Blnk","Bmi1","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd37","Cd74","Cd79b","Cxcl13","Fas","Fcgr2b","Gpr183","H2-Q10","Ikbkb","Ikbkg","Ikzf1","Il13","Il13ra1","Il1r2","Il21","Il2rg","Il5","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Nt5e","Prdm1","Sh2b2","Stat5b","Syk","Ticam1","Tirap","Tnfaip3","Tnfrsf4","Tnfsf4"),
  "Mature_T_function"=c("Anxa1","Bcl10","Card11","Ccl19","Ccl2","Ccl3","Ccr2","Ccr7","Cd1d2","Cd247","Cd274","Cd276","Cd48","Cd5","Cd59b","Cd80","Cd83","Cd86","Ctsh","Cxcl13","Eomes","Fasl","Fcgr4","Foxj1","Fut7","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Icam1","Ido1","Ifnar1","Ifnb1","Ikzf1","Il12a","Il12rb1","Il23a","Il2rg","Il6st","Il7r","Itgal","Itgam","Itgax","Itgb2","Itk","Lcp1","Mapk1","Mill2","Psen1","Psen2","Psmb10","Pvr","Pvrl2","Rag1","Rorc","Rps6","Spn","Stat5b","Syk","Tap1","Tcf7","Tgfb1","Tgfb2","Tigit","Tnfsf18","Tnfsf8","Traf2","Txk","Vcam1","Xcl1","Zap70"),
  "MHC1_MHC2"=c("Cd160","Cd1d1","Cd40lg","Cd74","Ciita","Ctsh","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Il10","Klrk1","Lag3","Mr1","Nlrc5","Pml","Tap1","Tap2","Tapbp"),
  "Microglial_function"=c("CASP1","CX3CR1","NOD2","TLR6","TLR7"),
  "NK_function"=c("Axl","Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Cd2","Cd244","Cd247","Cd7","Cd96","Flt3l","H2-M3","H60a","Ifi27","Ikzf1","Il11ra1","Il12a","Il12b","Il12rb1","Il15","Il15ra","Il21","Il23a","Il2rb","Itgb2","Klra1","Klra15","Klra17","Klra2","Klra20","Klra21","Klra27","Klra3","Klra4","Klra5","Klra6","Klra7","Klrb1c","Klrc1","Klrd1","Klrk1","Lag3","Mertk","Mill2","Ncam1","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Slamf7","Stat5b","Ulbp1"),
  "Pathogen_response"=c("Ccl2","Cd14","Fos","Hmgb1","Ifnb1","Ifng","Il10","Il12a","Il1b","Il6","Irak1","Irf3","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Rela","Ripk2","Tbk1","Ticam1","Tlr3","Tlr4","Tlr6","Tlr7","Tlr8","Tnf","Tnfrsf1a"),
  "Phagocytosis"=c("Anxa1","C3","Cd14","Cd36","Cd44","Cd47","Clec7a","Colec12","Crp","Csf1","Csf2","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Ifng","Il1rl1","Itgam","Itgb2","Marco","Mbl2","Mfge8","Mif","Myd88","Nod1","Pecam1","Siglec1","Ticam1","Tlr3","Tlr9","Tnf","Tnfsf11"),
  "Regulation_inflammatory_response"=c("Bcl6","C3","C3ar1","C4b","Ccl1","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl3","Ccl4","Ccl7","Ccl8","Ccr2","Ccr3","Ccr4","Ccr7","Cd14","Cd40","Cd40lg","Cebpb","Crp","Csf1","Cxcl1","Cxcl10","Cxcl11","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Fos","Il10","Il1r1","Il1rap","Il22","Il9","Itgb2","Ly96","Myd88","Nfkb1","Nos2","Ripk2","Sele","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr9","Tollip"),
  "Treg_function"=c("Ccr6","Foxp3","Ikzf2","Il9","Irf4","Irf8","Pou2f2","Rel","Tnfsf11"),
  "Senescence"=c("Abl1","Atm","Bmi1","Cdkn1a","Egr1","Ets1","Hras","Ifng","Igf1r","Irf3","Irf5","Irf7","Map2k1","Mapk14","Myc","Nfkb1","Plau","Prkcd","Serpinb2","Tgfb1","Trp53","Twist1"),
  "T_cell_function"=c("Ada","Anxa1","Bcl10","Bcl2","Bcl6","Btla","Card11","Casp3","Ccl11","Ccl19","Ccl2","Ccl3","Ccl5","Ccl7","Ccnd3","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Cd1d1","Cd1d2","Cd2","Cd247","Cd27","Cd274","Cd276","Cd28","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd40lg","Cd47","Cd48","Cd5","Cd59b","Cd74","Cd80","Cd83","Cd86","Cd8a","Cd8b1","Cebpb","Cma1","Crebbp","Csf2","Ctla4","Ctsh","Cxcl12","Cxcl13","Cxcr3","Cxcr4","Dpp4","Egr1","Eomes","Fas","Fasl","Fcgr4","Flt3","Foxj1","Foxp3","Fut7","Gata3","Gfi1","Gpr44","Gzmb","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Havcr2","Icam1","Icos","Icosl","Ido1","Ifnar1","Ifnb1","Ifng","Ikzf1","Ikzf2","Il10","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il15","Il17a","Il18","Il18r1","Il18rap","Il1b","Il1r1","Il2","Il21","Il23a","Il25","Il27","Il2ra","Il2rg","Il4","Il4ra","Il5","Il6","Il6st","Il7","Il7r","Il9","Irf1","Irf4","Irf8","Itch","Itga1","Itgal","Itgam","Itgax","Itgb2","Itk","Jak1","Jak2","Jak3","Lag3","Lck","Lcp1","Lgals3","Maf","Map3k7","Mapk1","Mapk8","Mill2","Nfatc1","Nfatc2","Nfkb1","Nos2","Pdcd1","Pdcd1lg2","Pou2f2","Psen1","Psen2","Psmb10","Ptprc","Pvr","Pvrl2","Rag1","Rel","Relb","Ripk2","Rora","Rorc","Rps6","Sell","Socs1","Socs3","Spn","Spp1","Stat1","Stat4","Stat5b","Stat6","Syk","Tap1","Tbx21","Tcf7","Tgfb1","Tgfb2","Tgfb3","Thy1","Tigit","Tlr4","Tlr6","Tmed1","Tnf","Tnfrsf13c","Tnfrsf14","Tnfrsf4","Tnfrsf8","Tnfsf11","Tnfsf13b","Tnfsf14","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf6","Trp53","Txk","Tyk2","Vcam1","Vegfa","Xcl1","Yy1","Zap70"),
  "T_cell_anergy"=c("Cma1","Ctla4","Gzmb","Icos","Itch","Itga1","Jak1","Lgals3","Pdcd1","Sell"),
  "T_cell_differentiation"=c("Ada","Bcl2","Cd1d1","Cd3d","Cd4","Cd74","Egr1","Flt3","Il27","Il7","Irf4","Nos2","Socs1"),
  "T_cell_proliferation"=c("Casp3","Ccnd3","Cd3e","Cxcl12","Cxcr4","Icosl","Il15","Il1b","Pdcd1lg2","Ptprc","Ripk2","Spp1","Tnfrsf13c","Tnfsf13b","Tnfsf14","Traf6","Trp53"),
  "T_cell_regulation"=c("Btla","Cd27","Cd3g","Cd47","Cd8a","Cd8b1","Dpp4","Fas","Foxp3","Icam1","Il2","Lag3","Lck","Map3k7","Tgfb1","Thy1","Tnfrsf14"),
  "Th1_Th2_differentiation"=c("Cd2","Cd28","Cd40","Cd40lg","Ifng","Il12b","Il18","Il2ra","Il4","Il6","Relb"),
  "Th1_function"=c("Ccr5","Csf2","Cxcr3","Havcr2","Il12rb2","Il18r1","Il18rap","Il2","Il27","Irf1","Nfkb1","Socs1","Stat1","Stat4","Tbx21","Tlr4","Tlr6","Tnf","Vegfa"),
  "Th2_function"=c("Bcl6","Ccl11","Ccl5","Ccl7","Ccr3","Ccr4","Cebpb","Gata3","Gfi1","Gpr44","Icos","Il10","Il13","Il13ra1","Il1r1","Il25","Il4ra","Il5","Jak1","Jak3","Maf","Nfatc1","Nfatc2","Stat6","Tmed1"),
  "Th17_function"=c("Il17a","Il21","Rora","Rorc"),
  "TLR"=c("Cd86","Gfi1","Irak1","Irak2","Irf3","Irf4","Map3k7","Mapkapk2","Myd88","Nfkbia","Prkce","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Traf3","Traf6"),
  "TNF_superfamily"=c("Cd27","Cd40","Cd40lg","Cd70","Fas","Fasl","Lta","Ltb","Ltbr","Tnf","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_members"=c("Cd40lg","Cd70","Fasl","Lta","Ltb","Tnf","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_receptors"=c("Cd27","Cd40","Fas","Ltbr","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9"),
  "Transcription_factors"=c("Atf1","Atf2","Batf","Bcl10","Bcl2","Bcl6","Card11","Card9","Cd34","Cd36","Cd40","Cdh1","Cebpb","Ciita","Clu","Cmklr1","Creb1","Creb5","Crebbp","Cyld","Ddx58","Ecsit","Egr1","Egr2","Elk1","Eomes","Ep300","Ets1","Fos","Foxj1","Foxp3","Gata3","Gfi1","Gtf3c1","Hck","Hmgb1","Icam1","Ikbkb","Ikbkg","Ikzf1","Il10","Il4","Il5","Irak1","Irak2","Irak3","Irf1","Irf2","Irf3","Irf4","Irf5","Irf7","Irf8","Itch","Itgb2","Jun","Maf","Mapk1","Mapk3","Mavs","Mef2c","Mnx1","Myc","Nfatc2","Nfatc3","Nfatc4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pax5","Pou2f2","Pparg","Pycard","Rel","Rela","Relb","Rora","Rorc","Runx1","Runx3","Sigirr","Smad4","Snai1","Stat1","Stat3","Stat4","Stat6","Tal1","Tbx21","Tcf7","Tfe3","Tfeb","Ticam1","Tlr3","Tlr9","Tmem173","Tnfrsf11a","Tnfrsf4","Tnfrsf8","Tnfsf18","Traf2","Traf3","Traf6","Twist1","Xbp1","Yy1"),
  "Transmembrane_transporter"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Apoe","App","Atg10","Atg16l1","Atg7","Axl","Bax","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd36","Cd3g","Cmah","Col1a1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fez1","Fyn","Icam1","Ifit1","Igf2r","Il13","Il1b","Il4","Itgb3","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Mertk","Msr1","Nefl","Nfatc1","Nup107","Pparg","Prkcd","Prkce","Psen1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Tmed1"),
  "Transporter_function"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Anxa1","Apoe","App","Atg10","Atg12","Atg16l1","Atg5","Atg7","Axl","Bax","C3","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd14","Cd36","Cd3g","Cd44","Cd47","Clec7a","Cmah","Col1a1","Colec12","Crp","Csf1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Fez1","Fyn","Icam1","Ifit1","Ifng","Igf2r","Il13","Il1b","Il1rl1","Il4","Itgam","Itgb2","Itgb3","Lamp1","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Marco","Mbl2","Mertk","Mfge8","Mif","Msr1","Myd88","Nefl","Nfatc1","Nod1","Nup107","Pecam1","Pparg","Prkcd","Prkce","Psen1","Siglec1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Ticam1","Tlr3","Tlr9","Tmed1","Tnf","Tnfsf11")
)

seu<-AddModuleScore(agg, 
                    features =hu.geneset.list, 
                    name=c("Adhesion","Adaptive_immunity","Antigen_processing",
                           "Apoptosis","Autophagy","B_activation","B_differentiation",
                           "B_function","B_proliferation","B_regulation",
                           "Bacterial_response","Basic_cell_function",
                           "Cancer_progression","CD4_T_function","CD_molecules",
                           "Cell_cycle","Chemokine_receptors","Chemokines",
                           "Chemokines_chemokine_receptors","Chemotaxis",
                           "Complement_pathway","Cytokine_receptors","Cytokines",
                           "Cytokines_cytokine_receptors","Cytotoxicity","DC_function",
                           "Humoral_immune_rersponse","Inflammation",
                           "Innate_immune_response","Interferon_response","Interleukins",
                           "Leukocyte_function","Macrophage_function","Mast_cell_function",
                           "Mature_B_function","Mature_T_function","MHC1_MHC2",
                           "Microglial_function","NK_function","Pathogen_response",
                           "Phagocytosis","Regulation_inflammatory_response",
                           "Treg_function","Senescence","T_cell_function",
                           "T_cell_anergy","T_cell_differentiation","T_cell_proliferation",
                           "T_cell_regulation","Th1_Th2_differentiation","Th1_function",
                           "Th2_function","Th17_function","TLR","TNF_superfamily",
                           "TNF_superfamily_members","TNF_superfamily_receptors",
                           "Transcription_factors","Transmembrane_transporter",
                           "Transporter_function") 
)

mMarkers <- immune.markers
# Filter
mMarkers <- mMarkers[ mMarkers$p_val_adj<0.01, ]
# Sort
mMarkers <- mMarkers[ order(mMarkers$cluster, -mMarkers$avg_logFC), ]

#Set your number of clusters
mGroups <- c("0":"15")

#Entrez4Enrichment- convert genes to EntrezID for GO conversion
for (group in mGroups)
{
  GeneSymbol4Enrichment<-mMarkers[mMarkers$cluster==group,] %>% group_by(cluster) %>% top_n(n = NTop, wt = avg_logFC)
  GeneSymbol4Enrichment.group<-GeneSymbol4Enrichment$gene
  converted <- bitr(GeneSymbol4Enrichment.group, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  Entrez4Enrichment <- converted$ENTREZID
  assign(paste0("mEntrez4Enrichment", group), converted$ENTREZID)
}


# background gene selection
NTopBkgGenes<-2000
mExpressionVector<-rowMeans(as.matrix(GetAssayData(immune.combined, slot = "data")))
mExpressionVector<-sort(mExpressionVector, decreasing = T)
mGeneSymbol4Bkg<-names(mExpressionVector[1:NTopBkgGenes])
mconverted <- bitr(mGeneSymbol4Bkg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
mEntrez4Bkg <- mconverted$ENTREZID

#choose the gene ontologies you want to check
Onts <- c("CC","MF","BP")

for (group in mGroups)
{
  Entrez_counter <- paste("mEntrez4Enrichment", group, sep = "")
  Entrez_counter <- eval(parse(text = Entrez_counter))
  for (ont in Onts)
  {
    print("Working on")
    print(group)
    print(ont)
    enrich.res<-enrichGO(gene= Entrez_counter, universe = mEntrez4Bkg,
                         ont=ont,OrgDb='org.Mm.eg.db' ,readable=T,
                         minGSSize = 20, maxGSSize = 500, qvalueCutoff = 0.05)
    if(is.null(enrich.res) == T){
      
    }
    else{
      print("computed")
      assign(paste0("enrich.res", group, ont), enrich.res)
    }
  }
}


imp.neut1 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 0 & gemgroup == 1)
imp.neut2 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 0 & gemgroup == 2)
imp.neut3 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 0 & gemgroup == 3)

imp.t1 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 3 & gemgroup == 1)
imp.t2 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 3 & gemgroup == 2)
imp.t3 <- subset(immune.combined, stim == "IMP" & seurat_clusters == 3 & gemgroup == 3)

spo.neut1 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 0 & gemgroup == 1)
spo.neut2 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 0 & gemgroup == 2)
spo.neut3 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 0 & gemgroup == 3)

spo.t1 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 3 & gemgroup == 1)
spo.t2 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 3 & gemgroup == 2)
spo.t3 <- subset(immune.combined, stim == "SPO" & seurat_clusters == 3 & gemgroup == 3)


dotplot(enrich.res0MF, showCategory=20) +ggtitle("Significantly Enriched Pathways")

enrich.res<-enrichGO(gene= mEntrez4Enrichment1, universe = mEntrez4Bkg,
                     ont="MF",OrgDb='org.Mm.eg.db' ,readable=T,
                     minGSSize = 20, maxGSSize = 500, qvalueCutoff = 0.05)



imp1 <- subset(imp, seurat_clusters == 1)
imp2 <- subset(imp, seurat_clusters == 2)
shan.imp <- subset(immune.combined, stim == "IMP")
shan.spo <- subset(immune.combined, stim == "SPO")
diversity(shan.imp$nFeature_RNA, index = "simpson" )
diversity(shan.spo$nFeature_RNA, index = "simpson")
diversity(shan.imp$nCount_RNA)
diversity(shan.spo$nCount_RNA)
specnumber(immune.combined$nFeature_SCT, groups = immune.combined$seurat_clusters)
imp.0 <- subset(immune.combined, immune.combined$seurat_clusters == 0 & stim == "IMP")
imp.1 <- subset(immune.combined, immune.combined$seurat_clusters == 1 & stim == "IMP")
imp.2 <- subset(immune.combined, immune.combined$seurat_clusters == 2 & stim == "IMP")
imp.3 <- subset(immune.combined, immune.combined$seurat_clusters == 3 & stim == "IMP")
imp.4 <- subset(immune.combined, immune.combined$seurat_clusters == 4 & stim == "IMP")

spo.0 <- subset(immune.combined, immune.combined$seurat_clusters == 0 & stim == "SPO")
spo.1 <- subset(immune.combined, immune.combined$seurat_clusters == 1 & stim == "SPO")
spo.2 <- subset(immune.combined, immune.combined$seurat_clusters == 2 & stim == "SPO")
spo.3 <- subset(immune.combined, immune.combined$seurat_clusters == 3 & stim == "SPO")
spo.4 <- subset(immune.combined, immune.combined$seurat_clusters == 4 & stim == "SPO")


length(imp.0$nFeature_RNA)



M <- as.table(rbind(c(2843, 1039, 686, 805, 120), c(200, 816, 289, 63, 19)))
dimnames(M) <- list(model = c("Implanted", "Spontaneous"),
                    cluster = c("0","1", "2", "3", "4"))
(Xsq <- chisq.test(M))  # Prints test summary


p1.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 1)
p1.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 1)
p1.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 1)
p1.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 1)
p1.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 1)
p1.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 1)
p1.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 1)

p2.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 2)
p2.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 2)
p2.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 2)
p2.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 2)
p2.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 2)
p2.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 2)
p2.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 2)

p3.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 3)
p3.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 3)
p3.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 3)
p3.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 3)
p3.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 3)
p3.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 3)
p3.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 3)

p4.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 4)
p4.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 4)
p4.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 4)
p4.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 4)
p4.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 4)
p4.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 4)
p4.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 4)

p5.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 5)
p5.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 5)
p5.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 5)
p5.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 5)
p5.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 5)
p5.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 5)
p5.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 5)

p6.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 6)
p6.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 6)
p6.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 6)
p6.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 6)
p6.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 6)
p6.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 6)
p6.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 6)

p7.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 7)
p7.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 7)
p7.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 7)
p7.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 7)
p7.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 7)
p7.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 7)
p7.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 7)

p8.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 8)
p8.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 8)
p8.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 8)
p8.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 8)
p8.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 8)
p8.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 8)
p8.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 8)

p9.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 9)
p9.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 9)
p9.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 9)
p9.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 9)
p9.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 9)
p9.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 9)
p9.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 9)

p10.0 <- subset(agg, agg$seurat_clusters == 0 & gemgroup == 10)
p10.1 <- subset(agg, agg$seurat_clusters == 1 & gemgroup == 10)
p10.2 <- subset(agg, agg$seurat_clusters == 2 & gemgroup == 10)
p10.3 <- subset(agg, agg$seurat_clusters == 3 & gemgroup == 10)
p10.4 <- subset(agg, agg$seurat_clusters == 4 & gemgroup == 10)
p10.5 <- subset(agg, agg$seurat_clusters == 5 & gemgroup == 10)
p10.6 <- subset(agg, agg$seurat_clusters == 6 & gemgroup == 10)

m1.imp <- subset(imp, gemgroup == 1)
m2.imp <- subset(imp, gemgroup == 2)
m3.imp <- subset(imp, gemgroup == 3)

m1.spo <- subset(spo, gemgroup == 1)
m2.spo <- subset(spo, gemgroup == 2)
m3.spo <- subset(spo, gemgroup == 3)


