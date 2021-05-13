if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")
BiocManager::install("cytofkit")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('satijalab/seurat-wrappers')
install.packages("vegan")

library(Seurat)
library(dplyr)
library(cowplot)
library(reticulate)
library(ggplot2)
set.seed(42)
library(monocle3)
library(cytofkit)
library(SeuratWrappers)
library(vegan)


#cds <- load_mm_data(mat_path = "~/Documents/Thesis/IMPx/IMPx5/Matrices/ImpAgg/matrix.mtx", 
#                    feature_anno_path = "~/Documents/Thesis/IMPx/IMPx5/Matrices/ImpAgg/features.tsv", 
#                    cell_anno_path = "~/Documents/Thesis/IMPx/IMPx5/Matrices/ImpAgg/barcodes.tsv")
imp.mono <- load_cellranger_data("~/Documents/Thesis/IMPx/IMPx5/ImpAgg/")
spo.mono <- load_cellranger_data("~/Documents/Thesis/IMPx/IMPx5/SpoAgg/")
hu.mono <- load_cellranger_data("~/Documents/Thesis/IMPx/IMPx5/HuAGG3/")
comb.mono <- combine_cds(list(imp.mono,spo.mono))

imp.mono <- estimate_size_factors(imp.mono)
spo.mono <- estimate_size_factors(spo.mono)
hu.mono <- estimate_size_factors(hu.mono)
comb.mono <- estimate_size_factors(comb.mono)

imp.mono <- preprocess_cds(imp.mono, num_dim = 100)
spo.mono <- preprocess_cds(spo.mono, num_dim = 100)
hu.mono <- preprocess_cds(hu.mono, num_dim = 100)
comb.mono <- preprocess_cds(comb.mono, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
imp.mono <- align_cds(imp.mono, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
imp.mono <- reduce_dimension(imp.mono)
spo.mono <- reduce_dimension(spo.mono)
hu.mono <- reduce_dimension(hu.mono)
comb.mono <- reduce_dimension(comb.mono)


## Step 4: Cluster the cells
imp.mono <- cluster_cells(imp.mono)
plot_cells(imp.mono)
spo.mono <- cluster_cells(spo.mono)
plot_cells(spo.mono)
hu.mono <- cluster_cells(hu.mono)
plot_cells(hu.mono)
tru.hu.mono <- cluster_cells(tru.hu.mono)
plot_cells(tru.hu.mono, show_trajectory_graph = F)
comb.mono <- cluster_cells(comb.mono)
plot_cells(comb.mono)
tru.comb.mono <- cluster_cells(tru.comb.mono)
plot_cells(tru.comb.mono, show_trajectory_graph = F)

## Step 5: Learn a graph
imp.mono <- learn_graph(imp.mono)
spo.mono <- learn_graph(spo.mono)
hu.mono <- learn_graph(hu.mono)
comb.mono <- learn_graph(comb.mono)
tru.hu.mono <- learn_graph(tru.hu.mono)
tru.comb.mono <- learn_graph(tru.comb.mono)

## Step 6: Order cells
imp.mono <- order_cells(imp.mono)
spo.mono <- order_cells(spo.mono)
hu.mono <- order_cells(hu.mono)
comb.mono <- order_cells(comb.mono)
tru.hu.mono <- order_cells(tru.hu.mono)
tru.comb.mono <- order_cells(tru.comb.mono)

plot_cells(imp.mono)
plot_cells(spo.mono)
plot_cells(hu.mono)
plot_cells(comb.mono)
plot_cells(tru.hu.mono)
plot_cells(tru.comb.mono)


#DEG's
imp_marker <- top_markers(imp.mono, group_cells_by="cluster", cores=8)

imp_specific_markers <- imp_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

spo_marker <- top_markers(spo.mono, group_cells_by="cluster", cores=8)

spo_specific_markers <- spo_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

hu_marker <- top_markers(hu.mono, group_cells_by="cluster", cores=8)

hu_specific_markers <- hu_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

tru_hu_marker <- top_markers(tru.hu.mono, group_cells_by="cluster", cores=8)

tru_hu_specific_markers <- tru_hu_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

comb_marker <- top_markers(comb.mono, group_cells_by="cluster", cores=8)

comb_specific_markers <- comb_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

tru_comb_marker <- top_markers(tru.comb.mono, group_cells_by="cluster", cores=8)

tru_comb_specific_markers <- tru_comb_marker %>% filter(fraction_expressing >= 0.10) %>% 
  group_by(cell_group) %>% top_n(3, pseudo_R2)

HSCvsCPSup <- list("AATF", "ACTL7A", "ACTR6", "ADAM30", "ADGRL4", "AFF4", "AHSG",	"ALB", "ANGPTL8",
                   "APOH", "APOL1",	"ARF3",	"ARMC1", "ARNTL", "ATAD2", "AUNIP", "BARX2", "BNC1",
                   "BRIX1", "BZW1", "C15orf39", "C1GALT1C1", "C1QL1", "CACNA1H", "CASP1", "CD2",
                   "CD247", "CD40LG", "CD59", "CDC23", "CDC27",	"CDCA3", "CDK5R1", "CEBPD",	
                   "CERNA1", "CLN3", "CNPY4",	"COPE",	"CPNE6",	"CREB3L2",	"CRIM1",	"CRYL1",
                   "CSF1R",	"CSTF1",	"CXCL2",	"CXCR6",	"CYLC1",	"DCSTAMP",	"DDX31",	"DDX39A",
                   "DGAT1",	"DHCR7",	"DOC2A", "DOP1A",	"DPF3",	"DPPA4", "DUSP6",	"ELL",
                   "ENO1",	"ESR1",	"ETFDH",	"EXT1",	"EXTL1",	"FA2H",	"FAM8A1",	"FBN2",	
                   "FLT4",	"FOSL1",	"FSD1",	"FUT8",	"GIMAP5",	"GLRX",	"GORASP2",	"GP9",
                   "GUCY1A2",	"H3-3B", "HBE1",	"HBEGF",	"HOPX",	"HOXB3",	"HOXB6",	"HSD17B14",
                   "HSPA1A",	"HSPA5",	"HSPA9",	"IDH3G",	"IFT57",	"IL6R",	"IMP3",	"INAVA",
                   "INTS11", "ITFG1",	"KATNBL1",	"KCNC4",	"KLHL1",	"KLRG1",	"KRT75",	"LAMA5",
                   "LILRA4",	"LINC01558",	"LIPG",	"MAGEB3",	"MAPK13",	"MICA",	"MIR22HG",	
                   "MKNK1",	"MMRN2", "MON1B",	"MRPL49",	"MST1R",	"MTMR11",	"MTNR1B",	"MYBPC1",
                   "MYDGF",	"MYLK",	"NAA38",	"NCK2",	"NDNF",	"NELFCD",	"NEU1",	"NFYA",	"NGB",
                   "NINJ2",	"NOVA2",	"NSUN5",	"NTRK2",	"NUDT6",	"OSBPL7",	"P2RY14",	"P3H2",
                   "PADI1",	"PAIP1",	"PCID2",	"PDZD2",	"PEA15",	"PITPNA",	"PLCB1",	"PLEKHF1",
                   "PPIH",	"PRR7",	"PRSS2",	"PRSS8",	"PRTN3",	"PSG4",	"PSMD11",	"PTH2R",	
                   "PYGO1",	"RAB33B",	"RAD1",	"RANGRF",	"RMDN3",	"RNF126",	"ROBO1",	"RPH3AL",
                   "RTF2",	"SAMD4A",	"SEMG1",	"SERPINA4",	"SERPINE1",	"SERPING1",	"SH2D1A",	
                   "SIGMAR1",	"SLC22A4",	"SLC2A1",	"SLC5A4",	"SLC9A1",	"SLCO3A1",	"SNRNP40",
                   "TBC1D31",	"TCN2",	"TDO2",	"TEF",	"TEK",	"TEKT2",	"TENM1",	"TFE3",	"THRAP3",
                   "TLR3",	"TMCO1",	"TMEM30B",	"TMEM51",	"TRPV6",	"TSPAN15",	"TVP23B",	"UFD1",
                   "UGGT1",	"ULBP2",	"URM1",	"VKORC1",	"WASF3",	"WSB2",	"XRCC2",	"YES1",	"ZBBX",
                   "ZBTB16",	"ZNF205", "ADAR",	"ANKRD12",	"APOA2",	"ARHGAP18",	"ARMC5",	
                   "ATP10A",	"BCL9",	"BDKRB2",	"CAMKK2",	"CASP6",	"CBR3",	"CBX1",	"CCDC171",
                   "CCR8",	"CD200",	"CD38",	"CD5",	"CD6",	"CEBPZ",	"CEP152",	
                   "CEP97",	"CHST15",	"CNPY3",	"CROT",	"CTLA4",	"CYP51A1",	"DCTN4",	"DDX23",
                   "DGKZ",	"DHRS3",	"DNMT3A",	"DUSP1",	"E2F8",	"EMC3",	"EOMES",
                   "EPHX1",	"EWSR1",	"FAM53B",	"FBH1",	"FBXL19",	"FMNL1",	"FURIN",	"FZD10",
                   "GET4",	"GFI1",	"GIT2",	"GML",	"GPR83",	"GPR88",	"GRAMD1B",	"GRSF1",	
                   "GSTO1",	"GULOP", "H2AC15",	"HACL1",	"HIVEP3",	"HOOK1",	"HSPA4L",	"HTR2C",
                   "IBA57",	"ID2",	"IKZF2",	"IL2RB",	"ITGA9",	"ITIH5",	"ITM2A",	"IZUMO1R",
                   "JARID2",	"KCNF1",	"KCTD12",	"KDM4A",	"KLB",	"KRTAP6-2",	"LPCAT1",	"LRMP",
                   "LRRK2",	"LYSMD1",	"MAN1C1",	"MAP3K2",	"MAP4K4",	"MAP4K5",	"MAPK6",	"MAPKAPK2",
                   "MARCKS",	"MCOLN1",	"MIR15B",	"MIR592",	"MIRLET7G",	"MYB",	"NAA35",	"NOD1",
                   "NOTCH1",	"NR4A1",	"NUDT13",	"OBI1",	"ORAI2",	"PCED1A",	"PFKFB3",	"PFKFB4",
                   "PHACTR2",	"PIK3CD",	"PLCG1",	"PLSCR1",	"PTK2B",	"PTPN12",	"PTPN23",	
                   "PTTG1IP",	"PUS3",	"RAB13",	"RABL3",	"RANBP3L",	"REPS1",	"RFC1",	"RGS16",
                   "RGS3",	"RNF149",	"RNF2",	"RSF1",	"SBK2",	"SBNO2",	"SENP7",	"SERP1",	
                   "SGPP2",	"SH3YL1",	"SHISAL2B",	"SKI",	"SLC13A4",	"SLC9B2",	"SMC2",	"SMC6",
                   "SMPD3",	"SPAG1",	"SPAG17",	"SSBP3",	"ST3GAL4",	"ST8SIA1",	"STAT5B",	"SUZ12",
                   "SWAP70",	"SYPL1",	"SYT11",	"TAF5",	"TBX1",	"TGFBR3",	"TMEM248",	"TNFRSF4",
                   "TNFRSF9",	"TNFSF9",	"TNIK",	"TNIP2",	"TOP2A",	"TRIM9",	"TSC22D2",	"TTC3",
                   "URI1",	"UTF1",	"XCL1",	"XKRX",	"ZFP37",	"ZNF318",	"ZNRF1",	"ZRANB1")
group <- list("1")
HSC.up.df <- stack(setNames(HSCvsCPSup, group))
HSC.up.df$ind<- "1"

#top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_cells(tru.comb.mono, show_trajectory_graph = F, genes = list("Cd3d", "Cd3e", "Cd3g", "Klrd1", "Nkg7",
                                                              "Nktr", "Cd79a", "Cd79b", "Ms4a7"))
plot_cells(tru.comb.mono, show_trajectory_graph = F, genes = list("Cd14", "Grn", "C1qa", "Itgam", "Pirb",
                                                              "Il1b", "Mrc1", "Itgax", "S100a8"))
plot_cells(hu.mono, show_trajectory_graph = F, genes = list("CD3D", "CD3E", "CD3G", "KLRD1", "NKG7",
                                                              "NKTR", "CD79A", "CD79B", "MS4A7"))
plot_cells(tru.hu.mono, show_trajectory_graph = F, genes = list("CD14", "GRN", "C1QA", "ITGAM", "LILRB1",
                                                              "IL1B", "MRC1", "ITGAX", "S100A8"))

plot_cells(tru.hu.mono, show_trajectory_graph = F, genes = list("CD14", "CD2", "NKTR", "HCST"))
plot_cells(tru.comb.mono, show_trajectory_graph = F, genes = list("Cd14", "Cd2", "Nktr", "Hcst"))

tru.hu.mono <- choose_cells(hu.mono, reduction_method = "UMAP")
tru.comb.mono <- choose_cells(comb.mono, reduction_method = "UMAP")

plot_cells(hu.mono, show_trajectory_graph = F, genes = HSC.up.df)

plot_cells(tru.hu.mono)
HSC.up.df <- rename(HSC.up.df, gene_short_name =values , group = ind)

imp.cds <- as.cell_data_set(imp)

## Calculate size factors using built-in function in monocle3
imp.cds <- estimate_size_factors(imp.cds)

## Add gene names into CDS
imp.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(imp[["RNA"]])
colnames(imp.cds)[2] <- "gene_short_name"

imp.cds <- preprocess_cds(imp.cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
imp.cds <- reduce_dimension(imp.cds)

## Step 4: Cluster the cells
imp.cds <- cluster_cells(imp.cds)

# With regression:
gene_fits <- fit_models(cds, model_formula_str = "~embryo.time")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")
emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)

# With graph autocorrelation:
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))

#DEG's
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                              cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

colData(cds)

## Step 5: Learn a graph
imp.cds <- learn_graph(imp.cds)

## Step 6: Order cells
imp.cds <- order_cells(imp.cds)

plot_cells(imp.cds)

spo.cds <- as.cell_data_set(spo)

## Calculate size factors using built-in function in monocle3
spo.cds <- estimate_size_factors(spo.cds)

## Add gene names into CDS
spo.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(spo[["RNA"]])
colnames(spo.cds)[2] <- "gene_short_name"

spo.cds <- preprocess_cds(spo.cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
spo.cds <- reduce_dimension(spo.cds)

## Step 4: Cluster the cells
spo.cds <- cluster_cells(spo.cds)
plot_cells(spo.cds)

plot_cells(spo.cds, color_cells_by = "gemgroup")

## Step 5: Learn a graph
spo.cds <- learn_graph(spo.cds)

## Step 6: Order cells
spo.cds <- order_cells(spo.cds)
plot_cells(spo.cds)

imp.com <- as.cell_data_set(imp)
## Add gene names into CDS
imp.com@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(imp[["RNA"]])
colnames(imp.com)[2] <- "gene_short_name"
spo.com <- as.cell_data_set(spo)
## Add gene names into CDS
spo.com@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(spo[["RNA"]])
colnames(spo.com)[2] <- "gene_short_name"
mo_combined_cds <- combine_cds(list(imp.com, spo.com))

## Calculate size factors using built-in function in monocle3
mo_combined_cds <- estimate_size_factors(mo_combined_cds)


agg.cds <- as.cell_data_set(agg)

## Calculate size factors using built-in function in monocle3
agg.cds <- estimate_size_factors(agg.cds)

## Add gene names into CDS
agg.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(agg[["RNA"]])
colnames(agg.cds)[2] <- "gene_short_name"

agg.cds <- preprocess_cds(agg.cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
agg.cds <- reduce_dimension(agg.cds)

## Step 4: Cluster the cells
agg.cds <- cluster_cells(agg.cds)

plot_cells(agg.cds)

plot_cells(agg.cds, color_cells_by = "gemgroup")

## Step 5: Learn a graph
agg.cds <- learn_graph(agg.cds)

## Step 6: Order cells
agg.cds <- order_cells(agg.cds)
plot_cells(agg.cds)



