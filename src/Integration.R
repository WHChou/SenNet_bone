library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(repr)
library(circlize)
library(ComplexHeatmap)
library(grid)

obj_list <- readRDS("RDS/bm_merge_inhouse0805_TanGrimesRef_v4_SCT_reMerge_obj_list.rds")

# Try using the union of variable features as integration anchor
var_genes <- lapply(X = obj_list, FUN = VariableFeatures)
kt_var <- var_genes$`Kai Tan reference` 
ortho_var <- var_genes$`In-house orthopedic` 
pro_var <- var_genes$`Grimes reference`
a <- union(kt_var, ortho_var)
all_var <- union(a, pro_var)


# Integrate
#features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 4000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = all_var)
obj_list <- lapply(X = obj_list, FUN = RunPCA, features = all_var)
rna.anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = all_var, dims = 1:30, reduction = "rpca", k.anchor = 5)
combined <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30, k.weight=100)

# Load genes
gs <- readRDS('/diskmnt/Projects/myeloma_scRNA_analysis/MMY_IRD/analysis/resources/gene_sets.rds')
s.genes <- gs$s.genes
g2m.genes <- gs$g2m.genes
ig_genes <- gs$ig_genes #length=436

# remove Ig genes from var.features for PCA
var_genes <- combined@assays$integrated@var.features #4000
var_genes_nomatch <- !var_genes %in% ig_genes
var_genes_noIG <- var_genes[var_genes_nomatch] #3974
combined <- RunPCA(combined, npcs=50, features=var_genes_noIG, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5, algorithm=4, method='igraph')

# Plot RPCA integration UMAPs
pdf("scRNA/bm_merge_inhouse0805_TanGrimesRef_v4_2_SCT_RPCAIntegration_PC30_UMAPs.pdf", width = 8, height = 6)
print(DimPlot(combined, group.by = "seurat_clusters", label = TRUE))
print(DimPlot(combined, group.by = "subcluster_anno_RPCA_WC", label = TRUE) + NoLegend())
print(DimPlot(combined, group.by = "Sample_ID", label = TRUE))
print(DimPlot(combined, group.by = "data_source", label = TRUE))
print(DimPlot(combined, group.by = "data_type_rna", label = TRUE))
print(DimPlot(combined, group.by = "Site_of_origin", label = TRUE))
dev.off()
pdf("scRNA/bm_merge_inhouse0805_TanGrimesRef_v4_2_SCT_RPCAIntegration_PC30_UMAP_origIdent.pdf", width = 20, height = 8)
print(DimPlot(combined, group.by = "orig_anno", label = TRUE))
dev.off()

pdf("scRNA/bm_merge_inhouse0805_TanGrimesRef_v4_1_SCT_RPCAIntegration_PC30cellType_dotPlot.pdf", width = 12, height = 8)
print(DotPlot(combined, group.by = "subcluster_anno_RPCA_WC", features = c("CD34", "AVP", "CRHBP", "SMIM24", "SPINK2", 
                                                                "ELANE", "MPO", "CTSG", "AZU1", "CAMP", "LTF", "MMP9", "AQP9", 
                                                                 "CD14", 'S100A9', 'SELL', 'VCAN', #CD14 monocyte
                                                              'FCGR3A', 'MS4A7', 'TNFRSF1B', #CD16
                                                                 'CD68', 'AIF1', 'IFITM2', 'LST1',  #MAC
                                                             'C1QA', 'C1QB', 'C1QC', 
                                                                 'CLEC4C', 'GZMB', 'IL3RA', 'IRF8', #pDC
                                                             'CLEC10A', 'CD1C','FCER1A')) + RotatedAxis()) #cDC2
print(DotPlot(combined, group.by = "subcluster_anno_RPCA_WC", features = c("CD34", "AVP", "CRHBP", "SMIM24", "SPINK2", "DNTT", "VPREB1", "VPREB3",
                                                               "MS4A1", "CD19", "BANK1", "SDC1", "JCHAIN", "MDC1", "TNFRSF17",
                                                               "PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "CD8B", 
                                                               "CCR7", "SELL", "IL7R", "LEF1", "GATA3", "FOXP1", # naive?
                                                               "TNF", "GZMA", "GZMB", "GZMK", "PRF1", "LAMP1", "FASLG",# cytotoxic
                                                               "NCAM1", "FCGR3A", "CD69")) + RotatedAxis()) # CD56 NK
print(DotPlot(combined, group.by = "subcluster_anno_RPCA_WC", features = c("CD34", "AVP", "CRHBP", "SMIM24", "SPINK2", 
                                                               "FCER1A", "TIMP3", "SLC40A1", "PKIG",  #MEP from Reyka
                                                               "APOC1", "BLVRB", "CA1", "HBB", "KLF1", "TFR2", #ERP from Reyka
                                                               "TFRC", "GATA1", "GATA2", "AHSP", "HEMGN", "HBD", "PLEK", "PF4", "ITGA2B")) + RotatedAxis())
print(DotPlot(combined, group.by = "subcluster_anno_RPCA_WC", features = c("ACTA2", "MYH11", "RGS5", "NOTCH3",   ## vSMC/pericyte
                                 "DCN", "LUM", "COL1A1", "VIM", "PDGFRA", "PDGFRB",  ## Fibro
                                 "PECAM1", "VWF", "FLT1", "SEMA3A", "FLT4",  ## endothelial
                                 "ADIPOQ", "LPL", "PPARG", "PLIN4",
                                  "SPP1", "RUNX2", "BGLAP",
                                 "LEPR", "CXCL12", "KITLG", "NT5G")) + RotatedAxis())
dev.off()



saveRDS(combined, 'RDS/bm_merge_inhouse0805_TanGrimesRef_v4_SCT_reMerge_2_RPCAIntegration.rds', compress=F)
