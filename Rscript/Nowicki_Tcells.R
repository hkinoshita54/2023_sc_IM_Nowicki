# load packages----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(biomaRt)
library(msigdbr)
library(fgsea)
library(DESeq2)
library(speckle)
library(miloR)
library(SingleCellExperiment)
library(scater)

# Load data----
# https://cellxgene.cziscience.com/collections, UMAP of T cells, Seurat v4
seu <- readRDS("RDSfiles/Tcells.rds")
seu
View(seu@meta.data)
DefaultDimReduc(seu)
head(rownames(seu)) # ENSG instead of gene symbol

# convert ENSG to gene symbol, remove features without gene symbol----
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
lookup <- getBM(mart = mart, attributes = c('hgnc_symbol','ensembl_gene_id'),uniqueRows = TRUE)
rnames <- lookup$hgnc_symbol[match(rownames(seu), lookup$ensembl_gene_id)]
rnames[rnames == ""] <- NA
rownames(seu@assays$RNA@counts) <- rnames
rownames(seu@assays$RNA@data) <- rnames
features_keep <- rownames(seu)[!is.na(rownames(seu))]
seu <- subset(seu, features = features_keep)
rownames(seu@assays$RNA@meta.features) <- rownames(seu)
seu

# explore the dataset----
DimPlot(seu, group.by = "Celltypes_global") + NoAxes() 
DimPlot(seu, group.by = "cell_type") + NoAxes()
DimPlot(seu, group.by = "Detailed_Cell_Type") + NoAxes()
DimPlot(seu, group.by = "tissue") + NoAxes()
DimPlot(seu, group.by = "Tissue_in_paper") + NoAxes()
DimPlot(seu, group.by = "disease") + NoAxes() # normal includes various tissues
DimPlot(seu, group.by = "donor_id") + NoAxes()
table(seu$tissue) # body of stomach includes 10500 cells
table(seu$Tissue_in_paper) # CAG, GIM, NAG and NGB consist body of stomach
Idents(seu) <- "cell_type"

FeaturePlot(seu,features = "CD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend() 
FeaturePlot(seu,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "SELL", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# features are not well clustered... just use the authors annotation for DA analysis

saveRDS(seu, file = "RDSfiles/Tcells_annotated.RDS")

# subset body of stomach----
seu_bs <- subset(seu, tissue == "body of stomach")
Idents(seu_bs) <- "cell_type"
table(seu_bs$donor_id)

# Rank patients by epithelila ARID1A level----
# check average ARID1A in epithelial cells in "avg_bs" in excel
Idents(seu_bs) <- "donor_id"
ARID1A_H <- WhichCells(seu_bs, ident = c("Patient17",	"Patient14",	"Patient22",	"Patient16",	"Patient21",	"P5866",	"P5846",	"L05",	"L09",	"L07",	"L04"))
ARID1A_L <- WhichCells(seu_bs, ident = c("L08",	"L06",	"P6207",	"L03",	"L02",	"P5931",	"P6342",	"P6592",	"L01",	"P6709",	"P6649"))
seu_bs <- SetIdent(seu_bs, cells = ARID1A_H, value = "ARID1A_H")
seu_bs$ARID1A_level <- Idents(seu_bs)
seu_bs <- SetIdent(seu_bs, cells = ARID1A_L, value = "ARID1A_L")
seu_bs$ARID1A_level <- Idents(seu_bs)
DimPlot(seu_bs, group.by = "ARID1A_level") + NoAxes()

# DA analysis by speckle----
seu <- seu_bs
propres <- propeller(clusters = seu$cell_type, sample = seu$donor_id, group = seu$ARID1A_level)
write.table(propres, file = "results/DA//propeller_Tcells_ARID1AHvsL.txt", sep ="\t", col.names = T,row.names = F)
plotCellTypeProps(clusters = seu$cell_type, sample = seu$donor_id)

# remove samples which contain less than 100 cells
seu <- subset(seu, donor_id %in% c("L01", "L02", "L07", "L09", "P6342", "Patient17", "Patient22"), invert = TRUE)
table(seu$donor_id)
propres <- propeller(clusters = seu$cell_type, sample = seu$donor_id, group = seu$ARID1A_level)
write.table(propres, file = "results/DA//propeller_Tcells_ARID1AHvsL_2.txt", sep ="\t", col.names = T,row.names = F)
plotCellTypeProps(clusters = seu$cell_type, sample = seu$donor_id)

# DA analysis by miloR----
# seu <- JoinLayers(seu)
# options(Seurat.object.assay.version = "v3")
# seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(seu)
milo <- Milo(sce)
milo <- buildGraph(milo, k = 100, d = 30, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.2, k = 100, d=30, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo)

milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples="donor_id")
head(nhoodCounts(milo))

design <- data.frame(colData(milo))[,c("donor_id", "ARID1A_level")]
design <- distinct(design)
rownames(design) <- design$donor_id
design <- design[colnames(nhoodCounts(milo)), , drop=FALSE]
design

milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")
rownames(design) <- design$donor_id
contrast.1 <- c("ARID1A_levelARID1A_L - ARID1A_levelARID1A_H")
da_results <- testNhoods(milo, design = ~ 0 + ARID1A_level, design.df = design, model.contrasts = contrast.1)
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

milo <- buildNhoodGraph(milo)

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)


milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = ".umap_MinDist_0.2_N_Neighbors_15", colour_by="ARID1A_level", text_by = "cell_type", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout=".umap_MinDist_0.2_N_Neighbors_15",alpha=0.1) 

umap_pl + nh_graph_pl + plot_layout(guides="collect")

da_results <- annotateNhoods(milo, da_results, coldata_col = "cell_type")
head(da_results)
write.table(da_results, file = "results/DA//miloR_res_imm_t_DF_annotated_ARID1AHvsL.txt", sep ="\t", col.names = T,row.names = F)
plotDAbeeswarm(da_results, group.by = "cell_type")

