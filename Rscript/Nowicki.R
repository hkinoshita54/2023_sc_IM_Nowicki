# load packages----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(biomaRt)
library(msigdbr)
library(fgsea)
library(DESeq2)

# Load data----
# https://cellxgene.cziscience.com/collections, UMAP of Columnar cells, Seurat v4
seu <- readRDS("RDSfiles/all_columnar_cells.RDS")
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
DimPlot(seu, group.by = "Celltypes_global") + NoAxes() # identical to Figure 3A in the original paper
DimPlot(seu, group.by = "cell_type") + NoAxes()
DimPlot(seu, group.by = "Columnar_clusters") + NoAxes()
DimPlot(seu, group.by = "tissue") + NoAxes()
DimPlot(seu, group.by = "Tissue_in_paper") + NoAxes()
DimPlot(seu, group.by = "disease") + NoAxes() # normal includes various tissues
DimPlot(seu, group.by = "donor_id") + NoAxes()
table(seu$tissue) # body of stomach includes 41003 cells
table(seu$Tissue_in_paper) # CAG, GIM, NAG and NGB consist body of stomach
Idents(seu) <- "Celltypes_global"

FeaturePlot(seu,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "GHRL", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu, file = "RDSfiles/all_columnar_cells_annotated.RDS")

# subset body of stomach----
seu_bs <- subset(seu, tissue == "body of stomach")
Idents(seu_bs) <- "Celltypes_global"
table(seu_bs$donor_id)
FeaturePlot(seu_bs,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# Rank patients by MUC6 level----
avg_bs<-AverageExpression(seu_bs, group.by = "donor_id" ) 
write.table(as.matrix(avg_bs$RNA), "results/GeneExpression/avg_bs.txt", sep="\t",col.names = T, row.names = T)  # check it in excel
Idents(seu_bs) <- "donor_id"
MUC6_H <- WhichCells(seu_bs, ident = c("P5846",	"Patient22",	"Patient21",	"Patient14",	"L05",	"P5866",	"P6709",	"L07",	"L04",	"P5931",	"P6592"))
MUC6_L <- WhichCells(seu_bs, ident = c("L03",	"L06",	"L08",	"P6649",	"L02",	"L09",	"P6342",	"Patient16",	"P6207",	"L01",	"Patient17"))
seu_bs <- SetIdent(seu_bs, cells = MUC6_H, value = "MUC6_H")
seu_bs$MUC6_level <- Idents(seu_bs)
seu_bs <- SetIdent(seu_bs, cells = MUC6_L, value = "MUC6_L")
seu_bs$MUC6_level <- Idents(seu_bs)
DimPlot(seu_bs, group.by = "MUC6_level") + NoAxes()

# PB for DESeq2 following Seurat vignette----
bulk <- AggregateExpression(seu_bs, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("donor_id", "MUC6_level"))
tail(Cells(bulk))
bulk$donor_id <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
write.table(as.matrix(de_markers),"results/DEG//DESeq2_Seurat_FindMarkers_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

# DEseq2 from matrix, make rank for fgsea----
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.1, gene_name,"")), colour = "red", size = 3)
write.table(as.matrix(res),"results/DEG//DESeq2_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# prepare gene sets----
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})
# fgsea----
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
write.table(fgseaRes[,-8],"results/DEG/fgseaRes_C6_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
plotEnrichment(collections$C6[["EGFR_UP.V1_DN"]], ranks) + labs(title="C6_EGFR_UP.V1_DN")
plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
plotEnrichment(collections$C6[["MEK_UP.V1_DN"]], ranks) + labs(title="C6_MEK_UP.V1_DN")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_UP"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_UP")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

MUC6KO_UP_DN <- gmtPathways("gene_set/MUC6KO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = MUC6KO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_UP"]], ranks) + labs(title="MUC6KO_UP")
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_DN"]], ranks) + labs(title="MUC6KO_DN")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_MUC6KO_UP_DN_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

# Rank patients by ARID1A level----
avg_bs<-AverageExpression(seu_bs, group.by = "donor_id" ) 
write.table(as.matrix(avg_bs$RNA), "results/GeneExpression/avg_bs.txt", sep="\t",col.names = T, row.names = T)  # check it in excel
Idents(seu_bs) <- "donor_id"
table(seu_bs$donor_id)
ARID1A_H <- WhichCells(seu_bs, ident = c("Patient17",	"Patient14",	"Patient22",	"Patient16",	"Patient21",	"P5866",	"P5846",	"L05",	"L09",	"L07",	"L04"))
ARID1A_L <- WhichCells(seu_bs, ident = c("L08",	"L06",	"P6207",	"L03",	"L02",	"P5931",	"P6342",	"P6592",	"L01",	"P6709",	"P6649"))
seu_bs <- SetIdent(seu_bs, cells = ARID1A_H, value = "ARID1A_H")
seu_bs$ARID1A_level <- Idents(seu_bs)
seu_bs <- SetIdent(seu_bs, cells = ARID1A_L, value = "ARID1A_L")
seu_bs$ARID1A_level <- Idents(seu_bs)
DimPlot(seu_bs, group.by = "ARID1A_level") + NoAxes()

# PB for DESeq2 following Seurat vignette----
bulk <- AggregateExpression(seu_bs, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("donor_id", "ARID1A_level"))
tail(Cells(bulk))
bulk$donor_id <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$ARID1A_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$ARID1A_level <- factor(x = bulk$ARID1A_level, levels = c("L", "H"))
Idents(bulk) <- "ARID1A_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
write.table(as.matrix(de_markers),"results/DEG//DESeq2_Seurat_FindMarkers_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

# DEseq2 from matrix, make rank for fgsea----
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.1, gene_name,"")), colour = "red", size = 3)
write.table(as.matrix(res),"results/DEG//DESeq2_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# prepare gene sets----
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})
# fgsea----
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
write.table(fgseaRes[,-8],"results/DEG/fgseaRes_C6_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

# plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
# plotEnrichment(collections$C6[["EGFR_UP.V1_DN"]], ranks) + labs(title="C6_EGFR_UP.V1_DN")
# plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
# plotEnrichment(collections$C6[["MEK_UP.V1_DN"]], ranks) + labs(title="C6_MEK_UP.V1_DN")
# plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_UP"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_UP")
# plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
# plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
# plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

ARID1AKO_UP_DN <- gmtPathways("gene_set/ARID1AKO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = ARID1AKO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_UP"]], ranks) + labs(title="ARID1AKO_UP")
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_DN"]], ranks) + labs(title="ARID1AKO_DN")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_ARID1AKO_UP_DN_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)
