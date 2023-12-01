# subset only "10x 3' v2"----
seu_bs <- subset(seu_bs, assay =="10x 3' v2")

# Rank patients by MUC6 level----
avg_bs<-AverageExpression(seu_bs, group.by = "donor_id" ) 
write.table(as.matrix(avg_bs$RNA), "results/GeneExpression/avg_bs_v2.txt", sep="\t",col.names = T, row.names = T)  # check it in excel
Idents(seu_bs) <- "donor_id"
seu_bs <- subset(seu_bs, idents = "L03", invert = T)
MUC6_H <- WhichCells(seu_bs, ident = c("P5846",	"L05",	"P5866",	"P6709",	"L07",	"L04",	"P5931",	"P6592"))
MUC6_L <- WhichCells(seu_bs, ident = c("L06",	"L08",	"P6649",	"L02",	"L09",	"P6342",	"P6207",	"L01"))
seu_bs <- SetIdent(seu_bs, cells = MUC6_H, value = "H")
seu_bs$MUC6_level <- Idents(seu_bs)
seu_bs <- SetIdent(seu_bs, cells = MUC6_L, value = "L")
seu_bs$MUC6_level <- Idents(seu_bs)
DimPlot(seu_bs, group.by = "MUC6_level") + NoAxes()

# PB for DESeq2 following Seurat vignette----
bulk <- AggregateExpression(seu_bs, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("donor_id", "MUC6_level"))
tail(Cells(bulk))
bulk$donor_id <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
write.table(as.matrix(de_markers),"results/DEG//DESeq2_Seurat_FindMarkers_MUC6LvsH_v2_only.txt", sep ="\t", col.names = T,row.names = F)

# DEseq2 from matrix, make rank for fgsea----
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
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
write.table(as.matrix(res),"results/DEG//DESeq2_MUC6LvsH_v2_only.txt", sep ="\t", col.names = T,row.names = F)

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
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
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
