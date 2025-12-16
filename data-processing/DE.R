library(Seurat)
library(stats)
library(EnhancedVolcano)

final_integrated[["RNA"]] <- JoinLayers(final_integrated[["RNA"]])
DefaultAssay(final_integrated) <- "RNA"
de <- FindMarkers(final_integrated, ident.1 = 'Normal Epithelial', ident.2 = 'Tumor Epithelial' , test.use = "wilcox")
de$gene <- rownames(de)
if (!"p_val_adj" %in% names(de)) de$p_val_adj <- p.adjust(de$p_val, "BH")
  de$p_val_adj[is.na(de$p_val_adj)] <- 1

final_integrated@misc$de_epithelial <- de
de_negLFC <- de %>%
  dplyr::filter(avg_log2FC < 0)
topgenes =  c(de$gene[1:80], de_negLFC$gene[1:10])

volcano <- EnhancedVolcano(
  de, 
  lab = rownames(de),
  x ="avg_log2FC", 
  y ="p_val",
  pointSize = 0.5,
  #title = "Negative LFC (left)= UP in Tumor Epithelial",
  title = NULL,
  subtitle = NULL,
  labSize = 3.0,
  lab
  pCutoff = 1e-02,
  legendPosition = 'right',
  drawConnectors = FALSE,
  max.overlaps = 100,
  selectLab = topgenes,
  ylim = c(0,290)
)
