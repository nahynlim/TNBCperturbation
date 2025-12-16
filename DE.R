library(Seurat)
library(stats)
library(EnhancedVolcano)

final_integrated[["RNA"]] <- JoinLayers(final_integrated[["RNA"]])
de <- FindMarkers(final_integrated, ident.1 = 'Normal Epithelial', ident.2 = 'Tumor Epithelial' , test.use = "wilcox")
de$gene <- rownames(de)
if (!"p_val_adj" %in% names(de)) de$p_val_adj <- p.adjust(de$p_val, "BH")
de$p_val_adj[is.na(de$p_val_adj)] <- 1

final_integrated@misc$de_epithelial <- de
topgenes =  de$gene[1:15]
volcano <- EnhancedVolcano(de , 
                rownames(de),
                x ="avg_log2FC", 
                y ="p_val_adj",
                selectLab = topgenes,
                subtitle = "Negative LFC (left)= UP in Tumor Epithelial",
                title = NULL,
                legendPosition = "right",
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black'
                )
volcano
