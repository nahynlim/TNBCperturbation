library(Seurat)
library(SeuratObject)
library(SingleR)
library(celldex)
library(ggplot2)
library(patchwork)
library(Matrix)

## 161529--------------------------------------------------------------------------
raw_path <- "GSE161529_RAW"
sample_dirs <- list.dirs(raw_path, recursive = FALSE)
seurat_list <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  message("Processing: ", sample_name)
  barcodes <- file.path(sample_dir, list.files(sample_dir, pattern = "barcodes"))
  genes    <- file.path(sample_dir, list.files(sample_dir, pattern = "features|genes"))
  matrix   <- file.path(sample_dir, list.files(sample_dir, pattern = "matrix"))
  counts <- ReadMtx(mtx = matrix, cells = barcodes, features = genes)

  obj <- CreateSeuratObject(counts = counts)
  obj$orig.ident <- sample_name
  obj$cells <- rownames(obj@meta.data)
  obj$sample <- sample_name
  Idents(obj) <- "sample"

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = percent.mt < 20 & nCount_RNA < 25000 & nFeature_RNA > 200 & nFeature_RNA < 6000)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  seurat_list[[sample_name]] <- obj
}

pdf("QC_VlnPlots_161529.pdf", width = 10, height = 6)
for (i in 1:5) {
  print(
    VlnPlot(
      seurat_list[[i]],
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      pt.size = 0.1
    ) + ggplot2::ggtitle(paste0("Sample ", i, " QC VlnPlot"))
  )
}
dev.off()

anchors <- FindIntegrationAnchors(object.list = seurat_list)
integrated_161529 <- IntegrateData(anchorset = anchors)
saveRDS(integrated_161529, 'integrated_161529.RDS')

## 176078-------------------------------------------------------------------------------
raw_path <- "GSE176078_RAW"
sample_dirs <- list.dirs(raw_path, recursive = FALSE)
seurat_list <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  
  barcodes <- file.path(sample_dir, "count_matrix_barcodes.tsv")
  genes    <- file.path(sample_dir, "count_matrix_genes.tsv")
  matrix   <- file.path(sample_dir, "count_matrix_sparse.mtx")
  metadata <- file.path(sample_dir, "metadata.csv")

  counts <- ReadMtx(mtx = matrix,features = genes,cells = barcodes, feature.column = 1)
  obj <- CreateSeuratObject(counts = counts)
  obj$orig.ident <- sample_name
  bc <- colnames(obj)   # barcodes

  parts <- strsplit(bc, "_")
  sample_id <- sapply(parts, `[`, 1)          # "CID4465"
  well_id <- substr(barcode_full, 1, 2)       # "AA"
  
  obj$well_id <- well_id
  obj$barcode_full <- bc
  
  meta_df <- read.csv(metadata)
  obj <- AddMetaData(obj, meta_df)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = percent.mt < 20 & nCount_RNA < 25000 & nFeature_RNA > 200 & nFeature_RNA < 6000)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  
  seurat_list[[sample_name]] <- obj
}

anchors <- FindIntegrationAnchors(object.list = seurat_list)
integrated_176078 <- IntegrateData(anchorset = anchors)

integrated_176078$CID <- integrated_176078@active.ident

integrated_176078$GSM <- cid_to_gsm$GSM[match(integrated_176078$CID, cid_to_gsm$CID)]
head(integrated_176078@meta.data[, c("CID", "GSM")])
Idents(integrated_176078) <- integrated_176078$GSM

saveRDS(integrated_176078, 'integrated_176078.RDS')

## Integrate--------------------------------------------------------------
integrated_161529 <- readRDS('/restricted/projectnb/czshare/A00/nadinelim/integrated_161529.RDS')
integrated_176078 <- readRDS('/restricted/projectnb/czshare/A00/nadinelim/integrated_176078.RDS')
final_anchors <- FindIntegrationAnchors(list(integrated_161529,integrated_176078))
final_integrated <- IntegrateData(anchorset = final_anchors)
saveRDS(final_integrated, '/restricted/projectnb/czshare/A00/nadinelim/p/integrated.RDS')

final_integrated <- readRDS('/restricted/projectnb/czshare/A00/nadinelim/p/integrated.RDS')
DefaultAssay(final_integrated) <- "integrated"

final_integrated <- ScaleData(final_integrated, verbose = FALSE)
final_integrated <- RunPCA(final_integrated, verbose = FALSE)
ElbowPlot(final_integrated, ndims=50)
final_integrated <- RunUMAP(final_integrated, dims = 1:10, reduction = "pca")

final_integrated <- FindNeighbors(final_integrated, dims = 1:10)
final_integrated <- FindClusters(final_integrated, resolution = 0.2)

#cols_to_remove <- c("well_id", "barcode_full", "X", "percent.mito", "subtype", "celltype_subet", "integrated_snn_res.0.2", "integrated_snn_res.0.1", "integrated_snn_res.0.3", "integrated_snn_res.0.5", "seurat_clusters")
#final_integrated@meta.data[, cols_to_remove] <- NULL

## Annotate---------------------------------------------------------------------
#Auto-annotation using humanatlas reference
ref <- celldex::HumanPrimaryCellAtlasData()
excluded_celltypes <- c("Astrocyte", "Chondrocytes", "Embryonic_stem_cells", "Erythroblast", "Gametocytes", "Hepatocytes", 
                        "HSC_-G-CSF", "HSC_CD34+", "iPS_cells", "Keratinocytes", "MEP", "Neuroepithelial_cell", "Neurons", "Osteoblasts")
ref <- ref[, !(ref$label.main %in% excluded_celltypes)]
counts <- final_integrated@assays$integrated@data
pred <- SingleR(test=counts, ref=ref, label=ref$label.main)
final_integrated$singleR.labels <-pred$labels[match(rownames(final_integrated@meta.data), rownames(pred))]
tab <- table(Assigned=pred$labels, Clusters=final_integrated$seurat_clusters)
na_cells <- rownames(final_integrated@meta.data)[is.na(final_integrated$celltype_minor)]

FeaturePlot(final_integrated,
            features = "nCount_RNA",  # dummy feature
            cells = na_cells,
            cols = c("grey80", "red")) +
  ggtitle("Cells with NA in celltype_minor")

saveRDS(final_integrated, 'integrated.RDS')

#save featureplots
DefaultAssay(final_integrated) <- "RNA"

pdf("/restricted/projectnb/czshare/A00/nadinelim/p/FeaturePlots_all_panels.pdf", width = 10, height = 7)
  plot_features <- function(obj, feats, title_text) {
    print(
      FeaturePlot(obj, features = feats, label = TRUE) +
        ggtitle(title_text)
    )
  }
  plots <- list(
    plot_features(final_integrated,c("KRT5","KRT14","KRT17","EGFR","TP53"),"Basal-like TNBC"),
    plot_features(final_integrated,c("CD274","CD3E","CD8A","CD4"),"Immunomodulatory TNBC"),
    plot_features(final_integrated,c("AR","PIK3CA"),"LAR TNBC"),
    plot_features(final_integrated,c("CD79A","MS4A1","HLA-DOB","MZB1","POU2AF1","CD79B","FCRLA","FCRL5","RAB30","DERL3"),"B cell markers"),
    plot_features(final_integrated,c("AQP1","CDH5","CLEC14A","CXorf36","CYYR1","ECSCR","HYAL2","MALL","PLVAP","VWF"),"Endothelial markers"),
    plot_features(final_integrated,c("GPX1","TYROBP","HLA-DMA","CD74","HLA-DRA","FCGR3A","HLA-DPB1","PILRA","LAPTM5"),"Macrophage markers"),
    plot_features(final_integrated,c("CALD1","DCN","COL6A3","COL6A2","COL6A1","EGR1","C1S","SERPINF1"),"Stromal markers"),
    plot_features(final_integrated,c("CD3E","CCL5","LCK","CD2","CD3D","CLEC2D","IL32","CD98","LAT","CTSW"),"T cell markers"),
    plot_features(final_integrated,c("TFF1","TFF3","ESR1"),"Normal epithelial markers"),
    plot_features(final_integrated,c("SCGB2A2","SCGB1D2","S100A8","S100A9"),"Tumor epithelial markers"),
    plot_features(final_integrated,c("PDGFRA"),"Fibroblast (PDGFRA)"),
    plot_features(final_integrated,c("VWF"),"Endothelial (VWF only)"),
    plot_features(final_integrated,c("PTPRC","CD14"),"Monocyte markers"),
    plot_features(final_integrated,c("ARTN","L1CAM"),"Hypoxic tumor subpopulation")
  )
  n_per_page<-6
  n_pages<-ceiling(length(plots)/n_per_page)
  for(i in seq_len(n_pages)){
    start<-(i-1)*n_per_page+1;end<-min(i*n_per_page,length(plots))
    print(wrap_plots(plots[start:end],ncol=2))}
dev.off()

# Rename all identities
new_names <- c("0" = "Tumor Epithelial",
               "1" = "CD4+ T-cells",
               "2" = "Tumor Epithelial",
               "3" = "Tumor Epithelial",
               "4" = "Macrophage",
               "5" = "NK/T",
               "6" = "Tumor Epithelial",
               "7" = "Cycling Tumor Epithelial",
               "8" = "CAFs",
               "9" = "Plasmablasts",
               "10" = "B-cells", 
               "11" = "Monocyte",
               "12" = "Plasmablasts",
               "13" = "PVL",
               "14" = "Endothelial",
               "15" = "Cycling T-cells",
               "16" = "Normal Epithelial")

final_integrated <- RenameIdents(object = final_integrated, new_names)

final_integrated$celltype <- Idents(final_integrated)

markerlist <- list()

for (id in c("0","2","3","6","7","12","16")) {
  markers <- FindMarkers(annotated, ident.1 = id)
  markerlist[[id]] <- markers
}

for (i in 1:7) {
  cat("markers for cluster:", names(markerlist)[i], "\n")
  print(head(markerlist[[i]]))
  cat("\n")
}

FeaturePlot(final_integrated, c("SCGB2A2", "SCGB1D2", "S100A8", "S100A9"), ncol = 2, raster=FALSE)
VlnPlot(final_integrated, features = c("PECAM1", "VWF"), ncol = 3, pt.size=0, raster=FALSE) 

DimPlot(object = annotated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = FALSE)

saveRDS(final_integrated, "/restricted/projectnb/czshare/A00/nadinelim/p/integrated.rds")
