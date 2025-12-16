library(CellChat)
library(AUCell)
cc_db <- cellchat_obj@DB$interaction

cc_genesets <- split(
  unique(c(cc_db$ligand, cc_db$receptor)),
  cc_db$pathway_name
)

expr <- as.matrix(final_integrated@assays$integrated@data)
cells_rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
cc_auc <- AUCell_calcAUC(cc_genesets, cells_rankings, nCores = 8)
auc_mat <- t(as.data.frame(getAUC(cc_auc)))

common_cells <- intersect(rownames(auc_mat), rownames(final_integrated@meta.data))

final_integrated@meta.data <- cbind(
  final_integrated@meta.data,
  auc_mat[common_cells, ]
)

final_integrated@misc$pathway_scores <- auc_mat

set.seed(333)
par(mfrow=c(3,4)) 

tnbc_pathways_filtered <- c(
  "WNT",
  "PDGF",
  "IGF",
  "IFN-I",
  "COMPLEMENT",
  "CCL",
  "ANGPT",
  "ACTIVIN",
  "NODAL",
  "MSTN",
  "AVP",
  "THPO"
)

cells_assignment <- AUCell_exploreThresholds(cc_auc[tnbc_pathways_filtered, ], plotHist=TRUE, assign=TRUE, nCores = 8) 
