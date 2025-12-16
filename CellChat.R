library(CellChat)
library(Seuart)

DefaultAssay(final_integrated) <- "RNA"
final_integrated[["RNA"]] <- JoinLayers(final_integrated[["RNA"]])
data.input <- GetAssayData(final_integrated, assay = "RNA", slot = "data") # normalized data matrix
labels <- final_integrated$celltype
samples <- final_integrated@meta.data$orig.ident
unique(samples)
meta <- data.frame(group = labels, row.names = names(labels)) #create a dataframe of the cell labels
meta$samples <- samples

cellchat_obj <- createCellChat(object = final_integrated, group.by = "celltype", assay = "RNA")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
#CellChatDB.use <- subsetDB(CellChatDB, search = "all")
CellChatDB.use <- CellChatDB
unique(CellChatDB.use$interaction$annotation)
cellchat_obj@DB <- CellChatDB.use

# Subset the expression data of signaling genes for saving computation cost
cellchat_obj <- subsetData(cellchat_obj)

# Identify expressed ligand-receptor pairs
# Identify signaling genes that are significantly overexpressed in each cell group 
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
# Identify overexpressed ligand-receptor interaction pairs
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)

# Infer cell-cell communication at a ligand-receptor pair level
# parameter: type - computing the average gene expression per cell group
# - default type; triMean (fewer interactions but stronger interactions)
# - type "truncatedMean"; needs trim (0.1, 0.05) 
cellchat_obj <-  computeCommunProb(cellchat_obj, type = "triMean", nboot = 20)
cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)

# Infer cell-cell communication at a signaling pathway level
cellchat_obj <-  computeCommunProbPathway(cellchat_obj)
# Calculate aggregated cell-cell communication network
cellchat_obj <- aggregateNet(cellchat_obj)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat_obj@idents))
#par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of Interactions")
p2 <- netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction Strength")

png("interaction_number.png", width = 2400, height = 2400, res = 300)
p1
dev.off()

png("interaction_strength.png", width = 2400, height = 2400, res = 300)
p2
dev.off()

path_matrix <- cellchat_obj@netP
pathways <- path_matrix$pathways     
prob_arr <- path_matrix$prob      
ct_names <- dimnames(prob_arr)[[1]]   

sender_scores <- sapply(seq_along(pathways), function(i){
  mat <- prob_arr[,,i]  
  rowSums(mat)          
})
colnames(sender_scores) <- pathways

receiver_scores <- sapply(seq_along(pathways), function(i){
  mat <- prob_arr[,,i]
  colSums(mat)         
})
colnames(receiver_scores) <- pathways

cell_ct <- final_integrated$celltype

sender_per_cell   <- sender_scores[cell_ct, ]
receiver_per_cell <- receiver_scores[cell_ct, ]

colnames(sender_per_cell)   <- paste0("Sender_", colnames(sender_scores))
colnames(receiver_per_cell) <- paste0("Receiver_", colnames(receiver_scores))

final_integrated@meta.data <- cbind(
  final_integrated@meta.data,
  sender_per_cell,
  receiver_per_cell
)

pathways.use <- c("SPP1", "MIF", "CXCL", "TGF-B", "VEGF", "COLLAGEN")
pathways.use <- intersect(pathways.use, cellchat_obj@netP$pathways)
pathways.use

df.comm <- subsetCommunication(
  cellchat_obj,
  sources.use = c("Macrophage", "CAFs"),
  targets.use = "Tumor Epithelial",
  signaling   = pathways.use
)

df.heat <- df.comm %>%
  group_by(pathway_name, source) %>%
  summarise(prob = mean(prob), .groups = "drop") %>%  
  mutate(
    source       = factor(source, levels = c("Macrophage", "CAFs")),
    pathway_name = factor(pathway_name, levels = rev(pathways.use)) 
  )

cc_heatmap <- ggplot(df.heat, aes(x = source, y = pathway_name, fill = prob)) +
  geom_tile(color = "grey90", size = 0.2) +
  scale_fill_gradient(low = "white", high = "#d7301f",
                      name = "Communication\nstrength") +
  labs(
    x = "Sender cell type",
    y = "Signaling pathway",
    title = "Macrophage/CAF → Tumor epithelial signaling"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 9),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

senders.use <- c(
  "Macrophage", "CAFs", "Monocyte", "Endothelial",
  "B cells", "Plasmablasts",
  "CD4+ T cells", "Cycling T cells", "NK/T",
  "Perivascular-like cells"
)

netVisual_bubble(
  cellchat_obj,
  sources.use   = senders.use,
  targets.use   = "Tumor Epithelial",
  signaling     = pathways.use,
  remove.isolate = TRUE,
  max.dataset    = 1,
  title.name     = "Ligand–receptor interactions targeting Tumor epithelial cells"
)

netVisual_bubble(
  cellchat_obj,
  sources.use   = senders.use,
  targets.use   = "Tumor Epithelial",
  remove.isolate = TRUE,
  max.dataset    = 1,
  title.name     = "Global ligand–receptor network to Tumor epithelial cells"
)

