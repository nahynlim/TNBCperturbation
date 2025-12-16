setwd('~/Downloads')
integrated_processed <- readRDS('integrated_processed.RDS')

install.packages(c("dplyr", "bnlearn", "igraph"))
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
remotes::install_github("cmap/cmapR")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Rgraphviz")

library(dplyr)
library(bnlearn)
library(igraph)
library(cmapR)
library(Seurat)
library(Rgraphviz)
library(visNetwork)

genes <- rownames(integrated_processed@assays$integrated@data)
de <- FindMarkers(
  integrated_processed,
  ident.1 = "Tumor Epithelial",  
  ident.2 = "Normal Epithelial"  
)
up_genes <- rownames(de)[de$avg_log2FC > 1 & de$p_val_adj < 0.05]

gct_data <- read.delim("clue_be562.gct", skip = 2, header = TRUE, row.names = 1)
gct_data <- gct_data[-1,]

query_gct <- read.delim("clue_be562.gct", skip=2, stringsAsFactors = FALSE, check.names=FALSE)
gsea_gct  <- read.delim("clue_be562_2.gct", skip=2, stringsAsFactors = FALSE, check.names=FALSE)

head(query_gct)
head(gsea_gct)

names(query_gct)
names(gsea_gct)

if("NES" %in% names(gsea_gct)) {
  gsea_gct <- gsea_gct %>%
    group_by(id) %>%
    arrange(desc(NES)) %>%
    slice(1) %>%
    ungroup()
}

merged <- merge(query_gct, gsea_gct, by="id", all.x=TRUE, suffixes=c(".query",".gsea"))

cols_for_bn <- c("pert_iname", "cell_line", "pert_type", "moa", "NES", "pval", "pathway")
cols_for_bn <- cols_for_bn[cols_for_bn %in% names(merged)]
bn_data <- merged[, cols_for_bn, drop=FALSE]

if("NES" %in% names(bn_data)) {
  bn_data$NES <- cut(as.numeric(bn_data$NES), breaks=3, labels=c("low", "med", "high"))
}
if("pval" %in% names(bn_data)) {
  bn_data$pval <- cut(as.numeric(bn_data$pval), breaks=3, labels=c("low", "med", "high"))
}
factor_cols <- c("pert_iname", "cell_line", "pert_type", "moa", "pathway")
present_factor_cols <- intersect(factor_cols, names(bn_data))
bn_data[present_factor_cols] <- lapply(bn_data[present_factor_cols], as.factor)
bn_data <- na.omit(bn_data)

fit <- hc(bn_data)
bn_fit <- bn.fit(fit, bn_data)

if("NES" %in% names(bn_fit)) print(bn_fit$NES)
if("cell_line" %in% names(bn_data)) {
  cpquery(bn_fit, event = (NES == "high"), evidence = (cell_line == "MCF7"))
}

merged <- merged %>% filter(moa != "-666")

top_drugs <- names(sort(table(merged$pert_iname)))[1:80]
filtered <- merged %>% filter(pert_iname %in% top_drugs)

top_moas <- names(sort(table(filtered$moa), decreasing = TRUE))
top_moas <- top_moas[1:80]
filtered$moa <- ifelse(filtered$moa %in% top_moas, filtered$moa, "Other_moa")

filtered$cell_type <- factor(filtered$cell_iname.query)

drug_mat <- model.matrix(~ pert_iname - 1, data = filtered)
colnames(drug_mat) <- paste0("drug_", gsub("\\|", "_", gsub("^pert_iname", "", colnames(drug_mat))))
drug_mat <- as.data.frame(apply(drug_mat, 2, function(x) factor(ifelse(x==1, "yes", "no"))))

moa_mat <- model.matrix(~ moa - 1, data = filtered)
colnames(moa_mat) <- paste0("moa_", gsub("\\|", "_", gsub("^moa", "", colnames(moa_mat))))
moa_mat <- as.data.frame(apply(moa_mat, 2, function(x) factor(ifelse(x==1, "yes", "no"))))

cell_mat <- model.matrix(~ cell_type - 1, data = filtered)
colnames(cell_mat) <- paste0("cell_", gsub("\\|", "_", gsub("^cell_type", "", colnames(cell_mat))))
cell_mat <- as.data.frame(apply(cell_mat, 2, function(x) factor(ifelse(x==1, "yes", "no"))))

bn_data <- cbind(drug_mat, moa_mat, cell_mat, filtered[, setdiff(names(filtered), c("pert_iname", "moa", "cell_type"))])
bn_data[] <- lapply(bn_data, factor)

bad_vars <- sapply(bn_data, function(x) nlevels(x) < 2)
bn_data <- bn_data[, !bad_vars, drop = FALSE]

fit <- hc(bn_data)
bn_fit <- bn.fit(fit, bn_data, method = "bayes")

g <- as.graphNEL(fit)

edges_list <- lapply(nodes(g), function(n) {
  to <- edges(g)[[n]]
  if(length(to)==0) return(NULL)
  data.frame(from=n, to=to, stringsAsFactors=FALSE)
})
edges_df <- do.call(rbind, edges_list)

all_nodes <- nodes(g)

all_nodes <- grep("^(drug_|moa_|cell_)", all_nodes, value = TRUE)

edges_list <- lapply(all_nodes, function(n) {
  to <- edges(g)[[n]]
  to <- intersect(to, all_nodes) 
  if(length(to) == 0) return(NULL)
  data.frame(from = n, to = to, stringsAsFactors = FALSE)
})
edges_df <- do.call(rbind, edges_list)

drug_nodes <- grep("^drug_", all_nodes, value = TRUE)
moa_nodes  <- grep("^moa_",  all_nodes, value = TRUE)
cell_nodes <- grep("^cell_", all_nodes, value = TRUE)
cell_nodes <- cell_nodes[cell_nodes != "cell_iname.query"]

node_colors <- ifelse(all_nodes %in% drug_nodes, "lightskyblue",
                      ifelse(all_nodes %in% moa_nodes, "moccasin",
                             ifelse(all_nodes %in% cell_nodes, "plum", "white")))

nodes_df <- data.frame(
  id = all_nodes,
  label = gsub("^(drug_|moa_|cell_)", "", all_nodes), 
  color = node_colors,
  shape = "ellipse",
  size = 100,
  font = list(size = 150),
  stringsAsFactors = FALSE
)


visNetwork(nodes_df, edges_df) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

# Adding Layer of Probabilities 

cell_nodes <- grep("^cell_", names(bn_data), value = TRUE)

drug_nodes <- grep("^drug_", names(bn_data), value = TRUE)
prob_list <- list()

for (cell in cell_nodes) {
  for (drug in drug_nodes) {
    prob <- cpquery(
      bn_fit,
      event = (eval(as.name(cell)) == "yes"),
      evidence = (eval(as.name(drug)) == "yes"),
      n = 5000
    )
    prob_list[[paste(cell, drug, sep = "_")]] <- prob
  }
}

prob_df <- do.call(rbind, lapply(names(prob_list), function(x) {
  parts <- strsplit(x, "_")[[1]]
  data.frame(cell = paste(parts[1], parts[2], sep = "_"),
             drug = paste(parts[-c(1,2)], collapse = "_"),
             probability = prob_list[[x]],
             stringsAsFactors = FALSE)
}))

head(prob_df)
prob_df <- prob_df %>%
  mutate(
    from = drug,      
    to   = cell,    
    label = round(probability, 2)  
  )


edges_annotated <- edges_df %>%
  left_join(prob_df, by = c("from", "to")) %>%
  mutate(
    label = ifelse(is.na(label), "", as.character(label)), 
    font.size = ifelse(label == "", 0, 20)                
  )

visNetwork(nodes_df, edges_annotated) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

# Filtered Probabilities

prob_df <- prob_df %>%
  filter(probability > 0.40) %>% 
  mutate(
    from = drug,
    to   = cell,
    label = sprintf("%.2f", probability),
    font.size = 20
  )

edges_combined <- bind_rows(
  edges_df %>% mutate(label = "", font.size = 0),  
  prob_df %>% select(from, to, label, font.size)  
) %>%
  distinct(from, to, .keep_all = TRUE)  

visNetwork(nodes_df, edges_combined) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

# Inferencing of Bayesian Network for MOAs for Drug X 

colnames(prob_moa_cell)
colnames(prob_moa_drug)

prob_moa_cell_clean <- prob_moa_cell %>%
  rename(cell = to, moa = from) %>%   
  filter(!is.na(probability)) %>%
  mutate(probability = as.numeric(probability)) %>%
  group_by(cell, moa) %>%
  summarise(probability = mean(probability), .groups = "drop")

prob_moa_drug_clean <- prob_moa_drug %>%
  rename(drug = from, moa = to) %>%   
  filter(!is.na(probability)) %>%
  mutate(probability = as.numeric(probability)) %>%
  group_by(drug, moa) %>%
  summarise(probability = mean(probability), .groups = "drop")

moa_cell_mat <- prob_moa_cell_clean %>%
  pivot_wider(names_from = moa, values_from = probability, values_fill = 0)
rownames(moa_cell_mat) <- moa_cell_mat$cell
moa_cell_mat$cell <- NULL

drug_moa_mat <- prob_moa_drug_clean %>%
  pivot_wider(names_from = moa, values_from = probability, values_fill = 0)
rownames(drug_moa_mat) <- drug_moa_mat$drug
drug_moa_mat$drug <- NULL

common_moas <- intersect(colnames(drug_moa_mat), colnames(moa_cell_mat))
drug_moa_mat <- drug_moa_mat[, common_moas, drop = FALSE]
moa_cell_mat <- moa_cell_mat[, common_moas, drop = FALSE]

indirect_prob <- as.matrix(drug_moa_mat) %*% t(as.matrix(moa_cell_mat))
indirect_prob_df <- as.data.frame(indirect_prob)
indirect_prob_df$drug <- rownames(indirect_prob_df)

indirect_prob_long <- indirect_prob_df %>%
  pivot_longer(-drug, names_to = "cell", values_to = "probability") %>%
  filter(probability > 0.2)

head(indirect_prob_long)
g <- graph_from_data_frame(edges_combined, directed = TRUE, vertices = nodes_df)
degree(g, mode = "out")
betweenness(g)

DrugX_name <- "fomocaine"
DrugX_id <- paste0("drug_", DrugX_name)

required_node_cols <- c("id", "label", "value", "color")
for(col in required_node_cols){
  if(!col %in% colnames(nodes_df)){
    if(col == "value") nodes_df[[col]] <- 0
    else if(col == "color") nodes_df[[col]] <- "lightblue"
    else nodes_df[[col]] <- nodes_df$id  
  }
}

if(!DrugX_id %in% nodes_df$id){
  nodes_df <- nodes_df %>%
    add_row(
      id = DrugX_id,
      label = DrugX_name,
      value = 200,
      color = "red"
    )
}

DrugX_MOA <- "Voltage gated channel inhibitor"
MOA_id <- paste0("moa_", DrugX_MOA)

if(!MOA_id %in% nodes_df$id){
  nodes_df <- nodes_df %>%
    add_row(
      id = MOA_id,
      label = DrugX_MOA,
      value = 50,
      color = "moccasin"
    )
}

required_edge_cols <- c("from", "to", "label", "font.size", "arrows", "width", "color")
for(col in required_edge_cols){
  if(!col %in% colnames(edges_annotated)){
    edges_annotated[[col]] <- if(col %in% c("arrows","color")) "black" else 0
  }
}

edge_to_add <- data.frame(
  from = DrugX_id,
  to   = MOA_id,
  label = "",
  font.size = 20,
  arrows = "to",
  width = 2,
  color = "red"
)

if(!any(edges_annotated$from == edge_to_add$from & edges_annotated$to == edge_to_add$to)){
  edges_annotated <- bind_rows(edges_annotated, edge_to_add)
}

visNetwork(nodes_df, edges_annotated) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

existing_nodes <- nodes_df$id
edges_net <- edges_annotated %>%
  filter(from %in% existing_nodes & to %in% existing_nodes) %>%
  select(from, to)
g <- graph_from_data_frame(edges_net, directed = TRUE, vertices = nodes_df)

DrugX_name <- "NMS-1286937"
DrugX_id <- paste0("drug_", DrugX_name)

if(!DrugX_id %in% nodes_df$id){
  nodes_df <- nodes_df %>%
    add_row(
      id = DrugX_id,
      label = DrugX_name,
      value = 200,
      color = "red"
    )
}

drug_nodes <- nodes_df$id[grepl("^drug_", nodes_df$id) & nodes_df$id != DrugX_id]
moa_nodes <- nodes_df$id[grepl("^moa_", nodes_df$id)]

drug_similarity <- sapply(drug_nodes, function(d){
  neighbors(g, d, mode = "out")$name %>% intersect(moa_nodes) %>% length()
})

top_sim_drug <- names(which.max(drug_similarity))
message("Most similar drug to DrugX based on network: ", top_sim_drug)

predicted_moas <- neighbors(g, top_sim_drug, mode = "out")$name
predicted_moas <- predicted_moas[predicted_moas %in% moa_nodes]

for(moa in predicted_moas){
  if(!moa %in% nodes_df$id){
    nodes_df <- nodes_df %>%
      add_row(
        id = moa,
        label = gsub("^moa_", "", moa),
        value = 50,
        color = "moccasin"
      )
  }
  
  if(!any(edges_annotated$from == DrugX_id & edges_annotated$to == moa)){
    edges_annotated <- edges_annotated %>%
      add_row(
        from = DrugX_id,
        to = moa,
        label = "",
        font.size = 20,
        arrows = "to",
        width = 2,
        color = "red"
      )
  }
}

visNetwork(nodes_df, edges_annotated) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

visNetwork(nodes_df, edges_annotated) %>%
  visIgraphLayout(layout = "layout_with_fr") %>%   # <-- spreads nodes nicely
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

infer_drug_moas <- function(drug_name, nodes_df, edges_df, top_n = 3) {
  
  DrugX_id <- paste0("drug_", drug_name)
  
  if(!DrugX_id %in% nodes_df$id){
    nodes_df <- nodes_df %>%
      add_row(
        id = DrugX_id,
        label = drug_name,
        value = 200,
        color = "red"
      )
  }
  
  g <- graph_from_data_frame(
    edges_df %>% select(from, to),
    directed = TRUE,
    vertices = nodes_df
  )
  
  drug_nodes <- nodes_df$id[grepl("^drug_", nodes_df$id) & nodes_df$id != DrugX_id]
  cell_nodes <- nodes_df$id[grepl("^cell_", nodes_df$id)]
  moa_nodes  <- nodes_df$id[grepl("^moa_", nodes_df$id)]
  
  drug_cell_mat <- sapply(drug_nodes, function(d){
    as.numeric(cell_nodes %in% neighbors(g, d, mode = "out")$name)
  })
  rownames(drug_cell_mat) <- cell_nodes
  drug_cell_mat <- t(drug_cell_mat) 

  DrugX_cell_vec <- rep(0, length(cell_nodes))
  names(DrugX_cell_vec) <- cell_nodes
  DrugX_cell_vec[sample(length(cell_nodes), 2)] <- 1
  
  similarity <- apply(drug_cell_mat, 1, function(x) sum(x & DrugX_cell_vec))
  top_drugs <- names(sort(similarity, decreasing = TRUE))[1:top_n]
  
  predicted_moas <- unique(unlist(lapply(top_drugs, function(d){
    neighbors(g, d, mode = "out")$name[neighbors(g, d, mode = "out")$name %in% moa_nodes]
  })))
  
  for(moa in predicted_moas){
    if(!moa %in% nodes_df$id){
      nodes_df <- nodes_df %>%
        add_row(
          id = moa,
          label = gsub("^moa_", "", moa),
          value = 50,
          color = "moccasin"
        )
    }
    if(!any(edges_df$from == DrugX_id & edges_df$to == moa)){
      edges_df <- edges_df %>%
        add_row(
          from = DrugX_id,
          to = moa,
          label = "",
          font.size = 20,
          arrows = "to",
          width = 2,
          color = "red"
        )
    }
  }
  
  return(list(nodes_df = nodes_df, edges_df = edges_df, predicted_moas = predicted_moas))
}

res <- infer_drug_moas("linoleic-acid", nodes_df, edges_annotated)
nodes_df <- res$nodes_df
edges_annotated <- res$edges_df
res$predicted_moas

required_node_cols <- c("id", "label", "value", "color")
for(col in required_node_cols){
  if(!col %in% colnames(nodes_df)){
    if(col == "value") nodes_df[[col]] <- 0
    else if(col == "color") nodes_df[[col]] <- "lightblue"
    else nodes_df[[col]] <- nodes_df$id
  }
}

required_edge_cols <- c("from", "to", "label", "font.size", "arrows", "width", "color")
for(col in required_edge_cols){
  if(!col %in% colnames(edges_annotated)){
    edges_annotated[[col]] <- if(col %in% c("arrows", "color")) "black" else 0
  }
}

DrugX_name <- "linoleic-acid"
DrugX_id <- paste0("drug_", DrugX_name)

edges_annotated <- edges_annotated %>%
  mutate(color = ifelse(from == DrugX_id, "red", color),
         width = ifelse(from == DrugX_id, 3, width),
         arrows = ifelse(from == DrugX_id, "to", arrows))

visNetwork(nodes_df, edges_annotated) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)

visNetwork(nodes_df, edges_annotated) %>%
  visIgraphLayout(layout = "layout_with_fr") %>%   # spreads nodes evenly
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE)

unique_counts_df <- data.frame(
  column = names(query_gct),
  n_unique = sapply(query_gct, function(col) length(unique(col))),
  row.names = NULL
)

unique_counts_df

required_node_cols <- c("id", "label", "value", "color")
for(col in required_node_cols){
  if(!col %in% colnames(nodes_df)){
    if(col == "value") nodes_df[[col]] <- 0
    else if(col == "color") nodes_df[[col]] <- "lightblue"
    else nodes_df[[col]] <- nodes_df$id
  }
}

required_edge_cols <- c("from", "to", "label", "font.size", "arrows", "width", "color")
for(col in required_edge_cols){
  if(!col %in% colnames(edges_annotated)){
    edges_annotated[[col]] <- if(col %in% c("arrows","color")) "black" else 0
  }
}

infer_drug_moas <- function(drug_name, nodes_df, edges_df, top_n = 3) {
  
  DrugX_id <- paste0("drug_", drug_name)
  
  if(!DrugX_id %in% nodes_df$id){
    nodes_df <- nodes_df %>%
      add_row(
        id = DrugX_id,
        label = drug_name,
        value = 200,
        color = "red"
      )
  }
  
  g <- graph_from_data_frame(edges_df %>% select(from, to), directed = TRUE, vertices = nodes_df)
  
  drug_nodes <- nodes_df$id[grepl("^drug_", nodes_df$id) & nodes_df$id != DrugX_id]
  cell_nodes <- nodes_df$id[grepl("^cell_", nodes_df$id)]
  moa_nodes  <- nodes_df$id[grepl("^moa_", nodes_df$id)]

  drug_cell_mat <- sapply(drug_nodes, function(d){
    as.numeric(cell_nodes %in% neighbors(g, d, mode = "out")$name)
  })
  rownames(drug_cell_mat) <- cell_nodes
  drug_cell_mat <- t(drug_cell_mat)
  
  DrugX_cell_vec <- rep(0, length(cell_nodes))
  names(DrugX_cell_vec) <- cell_nodes
  DrugX_cell_vec[sample(length(cell_nodes), 2)] <- 1
  
  similarity <- apply(drug_cell_mat, 1, function(x) sum(x & DrugX_cell_vec))
  top_drugs <- names(sort(similarity, decreasing = TRUE))[1:top_n]
  
  predicted_moas <- unique(unlist(lapply(top_drugs, function(d){
    neighbors(g, d, mode = "out")$name[neighbors(g, d, mode = "out")$name %in% moa_nodes]
  })))
  
  for(moa in predicted_moas){
    if(!moa %in% nodes_df$id){
      nodes_df <- nodes_df %>%
        add_row(
          id = moa,
          label = gsub("^moa_", "", moa),
          value = 50,
          color = "orange"
        )
    }
    if(!any(edges_df$from == DrugX_id & edges_df$to == moa)){
      edges_df <- edges_df %>%
        add_row(
          from = DrugX_id,
          to = moa,
          label = "",
          font.size = 20,
          arrows = "to",
          width = 2,
          color = "red"
        )
    }
  }
  
  return(list(nodes_df = nodes_df, edges_df = edges_df, predicted_moas = predicted_moas))
}

res <- infer_drug_moas("linoleic-acid", nodes_df, edges_annotated)
nodes_df <- res$nodes_df
edges_annotated <- res$edges_df
res$predicted_moas

visNetwork(nodes_df, edges_combined) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(navigationButtons = TRUE, zoomView = TRUE) %>%
  visPhysics(enabled = FALSE)







