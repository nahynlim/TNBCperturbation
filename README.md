# Bayesian Network Modeling of Cell–Cell Communication for Drug Repositioning in Triple-Negative Breast Cancer (TNBC)  
**BE562 Course Project**  
Simran Kaur, Aiden Ly, Nadine Lim, Haeun Oh

## Overview
Triple-negative breast cancer (TNBC) is an aggressive breast cancer subtype lacking ER, PR, and HER2 receptors, which limits targeted treatment options. TNBC tumors are highly heterogeneous and driven by complex interactions between malignant and microenvironmental cell populations. This project develops a computational framework that integrates single-cell RNA-seq analysis, cell–cell communication inference, and Bayesian network modeling to identify candidate drugs and mechanisms of action for potential repositioning in TNBC.

## Approach
- Single-cell RNA-seq preprocessing, clustering, and annotation using *Seurat*
- Ligand–receptor–based cell–cell communication analysis with *CellChat*
- Pathway enrichment analysis using *KEGG* and *Gene Ontology*
- Integration of drug perturbation signatures into a Bayesian network–based knowledge graph

## Data Availability
Processed RDS files for the integrated scRNA-seq dataset and the corresponding CellChat object are available via Google Drive:  
https://drive.google.com/drive/folders/1Tf3r8dciJi7BK7YkJPbXR1ktoFJSVPJK
