## Code to generate `cycle.small` PCA embedding and clusters

# We are going to use single-cell expression data from "A single-cell resolution map of mouse hematopoietic stem and progenitorc cell differentation" by Nestorowa et al, Blood 2016 available here: https://pubmed.ncbi.nlm.nih.gov/27365425/

# This dataset is also utilized in the Seurat vignette (located here: https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html) for regressing out cell cycle effects, which is how I came across it

# We are only going to save PC embeddings from this data as a minimal example. Full analysis will be in a vignette on pkgdown

library(Seurat)

# The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "~/Desktop/cyclical_pseudotime/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    as.is = TRUE, row.names = 1)

# Define cell cycle genes - loaded with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- ScaleData(marrow, features = rownames(marrow))

# PCA by cell cycle genes
marrow <- RunPCA(marrow,reduction.key="PC_",features=c(s.genes,g2m.genes),reduction.name="pca_cell_cycle")

marrow <- FindNeighbors(marrow,reduction="pca_cell_cycle",dims=1:5)
marrow <- FindClusters(marrow,resolution=1)

marrow <- BuildClusterTree(marrow,reorder=T,reorder.numeric=T)

# Extract PCA embeddings
pca.x <- marrow@reductions$pca_cell_cycle@cell.embeddings[,c(1:3)]

cycle.small <- data.frame(pca.x,clusters=as.factor(marrow@meta.data$tree.ident))

# Save dataset
usethis::use_data(cycle.small)
