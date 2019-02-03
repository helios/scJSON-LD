library(Seurat)
library(RJSONIO)
library(dplyr)

if (file.exists("dataSeurat.Robj")){
  load("dataSeurat.Robj")
} else {
data10x <- Read10X("../../til_single_cell/datasets/10x/K_CRC/outs/filtered_gene_bc_matrices_mex/GRCh38/")
dataSeurat <- CreateSeuratObject(raw.data = data10x, min.cells = 3, min.genes = 200,
                           project = "10xMyExp")

#dataJson <- toJSON(dataSeurat)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataSeurat@data), value = TRUE)
percent.mito <- colSums(dataSeurat@raw.data[mito.genes, ])/colSums(dataSeurat@raw.data)

# AddMetaData adds columns to object@data.info, and is a great place to
# stash QC stats
dataSeurat <- AddMetaData(object = dataSeurat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = dataSeurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = dataSeurat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = dataSeurat, gene1 = "nUMI", gene2 = "nGene")


dataSeurat <- FilterCells(object = dataSeurat, subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))


dataSeurat <- NormalizeData(object = dataSeurat, normalization.method = "LogNormalize",
                      scale.factor = 10000)


dataSeurat <- FindVariableGenes(object = dataSeurat, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = dataSeurat@var.genes)

dataSeurat <- ScaleData(object = dataSeurat, vars.to.regress = c("nUMI", "percent.mito"))


dataSeurat <- RunPCA(object = dataSeurat, pc.genes = dataSeurat@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)


PrintPCA(object = dataSeurat, pcs.print = 1:5, genes.print = 5, use.full = FALSE)


VizPCA(object = dataSeurat, pcs.use = 1:2)


PCAPlot(object = dataSeurat, dim.1 = 1, dim.2 = 2)


dataSeurat <- ProjectPCA(object = dataSeurat, do.print = FALSE)


PCHeatmap(object = dataSeurat, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)


PCHeatmap(object = dataSeurat, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = dataSeurat)


dataSeurat <- FindClusters(object = dataSeurat, reduction.type = "pca", dims.use = 1:15,
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = dataSeurat)


dataSeurat <- RunTSNE(object = dataSeurat, dims.use = 1:15, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = dataSeurat)

save(dataSeurat, file = "dataSeurat.Robj")

}

dataSeurat.markers <- FindAllMarkers(object = dataSeurat, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
dataSeurat.markers %>% group_by(cluster) %>% top_n(2, avg_diff)


cluster1.markers <- FindMarkers(object = dataSeurat, ident.1 = 0, thresh.use = 0.25,
                                test.use = "roc", only.pos = TRUE)

VlnPlot(object = dataSeurat, features.plot = c("MS4A1", "CD79A"))
JoyPlot(object = dataSeurat, features.plot = c("FOXP3", "MAGEH1"))
CellPlot(object = dataSeurat, cell1 = 10, cell2 = 2000, features.plot = c("FOXP3", "MAGEH1"))

FeaturePlot(object = dataSeurat, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14",
                                              "LYZ", "CD8A"), cols.use = c("grey", "blue"),
            reduction.use = "tsne")


top10 <- dataSeurat.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = dataSeurat, genes.use = top10$gene,  slim.col.label = TRUE, remove.key = TRUE)


# First lets stash our identities for later
dataSeurat <- StashIdent(object = dataSeurat, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
dataSeurat <- FindClusters(object = dataSeurat, reduction.type = "pca", dims.use = 1:15,
                     resolution = 0.8, print.output = FALSE)
dataSeurat <- StashIdent(object = dataSeurat, save.name = "ClusterNames_0.8")

dataSeurat <- FindClusters(object = dataSeurat, reduction.type = "pca", dims.use = 1:15,
                           resolution = 0.4, print.output = FALSE)
dataSeurat <- StashIdent(object = dataSeurat, save.name = "ClusterNames_0.4")
dataSeurat <- FindClusters(object = dataSeurat, reduction.type = "pca", dims.use = 1:15,
                           resolution = 0.95, print.output = FALSE)
dataSeurat <- StashIdent(object = dataSeurat, save.name = "ClusterNames_0.95")


# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- TSNEPlot(object = dataSeurat, do.return = TRUE, group.by = "ClusterNames_0.4", no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = dataSeurat, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = TRUE)
plot3 <- TSNEPlot(object = dataSeurat, do.return = TRUE, group.by = "ClusterNames_0.8", no.legend = TRUE, do.label = TRUE)
plot4 <- TSNEPlot(object = dataSeurat, do.return = TRUE, group.by = "ClusterNames_0.95", no.legend = TRUE, do.label = TRUE)

plot_grid(plot1, plot2, plot3, plot4)

SetAllIdent(object = dataSeurat, id = 'ClusterName_0.6')
tcell.markers <- FindMarkers(object = dataSeurat, ident.1 = 0, ident.2 = 4)


SetAllIdent(object = dataSeurat, id = 'ClusterName_0.8')
tcell.markers.0.5 <- FindMarkers(object = dataSeurat, ident.1 = 0, ident.2 = 5)
tcell.markers.0.7 <- FindMarkers(object = dataSeurat, ident.1 = 0, ident.2 = 7)
tcell.markers.0.3 <- FindMarkers(object = dataSeurat, ident.1 = 0, ident.2 = 3)
tcell.markers.2.6 <- FindMarkers(object = dataSeurat, ident.1 = 2, ident.2 = 6)

SetAllIdent(object = dataSeurat, id = 'ClusterName_0.95')
tcell.markers.2.7 <- FindMarkers(object = dataSeurat, ident.1 = 2, ident.2 = 7)

SetAllIdent(object = dataSeurat, id = 'ClusterName_0.95')

dataSeurat.markers.0.95 <- FindAllMarkers(object = dataSeurat, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)

(dataSeurat.markers.0.95 %>% group_by(cluster) %>% top_n(2, avg_diff))$gene



DoHeatmap(object = dataSeurat, genes.use = (dataSeurat.markers.0.95 %>% group_by(cluster) %>% top_n(2, avg_diff))$gene,  slim.col.label = TRUE, remove.key = TRUE)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
top2 = (dataSeurat.markers.0.95 %>% group_by(cluster) %>% top_n(2, avg_diff))
FeaturePlot(object = dataSeurat, features.plot = top2$gene, cols.use = c("green", "blue"), pt.size = 0.2)

VlnPlot(object = dataSeurat, features.plot = top2$gene)
