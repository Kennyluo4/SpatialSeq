#Ziliang
#2022-4

#####load libraries
library(Seurat)
# library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
#####load the data
d1 <- Load10X_Spatial("/blue/wang/luoziliang/nodulation/SpatialSeq/ahy_spt_d1_mkfastq1/outs",
       filename = "filtered_feature_bc_matrix.h5",
       assay = "Spatial",
       slice = "D1",
       filter.matrix = TRUE,
       to.upper = FALSE,
       image = NULL)
# using sctransform (Hafemeister and Satija, Genome Biology 2019), which which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance. 
d1 <- SCTransform(d1, assay = "Spatial", verbose = FALSE)

######plot the molecular(RNA) counts 
pdf("molecule_count.pdf")
plot1 <- VlnPlot(d1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(d1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()

###### Gene expression visualization
pdf('gene_expression.pdf')
SpatialFeaturePlot(d1, features = c("LOC112703831", "LOC112775219"), alpha = c(0.1, 1),pt.size.factor = 1.6)
dev.off()
# features: Name of the feature to visualize. Provide either group.by OR features, not both. (just put gene name here)
# group.by: Name of meta.data column to group the data by
# pt.size.factor- This will scale the size of the spots. Default is 1.6
# alpha - minimum and maximum transparency. Default is c(1, 1).
# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression

###### Dimensionality reduction, clustering, and visualization

d1 <- RunPCA(d1, assay = "SCT", verbose = FALSE)
d1 <- FindNeighbors(d1, reduction = "pca", dims = 1:30)
d1 <- FindClusters(d1, verbose = FALSE)
d1 <- RunUMAP(d1, reduction = "pca", dims = 1:30)

# visualize the results of the clustering either in UMAP 
pdf('umap.pdf')
p1 <- DimPlot(d1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(d1, label = TRUE, label.size = 3)
# label parameter places a colored box at the median of each cluste
p1 + p2
dev.off()

pdf('umap_cell_types.pdf')
SpatialDimPlot(d1, cells.highlight = CellsByIdentities(object = d1, idents = c(0,1,2)), 
	facet.highlight = TRUE,
	alpha = c(0.5, 1))
dev.off()
#use cells.highlight parameter to demarcate particular cells of interest, idents are the cluster #

## Interactive plotting
SpatialDimPlot(d1, interactive = TRUE)


# Identification of Spatially Variable Features
# To identify molecular features that correlate with spatial location. 
# 1.perform differential expression based on pre-annotated anatomical regions.(clusters)
de_markers <- FindMarkers(d1, ident.1 = 0, ident.2 = 1)
#take top 3 marker to differentiate the tissue types
pdf('top_marker_gene_predefined_expression.pdf')
SpatialFeaturePlot(object = d1, 
				features = rownames(de_markers)[1:6],
				alpha = c(0.1, 1), 
				# ncol = 3,
				)
dev.off()
# 2. alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation
d1 <- FindSpatiallyVariableFeatures(d1, assay = "SCT", features = VariableFeatures(d1)[1:1000],
    selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(d1, selection.method = "markvariogram"), 6)
pdf('top_marker_gene_spatial_variable_expression.pdf')
SpatialFeaturePlot(d1, features = top.features, ncol = 3, alpha = c(0.1, 1))
dev.off()


###Subset out anatomical regions
cortex <- subset(d1, idents = c(0))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
