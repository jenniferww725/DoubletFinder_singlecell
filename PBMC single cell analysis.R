library(Seurat)
library(tidyverse)
library(ggplot2)

PBMC.sparse.m <-Read10X_h5("/Users/jenniferwang/Desktop/data/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5")

str(PBMC.sparse.m)

PBMC.sparse.obj<-CreateSeuratObject(counts=PBMC.sparse.m, project="PBMC", min.cells = 3, min.features = 200)


View(PBMC.sparse.obj@meta.data)

PBMC.sparse.obj[["percent.mt"]]<- PercentageFeatureSet(PBMC.sparse.obj, pattern = "^MT-")
View(PBMC.sparse.obj@meta.data)
VlnPlot(PBMC.sparse.obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)

FeatureScatter(PBMC.sparse.obj, feature1= "nCount_RNA", feature2 = "nFeature_RNA")+ geom_smooth(method="lm")

PBMC.sparse.obj<-subset(PBMC.sparse.obj, subset= nFeature_RNA>200 & nFeature_RNA< 5000 & percent.mt <5)

PBMC.sparse.obj <- NormalizeData(PBMC.sparse.obj)

PBMC.sparse.obj <- FindVariableFeatures(PBMC.sparse.obj, selection.method = "vst", nfeatures = 2000)

top10<- head(VariableFeatures(PBMC.sparse.obj),10)

plot1<- VariableFeaturePlot(PBMC.sparse.obj)
LabelPoints(plot=plot1, points=top10, repel=TRUE)

all.genes<- rownames(PBMC.sparse.obj)
PBMC.sparse.obj <- ScaleData(PBMC.sparse.obj, features=all.genes)
str(PBMC.sparse.obj)

PBMC.sparse.obj<-RunPCA(PBMC.sparse.obj,features = VariableFeatures(object=PBMC.sparse.obj))

print(PBMC.sparse.obj [["pca"]],dims=1:5, nfeatures=5)
DimHeatmap(PBMC.sparse.obj, dims=1, cells=500, balanced = TRUE)

ElbowPlot(PBMC.sparse.obj)

PBMC.sparse.obj <-FindNeighbors(PBMC.sparse.obj, dims=1:10)
PBMC.sparse.obj <- FindClusters(PBMC.sparse.obj, resolution = 0.3)


DimPlot(PBMC.sparse.obj, group.by = "RNA_snn_res.0.3", label=TRUE)

Idents(PBMC.sparse.obj) <- "RNA_snn_res.0.3"

head(Idents(PBMC.sparse.obj),5)

PBMC.sparse.obj<-RunUMAP(PBMC.sparse.obj, dims=1:10)
DimPlot(PBMC.sparse.obj, reduction="umap")


cluster2.marker <- FindMarkers(PBMC.sparse.obj, ident.1 = 2, min.pct = 0.25)
head(cluster2.marker, n=5)

PBMC.sparse.obj.markers <- FindAllMarkers(PBMC.sparse.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

PBMC.sparse.obj.markers %>%
  group_by(cluster)%>%
  slice_max(n=2, order_by = avg_log2FC)


VlnPlot(PBMC.sparse.obj, features = c("CCR7","LEF1"))
FeaturePlot(PBMC.sparse.obj, features = c("CCR7", "LEF1","CCR7","GAS5","ITGB1","KLRB1","ITGB1","CRIP1","S100A9" ))

top10 <- PBMC.sparse.obj.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) 

DoHeatmap(PBMC.sparse.obj, features=top10$gene) + NoLegend()



new.clusters.ids<- c("CCR7,LEF1", "ITGB1,MAF", "S100A9, LYZ", "GNLY, NKG7", "PTGDS,GZMB", "IGKC, IGLC3", "S100A9,S100A8", "HLA-DPB1, HLA-DPA1")
names(new.clusters.ids) <- levels(PBMC.sparse.obj)
PBMC.sparse.obj <- RenameIdents(PBMC.sparse.obj, new.clusters.ids)
DimPlot(PBMC.sparse.obj, reduction="umap", label= TRUE, pt.size = 0.5) + NoLegend()
