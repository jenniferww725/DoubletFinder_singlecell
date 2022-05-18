install.packages("DoubletFinder")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(tidyverse)

# data: 10k Human PBMCs, 3' v3.1, Chromium Controller
# Human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25-30 were obtained by 10x Genomics from AllCells.

#Libraries were generated from ~16,000 cells (11,485 cells recovered) as described in the Chromium Single Cell 3' Reagent Kits User Guide (v3.1 Chemistry Dual Index) (CG000315 Rev C) using the Chromium Controller and sequenced on an Illumina NovaSeq 6000 to a read depth of approximately 30,000 mean reads per cell.

#Paired-end, dual indexing Read 1: 28 cycles (16 bp barcode, 12 bp UMI) i5 index: 10 cycles (sample index) i7 index: 10 cycles (sample index) Read 2: 90 cycles (transcript)


cts <- ReadMtx(mtx="/Users/jenniferwang/Desktop/data/raw_feature_bc_matrix/matrix.mtx.gz",
        features="/Users/jenniferwang/Desktop/data/raw_feature_bc_matrix/features.tsv.gz",
        cells="/Users/jenniferwang/Desktop/data/raw_feature_bc_matrix/barcodes.tsv.gz")

pbmc.seurat <- CreateSeuratObject(counts=cts)
str(pbmc.seurat)

pbmc.seurat$mitopercent <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
View(pbmc.seurat@meta.data)
VlnPlot(pbmc.seurat,features = c("nCount_RNA", "nFeature_RNA"), ncol=2)

pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA >200 & nFeature_RNA>1000 & mitopercent<10)

pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered) %>% 
                        FindVariableFeatures() %>%
                        ScaleData() %>%
                        RunPCA()
ElbowPlot(pbmc.seurat.filtered)

pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:15)
pbmc.seurat.filtered <- FindClusters(pbmc.seurat.filtered, resolution = 0.5)
pbmc.seurat.filtered <- RunUMAP(pbmc.seurat.filtered, dims= 1:15)

##pK Identification (no ground-truth)
sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:15, sct= FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT= FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(x= pK, y= BCmetric, group=1)) +
  geom_point()+ geom_line()+theme_classic()

##

#Select the pK that corresponds to max bcmvn to optimize doublet detection
pK <- bcmvn_pbmc %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)

pK <- as.numeric(as.character(pK))

## Doublet Proportion Estimate

annotations <- pbmc.seurat.filtered@meta.data $seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(pbmc.seurat.filtered@meta.data)) # Assuming 7.5% doublet formation rate- should be tailered for dataset
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, PCs = 1:15, pN=0.25, pK=0.21, nExp=nExp_poi, reuse.pANN = FALSE, sct= FALSE)

pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, PCs= 1:15, pN=0.25, pK=0.21, nExp=nExp_poi.adj, reuse.pANN = FALSE, sct= FALSE)

# Visualize doublets
View(pbmc.seurat.filtered@meta.data)
names(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction="umap", group.by = "DF.classifications_0.25_0.21_737")

# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_737)
