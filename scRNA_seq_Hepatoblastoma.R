library(Seurat)
library(ggplot2)
library(plyr)
library(harmony)
library(tidyr)
library(dplyr)

setwd("/Users/jenniferwang/Desktop/data/SingleCellHBPilot/")
path="/Users/jenniferwang/Desktop/data/SingleCellHBPilot/"
dirs<- list.dirs(path, recursive = F, full.names = F)
dirs[1:7]

for (x in dirs) {
  print(x)
}
for (x in dirs) {
  name <- gsub("_filtered_feature_bc_matrix", "", x)
  print(name)
  cts <- ReadMtx(mtx=paste0(path, x, '/matrix.mtx.gz'),
                 features=paste0(path, x, '/features.tsv.gz'),
                 cells=paste0(path, x, '/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts=cts))
  print(name)
}

#merge dataset
# use ls() to check the list number as it will changed all the time
mergedhb.seurat <- merge(HB17_background, y=c(HB17_PDX, HB17_tumor,HB30_PDX, HB30_tumor, HB53_background, HB53_tumor),
                 add.cell.ids= ls() [9:15],
                 project= "HB")

View(mergedhb.seurat@meta.data)

# create a sample column
 mergedhb.seurat$sample <- rownames(mergedhb.seurat@meta.data)

#split sample column
 mergedhb.seurat@meta.data <- separate(mergedhb.seurat@meta.data, col="sample", into=c("Patient","Type","Barcode"), sep = "_")
 View(mergedhb.seurat@meta.data)
 
# QC and filtering
 mergedhb.seurat$mitopercent <- PercentageFeatureSet(mergedhb.seurat, pattern= "^MT-")

 mergedhb.seurat.filtered <- subset(mergedhb.seurat, subset= nCount_RNA>800 & nFeature_RNA > 500 & mitopercent < 10)
 
# perform standard workflow
 mergedhb.seurat.filtered <- NormalizeData(mergedhb.seurat.filtered) %>% 
                             FindVariableFeatures() %>%
                             ScaleData() %>%
                             RunPCA()
 # determine the "dimensionality" of the dataset(JackStraw (this is too slow), Elbow Plot)
 
  ElbowPlot(mergedhb.seurat.filtered)
 
 
 #harmony integration 
 
mergedhb.seurat.filtered = RunHarmony(mergedhb.seurat.filtered, group.by.vars = "Patient") # remove batch from each sample

mergedhb.seurat.filtered = RunUMAP(mergedhb.seurat.filtered, reduction="harmony", dims=1:30)
mergedhb.seurat.filtered = FindNeighbors(mergedhb.seurat.filtered, reduction="harmony",dims=1:30)
#use 3 resolutions

mergedhb.seurat.filtered = FindClusters(mergedhb.seurat.filtered, resolution= c(0.05, 0.1, 0.5, 1, 2, 3))
View (mergedhb.seurat.filtered@meta.data)
DimPlot(mergedhb.seurat.filtered, group.by = "RNA_snn_res.0.05", label = TRUE)

# Setting identity of clusters
ident <-Idents(mergedhb.seurat.filtered) <- "RNA_snn_res.0.05"
DimPlot(mergedhb.seurat.filtered, reduction= "umap", label=TRUE)

saveRDS(mergedhb.seurat.filtered, file = "/Users/jenniferwang/Desktop/data/SingleCellHBPilot/HB_integrate_reciPCA_7samples_harmony_perMito10.RDS")

#umap of cell classes

DimPlot(mergedhb.seurat.filtered, reduction= "umap",  group.by = "Patient", label=TRUE)+
  theme(legend.position = "bottom",
        legend.text= element_text(family="Helvetica", size=10),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size = 1)) + ggtitle("") 
 

# umap of sample groups 


DimPlot(mergedhb.seurat.filtered, reduction= "umap", group.by = "Type", cols = c("#3399FF","#FFCC99","#FF6666")) +
  theme(legend.position = "bottom",
        legend.text = element_text(family = "Helvetica", size = 10),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=1)) + ggtitle("")

# feature plots (figure 3b)




DefaultAssay(mergedhb.seurat.filtered) = "RNA"

FeaturePlot(mergedhb.seurat.filtered, features = c("GPC3","CYP3A4","COL6A3",
                                      "CD163","FLT1","PTPRC"), cols = c("#F5F5F5","red"), ncol = 3)





