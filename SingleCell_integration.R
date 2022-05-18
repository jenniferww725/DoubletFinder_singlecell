# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

path='/Users/jenniferwang/Desktop/data/GSE180665'
# get data location
dirs <- list.dirs(path, recursive = F, full.names = F)
dirs [1:7]

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  name
  url <- paste0 (path,x)
  cts <- Read10X(data.dir=url)
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}


