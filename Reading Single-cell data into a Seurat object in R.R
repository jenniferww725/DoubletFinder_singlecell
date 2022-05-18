library(Seurat)
library(SeuratDisk)

#.RDS format
rds_obj <- readRDS("/Users/jenniferwang/Desktop/data/ependymal_cells.rds")
str(rds_obj)

# 10X CellRanger. HDF5 format
hdf5_obj <- Read10X_h5("/Users/jenniferwang/Desktop/data/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

hdf5_obj [1:10, 1:10]

seuract_hdf5 <- CreateSeuratObject(counts = hdf5_obj)
str(seuract_hdf5)


# .mtx file
mtx_obj <- ReadMtx(mtx="/Users/jenniferwang/Desktop/data/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.tar.gz",       features = "/Users/jenniferwang/Desktop/data/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.tar.gz", cells = "/Users/jenniferwang/Desktop/data/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.tar.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

#.loom files
loom_oj <- Connect(filename="", mode= "r")
seuract_loom <- as. Seurat(loom_oj)

#.h5ad format
# step1 :covert AnnData object to an h5Seurat file
Convert("/Users/jenniferwang/Desktop/data/adata_SS2_for_download.h5ad", dest="h5seurat", overwrite=TRUE)

# step2: load h5seurat file into a Seurat object
seurat_anndata <- LoadH5Seurat("/Users/jenniferwang/Desktop/data/adata_SS2_for_download.h5seurat")





