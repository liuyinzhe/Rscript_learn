library(Seurat)

setwd("Single_Cell_transcriptome/data_fomart")
setdata <- getwd() 

file <- "xx.h5" 
a <- Read10X_h5(file , use.names = T)
b=gsub(".h5","",file)
mmRNA<-CreateSeuratObject(counts =a,project = b,min.cells = 3,min.features = 200)

saveRDS(mmRNA,"mmRNA.rds")


#############  multiple files #############

files <- list.files(setdata, pattern = "*.h5") 

seurat_list <- list()  

for (file in files) {
  a <- Read10X_h5(file, use.names = T)
  a_matrix <- as.matrix(a)
  filename <- gsub(".h5", "", file)  # Replace the suffix to obtain the file name
  

  seurat_obj <- CreateSeuratObject(counts = a_matrix, project = filename, min.cells = 3, min.features = 200)
  
  # Set orig.ident using file name
  seurat_obj$orig.ident <- filename
  
  # Store Seurat objects in a list
  seurat_list[[filename]] <- seurat_obj
}


# str(seurat_list)
# merge
mmRNA <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))


saveRDS(seurat_obj,"mmRNA.rds")


#############  multiple files #############

files <- list.files(pattern = "*.h5") 
seurat_list <- list()

for (file in files) {
  a <- Read10X_h5(file, use.names = T)
  b <- gsub(".h5", "", file)  
  a <- as.matrix(a)

  seurat_obj <- CreateSeuratObject(counts = a, project = b, min.cells = 3, min.features = 200)
  
  # Set orig.ident using file name
  seurat_obj$orig.ident <- b
  
  # Store Seurat objects in a list
  seurat_list[[b]] <- seurat_obj
}

# merge
mmRNA <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
print(mmRNA )

saveRDS(seurat_obj,"mmRNA.rds")