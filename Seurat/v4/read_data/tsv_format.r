library(Seurat)

setwd("Single_Cell_transcriptome/data_fomart")
setdata <- getwd() 


file_list <- list.files(pattern = "\\.tsv\\.gz$")


seurat_objects <- list()


for (file_name in file_list) {
  data <- read.table(file = gzfile(file_name), header=TRUE, check.names=FALSE)
  data_matrix <- as.matrix(data)
  # Extract the second column of gene names and set them as row names
  rownames(data_matrix) <- data_matrix[,1]

  # Avereps function implements taking the average of the same gene names by column
  #  limma avereps
  b <- data_matrix[,2:ncol(data_matrix)]
  dimnames <- list(rownames(b), colnames(b))
  numeric_matrix <- matrix(as.numeric(as.matrix(b)), nrow=nrow(b), dimnames=dimnames)
  
  seurat_obj <- CreateSeuratObject(counts = numeric_matrix, min.cells = 3, min.features = 200)
  
  
  # Set orig.ident using file name
  seurat_obj$orig.ident <- gsub("\\.tsv\\.gz$", "", file_name)
  
  # Store Seurat objects in a list
  seurat_objects[[file_name]] <- seurat_obj
}

seurat_obj 


mmRNA <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = names(seurat_objects))

saveRDS(mmRNA,"mmRNA.rds")