library(Seurat)
library(limma)
setwd("Single_Cell_transcriptome/data_fomart")
files <- "xx.csv.gz" 

df=read.csv(file = gzfile("xx.csv.gz"), header=T, check.names=F)
a_matrix=as.matrix(df)
rownames(a_matrix)=a_matrix[,1]
b_matrix=a_matrix[,2:ncol(a_matrix)]
dimnames=list(rownames(b_matrix),colnames(b_matrix))
c_matrix=matrix(as.numeric(as.matrix(b_matrix)),nrow=nrow(b_matrix),dimnames=dimnames)
# Avereps function implements taking the average of the same gene names by column
# Calculate the average value
c_matrix=avereps(c_matrix) #  limma avereps
mmRNA<-CreateSeuratObject(counts =c_matrix,min.cells = 3,min.features = 200)



saveRDS(mmRNA,"mmRNA.rds")


#############  multiple files #############


file_list <- list.files(pattern = "\\.csv\\.gz$")


seurat_objects <- list()


for (file_name in file_list) {

  data <- read.csv(file = gzfile(file_name), header=TRUE, check.names=FALSE)
  
  data_matrix <- as.matrix(data)
  
  # Set the second column(gene) as row index
  rownames(data_matrix) <- data_matrix[,1]
  
  # Convert data into a numerical matrix and remove duplicate values
  b <- data_matrix[,2:ncol(data_matrix)]
  dimnames <- list(rownames(b), colnames(b))
  numeric_matrix <- matrix(as.numeric(as.matrix(b)), nrow=nrow(b), dimnames=dimnames)
  
  # Avereps function implements taking the average of the same gene names by column
  #  limma avereps
  numeric_matrix <- avereps(numeric_matrix) 
  
  seurat_obj <- CreateSeuratObject(counts = numeric_matrix, min.cells = 3, min.features = 200)
  
  
  # Set orig.ident using file name
  seurat_obj$orig.ident <- gsub("\\.csv\\.gz$", "", file_name)
  
  # Store Seurat objects in a list
  seurat_objects[[file_name]] <- seurat_obj
}


# merge
mmRNA <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = names(seurat_objects))

saveRDS(mmRNA,"mmRNA.rds")
