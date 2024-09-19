library(Seurat)

setwd("Single_Cell_transcriptome/data_fomart")
setdata <- getwd() 
files <- "GSE72056.txt" 

a <- read.delim(files,header=TRUE,sep="\t")
b=gsub(".txt","",files)

# gene name -> Cell
Cell<- a[["Cell"]]
rownames(a)= Cell
a[["Cell"]]<-NULL
#class(a)
a_matrix <-as.matrix (a)
#class(a)
mmRNA <-CreateSeuratObject(counts=a_matrix,project =b,min.cells = 3,min.features = 200)


saveRDS(mmRNA,"mmRNA.rds")

#############  multiple files #############



files <- list.files(pattern = "*.txt") 

seurat_list <- list()  

for (file in files) {
  a <- read.delim(file, header = TRUE, sep = "\t")
  pj_name <- gsub(".txt", "", file)  
  rownames(a) <- a[,1]
  a_matrix <- as.matrix(a)
  
  seurat_obj <- CreateSeuratObject(counts = a_matrix, project = pj_name, min.cells = 3, min.features = 200)
  
  # Set orig.ident using file name
  seurat_obj$orig.ident <- pj_name
  
  # Store Seurat objects in a list
  seurat_list[[pj_name]] <- seurat_obj
}


str(seurat_list)
# merge
mmRNA <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))


saveRDS(mmRNA,"mmRNA.rds")