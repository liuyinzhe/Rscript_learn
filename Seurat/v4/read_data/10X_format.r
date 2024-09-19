library(Seurat)

setwd("Single_Cell_transcriptome/data_fomart")



dir <- "D:/Rscript/Single_Cell_transcriptome/data_fomart/10x/sample/"
pro <- "BC"
samples <- list.files(dir, pattern = "^BC", recursive = FALSE)

mmRNAList <- list()

for (i in seq_along(samples)) {
  mmRNA <- CreateSeuratObject(
    counts = Read10X(paste0(dir, samples[[i]])),
    project = pro,
    min.cells = 5,
    min.features = 200
  )
  # Assign a unique orig.ident identifier to each sample
  mmRNA@meta.data$orig.ident <- samples[[i]]  # paste0('p',samples[[i]])
  mmRNAList[[i]] <- mmRNA
}


# merge
mmRNA <- merge(x = mmRNAList[[1]],y= mmRNAList[ -1 ],add.cell.ids = samples)



saveRDS(mmRNA,"mmRNA.rds")
