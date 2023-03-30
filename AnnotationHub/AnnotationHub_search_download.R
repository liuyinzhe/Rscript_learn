
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\AnnotationHub')

# 除非模式如$title: org.Hs.eg.db.sqlite 这种，其它基本没用，找到和使用基因组使用上有问题，clusterProfiler(4.6.1) GO富集会是空的，应该是格式有问题,还是用emapper_annotations 注释的方法更好

library(AnnotationHub)
ah <- AnnotationHub()
org <- ah[ah$rdataclass == "OrgDb",]
# search 'Homo sapiens'
# 查看$title: org.Hs.eg.db.sqlite -> org.Hs.eg.db

# AnnotationHub with 1 record
# snapshotDate(): 2022-10-31
# names(): AH107059
# $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
# $species: Homo sapiens
# $rdataclass: OrgDb
# $rdatadateadded: 2022-09-29
# $title: org.Hs.eg.db.sqlite
# $description: NCBI gene ID based annotations about Homo sapiens
# $taxonomyid: 9606
# $genome: NCBI genomes
# $sourcetype: NCBI/ensembl
# $sourceurl: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, ftp://ftp.ensembl.org/pub/current_fasta
# $sourcesize: NA
# $tags: c("NCBI", "Gene", "Annotation") 
# retrieve record with 'object[["AH107059"]]' 

# 或者记录 AH107059 ID


#https://www.jianshu.com/p/47b5ea646932

if(FALSE){
  query(org, "Zea mays")
  hub <- AnnotationHub::AnnotationHub()
  Zea.OrgDb <- hub[["AH107469"]]
  columns(Zea.OrgDb)
  head(keys(Zea.OrgDb, keytype = "SYMBOL"))
  saveDb(Zea.OrgDb, file = "Zea.OrgDb")
  #loadDb(file = "Zea.OrgDb")
}

if(FALSE){
  query(org, "Oryza sativa")
  hub <- AnnotationHub::AnnotationHub()
  Oryza.OrgDb <- hub[["AH107685"]]
  columns(Oryza.OrgDb)
  head(keys(Oryza.OrgDb, keytype = "ACCNUM"))
  saveDb(Oryza.OrgDb, file = "Oryza.OrgDb")
  #loadDb(file = "Oryza.OrgDb")
}


if(FALSE){
  query(org, "Triticum aestivum")
  hub <- AnnotationHub::AnnotationHub()
  Triticum.OrgDb <- hub[["AH107386"]]
  columns(Triticum.OrgDb)
  keys(Triticum.OrgDb)
  saveDb(Triticum.OrgDb, file = "Triticum.OrgDb")
  #loadDb(file = "Triticum.OrgDb")
}




if(FALSE){
query(org, "Solanum lycopersicum")
hub <- AnnotationHub::AnnotationHub()
Solanum.OrgDb <- hub[["AH107912"]]
keys(Solanum.OrgDb)
saveDb(Solanum.OrgDb, file = "Solanum.OrgDb")
#loadDb(file = "Solanum.OrgDb")
}
#用AnnotationHub获取非模式物种注释信息
#https://www.bioinfo-scrounger.com/archives/512/

