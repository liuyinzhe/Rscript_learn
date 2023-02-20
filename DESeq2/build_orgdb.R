rm(list = ls())
options(stringsAsFactors = F)
setwd('/data/home/build_db')

library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
options(stringsAsFactors = F)#设置一下options
#https://www.jianshu.com/p/9d38cb9f1d02

# run.01.emapper.sh
# cd /data/home/06.emapper_annotations
# mkdir -p output
# unset PYTHONPATH
# /data/home/liuyinzhe/envs/rna_seq/bin/emapper.py \
#  -m diamond \
#  -i GCF_000014225.1_ASM132v1_protein.faa  \
#  --itype proteins \
#  --cpu 30 \
#  --output_dir /data/home/liuyinzhe/project/RNA_seq/RS23IMWY02/06.emapper_annotations/output \
#  --excel \
#  -o Sao  >log 2>err

# 对 emapper 注释结果 Sao.emapper.annotations 处理，选取特定列
# sed '/^##/d' *.emapper.annotations| sed 's/#//g'| awk -vFS="\t" -vOFS="\t" '{print $1,$9,$10,$12}' > Sao.annotations

emapper <- read.table("Sao.annotations",header = TRUE,sep = "\t",quote = "",stringsAsFactor = FALSE)
#将空值替换为NA，方便后续使用na.omit()函数提出没有注释到的行
emapper[emapper==""]<-NA
emapper[emapper=="-"]<-NA

# library(dplyr)
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()

gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())


gos_list <- function(x){
  the_gos <- str_split(x[2], ",", simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA", length(the_gos)))
  return(df_temp)
}
gene2gol <- apply(as.matrix(gos),1,gos_list)
gene2gol_df <- do.call(rbind.data.frame, gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)

gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2ko<-na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko),1,gos_list)
gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)

# 从以下网址，点击 Download json 下载
# https://www.genome.jp/kegg-bin/get_htext?ko00001

# 或者 以下命令直接下载
# curl 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=' \
#   -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' \
#   -H 'Accept-Language: zh-CN,zh;q=0.9,en;q=0.8,en-CA;q=0.7,en-AU;q=0.6,en-US;q=0.5' \
#   -H 'Connection: keep-alive' \
#   -H 'Referer: https://www.genome.jp/kegg-bin/get_htext?ko00001' \
#   -H 'Sec-Fetch-Dest: document' \
#   -H 'Sec-Fetch-Mode: navigate' \
#   -H 'Sec-Fetch-Site: same-origin' \
#   -H 'Sec-Fetch-User: ?1' \
#   -H 'Upgrade-Insecure-Requests: 1' \
#   -H 'User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36 Edg/110.0.1587.50' \
#   -H 'sec-ch-ua: "Chromium";v="110", "Not A(Brand";v="24", "Microsoft Edge";v="110"' \
#   -H 'sec-ch-ua-mobile: ?0' \
#   -H 'sec-ch-ua-platform: "Windows"' \
#   --compressed


# library(jsonlite)
# 下面的json = "ko00001.json"，如果你下载到其他地方，记得加上路径
update_kegg <- function(json = "ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  kegg <- fromJSON(json)
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}}
  # 删除 NA 行
  ko2pathway <- na.omit(ko2pathway)
  pathway2name <- na.omit(pathway2name)
  write.table(pathway2name,file = "./pathway2name.txt", sep = "\t", quote = F, row.names = F)
  write.table(ko2pathway,file = "./ko2pathway.txt", sep = "\t", quote = F, row.names = F)
  save(pathway2name, ko2pathway, file = "kegg_info.RData")}
# 调用函数后在本地创建kegg_info.RData文件，以后只需要载入 "kegg_info.RData"即可
update_kegg()
# 载入kegg_info.RData文件
load(file = "kegg_info.RData")



gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>% na.omit()



################### STEP6： 制作自己的Orgdb ##############################

# 查询物种的Taxonomy，例如要查sesame
# https://www.ncbi.nlm.nih.gov/taxonomy/?term=sesame
# https://www.ncbi.nlm.nih.gov/taxonomy 这里搜拉丁名
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=40324

tax_id = 1234
genus = "Escherichia " 
species = "coli"
#################### 


#gene2go <- unique(gene2go)
#gene2go <- gene2go[!duplicated(gene2go),]
#gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
#gene_info <- gene_info[!duplicated(gene_info),]

#colnames(gene_info)[1]="GID"
#colnames(gene2go)[1]="GID"
#colnames(gene2ko)[1]="GID"

#https://www.thinbug.com/q/28526860
#gene_info <- t(do.call(rbind, gene_info))
#x<-simplify2array(gene_info)
#type(x)
#x
#Reduce( rbind, lapply(gene_info, t) )

type(gene_info)
#gene_info <-simplify2array(gene_info)
#gene2go <- simplify2array(gene2go) 
#gene2ko <- simplify2array(gene2ko)
#gene2pathway <- simplify2array(gene2pathway)


#type(gene_info)
#gene_info
#write.table(gene_info, file ="gene_info.tsv", sep="\t",row.names = FALSE,col.names =TRUE, quote =TRUE)
#write.table(gene2go, file ="gene2go.tsv", sep="\t",row.names = FALSE,col.names =TRUE, quote =TRUE)
#write.table(gene2ko, file ="gene2ko.tsv", sep="\t",row.names = FALSE,col.names =TRUE, quote =TRUE)
#write.table(gene2pathway, file ="gene2pathway.tsv", sep="\t",row.names = FALSE,col.names =TRUE, quote =TRUE)

#gene_info <- read.table("gene_info.tsv",header = TRUE,sep = "\t",quote = "",stringsAsFactor = FALSE)
#gene2go <- read.table("gene2go.tsv",header = TRUE,sep = "\t",quote = "",stringsAsFactor = FALSE)
#gene2ko <- read.table("gene2ko.tsv",header = TRUE,sep = "\t",quote = "",stringsAsFactor = FALSE)
#gene2pathway <- read.table("gene2pathway.tsv",header = TRUE,sep = "\t",quote = "",stringsAsFactor = FALSE)

#x <- colnames(gene_info)
#type(x)
#x[1]
#colnameGIDs <- sapply(gene_info, function(x){colnames(x)[1]})
#colnameGIDs[1]=="GID"
#any(colnameGIDs == "GID")

#if(any(colnameGIDs != "GID"))
#  stop("The 1st column must always be the gene ID 'GID'")

# 去重
# https://www.jianshu.com/p/9d38cb9f1d02

gene2go <- dplyr::distinct(gene2go)
gene2ko <- dplyr::distinct(gene2ko)
gene_info <- dplyr::distinct(gene_info)
gene2pathway <- dplyr::distinct(gene2pathway)


makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="1.0", 
               maintainer="smt <somethin@github.com>", 
               author="smt <something@github.com>", 
               outputDir=".",
               tax_id=tax_id,  
               genus=genus, 
               species=species,
               goTable="go")
