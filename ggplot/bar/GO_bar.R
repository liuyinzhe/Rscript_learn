library(ggplot2)

# 代码修改自 https://www.jianshu.com/p/9a9655dd83dc

setwd("Rscript/GO_bar")
go = read.csv("GO_enrich.tsv",header=T,stringsAsFactors = F,sep="\t")
# 1.按照qvalue升序排序，分别选出前10个BP,CC,MF的条目，由于enrichGO函数生成的数据框默认是按照qvalue升序排序，
# 所以这里我们只用选取前十个就行了

go_MF<-go[go$ONTOLOGY=="MF",][1:10,]
go_CC<-go[go$ONTOLOGY=="CC",][1:10,]
go_BP<-go[go$ONTOLOGY=="BP",][1:10,]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                         Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                         GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                         type=factor(
                            c(rep("biological process", 10),
                            rep("cellular component", 10),
                            rep("molecular function",10)),
                            levels=c("molecular function", "cellular component", "biological process")
                            )
                        )


## numbers as data on x axis

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

## shorten the names of GO terms

shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
## shorten the names of GO terms
# shorten_names <- function(x, n_word=4) {
#     if (length(strsplit(x, " ")[[1]]) > n_word){
#         return(paste(paste(strsplit(x, " ")[[1]][1:n_word], collapse=" "), "...", sep=""))
#     }else{
#         return(x)
#     }
# }

labels=(sapply(
  levels(as.factor(go_enrich_df$Description))[as.numeric(as.factor(go_enrich_df$Description))],
  shorten_names))
  
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange

CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")+
  theme(
    legend.position = "top",
    legend.title=element_blank()
  )
#coord_flip(...)横向转换坐标：把x轴和y轴互换，没有特殊参数

p
