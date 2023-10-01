#!/usr/bin/Rscript
rm(list = ls())
options(stringsAsFactors = F)
setwd('/data/Rscript/bar')

library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(dplyr)
library(scales)
# colors
#define custom color scale
#a<-pal_igv("default", 1)

#show_col(pal_igv(palette = c("default"), alpha = 1)(51))
#pal_igv(palette = c("default"), alpha = 1)(51)


#myColors_map <- brewer.pal(12, "Paired")

#自定义颜色，name对应颜色
#order	name  color
name_lst<- read.table("target_name.tsv",comment.char = "$",sep='\t',header = TRUE)

# order name #根据name列排序，因为最后按照name顺序对应颜色分配
name_lst<- name_lst %>% arrange(name)


md<- read.table("input.tsv",header = TRUE,comment.char = "",stringsAsFactors = FALSE,sep="\t")
# Sample	Taxa	Percent	group

# order Taxa #Taxa 对应name，按照Taxa排序
md<- md %>% arrange(Taxa)
# 构造dataframe
df <- data.frame( color=name_lst$color, name=name_lst$name)
# 根据 taxa name 筛选 dataframe 行
new_df <- df %>% filter(.,name %in% unique(md$Taxa))
# 获得需要的颜色
myColors <- new_df$color


#custom_colors <- scale_colour_manual(values = c("#EA4732","#F1784B","#FDA94A","#FDA039","#E29A72","#DAA192","#C88B4A","#BC7239"))
# 指定填充颜色 ，custom_colors 作为函数封装
custom_colors <- scale_fill_manual(values = myColors)



# 读取分组，列重命名
group<- read.table('group.txt',comment.char = "",header = TRUE,stringsAsFactors = FALSE)
colnames(group) <- c("sampleID","group")


list = name_lst$name



xt = 8
p<- ggplot(data = md,aes(Sample,weight=Percent,fill=factor(Taxa,levels = list)))+
    geom_bar()+theme_gray()+custom_colors+
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=xt,angle = 90,hjust = 0.5))+
    ylab("Absolute Abundance")+#coord_cartesian(ylim = c(0.04,50))+
    facet_grid(.~group,scales = 'free',space = 'free')

  

print(p)
png(paste("out",".png",sep=""),width =900,height =600,res=100,units = "px")
print(p)
dev.off()
pdf(paste("out",".pdf",sep=""),width =10,height =7)
print(p)
dev.off()



