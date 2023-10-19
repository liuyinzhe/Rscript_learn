#!/usr/bin/Rscript
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\pheatmap')

library(gplots)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(ggtext)
library(glue)
df <- read.csv('phylum_abundance.sorted.new.txt',sep="\t", row.name = 1)

# 行注释信息
annotation_row <- read.csv('phy_class_group_row.txt',sep="\t", row.name = 1)
# 列注释信息
annotation_col <- read.csv('group_sample_col.txt',sep="\t", row.name = 1)

ann_colors = list(
  Group= c(A='#3B8841',B='#EE316D'),
  environment = c(X1 = "slateblue3", X2 = "red2"),
  Phylum = c(Phylum = "#4DB858")
)

# https://zhuanlan.zhihu.com/p/366674882

p <- pheatmap(df, scale = "row", ##AF3160
              color = colorRampPalette(c("#AF3160", "white", "#278C43"))(50),
              # 聚类
              cluster_rows = TRUE,# 行不聚类
              cluster_cols = FALSE,
              # 注释列
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors=ann_colors,
              # 单元格
              cellwidth = 50,
              cellheight = 25, # 单元格指定宽高
              
              # 图例
              legend = TRUE, # 隐藏图例
              annotation_names_col = TRUE,
              # 标签
              
              #fontsize_row = 8,
              #fontsize_col = 12,
              #fontsize = 8, #统一字体大小
              
              angle_col = 0,
              # 边框
              #border_color = "white",#边框颜色，白色
              border = FALSE, # 是否保留边框，FALSE 则是紧凑型
              ) 
# 转换为ggplot对象，调整右边距
# https://www.jianshu.com/p/0ddd67343bdb
pobj <- as.ggplot(p) + theme(plot.margin = margin(l = 5, r = 5, t= 5, b= 5)) #+theme(axis.text.y = element_markdown())
pobj

png(paste("species.pheatmap",".png",sep=""),width =5500,height =5500,res=600,units = "px") 
print(pobj)
dev.off()
pdf(paste("species.pheatmap",".pdf",sep=""),width =10,height =10)
print(pobj)
dev.off()
