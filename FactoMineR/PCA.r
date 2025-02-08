rm(list = ls())
options(stringsAsFactors = F)

setwd('D:\\Rscript\\PCA')
#标识样本名称，使用 ggplot2 的拓展包 ggrepel 来完成
library(ggrepel)

#读取基因表达值矩阵
#推荐使用 log 转化后的基因表达值，降低不同基因表达水平数量级相差过大的问题
gene <- read.delim('all.counts.txt', row.names = 1, sep = '\t')
# gene <- log2(gene + 1)
#将基因表达值矩阵作个转置，使行为样本，列为基因
gene <- t(gene)

#我们使用 FactoMineR 包中的方法，实现 PCA 分析和聚类添加
library(FactoMineR)

#样本中基因表达值的 PCA 分析
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
# plot(gene.pca)  #PCA 简图


#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )



#读取并合并样本分组信息
group <- read.delim('group.txt', row.names = 1, sep = '\t', check.names = FALSE,colClasses=c("character","character"))
group <- group[rownames(pca_sample), ]
pca_sample <- cbind(pca_sample, group)

pca_sample$samples <- rownames(pca_sample)
head(pca_sample)  #作图数据中包含了样本坐标和分组信息

#ggplot2 绘制二维散点图
library(ggplot2)

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
   geom_point(aes(color = group), size = 3) +  #根据样本坐标绘制二维散点图
   scale_color_manual(values = c('#91AEB0','#776C6B','#E3B386','#DC8E75')) +  #自定义颜色
  geom_text_repel(aes(label = samples), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(color = 'white', fill = 'transparent'),
    legend.text = element_text(size = 18)
    
    ) + 
   labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

# 95%置信椭圆
# p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
#   scale_fill_manual(values = c('#91AEB0','#776C6B','#E3B386','#DC8E75'))


#多边形连接同类别对象边界的样式，适用于各组样本数大于 3 个的情况
library(plyr)
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])

#p+geom_polygon(data = cluster_border, aes(color = group), fill = NA, show.legend = FALSE)

P<-p+geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE)
P
ggsave('PCA.group.png',plot=P, height=12, width=12, units="in",dpi=300)
ggsave('PCA.group.pdf',plot=P, height=12, width=12, units="in",dpi=300)
