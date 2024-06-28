setwd('.')
rm(list = ls())
# Load package
library(dplyr)
library(vegan)
library(ggplot2)
library(ggthemes)
# Load data
#raw_otu<- read.table("format_data.txt",row.names = 1,header = T,sep="\t")


#otu<-read.table('otu.txt',row.names = 1,header = T)
# otu <- raw_otu %>% select(-group)

#data_comm <- md[,3:length(colnames(md))]


md<- read.table("format_data.txt",header = T,sep="\t")
data_comm <- md[,3:length(colnames(md))]

#ADONIS
adonis_result<-adonis2(data_comm ~ group, md, permutations = 999, method = "bray")
data_adonis <- paste0("adonis R2: ",round(adonis_result$R2,2), "; P-value: ", adonis_result$`Pr(>F)`)
# data_adonis 用于标注

#pcoa
# vegdist函数，计算距离；method参数，选择距离类型
distance <- vegdist(data_comm, method = 'bray')
# 对加权距离进行PCoA分析
pcoa <- cmdscale(distance, k = (nrow(data_comm) - 1), eig = TRUE)

## plot data
# 提取样本点坐标
plot_data <- data.frame({pcoa$point})[1:2]

# 提取列名，便于后面操作。
plot_data$ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')


## plot data
# 提取样本点坐标
plot_data <- data.frame({pcoa$point})[1:2]

# 提取列名，便于后面操作。
plot_data$ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig

#为样本点坐标添加分组信息
plot_data$group <-md$group
#plot_data <- merge(plot_data, group, by = 'ID', all.x = TRUE)
#head(plot_data)

# figure1
# ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, fill=group)) + # color = group
#   geom_point(shape = 21,color = 'black',size=4) +
#   scale_fill_manual(values = c('#73bbaf','#d15b64','#592c93','blue'))+
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
#   geom_hline(yintercept=0, linetype=4) +    
#   geom_vline(xintercept=0 ,linetype=4)+   
#   stat_ellipse(level = 0.95,alpha =0.1)+
#   theme_few()+
#   theme(legend.position = c(0.9, 0.2),
#         legend.title = element_blank(),
#         legend.background = element_rect(colour ="black"))

# https://zhuanlan.zhihu.com/p/442649228

# 圆圈图

ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, color=group)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=data_adonis) +
  geom_point(size=4) + 
  stat_ellipse(level=0.6) +
  theme_classic()

ggsave('pcoa_circle.pdf',width = 4,height = 4)

# 覆盖图
library(ggalt)
ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, color=group)) +
  labs(
    x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
    y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
    title=data_adonis) +
  geom_point(size=5) + 
  geom_encircle(aes(fill=group), alpha = 0.1, show.legend = F) +
  theme_classic() + 
  coord_fixed(1)

ggsave('pcoa_shadow.pdf',width = 4,height = 4)
