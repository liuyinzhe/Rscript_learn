rm(list=ls())#clear Global Environment
options(stringsAsFactors = F)
setwd("/data/draw/Rscript/ggplot2/")
#加载包
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(aplot) # Decorate a 'ggplot' with Associated Information

#加载数据
df <- read.table("data.txt",header = 1,check.names = F)
df$group <- factor(df$group,levels = c("F","CK")) #若不指定参数levels，则因子水平默认按字母顺序。
df$species <- factor(df$species,levels = c("Pedobacter","Aridibacter","Devosia","Rhizobium",
                                           "Phenylobacterium","Arthrobacter","Bradvrhizobium",
                                           "Pseudomonas","Gemmatimonas","Sphingomonas"))

##绘制热图
#设置渐变色
col <- colorRampPalette(c("#0066b2","#fdbd10","#ec1c24"))(50) #设置渐变色
p1 <- ggplot(df,aes(group,species,fill=value))+
  #绘制热图
  geom_tile(color="black")+
  #添加数值
  geom_text(aes(label = value), color = 'white', size = 5)+
  #标题
  labs(x=NULL,y=NULL,fill=NULL)+
  #颜色
  scale_fill_gradientn(colours = col)+
  #主题
  theme_void()+
  theme(axis.text.x = element_text(color = "black",size=12),
        axis.text.y = element_text(color = "black",size=12,hjust = 1),
        legend.position = "right")#去除图例
p1

##绘制柱状堆积图
p2 <- ggplot(df,aes(species,value,fill=group))+
  #柱状堆积图
  geom_col()+
  coord_flip()+
  #标题
  labs(x=NULL,y=NULL,fill=NULL)+
  #颜色
  scale_fill_manual(values = c("#f0b240","#62d7f6"))+
  #主题
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",size=12),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = "black",linewidth = 0.8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(-0.15, "cm"),
        axis.ticks.x = element_line(color = "black",linewidth = 0.8),
        legend.position = "right")+#图例设置
  scale_y_continuous(expand = c(0,0),breaks = c(0,5,10,15,20,25))#设置刻度从0开始
p2

##拼图
draw_obj <- p1%>%insert_right(p2,width = 2)


png(paste("heatmap.bar_stacked",".png",sep=""),width =1350,height =900,res=150,units = "px")
print(draw_obj)
dev.off()
pdf(paste("heatmap.bar_stacked",".pdf",sep=""),width =10,height =7)
print(draw_obj)
dev.off()
