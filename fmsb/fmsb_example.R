rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\fmsb')

df <- read.csv('xx.tsv',sep='\t',row.names=1)
#加载包
library(fmsb)
#加入填充色
library(scales)
col<- rainbow(5)
library(scales)
colors_in <- alpha(col,0.3)
#绘图
png('radar_chart.png',width =1200,height =800,res=120,units = "px")

a<-dev.cur()#记录png device
pdf('radar_chart.pdf',width =12,height =8)
dev.control('enable')#打开图形设备控制

radarchart(df,#数据
           # Customize the polygon
           pcol=c("#00AFBB", "#E7B800", "#FC4E07") ,#rainbow(5),#多边形特征：线的颜色
           pfcol=scales::alpha(c("#00AFBB", "#E7B800", "#FC4E07"),0.5) ,#colors_in,#多边形特征：填充色
           #plwd=2,#多边形特征：线宽
           #plty=2,#多边形特征：线形
           # Customize the grid
           cglcol='grey',#网格特征:网格颜色
           cglty=1,#网格特征:网格线形
           cglwd=0.8,#网格特征:网格线宽
           axistype=0,#坐标轴类型
           #axislabcol='red',#网格特征:轴颜色
           #caxislabels=seq(0,50,5),#网格特征:轴范围
          
           vlcex=0.65#组标签大小
           )
#Add an horizontal legend
legend(x=1.2, y=1.3, legend = rownames(df[-c(1,2),]), 
       bty = "n", pch=20 , col=rainbow(5) , 
       text.col = "black", 
       cex=1.0, # 文字大小
       pt.cex=3)

dev.copy(which=a)#复制pdf的图形给png

dev.off()#保存pdf
dev.off()#保存png


# # Add an horizontal legend
# legend(
#   x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
#   bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
#   text.col = "black", cex = 1, pt.cex = 1.5
# )
