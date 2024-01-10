
library(UpSetR)


#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#library('BiocManager',quietly=TRUE)
#BiocManager::install('UpSetR')

setwd('D:\\Rscript')
#mydata <- read.table('flower.txt',sep='\t',header=T,row=)

head(mydata)

png(paste('upset_draw',".png",sep=""),bg = "white",width =1200,height =800,res=120,units = "px")
a<-dev.cur()#记录png device
pdf(paste('upset_draw',".pdf",sep=""),width =12,height =8)
dev.control('enable')#打开图形设备控制

upset(fromList(mydata),  # fromList一个函数，用于将列表转换为与UpSetR兼容的数据形式。
      nsets = 50,     # 绘制的最大集合个数
      nintersects = 40, #绘制的最大交集个数，NA则全部绘制
      order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree"根据集合的大小排序。
      keep.order = F, # 保持设置与使用sets参数输入的顺序一致。默认值是FALSE
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 1 # 文字标签的大小
)

dev.copy(which=a)#复制pdf的图形给png

dev.off()#保存pdf
dev.off()#保存png
