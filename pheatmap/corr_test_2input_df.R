
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\corr')
#library(Hmisc)
#res <- Hmisc::rcorr(as.matrix(df1),as.matrix(df2))


#library(showtext)
## 使用Windows自带字体#
#font_add("simhei", "C:\\Windows\\Fonts\\simhei.ttf")

# 为了方便 保存pdf 时保存中文
library(showtext)
showtext.auto(enable=TRUE)


library(psych)  #psych包用于计算相关性、p值等信息
library(pheatmap) #用于绘制热图
library(reshape2) #reshape2包用于输出数据的整合处理

 #T1
df1 <- read.csv("Clinical_information_T.tsv",header=T,sep = '\t',row.names = 1) #读取微生物丰度信息表 #row.names = 1, row.names='taxonomy'
df2 <- read.csv("otu_table_D5_t1.tsv",header=T,sep = '\t',row.names = 1) #读入化合物丰度列表 # row.names = 1,,check.names=F

#T2
#df1 <- read.csv("Clinical_information_T.tsv",header=T,sep = '\t',row.names = 1) #读取微生物丰度信息表 #row.names = 1, row.names='taxonomy'
#df2 <- read.csv("otu_table_D5_t2.tsv",header=T,sep = '\t',row.names = 1) #读入化合物丰度列表 # row.names = 1,,check.names=F


#No1
#df1 <- read.csv("Clinical_information_NO.tsv",header=T,sep = '\t',row.names = 1) #读取微生物丰度信息表 #row.names = 1, row.names='taxonomy'
#df2 <- read.csv("otu_table_D5_no1.tsv",header=T,sep = '\t',row.names = 1) #读入化合物丰度列表 # row.names = 1,,check.names=F


#No2
#df1 <- read.csv("Clinical_information_NO.tsv",header=T,sep = '\t',row.names = 1) #读取微生物丰度信息表 #row.names = 1, row.names='taxonomy'
#df2 <- read.csv("otu_table_D5_no2.tsv",header=T,sep = '\t',row.names = 1) #读入化合物丰度列表 # row.names = 1,,check.names=F


#is.numeric(df1[,2])
#res <- Hmisc::rcorr(as.matrix(df1),as.matrix(df2))

res <- corr.test(df1,df2,method = "pearson") #method可选“pearson”、“spearman”、“kendall” #,alpha = 0.05



result_p <- res$p #提取p值
result_r <- res$r #提取cor值
p.out<-cbind(rownames(result_p),result_p) #输出p值
r.out<-cbind(rownames(result_r),result_r) #输出cor值
write.csv(p.out,"pout.tsv",row.names = F,fileEncoding='GB2312')
write.csv(r.out,"rout.tsv",row.names = F,fileEncoding='GB2312')

result_p[is.na(result_p)] = 1
result_r[is.na(result_r)] = 0

df <-melt(result_r,value.name="cor")




df$pvalue <-as.vector(result_p)  #宽格式转长格式
df2 <- subset(df,abs(df$cor)>0.75&df$pvalue<0.05)#筛选
write.csv(df,"melt_p_r.csv",row.names = F,fileEncoding='GB2312')#输出长格式
write.csv(df2,"melt_p005_r075.csv",row.names = F,fileEncoding='GB2312')



if (!is.null(result_p)){
  ssmt <- result_p< 0.01
  result_p[ssmt] <-'**'
  smt <- result_p >0.01& result_p <0.05
  result_p[smt] <- '*'
  result_p[!ssmt&!smt]<- ''
} else {
  result_p <- F
}  #判断p值大小，若p<0.01，则'**'，若0.01<p<0.05,则'*'，否则无显著。


pheatmap(result_r,scale = "row",
         cluster_row = F, cluster_col = T,
         border=NA,display_numbers = result_p, 
         fontsize_row=5,
         filename = "相关性热图.jpeg",
         number_color = "black") #根据自己的喜好绘制热图,cluster_row = F,number_color = "white"


pdf("相关性热图.pdf")
showtext_begin()

pheatmap(result_r,scale = "row",
         cluster_row = F, cluster_col = T,  # 过多的NA  产生需要去掉 cluster_row的聚类 #报错：NA/NaN/Inf in foreign function call (arg 10)
         border=NA,display_numbers = result_p, 
         fontsize_row=5,
         number_color = "black") #根据自己的喜好绘制热图,cluster_row = F,number_color = "white"

showtext_end()
dev.off()

#通过R分析物种与代谢物的相关性及绘制热图
#https://www.jianshu.com/p/0cae60aa5d54

#原文链接：https://blog.csdn.net/weixin_41772125/article/details/112285185
