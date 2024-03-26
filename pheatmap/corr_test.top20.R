
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\draw')

#ggplot2 绘制带星号和 Pvalue 值的相关系数热图
#https://zhuanlan.zhihu.com/p/102284049

library(psych)  #psych包用于计算相关性、p值等信息
library(pheatmap) #用于绘制热图
library(reshape2) #reshape2包用于输出数据的整合处理
library(ggplotify) # 转为 ggplot2 对象

 #T1
df1 <- read.csv("species_abundance.new.txt",header=T,sep = '\t',row.names = 1) #读取微生物丰度信息表 #row.names = 1, row.names='taxonomy'
df2 <- read.csv("metabolome_BD.top20.tsv",header=T,sep = '\t',row.names = 1) #读入化合物丰度列表 # row.names = 1,,check.names=F

df1_t <- t(df1)
df2_t <- t(df2)


res <- corr.test(df1_t,df2_t,method = "spearman") #method可选“pearson”、“spearman”、“kendall” #,alpha = 0.05



result_p <- res$p #提取p值
result_r <- res$r #提取cor值

# NA空数据处理
result_p[is.na(result_p)] = 1
result_r[is.na(result_r)] = 0

p.out<-cbind(rownames(result_p),result_p) #输出p值
r.out<-cbind(rownames(result_r),result_r) #输出cor值

#write.csv(p.out,"pout.csv",row.names = F,fileEncoding='GB2312')
#write.csv(r.out,"rout.csv",row.names = F,fileEncoding='GB2312')

write.table(p.out, file ="pout.top20.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
write.table(r.out, file ="rout.top20.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# 宽转长
df <-melt(result_r,value.name="cor")
df$pvalue <-as.vector(result_p)  #宽格式转长格式


df3 <- subset(df,abs(df$cor)>0.75&df$pvalue<0.05)#筛选
#write.csv(df,"melt_p_r.csv",row.names = F,fileEncoding='GB2312')#输出长格式
#write.csv(df3,"melt_p005_r075.csv",row.names = F,fileEncoding='GB2312')

write.table(df, file ="melt_p_r.top20.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


write.table(df3, file ="melt_p005_r075.top20.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


#result_r_sig <- result_r

if (!is.null(result_p)){
  sssmt <- result_p< 0.001
  result_p[sssmt] <-'***'
  ssmt <- result_p< 0.01
  result_p[ssmt] <-'**'
  smt <- result_p >0.01& result_p <0.05
  result_p[smt] <- '*'
  result_p[!sssmt&!ssmt&!smt]<- ''
} else {
  result_p <- F
}  #判断p值大小，若p<0.01，则'**'，若0.01<p<0.05,则'*'，否则无显著。




#result_p<-read.csv("pout_new.csv",header=T,sep = '\t',row.names = 1)

#result_p_r<- paste(result_p,result_r,sep='\n')

P<- pheatmap(result_r,scale = "none",#"raw"
             # 聚类
             cluster_rows = TRUE,# 行不聚类
             cluster_cols = TRUE,
             # 注释列
             #annotation_row = annotation_row,
             #annotation_col = annotation_col,
             #annotation_colors=ann_colors,
             
             # 单元格
             #cellwidth = 50,
             #cellheight = 45, # 单元格指定宽高
             
             # 边框
             border=FALSE,#NA,
             border_color="white",
             

             
             # 文字大小
             fontsize = 10, # 设置所有字体的大小，包括图例文字的，其它不同字体根据特定参数设置
             fontsize_row = 10,
             fontsize_col = 10,
             
             fontsize_number = 12, #热图上数值的字体大小
             # 文字格式
             display_numbers = result_p,#TRUE,#result_p, 
             #number_format = "%.2f",
             
             # 角度
             angle_col = 45,
             #filename = "相关性热图.jpeg",
             number_color = "black") #根据自己的喜好绘制热图,cluster_row = F,number_color = "white"

P
# 转换为ggplot对象，调整右边距
# https://www.jianshu.com/p/0ddd67343bdb
pobj <- as.ggplot(P) + theme(plot.margin = margin(l = 5, r = 5,t=5,b=5 )) 
pobj

ggsave('all_heatmap.top20.png',plot=pobj, height=12, width=24, units="in",dpi=600, bg = "white")

ggsave('all_heatmap.top20.pdf',plot=pobj, height=12, width=24, units="in",dpi=600)



#通过R分析物种与代谢物的相关性及绘制热图
#https://www.jianshu.com/p/0cae60aa5d54

#原文链接：https://blog.csdn.net/weixin_41772125/article/details/112285185



