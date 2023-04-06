rm(list = ls())
library(ggplot2)
setwd("D:\\Rscript\\wilcox")


# 手动处理保证样品数量一致
A_data=read.csv('A.txt',header=T,sep="\t",check.names=F)
B_data=read.csv('B.txt',header=T,sep="\t",check.names=F)

#sample  group factor1 factor2 factor3
#x1  groupA  0.3 0.2 0.3
#x2  groupA  0.6 0.4 0.5
#x3  groupA  0.8 0.5 0.6

A_data<-A_data[,3:ncol(A_data)]
B_data<-B_data[,3:ncol(B_data)]


#循环参考自https://cloud.tencent.com/developer/ask/sof/884030

# 创建空frame,行数
result <- as.data.frame(matrix(as.numeric(), nrow=ncol(A_data), ncol=0))
#名字加进去
result$info<-colnames(A_data)
# pvalue
result$pvalue<-NA

val1 <- as.numeric() 
val2 <- as.numeric()


p_value=c() #设置一个空数组用于储存p值
# unlist () 函数用于将列表转换为向量
#wilcox.test(c(1,3,4,5,6),c(3,2,1,4,5), paired = F)

for (i in 1:ncol(A_data)) {
  head(i)
  val1 <- as.numeric(unlist(A_data[ , i]))  
  val2 <- as.numeric(unlist(B_data[ , i]))
  res <- wilcox.test(val1, val2, paired =F,alt='two.sided',exact=F,correct=T)
  # i 行,2列 加入
  result[i, 2] <- as.numeric(res$p.value)
  p_value[i]=res$p.value
}
# FDR 矫正
#https://zhuanlan.zhihu.com/p/462036454
q_value<-p.adjust(
  p_value,      # P值列表
  method ="BH"  # FDR校正的方法
)

result<-data.frame(result,q_value) #添加p值到数据最后一列


write.table(result, file ="wilcox_test.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


