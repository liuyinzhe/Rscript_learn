library(ggplot2)
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\LV2\\class')

data = read.table('data\\KEGG_bar_plot.tsv',header = T,sep='\t',check.names=F)

#条件创建label,0 填写NA
data$label[data$gene_count == 0 ] = NA
data$label[data$gene_count != 0 ] = data$gene_count

library(tidyverse)
#数据框转换成tibble
kegg_data<-as_tibble(data)
#R 语言中的 class() 函数用于返回用作参数的数据类
class(kegg_data)

if(FALSE){
#按照pvalue -log10(p_value)升序排序
  
#原始数据筛选（kegg_class,description，p_value)散列，按照kegg_class，-log10(p_value)排序
data<-kegg_data%>%select(kegg_class,description,p_value)%>%arrange(kegg_class,desc(-log10(p_value)))
#画图时改变geom_bar的自动排序
data$description<-factor(data$description,levels = unique(data$description),ordered = T)
}

if(TRUE){
  #按照gene_count 升序排序
  
  #原始数据筛选（kegg_class,description，gene_count)散列，按照kegg_class，gene_count排序
  data<-kegg_data%>%select(kegg_class,description,gene_count)%>%arrange(kegg_class,desc(gene_count))
  #画图时改变geom_bar的自动排序
  data$description<-factor(data$description,levels = unique(data$description),ordered = T)
}


#作图
P<- ggplot(data)+
  geom_bar(aes(x=description,y=gene_count,fill=kegg_class),stat = 'identity')+
  labs(x="Description",
       y="Number of Gene",
       fill="Classification",
       colour = "Number of cylinders",
  ) + ##title="KEGG pathway anotation",
  coord_flip() # 翻转




ggsave('kegg_leve2_class.png',plot=P, height=12, width=12, units="in",dpi=600)
ggsave('kegg_leve2_class.pdf',plot=P, height=12, width=12, units="in",dpi=600)


#https://zhuanlan.zhihu.com/p/497293882
#颜色按categroy分类、并且p_value由小到大排序
#https://www.jb51.net/article/208762.htm
