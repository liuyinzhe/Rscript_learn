

rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\boxplot')
library(tidyr) 
library(ggpubr)
library(ggplot2)
library(ggsci)
library(dplyr)
library(patchwork)

data <- read.csv("index.txt", header = T,sep='\t')
head(data)
# sample	observed_species	chao1	shannon	simpson	goods_coverage	group

# 转换为长格式
long_data <- data %>%
  pivot_longer(
    cols = c(observed_species, chao1, shannon, simpson, goods_coverage),
    names_to = "type",
    values_to = "value"
  )

long_data$type=as.factor(long_data$type)
long_data$group=as.factor(long_data$group)

df_1<-filter(long_data,type=="observed_species")
df_2<-filter(long_data,type=="chao1")
df_3<-filter(long_data,type=="shannon")
df_4<-filter(long_data,type=="simpson")
#df_5<-filter(long_data,type=="goods_coverage")



# png(filename = "igg.png",width=2600,height=2600,res=500)
# your_font_size <- 5 # 改变box上p值大小的
# coreggplot(data = dt2, aes(x = time, y = value,fill = group)) +
# geom_boxplot(outlier.colour = NA) +
# geom_dotplot(binaxis='y', stackdir='center', binwidth = 100, position = position_dodge(0.75))+
# labs(fill = "group", x = "time", y = "Value") +
# scale_y_continuous(lim = c(0, 8000)) +
# theme_classic(base_size = 20) +
# theme(legend.position = c(0.9, 0.7))+
# stat_compare_means(method = "wilcox.test",label = "p.format", size = your_font_size)
# dev.off()

your_font_size <- 5 # 改变box上p值大小的
my_comparisons = list( c("A", "C"))



p1 <- ggplot(df_1, aes(x=group,y=value,fill=group)) + 
  stat_boxplot(geom='errorbar',width=0.1)+ # 添加上下限
  geom_boxplot()+
  geom_jitter(width=0.1, size=1) + # 再叠加散点
  #stat_boxplot()+
  theme_bw()+ theme(axis.text.x=element_text(angle=0,vjust=1,hjust=1))+ #base_family="Arial"
  #scale_fill_npg()+
  scale_fill_manual( # 自定义散点颜色
    values = c("#F0A3FF","#0075DC")
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     method = "t.test", # "wilcox.test", #
                     hide.ns=F,
                     #label.y=92,
                     size = your_font_size)+
  #stat_compare_means(method = "t.test",label = "p.signif")#wilcox.test#label = "p.format", size = your_font_size
  labs(x="",y="",fill="Group",title="")+
  theme(plot.title = element_text(hjust = 0.5),legend.position ="none")+

  facet_wrap(~type) 
  #theme()+#legend.position ="none"
  #theme(axis.text.x=element_text(angle=45))#vjust=1,hjust=1,size=8#,text=element_text(family="Arial")
  
#print(p1)

p2 <- ggplot(df_2, aes(x=group,y=value,fill=group)) + 
  stat_boxplot(geom='errorbar',width=0.1)+ # 添加上下限
  geom_boxplot()+
  geom_jitter(width=0.1, size=1) + # 再叠加散点
  #stat_boxplot()+
  theme_bw()+ theme(axis.text.x=element_text(angle=0,vjust=1,hjust=1))+ #base_family="Arial"
  #scale_fill_npg()+
  scale_fill_manual( # 自定义散点颜色
    values = c("#F0A3FF","#0075DC")
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     method = "t.test", # "wilcox.test", #
                     hide.ns=F,
                     #label.y=92,
                     size = your_font_size)+
  #stat_compare_means(method = "t.test",label = "p.signif")#wilcox.test#label = "p.format", size = your_font_size
  labs(x="",y="",fill="Group",title="")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~type)



p3 <- ggplot(df_3, aes(x=group,y=value,fill=group)) + 
  stat_boxplot(geom='errorbar',width=0.1)+ # 添加上下限
  geom_boxplot()+
  geom_jitter(width=0.1, size=1) + # 再叠加散点
  #stat_boxplot()+
  theme_bw()+ theme(axis.text.x=element_text(angle=0,vjust=1,hjust=1))+ #base_family="Arial"
  #scale_fill_npg()+
  scale_fill_manual( # 自定义散点颜色
    values = c("#F0A3FF","#0075DC")
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     method = "t.test", # "wilcox.test", #
                     hide.ns=F,
                     #label.y=92,
                     size = your_font_size)+
  #stat_compare_means(method = "t.test",label = "p.signif")#wilcox.test#label = "p.format", size = your_font_size
  labs(x="",y="",fill="Group",title="")+
  theme(plot.title = element_text(hjust = 0.5),legend.position ="none")+
  facet_wrap(~type)
  


p4 <- ggplot(df_4, aes(x=group,y=value,fill=group)) + 
  stat_boxplot(geom='errorbar',width=0.1)+ # 添加上下限
  geom_boxplot()+
  geom_jitter(width=0.1, size=1) + # 再叠加散点
  #stat_boxplot()+
  theme_bw()+ theme(axis.text.x=element_text(angle=0,vjust=1,hjust=1))+ #base_family="Arial"
  #scale_fill_npg()+
  scale_fill_manual( # 自定义散点颜色
    values = c("#F0A3FF","#0075DC")
  ) +
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     method = "t.test", # "wilcox.test", #
                     hide.ns=F,
                     #label.y=92,
                     size = your_font_size)+
  #stat_compare_means(method = "t.test",label = "p.signif")#wilcox.test#label = "p.format", size = your_font_size
  labs(x="",y="",fill="Group",title="")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~type)


p<- p1 + p2 + p3 + p4  + plot_layout(ncol = 2) + 
  plot_annotation(title = "",
                  theme=theme(plot.title = element_text(size = 12,hjust = 0.5))
                  )#theme(plot.title = element_text(size = 16))))
  #plot_annotation(title = "A closer look at the effect of drive train in cars",caption = "Source: mpg dataset in ggplot2")
  ggsave('boxplot.t-test.pvalue.png',plot=p,height=12, width=15, units="in", dpi=400) # 
  ggsave('boxplot.t-test.pvalue.pdf',plot=p,height=12, width=15, units="in",dpi=600)#
p
  