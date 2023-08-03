
library(dplyr) 
rm(list = ls())
options(stringsAsFactors = F)
setwd('~/box_ttest')
library(ggpubr)
library(ggplot2)
library(ggsci)



mydata <- read.csv("kegg.input.all.tsv", header = T,sep='\t')
head(mydata)


mydata$type=as.factor(mydata$type)
mydata$group=as.factor(mydata$group)



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

your_font_size <- 4.5 # 改变box上p值大小的
my_comparisons = list( c("Negative", "Positive"))

p <- ggplot(mydata, aes(x=type,y=value,fill=group)) + 
  geom_boxplot()+
  stat_boxplot()+
  theme_bw()+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+ #base_family="Arial"
  #scale_fill_nejm()+
  stat_compare_means(#comparisons = my_comparisons,
                     label = "p.signif",
                     method = "t.test",
                     hide.ns=T,
                     size = your_font_size)+
  #stat_compare_means(method = "t.test",label = "p.signif")#wilcox.test#label = "p.format", size = your_font_size
  #labs(x="",y="Metabolism",fill="Pathway category",title="")+
  #labs(x="",y="Genetic Information Processing",fill="Pathway category",title="")+
  #labs(x="",y="Environmental Information Processing",fill="Pathway category",title="")+
  #labs(x="",y="Cellular Processes",fill="Pathway category",title="")+
  #labs(x="",y="Organismal Systems",fill="Pathway category",title="")+
  labs(x="",y="Abundance",fill="Group",title="")
  #theme()+#legend.position ="none"
  #theme(axis.text.x=element_text(angle=45))#vjust=1,hjust=1,size=8#,text=element_text(family="Arial")
print(p)



  ggsave('boxplot.t-test.png',plot=p, height=6, width=8, units="in",dpi=300)
  ggsave('boxplot.t-test.pdf',plot=p, height=6, width=8, units="in",dpi=600)


  
  p2 <- ggplot(mydata, aes(x=type,y=value,fill=group)) + 
    geom_boxplot()+
    stat_boxplot()+
    theme_bw()+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+ #base_family="Arial"
    #scale_fill_nejm()+
    stat_compare_means(#comparisons = my_comparisons,
      label = "p.signif",
      method = "t.test",
      hide.ns=F,
      size = your_font_size)+
    labs(x="",y="Abundance",fill="Group",title="")
  print(p2)
  
  ggsave('boxplot.t-test.ns.png',plot=p2, height=6, width=8, units="in",dpi=300)
  ggsave('boxplot.t-test.ns.pdf',plot=p2, height=6, width=8, units="in",dpi=600)
  
  
  p3 <- ggplot(mydata, aes(x=type,y=value,fill=group)) + 
    geom_boxplot()+
    stat_boxplot()+
    theme_bw()+ theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+ #base_family="Arial"
    #scale_fill_nejm()+
    stat_compare_means(#comparisons = my_comparisons,
      label = "p.format",
      method = "t.test",
      hide.ns=F,
      size = 2)+
    labs(x="",y="Abundance",fill="Group",title="")
  print(p3)
  
  ggsave('boxplot.t-test.pvalue.png',plot=p3, height=6, width=8, units="in",dpi=300)
  ggsave('boxplot.t-test.pvalue.pdf',plot=p3, height=6, width=8, units="in",dpi=600)
  
  
