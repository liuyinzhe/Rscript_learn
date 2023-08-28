library(ggplot2)
library(ggpubr)
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\draw\\boxplot_group')

md<- read.table("plotdata.txt",header = T)
my_comparisons <- list(c("A", "B"))
p<- ggplot(md,aes(group,gnum))+
    geom_violin(aes(fill=group))+
    geom_boxplot(width = 0.1)+
    theme_bw()+
    stat_compare_means(#comparisons = my_comparisons,
                    #label = "p.format",#"p.signif",
                    method = "wilcox.test", #"t.test", #"wilcox.test", #"anova", 
                    #hide.ns=T,
                    #size = 3
                    )+
    theme(panel.grid = element_blank())+
    ylab("Number of element")

p
png("boxplot_group_ttest_pvalue.png",width=800,height=600,res=100,units="px")
print(p)
dev.off()
pdf("boxplot_group_ttest_pvalue.pdf",width=8,height=6)
print(p)
dev.off()

# https://zhuanlan.zhihu.com/p/652379786
