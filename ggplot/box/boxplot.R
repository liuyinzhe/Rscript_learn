library(ggplot2)
setwd("D:\\Rscript\\example")


# 指定第一列Geneid作为行index
fpkm = read.table('all.gene_FPKM.tsv',header = T,row.names='Geneid',check.names=F)

# 同时转为 log2(x+1)
n_fpkm = gather(log2(fpkm+1))

# FPKM
condition <- factor(c(rep("A",3),rep("B",3)),levels = c("A","B"))
df1 <- data.frame(key = colnames(countData),
                  group = condition
)
new_fpkm = merge(n_fpkm,df1,type=left)

box<- ggplot(new_fpkm, aes(x=key,y=value,fill=group)) +
  geom_boxplot()+
  stat_boxplot(geom = "errorbar",
               lwd=0.5,
               width=0.5)+
  #ggtitle("分组绘制箱式图")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = rel(0.8)))+
  labs(x="Sample",y = "log2(FPKM+1)",title = "")#+ylim(0,1)
#theme_minimal()+
#scale_fill_nejm()
#theme(legend.position ="none")
#box
#ggsave('box.png',box)
#box
ggsave(paste(group_str, 'box.png', sep = "."),plot=box, height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'box.pdf', sep = "."),plot=box, height=8, width=12, units="in",dpi=300)
