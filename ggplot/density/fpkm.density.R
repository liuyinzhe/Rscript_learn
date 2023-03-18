library(ggplot2)
setwd("D:\\Rscript\\example")

# 指定第一列Geneid作为行index
fpkm = read.table('all.gene_FPKM.tsv',header = T,row.names='Geneid',check.names=F)

# 多列转为适合箱式图的数据
n_fpkm = gather(fpkm)
if(FALSE){


colnames(n_fpkm) <-c("sample","value")
n_fpkm$class[0 <= n_fpkm$value & n_fpkm$value<= 0.5] = 'FPKM 0-0.5'
n_fpkm$class[0.5 < n_fpkm$value & n_fpkm$value<= 1] = 'FPKM 0.5-1'
n_fpkm$class[1 < n_fpkm$value & n_fpkm$value<= 10] = 'FPKM 1-10'
n_fpkm$class[10 < n_fpkm$value  ] = 'FPKM >=10'

n_fpkm$ord[0 <= n_fpkm$value & n_fpkm$value<= 0.5] = 0
n_fpkm$ord[0.5 < n_fpkm$value & n_fpkm$value<= 1] = 1
n_fpkm$ord[1 < n_fpkm$value & n_fpkm$value<= 10] = 2
n_fpkm$ord[10 < n_fpkm$value  ] = 3
write.table(n_fpkm, file =paste(group_str, "fpkm_plot.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
}
##################### fpkm  密度分布图 #########################
### 密度分布图 版本1 
library(reshape2)

density_input = n_fpkm
colnames(density_input) <-c("Sample","value")
#density_input['value'] = log10(density_input['value']+1)
#density_input['value']<-log2(density_input['value']+1)
density_plot <- ggplot(density_input,aes(log10(value), color=Sample)) + # fill=Sample,
  xlab("log10(FPKM)") + #"Expression level"
  geom_density(alpha = 0.6) + # coord_cartesian(xlim = c(-10, NA), ylim = c(0, NA))
   theme_bw() # geom_rug() +
ggsave(paste(group_str, 'fpkm.density.png', sep = "."),plot=density_plot, height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'fpkm.density.pdf', sep = "."),plot=density_plot, height=8, width=12, units="in",dpi=300)

# 使用geom_line函数绘制密度分布曲线
density_plot <-  ggplot(density_input,aes(log10(value),after_stat(density), color=Sample))  +
  geom_line(stat="density") +
  theme_bw() + facet_wrap(.~Sample) +
  theme(axis.title = element_text(size=16),
        axis.text=element_text(size=16))

ggsave(paste(group_str, 'fpkm.density.facet.png', sep = "."),plot=density_plot, height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'fpkm.density.facet.pdf', sep = "."),plot=density_plot, height=8, width=12, units="in",dpi=300)

#https://www.modb.pro/db/171989

#https://m.yunbios.net/cn/h-nd-1082.html
##################### FPKM  密度分布图 #########################
