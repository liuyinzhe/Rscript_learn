
remove(list = ls()) #清除 Global Environment
# rm(list = ls())
options(stringsAsFactors = F)
# par(ask=T) # 若为TRUE（且当前的R会话是可交互状态），则在绘制新图像之前会要求用户输入确认信息。
par(ask=F)
setwd('D:\\Rscript\\volcano_gene\\draw\\sample')

# 加载包
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)   # 用于作图的包
  library(ggrepel)   # geom_text_repel
})

data_pre <- read.csv('sample.genename.xls',sep='\t',comment.char = "",header = TRUE,stringsAsFactors = FALSE)
# 删除 pval列 为Na的行
data <- data_pre %>% drop_na(pval)

#colnames(data) <- c("gene_name", "FC", "log2FC", "PValue", "FDR", "log10FDR","log10PValue") #重新命名行名
# GeneName，FoldChange，Log2FoldChange，pval，padj，log10(padj),log10(pval)
data$log10pval<- -log10(data$pval)


###数据处理——根据FDR进行差异基因筛选###
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =1 #设置log2FC的阈值
# data$Sig = ifelse(data$pval < cut_off_FDR &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
#                         abs(data$Log2FoldChange) >= cut_off_log2FC,  #abs绝对值
#                       ifelse(data$Log2FoldChange > cut_off_log2FC ,'Up','Down'),'no')

# pval<0.05 并且log2 (fold change) 绝对值大于=1,红色；log2FoldChange 大于1 取蓝色
data$Sig <- ifelse(data$pval<cut_off_FDR & abs(data$Log2FoldChange)>= cut_off_log2FC,ifelse(data$Log2FoldChange > cut_off_log2FC,'Up','Down'),'Stable')


# 指定为有顺序的因子
data$Sig <- factor(data$Sig,levels = c("Up","Stable","Down"))

data = data.frame(data)
table(data$Sig) #查看数据统计情况

###绘图——基础火山图###
p1 <- ggplot(data, aes(x =Log2FoldChange, y=log10pval, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  #scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values = c("Up" = "#ff4757", "Down" = "#546de5", "Stable" = "#d2dae2")) + 
  xlim(c(-5, 5)) +
  # + coord_cartesian(xlim =c(-10, 10))+

  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="grey",lwd=0.6) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="grey",lwd=0.6) +  #添加y轴辅助线
  labs(x="log2(Fold Change)", y="-log10(P-value)") +  #x、y轴标签
  ggtitle("sample") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        #plot.title = element_text(size = 16) # 图表标题文字大小
        legend.position="right", 
        legend.title = element_blank(),
        #legend.title = element_text(size = 12), # 图例标题文字大小
        # 文字大小
        axis.title = element_text(size = 12), # 轴标题文字大小
        axis.text = element_text(size = 11), # 轴刻度文字大小
        legend.text = element_text(size = 12), # 图例项文字大小
        
  ) 


p1 #出图


# 使用 filter 函数筛选出包含 target_values 中的值的行
target_gene_df <- data %>% filter(GeneName %in% c("IRF5","PSMB9","TBXAS1","TNF"))


###添加基因名标记###
p2 <- p1 + geom_label_repel(
    #data = subset(data, data$pval < cut_off_FDR & abs(data$Log2FoldChange) >= cut_off_log2FC),# 可以设置跟上面不同的阈值，用数值替换即可
    data = target_gene_df,
    aes(label = GeneName,color = Sig), size = 2,
    seed = 123,
    #max.overlaps = Inf, # 控制重叠情况, inf 不限制重叠
    min.segment.length = 0, # 间隔 不限制
    label.padding = 0.25, # 调大则标签空白增大
    label.r = 0.15, # label 中心到边框半径，调高了就是椭圆
    label.size = 0.25, # 边框粗细
    segment.size = 0.6,# 连线粗细
    box.padding = unit(0.8, "lines"), # 连线距离
    point.padding = unit(0.8, "lines"),  # 连线与文字之间 留白
    #color="black",
    segment.color = "black", 
    #fill="white", # 标签填充
    show.legend = FALSE ) # 都是stable 则设置 color="black",
p2 #出图


p3 <- p2 + geom_point(data=target_gene_df,aes(x =Log2FoldChange, y=log10pval),
                      colour="black",size=2, shape = 1, stroke = 0.5,alpha=0.6) # stroke 圈的粗细


p3 #出图

ggsave('sample.volcano_gene.png',plot=p3,height=8, width=12, units="in",dpi=200)
ggsave('sample.volcano_gene.pdf',plot=p3,height=8, width=12, units="in",dpi=100)
