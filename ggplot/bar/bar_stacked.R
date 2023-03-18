library(ggplot2)
setwd("D:\\Rscript\\example")

data = read.table("bar_stacked_input.tsv",header=T,sep="\t",check.names=F)

#对样本和物种的因子排序，设置为表格导入时的顺序

#sorder = factor(data$sample,levels=unique(data$sample),order=TRUE)
#porder = factor(data$class,levels=unique(data$class),order=TRUE)

#以样本为横坐标，丰度值为纵坐标，分类fill 颜色，绘制堆叠图，包括数字显示
#P <- ggplot(data=data,aes(x=sorder,y=value,fill=porder))  + 
P <- ggplot(data=data,aes(x=sample,y=value,fill=class))  +
  geom_bar(stat="identity",position="stack") + scale_color_manual(values = c("#00468B","#ED0000","#42B540","#0099B4"),) +
  labs(x="Sample Names",y="Gene number",fill="class",title="") +
 # geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  #geom_text(aes(x = sample, label = value),position = "fill") +
  geom_text(aes(x = sample, y = value, label = value, group = class), position = position_stack(vjust = .05)) +
  #http://cn.voidcc.com/question/p-mywqjccf-h.html
  stat_summary(fun = sum, aes(label = after_stat(y), group = sample), geom = "text",position = position_stack(vjust = 1.015)) +
  # https://qa.1r1g.com/sf/ask/2718710711/
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,face="bold",size=15)) #设置横坐标的角度

ggsave('bar_stacked.png',,plot=P, height=12, width=12, units="in",dpi=600)
ggsave('bar_stacked.pdf',plot=P, height=12, width=12, units="in",dpi=600)
