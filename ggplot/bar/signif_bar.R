library(ggplot2)
library(ggsignif)
library(ggpubr)


# 代码修改自 https://blog.csdn.net/kanghua_du/article/details/138348688

setwd("Rscript/signif_bar")


data <- data <- read.table(text = "
group	sample	value
lncRNA-1	CK	4
lncRNA-1	CK	4
lncRNA-1	CK	5
lncRNA-1	CK	6
lncRNA-1	CK	3
lncRNA-1	CK	5
lncRNA-1	CK	5
lncRNA-1	Treat	1
lncRNA-1	Treat	2
lncRNA-1	Treat	3
lncRNA-1	Treat	3
lncRNA-1	Treat	3
lncRNA-1	Treat	1
lncRNA-1	Treat	5
mRNA	CK	6
mRNA	CK	7
mRNA	CK	4
mRNA	CK	4
mRNA	CK	5
mRNA	CK	7
mRNA	CK	4
mRNA	Treat	2
mRNA	Treat	6
mRNA	Treat	7
mRNA	Treat	4
", header = TRUE, row.names = NULL)

## 
# data[1:5,1:3] # 行，列


ggplot(data, aes(x = sample, y = value))+
  ## 绘制柱状图
  geom_bar(aes(fill = group), stat = "summary", position = position_dodge(1),
           color = "black",
           fun = mean, size = 0.5)+
  ## 添加误差线
  stat_summary(fun.data = "mean_sd", geom = "errorbar",
               width = 0.2, size = 1)+
  ## 添加显著性
  geom_signif(comparisons = list(c("CK","Treat")),
              map_signif_level= F,  ## T:显示*号，F显示数字
              tip_length=0, 
              size=1, 
              test = "t.test")+  ## t.test, wilcox.test 
   facet_wrap(~group)+
  ## X轴和Y轴坐标
  labs(x = "", y = "Expression levels",title = NULL)+
  ## 设置颜色
  scale_fill_manual(values = c("#4E659B","#8A8CBF", "#B8A8CF","#E7BCC6","#FDCF9E","#EFA484","#B6766C"))+
  theme_classic()+  # x y 轴风格
  theme(axis.line = element_line(size = 1),  ## 线的粗细
        text=element_text(family = "sans",colour ="black",size = 12), # 全局文字 字体 颜色 大小
        axis.text.x = element_text(color = "black", size = 12),  # x轴文字  颜色 大小
        axis.text.y = element_text(color = "black",size = 12), # y轴文字  颜色 大小
        axis.ticks = element_line(size = 1,colour = "black"), # 刻度文字 大小  颜色
        strip.text = element_text(color = "black",size = 16), # facet_wrap 的文字 颜色与大小
        axis.title = element_text(color = "black",size = 18), # 标题 的文字 颜色与大小
        legend.position = "none", # 关闭图注
        strip.background = element_blank() # 关闭facet_wrap 的背景
  )

#ggsave("signif_bar.pdf",width = 6, height = 4)

