# Rscript_learn


坐标轴位置
  scale_x_continuous(position = "top")#+
  #scale_y_continuous(position = "right")

请问ggplot2作图：如何把图形的title放置到图形的底下？
https://www.zhihu.com/question/53119447
  labs(caption = "Figure 1. GDP Growth Rate") +
  theme(plot.caption=element_text(colour = "blue", hjust=0.5))
