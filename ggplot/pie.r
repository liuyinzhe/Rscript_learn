
#library(this.path)
#setwd(this.path::this.dir())

library(ggplot2)
library(scales)

#getwd()
setwd("D:\\xx\\xx\\xx\\xxx")

df_res=read.table("input.tsv",header=T,sep="\t")
df_res$Precent = df_res$Number/sum(df_res$Number)


a<-ggplot(df_res, aes(x='', y=Precent, fill=Type)) +
  geom_bar(stat = 'identity',width = 1,position='stack') + #堆叠柱状图
  scale_y_continuous(expand = c(0,0))+
  coord_polar("y", start=0) + # stackpie ，# 堆叠图变为饼图的关键
  theme_bw()+
  
  theme(legend.position = 'right', 
        legend.text = element_text(colour = 'black',size = 10),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
  )+
  
 geom_text(label = scales::percent(df_res$Precent),
          hjust = "outward",
          position = position_stack(vjust = 0.5))


a

ggsave("pie.png",width = 10, height = 6, dpi = 300)
