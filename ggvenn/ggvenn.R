#!/usr/bin/Rscript
rm(list = ls())
#options(stringsAsFactors = F)
setwd('filter')

library(ggvenn)

data <- read.csv("venn_plot.tsv",sep='\t')


#data[data == "TRUE"] <- TRUE
#data[data == "FALSE"] <- FALSE

####################################  demo ########################################
# # 生成样例数据
# d <- tibble(value   = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13),
#             `Set 1` = c(T, F, T, T, F, T, F, T, F,  F,  F),
#             `Set 2` = c(T, F, F, T, F, F, F, T, F,  F,  T),
#             `Set 3` = c(T, T, F, F, F, F, T, T, F,  F,  F),
#             `Set 4` = c(F, F, F, F, T, T, F, F, T,  T,  F))

mypal= c("#69C2A5","#8DA0CC","#EE91C4","#FC8E65")
# 可视化绘制
ggvenn_4 <- ggplot(data, aes(A = B042D_B04ES2D, B = B043D_B042D, C = B043D_B04ES3D, D = B04ES3D_B04ES2D)) +
  geom_venn(
            data = data,
            # 颜色填充
            fill_color=mypal,
            fill_alpha = .5,  #透明度
            # 边线设置参数
            #stroke_color = mypal, 
            #stroke_alpha = 0.5,
            stroke_size = 0.8, #线宽
            stroke_linetype = "longdash", # "solid",
            # 集合名字设置
            set_name_color = mypal, # 文字颜色
            set_name_size = 4,
            # 图形中字体
            #text_color = 'black',
            text_size= 4,
            # 百分比/小数点
            show_percentage = FALSE,#TRUE, # 显示百分比
            digits = 2 # 百分比下属点后位数
            ) + 
  theme_void()+
  coord_fixed() +
  labs(title = "", #Example of ggvenn:: geom_venn function",
       subtitle = "", #processed charts with geom_venn()",
       caption = "")+  #Visualization by DataCharm") +
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                  size = 20, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
       plot.caption = element_text(face = 'bold',size = 12))

ggvenn_4

ggsave("venn.png",plot=ggvenn_4,height=10, width=15, units="cm")#,dpi=300) 

ggsave("venn.pdf", device='pdf',plot=ggvenn_4,height=10, width=15,units="in")#,dpi=300)




if(FALSE){
  # 代码参考自:
  #R语言绘图包05--韦恩图的绘制：ggvenn和VennDiagram
  #https://www.jianshu.com/p/6a1d5457f358
  #ggvenn绘制ven图
  #https://www.jianshu.com/p/06d180d9b3c8
}



