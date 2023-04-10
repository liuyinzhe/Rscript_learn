library(ggplot2)

# 读取富集气泡图数据文件
df= read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/bubble/data.txt")# 这里读取了网络上的demo数据，将此处换成你自己电脑里的文件

# 绘图
ggplot(df,aes(x = Ratio, 
              y = reorder(Term,Enrichment,sum), # 按照富集度大小排序
              size = Size,
              colour=Enrichment)) +
  geom_point(shape = 16) +                    # 设置点的形状
  labs(x = "Ratio", y = "Pathway")+           # 设置x，y轴的名称
  scale_colour_continuous(                    # 设置颜色图例
    name="Enrichment",                        # 图例名称
    low="green",                              # 设置颜色范围
    high="red")+
  scale_radius(                               # 设置点大小图例
    range=c(2,4),                             # 设置点大小的范围
    name="Size")+                             # 图例名称
  guides(   
    color = guide_colorbar(order = 1),        # 决定图例的位置顺序
    size = guide_legend(order = 2)
  )+
  theme_bw()                                  # 设置主题


# 代码来自 https://blog.csdn.net/qq_35294674/article/details/122129878
