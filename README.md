# Rscript_learn


坐标轴位置  
  scale_x_continuous(position = "top")#+  
  #scale_y_continuous(position = "right")  

请问ggplot2作图：如何把图形的title放置到图形的底下？  
https://www.zhihu.com/question/53119447   
  labs(caption = "Figure 1. GDP Growth Rate") +   
  theme(plot.caption=element_text(colour = "blue", hjust=0.5))  



#### 检查当前库路径
.libPaths()

#### 指定新库路径
new_path <- "/data/my_library"

#### 创建新库路径
dir.create(new_path, showWarnings = FALSE)  # 创建新目录（如果尚未存在）

#### 加如新库路径
.libPaths(c(new_path, .libPaths()))  # 将新路径添加到库路径

#### 安装到指定路径
install.packages("ggplot2", lib = new_path)

#### 从指定路径加载包
library(ggplot2, lib.loc = new_path)

#### 降级安装
```
#网站查看 https://cran.r-project.org/src/contrib/Archive 包的旧版本
#指定安装旧版本
install.packages("plotly",version="2.2",repos="http://cran.r-project.org")
```

