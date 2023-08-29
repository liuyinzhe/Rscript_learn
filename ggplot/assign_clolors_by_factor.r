library(ggplot2)
library(RColorBrewer)

#define custom color scale
myColors_map <- brewer.pal(12, "Paired")


colourCount <-  length(unique(iris$Species))
myColors <- colorRampPalette(myColors_map)(colourCount)
names(myColors) <- levels(iris$Species)
custom_colors <- scale_colour_manual(name = "Species Names", values = myColors)

ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) +
  geom_point() +
  custom_colors

# assign_clolors_by_factor

# How to Assign Colors by Factor in ggplot2 (With Examples)
#https://www.statology.org/color-by-factor-ggplot2/
# RColorBrewer 知乎
#https://zhuanlan.zhihu.com/p/546088806
# RColorBrewer | 再多的配色也能轻松搞定！~（二）
#https://zhuanlan.zhihu.com/p/569201131
