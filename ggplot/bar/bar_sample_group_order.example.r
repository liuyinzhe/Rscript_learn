
rm(list = ls())
options(stringsAsFactors = F)

setwd("D:\\Rscript\\meta")


# 加载必要的包
library(ggplot2)
library(dplyr)

# 1. 读取group.txt文件
# 假设文件以制表符或逗号分隔，列名包含"#SampleID"和"Group"
group_data <- read.delim("group.txt", header = TRUE, stringsAsFactors = FALSE) # 如果是制表符分隔
# 如果group.txt是CSV格式，使用：group_data <- read.csv("group.txt", header = TRUE, stringsAsFactors = FALSE)

# 处理列名：将"#SampleID"更改为合法的R列名
names(group_data) <- c("SampleID", "Group")

# 2. 设置因子水平顺序
# SampleID的顺序按照其在group.txt中出现的顺序
group_data$SampleID <- factor(group_data$SampleID, levels = unique(group_data$SampleID))
# Group的顺序也按照其在group.txt中首次出现的顺序
group_data$Group <- factor(group_data$Group, levels = unique(group_data$Group))

# 3. 创建示例绘图数据 (假设每个SampleID有一个对应的数值)
# 这里随机生成一些数值，您应替换为实际数据
set.seed(123) # 设置随机种子保证示例可重现
plot_data <- data.frame(
  SampleID = group_data$SampleID,
  Value = runif(nrow(group_data), min = 10, max = 100)
)
#plot_data
write.table(plot_data, file ="plot_data.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# 4. 将分组信息合并到绘图数据中，并确保顺序一致
plot_data <- plot_data %>%
  left_join(group_data, by = "SampleID") %>%
  mutate(
    SampleID = factor(SampleID, levels = levels(group_data$SampleID)), # 保持SampleID顺序
    Group = factor(Group, levels = levels(group_data$Group)) # 保持Group顺序
  )

# 5. 绘制柱状图
ggplot(plot_data, aes(x = SampleID, y = Value, fill = Group)) +
  geom_col(colour = "black", size = 0.2, width = 0.7) + # 绘制柱状图，设置边框、宽度等
  labs(title = "分组柱状图 (样品顺序和分组顺序已指定)", 
       x = "样品ID", 
       y = "测量值") +
  theme_bw() + # 使用黑白主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 将X轴标签旋转45度防止重叠
