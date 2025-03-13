# 加载必要的包
library(ggplot2)
library(vegan)
library(ggrepel)  # 用于防止标签重叠

# 创建示例数据
set.seed(123)
data <- data.frame(
  Sample = paste0("Sample", 1:12),
  Group = rep(c("Control", "TreatmentA", "TreatmentB"), each = 4),
  Var1 = c(rnorm(4, 10, 2), rnorm(4, 15, 3), rnorm(4, 8, 1)),
  Var2 = c(rnorm(4, 5, 1), rnorm(4, 7, 2), rnorm(4, 4, 0.5)),
  Var3 = c(rnorm(4, 20, 3), rnorm(4, 25, 4), rnorm(4, 18, 2))
)

# 查看数据结构
head(data)

# 计算Bray-Curtis距离矩阵
dist_matrix <- vegdist(data[, 3:5], method = "bray")

# 执行PCoA分析
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# 提取坐标和解释度
coordinates <- as.data.frame(pcoa_result$points)
colnames(coordinates) <- c("PCoA1", "PCoA2")
variance <- (pcoa_result$eig / sum(pcoa_result$eig)) * 100

# 添加样本信息和分组
coordinates$Sample <- data$Sample
coordinates$Group <- data$Group

# 绘制PCoA图
ggplot(coordinates, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
  # 隐藏散点但保留图例映射
  geom_point(size = 0, alpha = 0) +  # 双重保险隐藏散点（size=0 + alpha=0）
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    box.padding = 0.5, #散点 与 标签的距离
    show.legend = FALSE,
    min.segment.length=30,
    segment.size = 0,   # 点与标签连线的粗细,0隐藏
    #segment.color = NA # 点与标签连线的颜色,NA隐藏
  ) +
  # 强制显示图例并修正图例样式
  guides(
    color = guide_legend(
      override.aes = list(
        size = 2,        # 图例符号大小
        alpha = 1,        # 图例不透明度
        shape = c(16, 17, 15)  # 手动指定形状（需与分组顺序一致）
      )
    ),
    shape = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  
  scale_color_manual( # 自定义散点颜色
    values = c("Control" = "#1f77b4",
               "TreatmentA" = "#ff7f0e",
               "TreatmentB" = "#2ca02c")
  ) +
  scale_shape_manual(  # 自定义散点形状
    values = c("Control" = 16,
               "TreatmentA" = 17,
               "TreatmentB" = 15)
  ) +
  labs(x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
       y = paste0("PCoA2 (", round(variance[2], 1), "%)")) +
  theme_minimal()
