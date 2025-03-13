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
ggplot(coordinates, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Control" = "#1f77b4",
                               "TreatmentA" = "#ff7f0e",
                               "TreatmentB" = "#2ca02c")) +
  labs(x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
       y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
       title = "PCoA Plot with Sample Labels") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
