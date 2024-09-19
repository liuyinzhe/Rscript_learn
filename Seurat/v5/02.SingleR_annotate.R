if(F){

library(Seurat)
library(celldex)
library(tidydr)
library(ggalluvial)
library(ggplot2)
library(SingleR)

# 设置默认工作路径
setwd('E:/Rscript/seurat_v5')

# V4 分析内容,读取V5 读取结果进行
load('scRNA_harmony.Rdata') # scRNA_harmony
#FindAllMarkers 函数将返回一个包含标记基因信息的数据框。数据框的每一行代表一个标记基因，包含以下列：

#gene：基因名称。
#p_val：统计学显著性检验的 p 值。
#avg_log2FC：平均 log2 倍变化。
#pct.1 和 pct.2：标记基因在两个组中的表达百分比。
#cluster：标记基因所属的聚类（如果进行了聚类分析）。

mmRNA_harmony.markers <- FindAllMarkers(scRNA_harmony,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

# only.pos = TRUE：表示只寻找上调的标记基因（即表达增加的基因）。如果将其设置为 FALSE，则会同时寻找上调和下调的标记基因
# min.pct = 0.25：指定标记基因在至少多少百分比的细胞中表达。这里设置为 0.25，表示标记基因必须在至少 25%的细胞中表达
# logfc.threshold = 0.25：指定 log2 倍变化阈值。只有当基因的表达在不同组之间的差异达到或超过这个阈值时，才会被认为是标记基因

# 这段代码使用了 dplyr 包中的函数对 pbmc.markers 数据框进行分组，并在每个组中选择前 n 个具有最大 avg_log2FC 值的行
mmRNA_harmony.markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)


###展示前5个差异marker的热图
top5 <- mmRNA_harmony.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 展示每个类群前10个基因
####这里的长宽需要按照自己的cluster数量来调整，cluster越多，width的数值就越大
pdf(file = "heatmap.pdf", width = 45, height = 15)
pdf("heatmap.pdf")
DoHeatmap(scRNA_harmony, features = top5$gene) + NoLegend()
dev.off()


###导出marker基因的列表
write.csv(mmRNA_harmony.markers, file = "cluster_markers.csv")


##############################  singleR自动注释  ##############################
##########singleR自动注释
# 获取单细胞数据
SingleRdata <- HumanPrimaryCellAtlasData()
###运行完可以保存
saveRDS(SingleRdata, "singleRdata.rds")
# 获取细胞类型的唯一标签
unique_cell_types <- unique(SingleRdata@colData@listData[["label.main"]])

# 提取测试数据集的标准化表达矩阵
#Stand_data <- scRNA_harmony@assays$RNA@data # V4
Stand_data <- scRNA_harmony@assays$RNA$data  # V4 @layers$data
dim(Stand_data)

# 获取聚类结果作为参考
clusters <- scRNA_harmony@meta.data$seurat_clusters

# 使用SingleR进行细胞类型预测
cellpred <- SingleR(test = Stand_data, 
                    ref = SingleRdata, 
                    labels = SingleRdata$label.main,
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

# 展示预测结果的结构
str(cellpred, max.level = 3)

# 提取元数据信息
metadata <- cellpred@metadata

# 将预测结果与细胞类型合并
celltype <- data.frame(ClusterID = rownames(cellpred), 
                       celltype = cellpred$labels, 
                       stringsAsFactors = FALSE)

# 导出细胞类型信息
write.csv(celltype, "Auto_anno_SingleR.csv", row.names = FALSE)

# 绘制打分热图
pdf("singleR.pdf", width = 7.5, height = 7.5)
p <- plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")
p
dev.off()

# 应用SingleR注释后的结果进行可视化，并将注释结果整合到meta.data中
newLabels <- cellpred$labels
names(newLabels) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, newLabels)

# 进行UMAP和tSNE可视化
scRNA_harmony  <- Seurat::RunUMAP(scRNA_harmony,dims = 1:20)
###可视化
pdf("singleR_UMAP.pdf", width = 7, height = 5.5)
DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 1) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right")
dev.off()

pdf("singleR_TSEN.pdf", width = 7, height = 5.5)
DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 1) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right")
dev.off()

# 为mmRNA_harmony添加细胞类型列
scRNA_harmony$Celltype <- Idents(scRNA_harmony)

####保存
saveRDS(scRNA_harmony,"After_auto_anno_mmRNA_harmony.rds")

###
#绘制堆叠图##########
cluster_colors <- scales::hue_pal()(length(levels(scRNA_harmony$seurat_clusters)))

# 将颜色与簇名称对应起来
cluster_color_map <- setNames(cluster_colors, levels(scRNA_harmony$seurat_clusters))

# 打印UMAP图
umap_integrated_split <- DimPlot(scRNA_harmony, 
                                 reduction = "umap", 
                                 label = TRUE, 
                                 split.by = "orig.ident",  # 按样本拆分
                                 group.by = "seurat_clusters",  # 在每个样本内按cluster分组
                                 combine = FALSE,  # 防止自动组合和显示标题
                                 ncol = 3)
umap_integrated_split
combined_plot <- wrap_plots(umap_integrated_split, ncol = 3)  # 假设您希望以3列的方式排列
ggsave("umap_integrated_split.pdf", plot = combined_plot, width = 20, height = 5)

# 输出颜色映射，以便在第二段代码中使用
cat("Cluster color mapping:\n")
print(cluster_color_map)
# 计算每个样本中每个细胞类型的占比
head(scRNA_harmony@meta.data)
df <- scRNA_harmony@meta.data %>%
  group_by(orig.ident, seurat_clusters)

# 重命名列名
colnames(df) <- c("Sample", "Cluster", "Number", "ratio")

# 创建堆叠柱状图
# 这里我们假设 cluster_color_map 是在第一段代码中提取的颜色映射


###查看某个基因的表达情况
VlnPlot(scRNA_harmony,features = c("KIRREL3"))
VlnPlot(scRNA_harmony,features = c("EPHA2"))
#VlnPlot(scRNA_harmony,features = c("KRT5"))
###umap图展示某个基因的表达情况
FeaturePlot(scRNA_harmony, 
            features = c("KRT5"))

}
