rm(list = ls())
options(stringsAsFactors = F)
setwd('./')


library(igraph)
library(WGCNA)
#library(BiocManager)
#BiocManager::install('WGCNA')

CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat) # 由于相关性矩阵是对称的，取上三角矩阵计算即可
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}


# 读取otu-sample矩阵，行为sample，列为otu
otu = read.csv("otu_table.tsv",sep='\t',header=TRUE, check.names=FALSE,row.names=1)#,head=T,row.names=1)


occor <- corAndPvalue(otu, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系



cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边


igraph <- graph_from_data_frame(cor_df, directed = FALSE ) #,direct=F
length(V(igraph)) # 查看节点数
length(E(igraph)) # 查看边数

V(igraph)$size <- degree(igraph)*0.8 # 节点的大小与节点的度成正比，进行0.8倍放缩 #2.8


# 定制颜色
#color37 
cols <- c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

#cols <- c('#00A6FB', '#0582CA', '#fff0f3', '#006494', '#c9e4ca', '#31572c', '#90a955', '#ecf39e', '#a4133c', '#c9184a', '#ff4d6d')
# V(igraph)$color <- sample(cols, length(V(igraph)), replace = T) # 从随机颜色中采样给节点上色

# 第一列作为行索引
otu_pro = read.csv("otu_pro.tsv",sep='\t',header=TRUE, check.names=FALSE,row.names=1)


# 从otu名字提取对应 dataframe
otu_name_df = otu_pro[V(igraph)$name,] 
# 指定的一个 taxonomy 类别如 phylum,对应颜色
V(igraph)$color<-cols[factor(otu_name_df$phylum)]

# 添加其它 taxonomy
V(igraph)$kingdom <- otu_name_df$kingdom
V(igraph)$phylum <- otu_name_df$phylum
V(igraph)$class <- otu_name_df$class
V(igraph)$order <- otu_name_df$order
V(igraph)$family <- otu_name_df$family
V(igraph)$genus <- otu_name_df$genus
V(igraph)$species <- otu_name_df$species


# 获得lenged lable/color
legend_lable <- unique(factor(otu_name_df$phylum))
legend_color <- unique(cols[factor(otu_name_df$phylum)])

# 属性调节 
E(igraph)$color[E(igraph)$cor >= 0.6] <- "red"	 # 正相关则边为红色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "blue"	 # 负相关则边为蓝色

#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
#E(igraph)$correlation <- E(igraph)$cor#将weight赋值为correlation，正负代表相关性.
#E(igraph)$weight <- abs(E(igraph)$weight)#将weight取绝对值
E(igraph)$width <-  abs(E(igraph)$cor)*0.25 # 边的粗细与相关系数成正比，进行0.5倍放缩
#cor 就是 weight 绝对值，0.5倍缩放

# 固定随机数，保证出图一致性
set.seed(567)
coords <- layout_with_fr(igraph, niter=9999,grid="nogrid") # 生成布局
#coords <- layout_randomly(igraph)#, niter=9999)#,grid="nogrid") # 生成布局
#coords <- layout_with_kk(igraph)
#coords <- layout_components(igraph)
#coords <- layout_nicely(igraph)
#coords <- layout_with_graphopt(igraph)
#coords <- layout_as_star(igraph)
#coords <- layout_with_mds(igraph) #
#coords <- layout_on_sphere(igraph)


# pdf("Figure1.pdf", height = 100,width = 100) # 保存图为PDF，指定宽和高
# plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) # 画图
# dev.off()


pdf("Figure.pdf", height = 10,width = 14) # 保存图为PDF，指定宽和高
plot(
  igraph, 
  layout=coords, 
  #main="Co-occurrence network",
  cex.main = 6,
  vertex.label = NA, 
  vertex.frame.color=NA,
  margin=c(0,0,0,1) # 上下左右边框空余
) # 画图
legend(x=1.2,y=0.35,#"right", 
       title = "phylum",
       title.font = 2, # 黑体加粗
       title.adj = 0.2,
       title.cex = 2, # 文字大小
       legend = legend_lable, 
       pch=21,  # 散点类型
       col=legend_color, 
       pt.bg=legend_color, 
       pt.cex=2, # 图注中图形的大小
       cex=1.5,  # 注释整体大小调节
       bty="n", # 关闭背景框
       ncol=1   # 注释拆分多列
)

legend( x=1.2,y=0.75,
        title = "correlation",
        title.adj = 0.8,
        title.cex = 2, # 文字大小
        legend = c("positive","negative"), 
        col = c("red","blue"), 
        title.font = 2, # 黑体加粗
        lty = 1,  # 线条类型
        lwd = 1,  # 线条宽度
        bty="n",  # 关闭背景框
        cex=1.5,  # 注释整体大小调节
        pt.cex=2, # 图注中图形的大小
)
title(main = "Co-occurrence network", cex.main = 4)

dev.off()

png("Figure.png",height = 1000, width = 1400, units = "px",res=100) # 保存图为PDF，指定宽和高
plot(
  igraph, 
  layout=coords, 
  #main="Co-occurrence network",
  cex.main = 6,
  vertex.label = NA, 
  vertex.frame.color=NA,
  margin=c(0,0,0,1) # 上下左右边框空余
) # 画图
legend(x=1.2,y=0.35,#"right", 
       title = "phylum",
       title.font = 2, # 黑体加粗
       title.adj = 0.2,
       title.cex = 2, # 文字大小
       legend = legend_lable, 
       pch=21,  # 散点类型
       col=legend_color, 
       pt.bg=legend_color, 
       pt.cex=2, # 图注中图形的大小
       cex=1.5,  # 注释整体大小调节
       bty="n", # 关闭背景框
       ncol=1   # 注释拆分多列
)

legend( x=1.2,y=0.75,
        title = "correlation",
        title.adj = 0.8,
        title.cex = 2, # 文字大小
        legend = c("positive","negative"), 
        col = c("red","blue"), 
        title.font = 2, # 黑体加粗
        lty = 1,  # 线条类型
        lwd = 1,  # 线条宽度
        bty="n",  # 关闭背景框
        cex=1.5,  # 注释整体大小调节
        pt.cex=2, # 图注中图形的大小
)
title(main = "Co-occurrence network", cex.main = 3)

dev.off()


####保存为其他软件格式，
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(igraph, 'network.graphml', format = 'graphml')
#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write_graph(igraph, 'network.gml', format = 'gml')



#边的输出样式
edge <- data.frame(as_edgelist(igraph, names = TRUE))    #igraph 的邻接列表转为边列表

df <- as.data.frame(E(igraph)$cor)
df[df>0] <- 1
df[df<0] <- -1
colnames(df) <- c('cor')
edge <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$width, # 权重只是0.5倍的correlation绝对值
  correlation = E(igraph)$cor,
  cor = df)

#write.table(edge, 'network.edge.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(edge, 'network.edge.csv',row.names = F)


#输出点的表格
#节点属性列表
node <- data.frame(
  id = names(V(igraph)),
  kingdom = V(igraph)$kingdom,
  phylum = V(igraph)$phylum,
  class = V(igraph)$class,
  order = V(igraph)$order,
  family = V(igraph)$family,
  genus = V(igraph)$genus,
  species = V(igraph)$species)

#write.table(node, 'network.node.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(node, 'network.node.csv',row.names = F)
