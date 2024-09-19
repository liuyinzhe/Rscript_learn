# 生信超详细seurat_v5标准流程全代码
# 代码来自 https://blog.csdn.net/ZxVSaccount/article/details/141062694
if(T){
  # 主要用于做单细胞分析的各种函数
  library(Seurat)
  # 主要与批次效应有关
  library(harmony)
  # 广泛使用的数据分析和可视化的工具集
  library(tidyverse)
  # 是一个功能强大的R绘图包
  library(ggplot2)
  # 增强绘图和图形组合的包
  library(cowplot)
  # 用于数据的操作和处理,管道,过滤,排序等
  library(dplyr)
  # 用于组合多个图形
  library(patchwork)
  # 用于文件的处理
  library(R.utils)
  # 自动注释
  library(celldex)
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  ############
  library('BiocManager',quietly=TRUE)
  #BiocManager::install('celldex')
  
  #BiocManager::install('SingleR')
  library(SingleR) 
}

# 设置默认工作路径
setwd('E:/Rscript/seurat_v5')

# 
#################################################### 1.批量读取单细胞数据 #######################################################
# 使用list.files函数列出文件夹中的所有文件
dir_name <- list.files('data/')
dir_name

# 创建一个空的list对象来存储转化之后的seurat对象
scRNAlist <- list()
for (i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste('E:/Rscript/seurat_v5/data/',dir_name[i],sep = ''))
  # 使用createseuratobject来创建seurat对象,使用counts矩阵,设置样本名为目录名
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,
                                       project = strsplit(dir_name[i],'_')[[1]][2],
                                       min.cells = 3,
                                       min.features = 200)  
}
scRNAlist # 查看样本信息
#################################################### 2.数据质控 #######################################################
############## 2.1.批量计算每个样本的线粒体和红细胞的比例(质控 --- 线粒体,红细胞,核糖体) #################
# 注意如果是小鼠为 '^MT-'  人为 '^mt-'
# 如果是小鼠的话,将49-53行改为sc[['HB_percent']] <- PercentageFeatureSet(sc,pattern = "^Hb[^(p)]")
for (i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]] # 获取scRNAlist中的第i个seurat对象
  # 计算线粒体的比例
  sc[['mt_percent']] <- PercentageFeatureSet(sc,pattern = '^MT-')
  # 计算红细胞的比例
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 定义红细胞的基因列表
  HB_m <- match(HB_genes,rownames(sc@assays$RNA)) # 在seurat对象的RNA数据中有红细胞基因的索引位置
  HB_genes <- rownames(sc@assays$RNA)[HB_m] # 获取匹配到的红细胞基因的行名
  HB_genes <- HB_genes[!is.na(HB_genes)] # 删除na值(未匹配到的基因)
  sc[['HB_percent']] <- PercentageFeatureSet(sc,features = HB_genes) # 计算红细胞基因的比例,并将其储存在新列HB_percent中
  # 将sc的值赋给scRNAlist
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}

#################################################### 2.2 批量绘制质控前的小提琴图 ####################################################
violin_before <- list()
for (i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c('nFeature_RNA','nCount_RNA','mt_percent','HB_percent'),
                                pt.size = 0.01,
                                ncol = 4)
}
violin_before # 输出质控前的图片,有多少个样本就会输出多少张图片
violin_before[[2]]

#################################################### 2.3 批量过滤细胞,MT,HB基因 ####################################################
scRNAlist <- lapply(scRNAlist, FUN = function(x){
  x <- subset(x, # 怀疑有双细胞的nFeature_RNA可以设置低一点,依据这个样本的小提琴图最低设置4000
              subset=nFeature_RNA > 200 & nFeature_RNA < 3500 &
                mt_percent < 7 &
                HB_percent < 3 & # 红细胞基本没有,随便过滤1-5%都可以
                nCount_RNA < quantile(nCount_RNA,0.97) & # 过滤最高的 3% 的nCount_RNA
                nCount_RNA > 1000)
})
# 一般默认线粒体的含量要少于20% ,红细胞的数据至少要少于5%
view(scRNAlist[[1]]@meta.data)

#################################################### 3. merge合并样本 ####################################################
# 手动设定一个样本名称
dir_name <- c('26','27','28')
scRNAlist_merge <- merge(x = scRNAlist[[1]],y = scRNAlist[-1],add.cell.ids = dir_name)
# 质控后的数据(小提琴图)
violin_after <- VlnPlot(scRNAlist_merge,
                          features = c('nFeature_RNA','nCount_RNA','mt_percent','HB_percent'),
                          split.by = 'orig.ident',
                          layer = 'counts',
                          pt.size = 0.01,
                          ncol = 4)
# 统计细胞数量
table(scRNAlist_merge$orig.ident)

# tips:seuratv5中比seuratv4中多了一个layers,并且layers层下的每个样本都是分开的,这样做的主要目的就是为了优化内存

#################################################### 3.1 meatadata增加分组信息 ####################################################
# 方法很多,也可以不增加分组,因为orig.ident就已经是一个分组了
# 这里简单介绍依照分组,但是这个案例中的数据不会使用,只有normal组
# split 使用"-" 拆分 scRNAlist_merge@assays$RNA 行名(行索引)，取第一列
group <- str_split(colnames(scRNAlist_merge@assays$RNA),'-',simplify = T)[,1] # 这里的simplify表示将结果简化为一个矩阵形式,并且取出第一列
group <- ifelse(str_detect(group,"N$"),'normal','tumor')
table(group)
# 增加一列group
scRNAlist_merge$group <- group

#################################################### 4.标准化&高变基因&归一化&PCA ####################################################
scRNAlist_merge <- NormalizeData(scRNAlist_merge)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge,nFeature_RNA = 2000) # 这里指定最后想要得到的基因数量
scRNAlist_merge <- ScaleData(scRNAlist_merge,vars.to.regress = c('mt_percent')) # 这里去除线粒体的影响,当然这里如果有细胞周期的影响,我们要把细胞周期也回归一下,去除影响
scRNAlist_merge <- RunPCA(scRNAlist_merge)

# 可以SCT来替换此步骤(初学者不建议使用) --- 109-111
# scRNAlist_merge <- SCTransform(scRNAlist_merge,vars.to.regress = 'mt_percent')

#################################################### 5.选择harmony,CCA,scVI,RPCA,FastMNN等整合去批次 ####################################################
# 以下为官网的五种去批次的方法,运算速率相比seuratv4快了很多,特别是CCA
# Anchor-based CCA integration(method = CCAIntegration)
# Anchor-based RPCA integration(method = RPCAIntegration)
# harmony (method = HarmonyIntegration)
# FastMNN (method = FastMNNIntegration)
# scVI (method = scVIInteration) scVI得创建python环境,新手不建议使用

# Integrated with CCA
scRNA_cca <- IntegrateLayers(object = scRNAlist_merge,
                             method = CCAIntegration, # 需要什么方法改变这里即可
                             orig.reduction = 'pca', # 这里降维必须选择pca
                             new.reduction = 'integrated.cca') # 储存在新的降维结果integrated.cca中

# Integrated with harmony
scRNA_harmony <- IntegrateLayers(object = scRNAlist_merge,
                             method = HarmonyIntegration, # 需要什么方法改变这里即可
                             orig.reduction = 'pca', # 这里降维必须选择pca
                             new.reduction = 'harmony') # 储存在新的降维结果integrated.cca中

#################################################### 6.joinLayers合并样本counts和data ####################################################
scRNA_harmony[['RNA']] <- JoinLayers(scRNA_harmony[['RNA']])

#################################################### 7.降维聚类umap ####################################################
# 注意: redution得选择harmony,如果是CCA就选择cca
ElbowPlot(scRNA_harmony,ndims = 50)
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = 'harmony',dims = 1:20)
# 这里的分辨率resolution选的是0.1 - 1 ,之后根据样本直接选择一个就行;这里seq 步长0.1设置了0.1到1 的
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
# 后面 看聚类树 没判断 0.3之后分叉变多
#scRNA_harmony <- FindClusters(scRNA_harmony,resolution = 0.3)

# 降维
scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:20,reduction = 'harmony')
scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:20,reduction = 'harmony')

# 保存降维聚类之后的数据
save(scRNA_harmony,file = 'scRNA_harmony.Rdata')
load('scRNA_harmony.Rdata')

# 查看聚类树，判断 分辨率，resolution ——>FindClusters(scRNA_harmony,resolution)
library(clustree)
clustree(scRNA_harmony) # 如何选择? 首先选择中间的,之后选择没有连线的,连线越多,被打乱就越多,所以我们选择0.3

# 用umap图查看harmony的整合情况
# 我们这三个样本可以看出融合的还是比较好的
# 如果融合的情况不好的话,比如说肿瘤细胞有很强的抑制性,很难与其他细胞融合,这样我们也可以将它单独拿出来,研究它与其他细胞的特异性等等
DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'orig.ident') +
  ggtitle('harmony')

# 将分辨率改为0.3,分成12个cluster
scRNA_harmony$RNA_snn_res.0.3
Idents(scRNA_harmony) <- 'RNA_snn_res.0.3'
DimPlot(scRNA_harmony,reduction = 'umap')

# 绘图
umap_integrated_1 <- DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'orig.ident',label = T) 
umap_integrated_2 <- DimPlot(scRNA_harmony,reduction = 'tsne', label = T) 
umap_integrated_3 <- DimPlot(scRNA_harmony,reduction = 'umap', label = T) 

# 合并图片
integrated_plot <- CombinePlots(list(umap_integrated_1,umap_integrated_2,umap_integrated_3))
# 输出到画板
integrated_plot
# 保存
ggsave('integrated_plot.png',integrated_plot,width = 10,height = 10)
# 保存数据 # 降维步骤后保存
# save(scRNA_harmony,file = 'scRNA_harmony_2.Rdata')
# load('scRNA_harmony_2.Rdata')

#################################################### 8. 双细胞去除 ####################################################
############## 关于双细胞是否去除 #################
# 在质控中,去除双细胞有没有必要
# 这个操作难度就有点大,因为你不知道哪些细胞是双细胞
# 双细胞与测序和样本组织有关,特别是癌组织如果想要去除双细胞的话很容易就去除正常细胞
# 所以这个操作要谨慎处理,或者之后了解清楚后再做处理
# 去除双细胞的方法: DoubleFinder

############## 关于细胞周期的影响 #################
# 某些样本,有的细胞处于休眠期,有的细胞增殖期
# 如果把所有细胞周期的影响都去除掉,则会影响这两类细胞的鉴定
# 可以通过pca看细胞周期的影响,影响不大就不用回归掉

#################################################### 9. 细胞周期计算 ####################################################
# CellCycleScoring函数计算每个细胞的细胞周期得分，并将计算出的S期和G2/M期的评分保存在metadata中，
#以及细胞处于G2M，S或G1期的预测分类。
#通过设置set.ident = TRUE，则CellCycleScoring将Seurat对象中每个细胞的分组信息设置为其所处的细胞周期阶段。

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA_harmony <- CellCycleScoring(scRNA_harmony,s.features = s.genes,g2m.features = g2m.genes,set.ident = T)
# plot 
DimPlot(scRNA_harmony,group.by = 'Phase') # 这里可以看出细胞周期的影响很小,几乎都是均匀分布的
# 这里是正常的,不用回归细胞周期
DimPlot(scRNA_harmony,group.by = 'Phase',reduction = 'pca')

scRNA_harmony <- ScaleData(scRNA_harmony,vars.to.regress = c('S.Score','G2M.Score'),
                           features = rownames(sc_example))

############## 提取count data #################
# data 通常包含经过一些预处理（如标准化、缩放等）后的数据，它反映了经过一定转换后每个细胞中基因的表达量
# counts 则是原始的基因计数数据，也就是直接统计得到的每个细胞中各个基因的表达计数
# 1.counts 
# function1
scRNA_counts <- LayerData(scRNA_harmony,assay = 'RNA',layer = 'counts')
scRNA_counts <- as.data.frame(scRNA_counts)
# function2
scRNA_counts <- GetAssayData(scRNA_harmony,assay = 'RNA',layer = 'counts')
scRNA_counts <- as.data.frame(scRNA_counts)
#write.table(scRNA_counts,file = "./scRNA_counts.tsv", sep = "\t", quote = F, row.names = F)
# 2.data 
scRNA_data <- LayerData(scRNA_harmony,assay = 'RNA',layer = 'data')
scRNA_data <- as.data.frame(scRNA_data)
# 想要保存的话用write.table()函数
#write.table(scRNA_data,file = "./scRNA_data.tsv", sep = "\t", quote = F, row.names = F)
############## 求每个cluster的平均表达量AverageExpression #################

#AverageExpression 函数的计算逻辑
#data是 apply(count, 2, function(x){ log(x/sum(x)*1e4 +1) })

av_seurat <- AverageExpression(scRNA_harmony,assays = 'RNA',group.by = 'RNA_snn_res.0.3')
av_seurat <- as.data.frame(av_seurat$RNA)
av_seurat

############## 伪bulk数据AggregateExpression #################
av_seurat <- AggregateExpression(scRNA_harmony,assays = 'RNA',group.by = 'RNA_snn_res.0.3')
av_seurat <- as.data.frame(av_seurat)

av_seurat <- AggregateExpression(scRNA_harmony,assays = 'RNA',group.by = 'orig.ident')
av_seurat <- as.data.frame(av_seurat)

av_seurat <- AggregateExpression(scRNA_harmony,assays = 'RNA',group.by = c('RNA_snn_res.0.3','orig.ident'))
av_seurat <- as.data.frame(av_seurat)
write.table(av_seurat,file = "./av_seurat.tsv", sep = "\t", quote = F, row.names = F)
# AggregateExpression
# 伪bulk数据提供了一种将单细胞数据转化为表观上的批次数据的方法，它可以用于各种单细胞数据分析，包括细胞类型鉴定、差异表达分析、细胞群体分析和数据集集成等。
