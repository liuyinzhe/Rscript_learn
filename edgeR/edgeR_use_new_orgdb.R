
rm(list = ls())
options(stringsAsFactors = F)
setwd('/A_vs_WT')
group_str <- 'A_vs_WT'
#paste("group_str", "xiaoming", sep = ".")

par(ask=T)
##################需要修改的 ###########
#[0] 分组信息condition <- factor(c(rep("A",3),rep("WT",3))) #condition <- factor(c(rep("A",3),rep("WT",3)),levels = c("A","WT"))
#[1] 分组最重要的地方，有时候上面分组会混乱，这里一锤定音 res <- results(dds1,contrast = c('condition', 'A', 'WT'))
#[2] 物种数据库 org.Rn.eg.db替换org.Target.eg.db  # 只能对应两个字母，并且第一个字母大写
#[3] KEGG 注释enrichKEGG 函数 organism   = "mmu", #Rno #三个字母#需要小写，https://www.genome.jp/kegg/catalog/org_list.html
#
#####################

############################### 创建目录   ############################


mkdir <- function(new_dir=new_dir){
  if (!dir.exists(new_dir)){
    dir.create(new_dir,recursive = TRUE)
  } 
}  

KEGG_dir="./KEGG/"
GO_dir="./GO/"

KEGG_dir_up="./KEGG/Up"
KEGG_dir_down="./KEGG/Down"

GO_dir_up="./GO/Up"
GO_dir_down="./GO/Down"


mkdir(KEGG_dir_up)
mkdir(KEGG_dir_down)
mkdir(GO_dir_up)
mkdir(GO_dir_down)
############################## 读取处理 ##################################

#差异分析
#https://www.zhihu.com/column/c_1522546675062800384

# 加载包
suppressMessages({
  library(edgeR)  
  library(pheatmap)  # 用于作热图的包
  library(ggplot2)   # 用于作图的包
  library(ggrepel)   # geom_text_repel
  library(tidyr)     # 数据处理 gather 函数
})

# 多个hisat2 比对bam 按照 处理组3个样品与对照组3个样品bam 顺序作为featurecounts输入 后的结果是 all.counts.txt
#样品整体，相关系数
raw_data = read.table('all.counts.txt',header = T,row.names='Geneid',check.names=F)
countData <- as.matrix(raw_data[,6:ncol(raw_data)])


############################### Count 数据过滤 #################################
# 去除表达量过低的基因
countData <- countData[rowMeans(countData)>0,]  

############################## 聚类热图 ##################################

# 使用log2(count+1) 计算spearman相关系数，并进行归一化，绘制热图，查看样品之间的相关性
cor_heat<- pheatmap(scale(cor(log2(countData+1))))
#ggsave('cor_heatmap.png',p)
ggsave(paste(group_str, 'cor_heatmap.png', sep = "."),plot=cor_heat, height=12, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'cor_heatmap.pdf', sep = "."),plot=cor_heat, height=12, width=12, units="in",dpi=300)


########################################## 转换函数 ################################################
#Counts FPKM RPKM TPM CPM 的转化
#http://events.jianshu.io/p/2474dec2ab2f


#FPKM/RPKM (Fragments/Reads Per Kilobase Million )  每千个碱基的转录每百万映射读取的Fragments/reads
#RPKM与FPKM分别针对单端与双端测序而言，计算公式是一样的
counts2FPKM <- function(count=count, efflength=efflen){ 
  PMSC_counts <- sum(count)/1e6   #counts的每百万缩放因子 (“per million” scaling factor) 深度标准化
  FPM <- count/PMSC_counts        #每百万reads/Fragments (Reads/Fragments Per Million) 长度标准化
  FPM/(efflength/1000)                                    
}


#TPM (Transcripts Per Kilobase Million)  每千个碱基的转录每百万映射读取的Transcripts
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)   #每千碱基reads (reads per kilobase) 长度标准化
  PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
  RPK/PMSC_rpk                    
}                   

#FPKM与TPM的转化
FPKM2TPM <- function(fpkm){
  fpkm/sum(fpkm)*1e6
}

########################################## 转换函数 ################################################


########################################## count 转 fpkm/tpm  ################################################

# 去掉染色体和坐标列，保留外显子长度
x=raw_data[,5:ncol(raw_data)]
# 去掉染色体坐标与外显子长度列
xcout=raw_data[,6:ncol(raw_data)]
#x$M2 <- with(x, counts2FPKM(x$M2, x$Length))
#colSums(x)


fpkm <- as.data.frame(apply(xcout,2,counts2FPKM,efflength=x$Length))
# 过滤掉全为0的
fpkm <- fpkm[rowMeans(fpkm)>0,]

tpm <- as.data.frame(apply(xcout,2,counts2TPM,efflength=x$Length))
tpm <- tpm[rowMeans(tpm)>0,]

# 写 FPKM 
# cbind 行index 加入新列变量名GeneID，就是新名字
gene_FPKM=cbind(GeneID=row.names(fpkm), fpkm)
write.table(gene_FPKM, file =paste(group_str, "gene_FPKM.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# 写 TPM 
gene_TPM=cbind(GeneID=row.names(tpm), tpm)
write.table(gene_TPM, file =paste(group_str, "gene_TPM.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


##########################################  分组信息 ################################################
#绘制 箱式图
# 列分组信息
condition <- factor(c(rep("A",1),rep("WT",1)),levels = c("A","WT"))
# 合成列分组信息
colData <- data.frame(row.names=colnames(countData), condition)

########################################## FPKM/TPM 箱式图 log2(x+1) ################################################

# 多列转为适合箱式图的数据
# 同时转为 log2(x+1)
n_fpkm = gather(log2(fpkm+1))
n_tpm = gather(log2(tpm+1))

# FPKM
df1 <- data.frame(key = colnames(countData),
                  group = condition
)
new_fpkm = merge(n_fpkm,df1,type=left)

fpkm_box<- ggplot(new_fpkm, aes(x=key,y=value,fill=group)) + 
  geom_boxplot()+
  stat_boxplot(geom = "errorbar",
               lwd=0.5,
               width=0.5)+
  #ggtitle("分组绘制箱式图")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = rel(0.8)))+
  labs(x="Sample",y = "log2(FPKM+1)",title = "")#+ylim(0,1)
#theme_minimal()+
#scale_fill_nejm()
#theme(legend.position ="none")
#fpkm_box
#ggsave('fpkm_box.png',fpkm_box)
#fpkm_box
ggsave(paste(group_str, 'fpkm_box.png', sep = "."),plot=fpkm_box, height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'fpkm_box.pdf', sep = "."),plot=fpkm_box, height=8, width=12, units="in",dpi=300)

# TPM

new_tpm = merge(n_tpm,df1,type=left)

tpm_box<- ggplot(new_tpm, aes(x=key,y=value,fill=group)) + 
  geom_boxplot()+
  stat_boxplot(geom = "errorbar",
               lwd=0.5,
               width=0.5)+
  #ggtitle("分组绘制箱式图")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = rel(0.8)))+
  labs(x="Sample",y = "log2(TPM+1)",title = "")#+ylim(0,1)
#theme_minimal()+
#scale_fill_nejm()
#theme(legend.position ="none")
#fpkm_box
#ggsave('fpkm_box.png',fpkm_box)
#tpm_box
ggsave(paste(group_str, 'tpm_box.png', sep = "."),plot=tpm_box, height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'tpm_box.pdf', sep = "."),plot=tpm_box, height=8, width=12, units="in",dpi=300)


########################################## FPKM/TPM 箱式图 ################################################




################################# 差异分析 #####################################
# 使用edgeR进行无重复差异表达分析
# https://blog.csdn.net/u012110870/article/details/102804557
# edgeR进行无重复的基因差异表达分析
# https://www.bilibili.com/read/cv19693692/


# 指定分组
group <- rep(c('A','WT'), each = 1)
# 构建 DGEList 对象：
dgelist <- DGEList(counts = countData, group = group)

## 数据过滤，保留在至少在一个样本里有表达的基因(CPM > 1)
keep = rowSums(cpm(dgelist)>1) >=1 

dgelist= dgelist[keep, , keep.lib.sizes = FALSE] 

## 标准化
dgelist = calcNormFactors(dgelist) 


## 表达差异分析
y_bcv = dgelist
bcv = 0.01  # 注意调整BCV
et = exactTest(y_bcv, dispersion = bcv ^ 2)


## 查看上调与下调的基因数量，太少可以适当调小BCV
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)

## 提取结果并筛选
genes = et$table

res_tmp = genes[which(genes$PValue<0.05),]
res_diff = res_tmp[which(res_tmp$logFC>1 | -1>res_tmp$logFC),]
head(res_diff)

# 根据上调、下调、不变为基因添加颜色信息
# padj<0.05 并且log2 (fold change) 绝对值大于=1,红色；logFC 大于1 取蓝色
genes$color <- NA
genes$color <- ifelse(genes$PValue<0.05 & abs(genes$logFC)>= 1,ifelse(genes$logFC > 1,'Up','Down'),'Stable')


# 替换变量修改名字
res = genes
colnames(res) <- c("log2FoldChange","logCPM","pvalue",'color')

#将结果储存为表格形式
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1_w = gene_diff_data=cbind(GeneId=row.names(res1), res1)
#输出
write.table(res1_w, file ="all_gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)

#筛选结果
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)
#输出
#write.table(res1_total, file ="diff_gene_pvalue.tsv", sep="\t",
#            row.names = TRUE,col.names =TRUE, quote =TRUE)

res1_up_w=cbind(GeneId=row.names(res1_up), res1_up)
res1_down_w=cbind(GeneId=row.names(res1_down), res1_down)

write.table(res1_up_w, file ="up_gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)
write.table(res1_down_w, file ="down_gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)

################################# 差异基因热图 #####################################
df <- countData[intersect(rownames(countData),rownames(res1_total)),]    
# 在原表达矩阵中找到差异表达基因
df2<- as.matrix(df)                                                 
heatmap <- pheatmap(df2,
                    show_rownames = F,
                    show_colnames = T,
                    cluster_cols = F,
                    cluster_rows=T,
                    height=10,  
                    scale = "row",
                    frontsize = 10,
                    angle_col=45, 
                    color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
                    clustering_method = 'single',
) 

#ggsave('heatmap_gene.png', heatmap)
ggsave(paste(group_str, 'heatmap_gene.png', sep = "."),plot=heatmap,height=12, width=10, units="in",dpi=300)
ggsave(paste(group_str, 'heatmap_gene.pdf', sep = "."),plot=heatmap,height=12, width=10, units="in",dpi=300)


################################# 火山图 #####################################
# 火山图
# fold change的意思是样本间表达量的差异倍数，log2 fold change的意思是取log2
# Padj是P值矫正之后的数值,一般选取小于等于0.05(显著差异)的基因
genes<- res1
# 根据上调、下调、不变为基因添加颜色信息


# padj<0.05 并且log2 (fold change) 绝对值大于=1,红色；log2FoldChange 大于1 取蓝色
genes$color <- ifelse(genes$pvalue<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'Up','Down'),'Stable')
color <- c(Up = "red",Stable = "gray",Down = "blue")

# 条件提取筛选出有意义的
genes_tmp <- subset(genes,color=="Up"|color=="Down")
# 指定排序
# https://www.bilibili.com/video/BV1q7411c79S
genes_tmp2<- genes_tmp[order(genes_tmp$pvalue,decreasing = FALSE),]

# 取top10,head 好处就是处理少于10行数据方便
top_10 <- head(genes_tmp2,5)
top_10$symbol <- rownames(top_10)

volcano <- ggplot(
  # 指定数据、映射、颜色
  genes, aes(log2FoldChange, -log10(pvalue), col = color)) +  
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  # 辅助线
  labs(x="log2 (fold change)",y="-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  # 图例
  theme(legend.position = "None",
        legend.title = element_blank(),
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
        
  ) + coord_cartesian(xlim =c(-10, 10))+
  # 添加 top10
  geom_label_repel(data = top_10,
                   aes(log2FoldChange, -log10(pvalue), label = symbol,color = color),
                   size = 2)
# 注释
#geom_text_repel(
#  data = subset(genes, pvalue < 1e-100 & abs(genes$log2FoldChange) >= 10),
#  aes(label = rownames(genes)),
#  size = 5,
#  box.padding = unit(0.35, "lines"),
#  point.padding = unit(0.3, "lines"))

#ggsave('volcano_gene.png', volcano)

ggsave(paste(group_str, 'volcano_gene.png', sep = "."),plot=volcano,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str, 'volcano_gene.pdf', sep = "."),plot=volcano,height=8, width=12, units="in",dpi=300)

#res1_total
gene_diff_data=cbind(GeneId=row.names(genes), genes)

# 保存表格
write.table(gene_diff_data, file ="differential_gene.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
###############################
library(dplyr)
load(file = "kegg_info.RData")
write.table(pathway2name,file = "./pathway2name.txt", sep = "\t", quote = F, row.names = F)

#remove.packages("org.Target.eg.db")
install.packages("./org.Target.eg.db", repos=NULL,type='source')

library(org.Target.eg.db)
columns(org.Target.eg.db)

# 查看所有基因
keys(org.Target.eg.db)
# 查看特定基因的信息
# library(dplyr)
#select(org.Target.eg.db, keys = "CW07G09620", columns = c("GO"))

pathway2gene <- AnnotationDbi::select(org.Target.eg.db, 
                                      keys = keys(org.Target.eg.db), 
                                      columns = c("Pathway","Ko")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

write.table(pathway2gene,file = "./pathway2gene.txt", sep = "\t", quote = F, row.names = F)

################################# GO注释 #####################################
# GO 分析
# 加载包
suppressMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Target.eg.db)
  library(enrichplot)
})
# ID 转换
ids<-rownames(res1_total)
Up_ids<-rownames(res1_up)
Down_ids<-rownames(res1_down)      
#ids<- bitr(rownames(res1_total), fromType = "SYMBOL",toType = c( "GID"),OrgDb = org.Target.eg.db ,drop = T)

#Up_ids<- bitr(rownames(res1_up), fromType = "SYMBOL",toType = c( "GID"),OrgDb = org.Target.eg.db ,drop = T)
#Down_ids<- bitr(rownames(res1_down), fromType = "SYMBOL",toType = c( "GID"),OrgDb = org.Target.eg.db ,drop = T)

############################总的GO############################
# GO富集分析

GO <- enrichGO(gene=ids,
               OrgDb=org.Target.eg.db,
               keyType="GID",
               ont="ALL",   #CC/BP/MF可选
               pvalueCutoff = 0.05, # 默认值0.05
               pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
               minGSSize = 10)  # 默认值，富集的最小基因数量



#barplot(GO)
#write.csv(GO,file="GO_enrich.csv")
# 写入文件
write.table(GO, file =paste(GO_dir, "GO_enrich.tsv", sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

goCC <- enrichGO(gene=ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量


goBP <- enrichGO(gene=ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量


goMF <- enrichGO(gene=ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量

################################  barplot #####################################
GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste(GO_dir, paste(group_str,'GO_bar.png', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_bar.pdf', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)

################################  dotplot #####################################
GO_dot<- dotplot(GO,showCategory=10)

ggsave(paste(GO_dir, paste(group_str,'GO_dot.png', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)


GO_dot<- dotplot(goCC,showCategory=10)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.CC.png', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.CC.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)


GO_dot<- dotplot(goMF,showCategory=10)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.MF.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.MF.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)


GO_dot<- dotplot(goBP,showCategory=10)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.BP.png', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_dot.BP.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)

################################  emapplot #####################################
GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)

ggsave(paste(GO_dir, paste(group_str,'GO_emap.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_emap.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)

GO_pt <- pairwise_termsim(goCC)
emap <- emapplot(GO_pt)

ggsave(paste(GO_dir, paste(group_str,'GO_emap.CC.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_emap.CC.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)

GO_pt <- pairwise_termsim(goMF)
emap <- emapplot(GO_pt)

ggsave(paste(GO_dir, paste(group_str,'GO_emap.MF.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_emap.MF.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)

GO_pt <- pairwise_termsim(goBP)
emap <- emapplot(GO_pt)

ggsave(paste(GO_dir, paste(group_str,'GO_emap.BP.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_emap.BP.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
################################  cnetplot #####################################
GO_cnet<- cnetplot(GO)

ggsave(paste(GO_dir, paste(group_str,'GO_cnet.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_cnet.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)

GO_cnet<- cnetplot(goCC)

ggsave(paste(GO_dir, paste(group_str,'GO_cnet.cc.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_cnet.cc.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
GO_cnet<- cnetplot(goMF)

ggsave(paste(GO_dir, paste(group_str,'GO_cnet.MF.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_cnet.MF.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
GO_cnet<- cnetplot(goBP)

ggsave(paste(GO_dir, paste(group_str,'GO_cnet.BP.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir, paste(group_str,'GO_cnet.BP.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)

################################  DGA  #########################################


pdf(file=paste(GO_dir, paste(group_str,'enrich.go.CC.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(GO_dir, paste(group_str,'enrich.go.BP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(GO_dir, paste(group_str,'enrich.go.MP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()
############################总的 GO ############################

############################UP GO############################

GO <- enrichGO(gene=Up_ids,
               OrgDb=org.Target.eg.db,
               keyType="GID",
               ont="ALL",   #CC/BP/MF可选
               pvalueCutoff = 0.05, # 默认值0.05
               pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
               minGSSize = 10)  # 默认值，富集的最小基因数量



#barplot(GO)
#write.csv(GO,file="GO.csv")
# 写入文件
#
write.table(GO, file =paste(GO_dir_up, paste(group_str,'Up_GO_enrich.tsv.tsv', sep = "."), sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

GO_dot<- dotplot(GO,showCategory=10)

ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_dot.png', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_dot.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)

GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)

ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_emap.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_emap.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)

goCC <- enrichGO(gene=Up_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量


goBP <- enrichGO(gene=Up_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量


goMF <- enrichGO(gene=Up_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量


GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")


ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_bar.png', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_bar.pdf', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)

GO_cnet<- cnetplot(GO)

ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_cnet.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir_up, paste(group_str,'Up_GO_cnet.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)


pdf(file=paste(GO_dir_up, paste(group_str,'Up_enrich.go.CC.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(GO_dir_up, paste(group_str,'Up_enrich.go.BP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(GO_dir_up, paste(group_str,'Up_enrich.go.MP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()
############################UP GO############################

############################down GO############################
GO<-enrichGO( gene=Down_ids,
              OrgDb = org.Target.eg.db,
              keyType = "GID",
              ont = "ALL",
              pvalueCutoff = 0.05, # 默认值0.05
              pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
              minGSSize = 10)  # 默认值，富集的最小基因数量


#barplot(GO)
#write.csv(GO,file="GO.csv")
# 写入文件

write.table(GO, file =paste(GO_dir_down, paste(group_str,'Down_GO_enrich.tsv', sep = "."), sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

GO_dot<- dotplot(GO,showCategory=10)

ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_dot.png', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_dot.pdf', sep = "."), sep = "/"),plot=GO_dot,height=8, width=12, units="in",dpi=300)


GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)
ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_emap.png', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_emap.pdf', sep = "."), sep = "/"),plot=emap,height=12, width=16, units="in",dpi=300)


goCC <- enrichGO(gene=Down_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量

goBP <- enrichGO(gene=Down_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量

goMF <- enrichGO(gene=Down_ids,
                 OrgDb = org.Target.eg.db,
                 keyType="GID",
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, # 默认值0.05
                 pvalueCutoff = 0.2,  # 默认值0.2，是pvalue 的校正值 更严格
                 minGSSize = 10)  # 默认值，富集的最小基因数量

GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")


ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_bar.png', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_bar.pdf', sep = "."), sep = "/"),plot=GO_bar,height=12, width=16, units="in",dpi=300)


GO_cnet<- cnetplot(GO)

ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_cnet.png', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(GO_dir_down, paste(group_str,'Down_GO_cnet.pdf', sep = "."), sep = "/"),plot=GO_cnet,height=12, width=18, units="in",dpi=300)


pdf(file=paste(GO_dir_down, paste(group_str,'Down_enrich.go.CC.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(GO_dir_down, paste(group_str,'Down_enrich.go.BP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(GO_dir_down, paste(group_str,'Down_enrich.go.MP.DGA.tree.pdf', sep = "."), sep = "/"),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()

############################down GO############################


################################# Kegg注释 #####################################
# kegg_enrich 富集分析
#kegg_enrich <- enrichKEGG( ids,
#              organism   = "mmu", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
#              pvalueCutoff = 0.05,
#              pvalueCutoff = 0.2)

kegg_enrich <- enricher(gene=ids, 
                        TERM2GENE = pathway2gene, 
                        TERM2NAME = pathway2name, 
                        pvalueCutoff = 0.05,  # 表示全部保留，可以设为0.05作为阈值
                        pvalueCutoff = 0.2, # 表示全部保留，可以设为0.2作为阈值
                        pAdjustMethod = "BH",
                        minGSSize = 10) # 默认值10，最小富集的基因数量

#kk_read <- DOSE::setReadable(kegg_enrich, 
#OrgDb="org.Target.eg.db", 
#keyType='ENTREZID')#ENTREZID to gene Symbol

write.table(kegg_enrich, file =paste(KEGG_dir, paste(group_str,'KEGG_enrich.tsv', sep = "."), sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(KEGG_dir, paste(group_str,'KEGG_dot.png', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(KEGG_dir, paste(group_str,'KEGG_dot.pdf', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)

ggsave(paste(KEGG_dir, paste(group_str,'KEGG_emap.png', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir, paste(group_str,'KEGG_emap.pdf', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich,showCategory=10)

ggsave(paste(KEGG_dir, paste(group_str,'KEGG_bar.png', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir, paste(group_str,'KEGG_bar.pdf', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kegg_enrich)

ggsave(paste(KEGG_dir, paste(group_str,'KEGG_cnet.png', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(KEGG_dir, paste(group_str,'KEGG_cnet.pdf', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)


############################Kegg############################

############################Up kegg############################
# kegg_enrich 富集分析
#kegg_enrich <- enrichKEGG( Up_ids,
#                           organism   = "mmu", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
#                           pvalueCutoff = 0.05,
#                           pvalueCutoff = 0.2)

kegg_enrich <- enricher(gene=Up_ids, 
                        TERM2GENE = pathway2gene, 
                        TERM2NAME = pathway2name, 
                        pvalueCutoff = 0.05,  # 表示全部保留，可以设为0.05作为阈值
                        pvalueCutoff = 0.2, # 表示全部保留，可以设为0.2作为阈值
                        pAdjustMethod = "BH",
                        minGSSize = 10) # 默认值10，最小富集的基因数量


write.table(kegg_enrich, file =paste(KEGG_dir_up, paste(group_str,'Up_KEGG_enrich.tsv', sep = "."), sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_dot.png', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_dot.pdf', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)

ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_emap.png', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_emap.pdf', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich,showCategory=10)

ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_bar.png', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_bar.pdf', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kegg_enrich)

ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_cnet.png', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(KEGG_dir_up, paste(group_str,'Up_KEGG_cnet.pdf', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)

############################Up kegg############################

############################down kegg############################
# kegg_enrich 富集分析
#kegg_enrich <- enrichKEGG( Down_ids,
#                           organism   = "mmu", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
#                           pvalueCutoff = 0.05,
#                           pvalueCutoff = 0.2)

kegg_enrich <- enricher(gene=Down_ids, 
                        TERM2GENE = pathway2gene, 
                        TERM2NAME = pathway2name, 
                        pvalueCutoff = 0.05,  # 表示全部保留，可以设为0.05作为阈值
                        pvalueCutoff = 0.2, # 表示全部保留，可以设为0.2作为阈值
                        pAdjustMethod = "BH",
                        minGSSize = 10) # 默认值10，最小富集的基因数量



write.table(kegg_enrich, file =paste(KEGG_dir_down, paste(group_str,'Down_KEGG_enrich.tsv', sep = "."), sep = "/"), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_dot.png', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_dot.pdf', sep = "."), sep = "/"),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)

ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_emap.png', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_emap.pdf', sep = "."), sep = "/"),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich,showCategory=10)

ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_bar.png', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_bar.pdf', sep = "."), sep = "/"),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kegg_enrich)

ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_cnet.png', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(KEGG_dir_down, paste(group_str,'Down_KEGG_cnet.pdf', sep = "."), sep = "/"),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)

############################down kegg############################
#https://zhuanlan.zhihu.com/p/518134934


