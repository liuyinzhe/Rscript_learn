rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript')
# 3个A，3个B group
group_str <- 'A_vs_B'
#paste("group_str", "xiaoming", sep = ".")

par(ask=T)
#差异分析
#https://www.zhihu.com/column/c_1522546675062800384

# 加载包
suppressMessages({
library(DESeq2)   
library(pheatmap)  # 用于作热图的包
library(ggplot2)   # 用于作图的包
library(ggrepel)   # geom_text_repel
library(tidyr)     # 数据处理 gather 函数
library(topGO)
})

#差异分析
raw_data = read.table('M_vs_K.counts.txt',header = T,row.names='Geneid')
countData <- as.matrix(raw_data[,6:ncol(raw_data)])

# 使用log2(count+1) 计算spearman相关系数，并进行归一化，绘制热图，查看样品之间的相关性
cor_heat<- pheatmap(scale(cor(log2(countData+1))))
#ggsave('cor_heatmap.png',p)
ggsave(paste(group_str, 'cor_heatmap.png', sep = "."),plot=cor_heat, units="in",dpi=300)
ggsave(paste(group_str, 'cor_heatmap.pdf', sep = "."),plot=cor_heat, units="in",dpi=300)

#Counts FPKM RPKM TPM CPM 的转化
#http://events.jianshu.io/p/2474dec2ab2f
#Counts FPKM RPKM TPM CPM 的转化
#http://events.jianshu.io/p/2474dec2ab2f

#FPKM/RPKM (Fragments/Reads Per Kilobase Million )  每千个碱基的转录每百万映射读取的Fragments/reads
#RPKM与FPKM分别针对单端与双端测序而言，计算公式是一样的
counts2FPKM <- function(count=count, efflength=efflen){ 
  PMSC_counts <- sum(count)/1e6   #counts的每百万缩放因子 (“per million” scaling factor) 深度标准化
  FPM <- count/PMSC_counts        #每百万reads/Fragments (Reads/Fragments Per Million) 长度标准化
  FPM/(efflength/1000)                                    
}
x=raw_data[,5:ncol(raw_data)]
xcout=raw_data[,6:ncol(raw_data)]
#x$M2 <- with(x, counts2FPKM(x$M2, x$Length))
#colSums(x)


#绘制 箱式图
# 列分组信息
condition <- factor(c(rep("A",3),rep("B",3)))
# 合成列分组信息
colData <- data.frame(row.names=colnames(countData), condition)

fpkm <- as.data.frame(apply(xcout,2,counts2FPKM,efflength=x$Length))

# 写 FPKM 
gene_FPKM=cbind(id=row.names(fpkm), fpkm)
write.table(gene_FPKM, file =paste(group_str, "gene_FPKM.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# 多列转为适合箱式图的数据
# https://zhuanlan.zhihu.com/p/451653391
n_fpkm = gather(fpkm)

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
  theme_bw()+
  labs(x="",y = "FPKM",title = "")+ylim(0,1)
#theme_minimal()+
#scale_fill_nejm()
#theme(legend.position ="none")
#fpkm_box
#ggsave('fpkm_box.png',fpkm_box)
fpkm_box
ggsave(paste(group_str, 'fpkm_box.png', sep = "."),plot=fpkm_box, units="in",dpi=300)
ggsave(paste(group_str, 'fpkm_box.pdf', sep = "."),plot=fpkm_box, units="in",dpi=300)







# 去除表达量过低的基因
countData <- countData[rowMeans(countData)>0,]  

# 列分组信息
condition <- factor(c(rep("A",3),rep("B",3)))
# 合成列分组信息
colData <- data.frame(row.names=colnames(countData), condition)
# colData ，行索引就是样品名，第一列对应分组；

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


################################# PCA检测 #####################################
rld <- rlog(dds, blind=FALSE)
#若用 rld 数据，还可使用DESeq2自带函数 ;condition 组名c('A','A','A','B','B','B')与列对应
pcaData <-  plotPCA(rld, ntop = 500, intgroup=c("condition"))#,returnData=TRUE)
#head(pcaData)
#plot(pcaData[,1:2],pch=19,col=c("red","red","red","blue","blue","blue"))
#text(pcaData[,1],pcaData[,2]+0.1,row.names(pcaData),cex=0.5)
ggsave(paste(group_str, 'PCA.png', sep = "."),plot=pcaData, units="in",dpi=300)
ggsave(paste(group_str, 'PCA.pdf', sep = "."),plot=pcaData, units="in",dpi=300)

################################# PCA检测 #####################################

#标准化
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
#将结果用result()函数来获取
res <- results(dds1)
#summary(res) 


#将结果储存为表格形式
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#筛选结果
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)
#输出
write.table(res1_total, file ="gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)
write.table(res1_up, file ="up_gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)
write.table(res1_down, file ="down_gene_pvalue.tsv", sep="\t",
            row.names = TRUE,col.names =TRUE, quote =TRUE)

#展示某个基因的表达量
#可以直接用DESeq2的plotCounts
dds <- makeExampleDESeqDataSet()

#
plotCounts(dds1, gene = "gene_name_ABCa",intgroup = "condition")   # 指定某个基因
#可用ggplot对图片样式进行修改，并用ggrepel进行标注
d <- data.frame(t(subset(countData,rownames(countData)=="gene_name_ABCa")))
gene1<-ggplot(d, aes(x = condition, y = gene_name_ABCa, color = condition))+
  geom_point(position=position_jitter(w=0.2,h=0))+
  geom_text_repel(aes(label=rownames(d)))+
  theme_bw()+
  #ggtitle("gene_name_ABCa")+
  theme(plot.title=element_text(hjust=0.5))

#ggsave('gene1.png', gene1)
ggsave(paste(group_str, 'gene1.png', sep = "."),plot=gene1, units="in",dpi=300)
ggsave(paste(group_str, 'gene1.pdf', sep = "."),plot=gene1, units="in",dpi=300)

d <- data.frame(t(subset(countData,rownames(countData)=="gene_name_ABCb")))
gene2<-ggplot(d, aes(x = condition, y = gene_name_ABCb, color = condition))+
  geom_point(position=position_jitter(w=0.2,h=0))+
  geom_text_repel(aes(label=rownames(d)))+
  theme_bw()+
  #ggtitle("gene_name_ABCb")+
  theme(plot.title=element_text(hjust=0.5))

#ggsave('gene2.png', gene1)
ggsave(paste(group_str, 'gene2.png', sep = "."),plot=gene2, units="in",dpi=300)
ggsave(paste(group_str, 'gene2.pdf', sep = "."),plot=gene2, units="in",dpi=300)


d <- data.frame(t(subset(countData,rownames(countData)=="gene_name_ABCc")))
gene3=ggplot(d, aes(x = condition, y = gene_name_ABCc, color = condition))+
  geom_point(position=position_jitter(w=0.2,h=0))+
  geom_text_repel(aes(label=rownames(d)))+
  theme_bw()+
  #ggtitle("gene_name_ABCc")+
  theme(plot.title=element_text(hjust=0.5))

#ggsave('gene3.png', gene3)
ggsave(paste(group_str, 'gene3.png', sep = "."),plot=gene3, units="in",dpi=300)
ggsave(paste(group_str, 'gene3.pdf', sep = "."),plot=gene3, units="in",dpi=300)


d <- data.frame(t(subset(countData,rownames(countData)=="gene_name_ABCd")))
gene4<-ggplot(d, aes(x = condition, y = gene_name_ABCc, color = condition))+
  geom_point(position=position_jitter(w=0.2,h=0))+
  geom_text_repel(aes(label=rownames(d)))+
  theme_bw()+
  #ggtitle("gene_name_ABCc")+
  theme(plot.title=element_text(hjust=0.5))

#ggsave('gene4.png', gene4)
ggsave(paste(group_str, 'gene4.png', sep = "."),plot=gene4, units="in",dpi=300)
ggsave(paste(group_str, 'gene4.pdf', sep = "."),plot=gene4, units="in",dpi=300)




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
ggsave(paste(group_str, 'heatmap_gene.png', sep = "."),plot=heatmap, units="in",dpi=300)
ggsave(paste(group_str, 'heatmap_gene.pdf', sep = "."),plot=heatmap, units="in",dpi=300)

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
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  # 图例
  theme(legend.position = "None",
        legend.title = element_blank(),
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
        
  ) +
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

ggsave(paste(group_str, 'volcano_gene.png', sep = "."),plot=volcano, units="in",dpi=300)
ggsave(paste(group_str, 'volcano_gene.pdf', sep = "."),plot=volcano, units="in",dpi=300)

gene_diff_data=cbind(id=row.names(genes), genes)

# 保存表格
#write.table(gene_diff_data, file ="gene_diff.tsv", sep="\t",
#            row.names = FALSE,col.names =TRUE, quote =TRUE)


# GO 分析
# 加载包
suppressMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(enrichplot)
})
# ID 转换
ids<- bitr(rownames(res1_total), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Rn.eg.db ,drop = T)

Up_ids<- bitr(rownames(res1_up), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Rn.eg.db ,drop = T)
Down_ids<- bitr(rownames(res1_down), fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Rn.eg.db ,drop = T)

############################总的GO############################
# GO富集分析
GO<-enrichGO( ids$ENTREZID,
              OrgDb = org.Rn.eg.db,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)

#barplot(GO)
#write.csv(GO,file="GO_enrich.csv")
# 写入文件
write.table(GO, file ="GO_enrich.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

GO_dot<- dotplot(GO,showCategory=10)
ggsave(paste(group_str,'GO_dot.png', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'GO_dot.pdf', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)

GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)
ggsave(paste(group_str,'GO_emap.png', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'GO_emap.pdf', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)

goCC <- enrichGO(ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goBP <- enrichGO(ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goMF <- enrichGO(ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste(group_str,'GO_bar.png', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'GO_bar.pdf', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)

GO_cnet<- cnetplot(GO)
ggsave(paste(group_str,'GO_cnet.png', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'GO_cnet.pdf', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)


pdf(file=paste(group_str,"enrich.go.CC.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(group_str,"enrich.go.BP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(group_str,"enrich.go.MP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()
############################总的 GO ############################

############################UP GO############################
GO<-enrichGO( Up_ids$ENTREZID,
              OrgDb = org.Rn.eg.db,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)

#barplot(GO)
#write.csv(GO,file="GO.csv")
# 写入文件
write.table(GO, file ="Up_GO_enrich.tsv.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

GO_dot<- dotplot(GO,showCategory=10)
ggsave(paste(group_str,'Up_GO_dot.png', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'Up_GO_dot.pdf', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)

GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)
ggsave(paste(group_str,'Up_GO_emap.png', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Up_GO_emap.pdf', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)

goCC <- enrichGO(Up_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goBP <- enrichGO(Up_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goMF <- enrichGO(Up_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste(group_str,'Up_GO_bar.png', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Up_GO_bar.pdf', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)

GO_cnet<- cnetplot(GO)
ggsave(paste(group_str,'Up_GO_cnet.png', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'Up_GO_cnet.pdf', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)


pdf(file=paste(group_str,"Up_enrich.go.CC.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(group_str,"Up_enrich.go.BP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(group_str,"Up_enrich.go.MP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()
############################UP GO############################

############################down GO############################
GO<-enrichGO( Down_ids$ENTREZID,
              OrgDb = org.Rn.eg.db,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)

#barplot(GO)
#write.csv(GO,file="GO.csv")
# 写入文件
write.table(GO, file ="Down_GO_enrich.tsv", sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

GO_dot<- dotplot(GO,showCategory=10)
ggsave(paste(group_str,'Down_GO_dot.png', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'Down_GO_dot.pdf', sep = "."),plot=GO_dot,height=8, width=12, units="in",dpi=300)

GO_pt <- pairwise_termsim(GO)
emap <- emapplot(GO_pt)
ggsave(paste(group_str,'Down_GO_emap.png', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Down_GO_emap.pdf', sep = "."),plot=emap,height=12, width=16, units="in",dpi=300)


goCC <- enrichGO(Down_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='CC',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goBP <- enrichGO(Down_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

goMF <- enrichGO(Down_ids$ENTREZID,
                 OrgDb = org.Rn.eg.db,
                 ont='MF',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 keyType = 'ENTREZID')

GO_bar<- barplot(GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste(group_str,'Down_GO_bar.png', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Down_GO_bar.pdf', sep = "."),plot=GO_bar,height=12, width=16, units="in",dpi=300)


GO_cnet<- cnetplot(GO)
ggsave(paste(group_str,'Down_GO_cnet.png', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'Down_GO_cnet.pdf', sep = "."),plot=GO_cnet,height=12, width=18, units="in",dpi=300)


pdf(file=paste(group_str,"Down_enrich.go.CC.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goCC)
dev.off()

pdf(file=paste(group_str,"Down_enrich.go.BP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goBP)
dev.off()

pdf(file=paste(group_str,"Down_enrich.go.MP.DGA.tree.pdf", sep = "."),width = 10,height = 15)
plotGOgraph(goMF)
dev.off()

############################down GO############################

############################Kegg############################

# kegg_enrich 富集分析
kegg_enrich <- enrichKEGG( ids$ENTREZID,
              organism   = "rno", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05)

kk_read <- DOSE::setReadable(kegg_enrich, 
                             OrgDb="org.Rn.eg.db", 
                             keyType='ENTREZID')#ENTREZID to gene Symbol

write.table(kk_read, file =paste(group_str, "KEGG_enrich.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(group_str,'KEGG_dot.png', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'KEGG_dot.pdf', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)
ggsave(paste(group_str,'KEGG_emap.png', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'KEGG_emap.pdf', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich)

ggsave(paste(group_str,'KEGG_bar.png', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'KEGG_bar.pdf', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kk_read)
ggsave(paste(group_str,'KEGG_cnet.png', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'KEGG_cnet.pdf', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)


############################Kegg############################

############################Up kegg############################
# kegg_enrich 富集分析
kegg_enrich <- enrichKEGG( Up_ids$ENTREZID,
                           organism   = "rno", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

kk_read <- DOSE::setReadable(kegg_enrich, 
                             OrgDb="org.Rn.eg.db", 
                             keyType='ENTREZID')#ENTREZID to gene Symbol

write.table(kk_read, file =paste(group_str, "Up_KEGG_enrich.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(group_str,'Up_KEGG_dot.png', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'Up_KEGG_dot.pdf', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)
ggsave(paste(group_str,'Up_KEGG_emap.png', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Up_KEGG_emap.pdf', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich)

ggsave(paste(group_str,'Up_KEGG_bar.png', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Up_KEGG_bar.pdf', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kk_read)
ggsave(paste(group_str,'Up_KEGG_cnet.png', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'Up_KEGG_cnet.pdf', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)

############################Up kegg############################

############################down kegg############################
# kegg_enrich 富集分析
kegg_enrich <- enrichKEGG( Down_ids$ENTREZID,
                           organism   = "rno", #需要小写，https://www.genome.jp/kegg/catalog/org_list.html
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

kk_read <- DOSE::setReadable(kegg_enrich, 
                             OrgDb="org.Rn.eg.db", 
                             keyType='ENTREZID')#ENTREZID to gene Symbol

write.table(kk_read, file =paste(group_str, "Down_KEGG_enrich.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


KEGG_dot<- dotplot(kegg_enrich,showCategory=10)
ggsave(paste(group_str,'Down_KEGG_dot.png', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)
ggsave(paste(group_str,'Down_KEGG_dot.pdf', sep = "."),plot=KEGG_dot,height=8, width=12, units="in",dpi=300)

KEGG_pt <- pairwise_termsim(kegg_enrich)
kegg_emap <- emapplot(KEGG_pt)
ggsave(paste(group_str,'Down_KEGG_emap.png', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Down_KEGG_emap.pdf', sep = "."),plot=kegg_emap,height=12, width=16, units="in",dpi=300)


kegg_bar<- barplot(kegg_enrich)

ggsave(paste(group_str,'Down_KEGG_bar.png', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)
ggsave(paste(group_str,'Down_KEGGO_bar.pdf', sep = "."),plot=kegg_bar,height=12, width=16, units="in",dpi=300)

kegg_cnet<- cnetplot(kk_read)
ggsave(paste(group_str,'Down_KEGG_cnet.png', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)
ggsave(paste(group_str,'Down_KEGG_cnet.pdf', sep = "."),plot=kegg_cnet,height=12, width=18, units="in",dpi=300)

############################down kegg############################
#https://zhuanlan.zhihu.com/p/518134934
