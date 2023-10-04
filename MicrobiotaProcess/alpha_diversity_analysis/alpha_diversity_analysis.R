r·m(list = ls())
options(stringsAsFactors = F)
setwd("/data/draw/Rscript/MicrobiotaProcess/")
#加载R包
library(MicrobiotaProcess) # A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(phyloseq) # Handling and analysis of high-throughput microbiome census data


# 01 加载数据
# 加载
#   OTU特征表
#   OTU物种注释表
#   样品信息表
sample <- read.table("sample.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
#samples	group
OTU<- read.table("otu.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
#OTU_id	A_1	A_2	A_3	A_4	B_1	B_2	B_3	B_4	C_1	C_2	C_3	C_4
#OTU1	421	411	505	424	56	55	75	50	50	20	35	25
Tax <- read.table("tax.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
#OTU_id	Kingdom	Phylum	Class	Order	Family	Genus	Species
#OTU1	Bacteria	Actinobacteria	Actinobacteria	Actinomycetales	Thermomonosporaceae	Unassigned	Unassigned


# 02 利用phyloseq包重新构造可转换为分析的数据格式
ps <- phyloseq(sample_data(sample),
               otu_table(as.matrix(OTU), taxa_are_rows=TRUE), 
               tax_table(as.matrix(Tax)))

#转换数据格式
df <- ps %>% as.MPSE()
df
#     OTU   Sample Abundance RareAbundance
#    <chr>   <chr>      <int>         <int>
#   1 OTU1    A_1          421           410
#   2 OTU2    A_1            7             5

# 03 计算α多样性指数，包括Observe, Chao1, ACE, Shannon, Simpson, Pielou等指数

df %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

# 04 可视化
  
#可视化
f1 <- df %>% 
  mp_plot_alpha(
    .group=group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou),
    test = "t.test",#默认为"wilcox.test"
    size = 0.5,
    textsize = 4,
    map_signif_level = T,
    tip_length = 0
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF","red"), guide="none") +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF","red"), guide="none")

##不添加分组信息
f2 <- df %>%
  mp_plot_alpha(
    .alpha=c(Chao1, ACE, Shannon, Simpson)
  )

# 拼图
f1 / f2

draw_obj <- f1 / f2



png(paste("alpha_diversity",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_obj)
dev.off()
pdf(paste("alpha_diversity",".pdf",sep=""),width =10,height =7)
print(draw_obj)
dev.off()

# 06 计算结果提取
tbl <- df %>% 
  mp_extract_sample()

# 写文件
write.table(tbl, file =paste("alpha_diversity", "result.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
##########



# 07 通过没有连线的柱状堆积图来显示丰度
draw_bar_wo_obj <- df %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = group,#
    taxa.class = Genus,
    relative = T,
    topn = 10,
    geom = 'bar',#指定采用柱状图来显示相对丰度
    #features.dist = 'euclidean',#聚类方法
    #features.hclust = 'average',
    #sample.dist = 'bray',#计算距离方法
    #sample.hclust = 'average'
  )+
  theme(legend.position = "right")

# without line
png(paste("plot_abundance_bar_stacked_wo_line",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_heatmap_obj)
dev.off()
pdf(paste("plot_abundance_bar_stacked_wo_line",".pdf",sep=""),width =10,height =7)
print(draw_heatmap_obj)
dev.off()


# 08 按照分组合并数据


##对各组中排名前10的门丰度进行可视化
p5 <- df %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=group,
    taxa.class = Phylum,
    topn = 10,
    plot.group = TRUE#显示各组的丰度
  )+
  theme(legend.position = "right")

# without line
png(paste("plot_abundance_bar_stacked_wol_group_Phylum",".png",sep=""),width =900,height =600,res=100,units = "px")
print(p5)
dev.off()
pdf(paste("plot_abundance_bar_stacked_wol_group_Phylum",".pdf",sep=""),width =10,height =7)
print(p5)
dev.off()


# 属水平
p6 <- df %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group= group,
    taxa.class = Genus,
    topn = 10,
    relative = T,
    plot.group = TRUE
  )+
  theme(legend.position = "right")

png(paste("plot_abundance_bar_stacked_group_Genus",".png",sep=""),width =900,height =600,res=100,units = "px")
print(p6)
dev.off()
pdf(paste("plot_abundance_bar_stacked_group_Genus",".pdf",sep=""),width =10,height =7)
print(p6)
dev.off()


#拼图
merge_obj <- p5 / p6


png(paste("merge_group_Genus",".png",sep=""),width =900,height =600,res=100,units = "px")
print(merge_obj)
dev.off()
pdf(paste("merge_group_Genus",".pdf",sep=""),width =10,height =7)
print(merge_obj)
dev.off()
