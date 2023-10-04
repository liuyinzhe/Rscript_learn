r·m(list = ls())
options(stringsAsFactors = F)
setwd("/data/draw/Rscript/MicrobiotaProcess/")
#加载R包
library(MicrobiotaProcess) # A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework # A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework
library(microeco) # Microbial Community Ecology Data Analysis
library(dplyr) # A Grammar of Data Manipulation
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics # Create Elegant Data Visualisations Using the Grammar of Graphics # Create Elegant Data Visualisations Using the Grammar of Graphics # Create Elegant Data Visualisations Using the Grammar of Graphics # Create Elegant Data Visualisations Using the Grammar of Graphics # Create Elegant Data Visualisations Using the Grammar of Graphics
library(phyloseq) # Handling and analysis of high-throughput microbiome census data

# 01 加载数据(以microeco自带数据集为例)

data(sample_info_16S)
data(otu_table_16S)
data(taxonomy_table_16S)
data(phylo_tree_16S)
sample_info_16S <- as.data.frame(sample_info_16S)
otu_table_16S <- as.data.frame(otu_table_16S)
taxonomy_table_16S <- as.data.frame(taxonomy_table_16S)
phylo_tree_16S <- as.data.frame(phylo_tree_16S)

# 写文件
write.table(sample_info_16S, file =paste("sample_info_16S", "data.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
# 写文件
write.table(otu_table_16S, file =paste("otu_table_16S", "data.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
# 写文件
write.table(taxonomy_table_16S, file =paste("taxonomy_table_16S", "data.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)
# 写文件
write.table(phylo_tree_16S, file =paste("phylo_tree_16S", "data.tsv", sep = "."), sep="\t",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# 02 数据清洗
#统一分类信息
taxonomy_table_16S %<>% tidy_taxonomy
## 构造microtable
df <- microtable$new(sample_table = sample_info_16S,
                     otu_table = otu_table_16S,
                     tax_table = taxonomy_table_16S,
                     phylo_tree = phylo_tree_16S,
                     auto_tidy = F)
##去除不属于非古菌和细菌的OTU
# df$tax_table %<>% subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
df$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
##去除“线粒体”和“叶绿体”污染
df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
##统一各数据的样本和OTU信息
df$tidy_dataset()
##检查序列号
df$sample_sums() %>% range
##重采样以减少测序深度对多样性测量的影响，使每个样本的序列号相等。
df$rarefy_samples(sample.size = 10000)
df$sample_sums() %>% range



# 03 利用phyloseq包重新构造可转换为分析的数据格式
ps <- phyloseq(sample_data(df$sample_table),
               otu_table(as.matrix(df$otu_table), taxa_are_rows=TRUE), 
               tax_table(as.matrix(df$tax_table)))

#转换数据格式
df <- ps %>% as.MPSE()
df


# 04 PCoA分析

###样本或组之间的距离计算
##数据标准化
df %<>% 
  mp_decostand(.abundance=Abundance,
               method = "hellinger" #数据标准化方法, 可选择'total', 'max', 'frequency', 'normalize', 'range', 'rank', 'rrank', 'standardize' 'pa', 'chi.square', 'hellinger' and 'log', 具体可见??decostand
               )

## 计算样本之间的距离，详细见??mp_cal_dist
df %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")

###PCoA分析,详见??mp_cal_pcoa
df %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray", action="add")

# 用adonis或anosim方法来检验各分组对群体的差异性是否显著
df %<>%
  mp_adonis(.abundance=hellinger, .formula=~Group, distmethod="bray", permutations=9999, action="add")
df %>% mp_extract_internal_attr(name=adonis)
#可视化PCoA结果
p1 <- df %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = Group, 
    .color = Group, 
    .size = 1.2,
    .alpha = 1,
    ellipse = T,
    show.legend = T # 是否显示置信圈的图例
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF","red")) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF","red")) 
# point的大小也可以映射到其他变量，如Observe或Shannon
# 然后alpha多样性和beta多样性将同时显示
df %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
df$Shannon
p2 <- df %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Group, 
    .color = Group, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE 
  ) +
  scale_fill_manual(
    values = c("#00A087FF", "#3C5488FF","red"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#00A087FF", "#3C5488FF","red"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
#更换其他分组
p3 <- df %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = Saline, 
    .color = Saline, 
    .size = 1.2,
    .alpha = 1,
    ellipse = T,
    show.legend = T # 是否显示置信圈的图例
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF")) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF")) 
p4 <- df %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = Type, 
    .color = Type, 
    .size = 1.2,
    .alpha = 1,
    ellipse = F,
    show.legend = T # 是否显示置信圈的图例
  ) +
  scale_fill_manual(values=c("#ff0000", "#fbb034","#ffdd00","#c1d82f","#00a4e4","#8a7967")) 
#拼图
p1 + p2 
p3 + p4

draw_shannon_obj <- p1+p2
draw_type_obj <- p3+p4



png(paste("PCoA.Shannon",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_shannon_obj)
dev.off()
pdf(paste("PCoA.Shannon",".pdf",sep=""),width =10,height =7)
print(draw_shannon_obj)
dev.off()


png(paste("PCoA.Type",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_type_obj)
dev.off()
pdf(paste("PCoA.Type",".pdf",sep=""),width =10,height =7)
print(draw_type_obj)
dev.off()
