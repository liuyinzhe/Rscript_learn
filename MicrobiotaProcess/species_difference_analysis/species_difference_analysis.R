r·m(list = ls())
options(stringsAsFactors = F)
setwd("/data/draw/Rscript/MicrobiotaProcess/")
#加载R包
library(MicrobiotaProcess) # A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework
library(dplyr) # A Grammar of Data Manipulation
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(phyloseq) # Handling and analysis of high-throughput microbiome census data
library(ggtree) # an R package for visualization of tree and annotation data
library(ggtreeExtra) # An R Package To Add Geometric Layers On Circular Or Other Layout Tree Of "ggtree"
library(ggstar) # Multiple Geometric Shape Point Layer for 'ggplot2'
library(forcats) # Tools for Working with Categorical Variables (Factors)
library(ggnewscale)

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

# 03 计算物种的相对丰度

df %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=group
  )

# 04 物种差异分析
  

df %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = group,
    tip.level = "OTU",
    force = FALSE,
    relative = TRUE,
    taxa.class = "all",
    first.test.method = "kruskal.test",
    first.test.alpha = 0.05,
    p.adjust = "fdr",
    filter.p = "fdr",
    strict = TRUE,
    fc.method = "generalizedFC",
    second.test.method = "wilcox.test",
    second.test.alpha = 0.05,
    cl.min = 4,
    cl.test = TRUE,
    subcl.min = 3,
    subcl.test = TRUE,
    ml.method = "lda",# 'lda' or 'rf'
    normalization = 1e+06,
    ldascore = 2,#LDA阈值
    bootnums = 30,
    sample.prop.boot = 0.7,
    ci = 0.95,
    seed = 123,
    type = "species"
  )

# 06 提取结果并基于ggtree等R包进行可视化

# 结果提取
taxa.tree <- df %>% 
  mp_extract_tree(type="taxatree")
taxa.tree

taxa.tree %>% 
  select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_group, pvalue, fdr) %>%
  dplyr::filter(!is.na(fdr))

##通过ggtree和ggtreeExtra可视化
draw_obj <- ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3) +
  geom_point(data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21)+
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label))+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, group, .fun=min),
      size = RelRareAbundanceBySample,
      fill = group,
      subset = RelRareAbundanceBySample > 0),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","blue"))+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_group,
      subset = !is.na(LDAmean)),
    orientation = "y",
    offset = 0.05,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2))+
  geom_tiplab(size=2, offset=11.2)+
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_group)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_group,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02","blue"))+
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9))


png(paste("species_difference_analysis.cycle_graph",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_obj)
dev.off()
pdf(paste("species_difference_analysis.cycle_graph",".pdf",sep=""),width =10,height =7)
print(draw_obj)
dev.off()



# 07 通过没有连线的柱状堆积图来显示丰度

draw_simplify_MP_obj<- df %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.1
  ) +
  scale_fill_manual(values=c("deepskyblue", "orange","red")) +
  scale_fill_manual(
    aesthetics = "fill_new", 
    values = c("deepskyblue", "orange","red")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", 
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999"
    )
  )

png(paste("species_difference_analysis.simplify_cycle_graph",".png",sep=""),width =900,height =600,res=100,units = "px")
print(draw_simplify_MP_obj)
dev.off()
pdf(paste("species_difference_analysis.simplify_cycle_graph",".pdf",sep=""),width =10,height =7)
print(draw_simplify_MP_obj)
dev.off()
