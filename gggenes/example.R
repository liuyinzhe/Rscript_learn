
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("BiocManager")
#library('BiocManager', quietly = TRUE)
#BiocManager::install(c("gggenes",'ggtree','ggplot2'),force =TRUE)

library(ggplot2)
library('gggenes')
#windows platform 
setwd('D:\\Rscript')

df <-read.csv('example.txt',header = T, sep = '\t')
#df = read.csv('example.txt',header = F, sep = '\t')
#name(df) <- c('molecule','gene','start','end','strand','direction')
#df$molecule<-factor(df$molecule)


#head(df)
ggplot(df, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

