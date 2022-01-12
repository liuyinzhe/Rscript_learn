#ref# https://guangchuangyu.github.io/2019/11/geom-motif/
#ref# https://yulab-smu.top/treedata-book/chapter13.html

library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)

get_genes <- function(data, genome) {
  filter(data, molecule == genome) %>% pull(gene)
}

g <- unique(example_genes[,1])
n <- length(g)
d <- matrix(nrow = n, ncol = n)
rownames(d) <- colnames(d) <- g
genes <- lapply(g, get_genes, data = example_genes)

for (i in 1:n) {
  for (j in 1:i) {
    jaccard_sim <- length(intersect(genes[[i]], genes[[j]])) / 
      length(union(genes[[i]], genes[[j]]))
    d[j, i] <- d[i, j] <- 1 - jaccard_sim
  }
}

tree <- ape::bionj(d) 

p <- ggtree(tree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = example_genes, geom = geom_motif, panel = 'Alignment',
             on = 'genE', label = 'gene', align = 'left') +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_continuous(expand=c(0,0)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'))

facet_widths(p, widths=c(1,2))
