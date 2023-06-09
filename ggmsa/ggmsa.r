rm(list = ls())
options(stringsAsFactors = F)

library(ggmsa)
library(ggplot2)


setwd('D:\\supplemental-ggmsa-main')

protein_sequences <-  "mafft_alignment_A.fasta"
pp<- ggmsa(protein_sequences,  char_width = 0.5,start = 0, end = 100, seq_name = TRUE,use_dot = TRUE) + geom_seqlogo() 
#show.legend = TRUE,

ggsave("Fig2.pdf", plot = pp, width = 12, height = 10,dpi = 100)
ggsave("Fig2.png", plot = pp, width = 12, height = 10,dpi = 100)
#ggsave("Fig2.tiff", plot = pp , width = 8, height = 10,dpi = 300)
