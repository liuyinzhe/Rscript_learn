
#library(this.path)
#setwd(this.path::this.dir())

library(ggplot2)


#getwd()
setwd("D:\\xx\\xx\\xx\\xx")

df_res=read.table("grouping_histogram.input.tsv",header=T,sep="\t")

#head(df_res$Sample1)

a<-ggplot(df_res, aes(x=SampleName, y=Freq, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  theme_bw()

a
ggsave("grouping_histogram.png",width = 10, height = 6, dpi = 300)
