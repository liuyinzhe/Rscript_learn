library(ggplot2)
#library(gcookbook)
setwd("D://Rscript/SSR_bar")
f1 <- read.table("SSR.test.txt")
#f1 <- read.table("SSR.test.txt")
#f <-f1[1:100, ]
# R 从1开始索引计数
f$V1 <- f$V6*f$V7

f$V2[5 <= f$V1 & f$V1<= 8] = '5-8'
f$V2[8 <= f$V1 & f$V1<= 13] = '9-12'
f$V2[13 <= f$V1 & f$V1<= 16] = '13-16'
f$V2[17 <= f$V1 & f$V1<= 20] = '17-20'
f$V2[21 <= f$V1 & f$V1<= 24] = '21-24'
f$V2[25 <= f$V1] = 'more biger'
f$V6[f$V6 == 1] = 'Mono-'#重编码
f$V6[f$V6 == 2] = 'Di-'
f$V6[f$V6 == 3] = 'Tri-'
f$V6[f$V6 == 4] = 'Tetra-'
f$V6[f$V6 == 5] = 'Penta-'
f$V6[f$V6 == 6] = 'Hexa-'

#pdf("SSR.stat.pdf",width = 16/2.54, height = 10/2.54)
ggplot(f,aes(x = f$V6, y = f$V1, fill=f$V2))+geom_bar(stat="identity") +labs(
  x="SSR motif unit", #从新命名
  y="Repeat counts", #从新命名
  title="Distribution of SSR motifs")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(title='repeat_type'))+
  scale_fill_discrete(limits=c("5-8", "9-12", "13-16","17-20","21-24","more biger"))
  #+scale_fill_discrete(guide = FALSE)#图例关闭
#dev.off()
ggsave("SSR.stat.pdf",width = 16, height = 10, units="cm")
