#!/usr/bin/Rscript
library(ggplot2)
library(RColorBrewer)
library(getopt)

library(extrafont)
loadfonts(device="all")

spec=matrix(c(
  "help","h",0,"logical",
  "infile","i",1,"character",
  "outfile","o",1,"character",
  "group","g",1,"character",
  "list","l",1,"character"
),byrow=TRUE, ncol=4)
opt=getopt(spec)
usage<- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
Options: 
--help        -h     NULL         get this help
--infile     -i     character     the input file [forced]
--outfile    -o   character  the output file [forced]
--group    -g   character  the group file [forced]
--list    -l   character  the list file [forced]
\n")
  q(status=1)
}
if (is.null(opt$infile)){usage(spec)}
md<- read.table(opt$infile,header = TRUE,comment.char = "",stringsAsFactors = FALSE,sep="\t")



#cl<- c("#56A36C","#5E8579","#77C34F","#2E68AA","#7E884F","#7C8489","#4fB3A4","#F5B977","#FDFC7F","#48A7C2","#00B38C","#B1C914","#54AB11","#79E8D0","#FFB8B8","#86E65A","#24D197") 
cl<- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")
group<- read.table(opt$group,comment.char = "",header = TRUE,stringsAsFactors = FALSE)
colnames(group) <- c("sampleID","group")
sp<- read.table(opt$list,header = T,sep="\t")
list = rev(sp$ID)


#=======================
tt <- table(group$group)
md2 <-md
for( i in levels(factor(md$group))){
  print(i)
  md2[which(md2$group==i),3] <- md2[which(md2$group==i),3]/tt[i]
}
#=======================



if (length(rownames(group)) < 50){
  xt = 8
}else if (length(rownames(group)) >= 50 & length(rownames(group)) < 100 ){
  xt = 5
}else if (length(rownames(group)) >= 100 & length(rownames(group)) < 200){
  xt = 3
}else{ xt = 0}
if (grepl("absolute",opt$infile)){
p<- ggplot(data = md,aes(Sample,weight=Percent,fill=factor(Taxa,levels = list)))+
    geom_bar()+theme_gray()+scale_fill_manual(values = cl[1:length(unique(md$Taxa))])+
    theme(
        legend.title = element_blank(),
        #panel.grid = element_blank(),
        #panel.grid 绘图区网格线
        #panel.grid.major 主网格线
        #panel.grid.minor 次网格线
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype ="dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # 分页标签
        strip.text.x = element_text(family = "Arial",size = 10),
        strip.text.y = element_text(family = "Arial",size = 10),
        #axis.line 主轴线
        axis.line = element_line(colour = "black",linetype = "solid"),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 14),
        #axis.text.x = element_text(size=xt,angle = 90,hjust = 0.5),
        #axis.text, axis.text.x, axis.text.y,
        axis.text.x = element_text(family = "Arial",size=xt,angle = 90,hjust = 0.5),
        axis.text.y = element_text(family = "Arial",size = 10),
        axis.title.y = element_text(family = "Arial",size = 10),
        legend.text = element_text(family = "Arial",size = 10),
        #legend.title = element_text(family = "Arial",size = 10),
        #panel.grid = element_text(family = "Arial",size = 20),
        plot.title =element_text(family = "Arial",size = 20),
        #plot.subtitle
        #plot.caption
        #plot.tag
        
    )+
    # y轴0贴近x主轴
    scale_y_continuous(expand = c(0,0))+
    ylab("Absolute Abundance")+
    coord_cartesian(ylim = c(0.04,1))+
    facet_grid(.~group,scales = 'free',space = 'free')
    #scale_x_discrete(limits=group$sampleID)
    
    p2 <- ggplot(data = md2,aes(group,weight=Percent,fill=factor(Taxa,levels = list)))+
    geom_bar()+theme_gray()+scale_fill_manual(values = cl[1:length(unique(md$Taxa))])+
    theme(
    legend.title = element_blank(),
    #panel.grid = element_blank(),
    #panel.grid 绘图区网格线
    #panel.grid.major 主网格线
    #panel.grid.minor 次网格线
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60",linetype ="dotted"),
    panel.grid.minor = element_blank(),
    #axis.line 主轴线
    axis.line = element_line(colour = "black",linetype = "solid"),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 14),
    #axis.text.x = element_text(size=xt,angle = 90,hjust = 0.5),
    #
    axis.text.x = element_text(family = "Arial",size=xt,angle = 90,hjust = 0.5),
    axis.text.y = element_text(family = "Arial",size = 10),
    axis.title.y = element_text(family = "Arial",size = 14),
    legend.text = element_text(family = "Arial",size = 10),
    #legend.title = element_text(family = "Arial",size = 10),
    #panel.grid = element_text(family = "Arial",size = 20),
    plot.title =element_text(family = "Arial",size = 20),
    )+scale_y_continuous(expand = c(0,0))+
    ylab("Absolute Abundance")
  
  
  
}else{
  p<- ggplot(data = md,aes(Sample,weight=Percent,fill=factor(Taxa,levels = list)))+geom_bar()+theme_gray()+scale_fill_manual(values = cl[1:length(unique(md$Taxa))])+
    theme(
    legend.title = element_blank(),
    #panel.grid = element_blank(),
    #panel.grid 绘图区网格线
    #panel.grid.major 主网格线
    #panel.grid.minor 次网格线
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60",linetype ="dotted"),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    # 分页标签
    strip.text.x = element_text(family = "Arial",size = 10),
    strip.text.y = element_text(family = "Arial",size = 10),
    #axis.line 主轴线
    axis.line = element_line(colour = "black",linetype = "solid"),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 14),
    #axis.text.x = element_text(size=xt,angle = 90,hjust = 0.5),
    #
    axis.text.x = element_text(family = "Arial",size=xt,angle = 90,hjust = 0.5),
    axis.text.y = element_text(family = "Arial",size = 10),
    axis.title.y = element_text(family = "Arial",size = 14),
    legend.text = element_text(family = "Arial",size = 10),    
    #legend.title = element_text(family = "Arial",size = 10),
    #panel.grid = element_text(family = "Arial",size = 20),
    plot.title =element_text(family = "Arial",size = 20),
    )+
    scale_y_continuous(expand = c(0,0))+
    # y轴0贴近x主轴
    ylab("Relative Abundance")+
    coord_cartesian(ylim = c(0.04,1))+
    facet_grid(.~group,scales = 'free',space = 'free')
    #scale_x_discrete(limits=group$sampleID)
  
  
  p2 <- ggplot(data = md2,aes(group,weight=Percent,fill=factor(Taxa,levels = list)))+
    geom_bar()+theme_gray()+scale_fill_manual(values = cl[1:length(unique(md$Taxa))])+
    theme(
    legend.title = element_blank(),
    #panel.grid = element_blank(),
    #panel.grid 绘图区网格线
    #panel.grid.major 主网格线
    #panel.grid.minor 次网格线
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60",linetype ="dotted"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.line 主轴线
    axis.line = element_line(colour = "black",linetype = "solid"),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 14),
    #axis.text.x = element_text(size=xt,angle = 90,hjust = 0.5),
    axis.text.x = element_text(family = "Arial",size=xt,angle = 90,hjust = 0.5),
    axis.text.y = element_text(family = "Arial",size = 10),
    axis.title.y = element_text(family = "Arial",size = 14),
    legend.text = element_text(family = "Arial",size = 10),
    #legend.title = element_text(family = "Arial",size = 10),
    #panel.grid = element_text(family = "Arial",size = 20),
    plot.title =element_text(family = "Arial",size = 20),
    )+scale_y_continuous(expand = c(0,0))+
     ylab("Relative Abundance")
  
  
}


ggsave(paste(opt$outfile,".pdf",sep=""),plot=p, height=900, width=600, units="px",dpi=100)
ggsave(paste(opt$outfile,".sum.pdf",sep=""),plot=p2, height=900, width=600, units="px",dpi=100)

if(FALSE){ # 注释掉
#png(paste(opt$outfile,".png",sep=""),width =900,height =600,res=100,units = "px")
#print(p)
#dev.off()
pdf(paste(opt$outfile,".pdf",sep=""),width =10,height =7)
print(p)
dev.off()

#png(paste(opt$outfile,".sum.png",sep=""),width =900,height =600,res=100,units = "px")
#print(p2)
#dev.off()
pdf(paste(opt$outfile,".sum.pdf",sep=""),width =10,height =7)
print(p2)
dev.off()

}
