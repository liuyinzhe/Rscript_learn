  library(ggplot2)
  rm(list = ls())
  options(stringsAsFactors = F)
  
  setwd('D:\\Rscript\\example')
  
  #myData2 = read.table('input.txt',header = T,sep='\t',check.names=F)
  #facet_grid(. ~ City, switch = 'x',scales = "free_x")  #加了个 ,scales = "free_x" x轴根据分页组拆分(x轴确实可拆分的情况下)
  # https://www.5axxw.com/questions/content/b5y5qz
  myData<-structure(list(Name = c("Rick", "Rick", "Rick", "Rick", "Rick", 
                                  "Rick", "Rick", "Rick", "Jane", "Jane", "Jane", "Jane", "Jane", 
                                  "Jane", "Jane", "Jane", "Ellen", "Ellen", "Ellen", "Ellen", "Ellen", 
                                  "Ellen", "Ellen", "Ellen"), 
                         City = c("Boston", "Boston", "Boston", 
                                  "Boston", "Seattle", "Seattle", "Seattle", "Seattle", "Boston", 
                                  "Boston", "Boston", "Boston", "Seattle", "Seattle", "Seattle", 
                                  "Seattle", "Boston", "Boston", "Boston", "Boston", "Seattle", 
                                  "Seattle", "Seattle", "Seattle"), 
                         Transport = c("Car", "Train", 
                                       "Bus ", "Plane", "Car", "Train", "Bus ", "Plane", "Car", "Train", 
                                       "Bus ", "Plane", "Car", "Train", "Bus ", "Plane", "Car", "Train", 
                                       "Bus ", "Plane", "Car", "Train", "Bus ", "Plane"), 
                         Time = c(0L, 
                                  1L, 0L, 9L, 0L, 0L, 3L, 7L, 1L, 3L, 2L, 0L, 0L, 2L, 3L, 0L, 1L, 
                                  3L, 3L, 4L, 8L, 4L, 7L, 7L)), 
                    class = "data.frame", row.names = c(NA, -24L))   

#https://www.jianshu.com/p/00c31558bdea


#myData = read.table('bar_plot.tsv',header = T,sep='\t',check.names=F)
#myData$significance <- factor(myData$significance,levels = c('up','down'))



P<-ggplot(myData, aes(Name, Time, fill = Transport)) +
  geom_col(position = "dodge") + #stack #dodge
  scale_fill_viridis_d() +
  scale_y_continuous() + #limits = c(0, 20), breaks = 0:20
  labs(x="City",y = "Time", title = "People") + 
  scale_fill_manual(values = c("#FF0000","#0000FF")) + #"#9BB5CE", "#F78282"
 
  #theme_minimal(base_size = 14) +
  theme_classic(base_size = 14)+
  #theme_grey()+
    facet_grid(. ~ kegg_class, scales = "free_x") + #switch = 'x',
  theme(#strip.placement = 'outside',  
        strip.background = element_blank(), #element_rect( fill="#F8F8FF", size=0.5), # 
        strip.text = element_text(size = 9, face = 'bold'),
        axis.text.x=element_text(angle=65,vjust=1,hjust=1,face="bold",size=10), #,family = 'sans'
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12, face = "bold"), # y轴
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(colour="#D3D3D3"),
        panel.grid.minor = element_line(colour="#EDEDED"),
        plot.title = element_text(size = 18, hjust = 0.5, face = 'bold')
        )+
  guides(fill=guide_legend(title=""))


ggsave('class_group.png',plot=P, height=12, width=12, units="in",dpi=600)
ggsave('class_group.pdf',plot=P, height=12, width=12, units="in",dpi=600)
