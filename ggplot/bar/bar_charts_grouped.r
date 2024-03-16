  library(ggplot2)
  rm(list = ls())
  options(stringsAsFactors = F)
  
  setwd('D:\\Rscript')
  
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
  
 P<-ggplot(myData, aes(Name, Time, fill = Transport)) +
    geom_col(position = "dodge") +  # stack # dodge
    scale_fill_viridis_d() +
    scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
    labs(x = "City", y = "Time", title = "People") + 
    facet_grid(. ~ City, switch = 'x',scales = "free_x") +
    #scale_fill_manual(values = c("#8ecfc9","#ffbe7a","#fa7f6f","#82b0d2",
    #                             "#beb8dc","#e7dad2","#63e398","#a9b8c6")) +#"#9BB5CE", "#F78282",
    #theme_minimal(base_size = 14) +
    theme_classic(base_size = 14)+
    #theme_test()+ #全框
    #theme_light()+ #全框浅线
    #theme_bw()+
    theme(strip.placement = 'outside',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = 'bold'),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5, face = 'bold')
          )
  
  ggsave('Bar_charts_grouped.png',plot=P, height=12, width=12, units="in",dpi=600)
  ggsave('Bar_charts_grouped.pdf',plot=P, height=12, width=12, units="in",dpi=600)

# https://www.jianshu.com/p/00c31558bdea
# 翻转
if (FALSE):
{
if(FALSE){
  
  library(ggplot2)
  rm(list = ls())
  options(stringsAsFactors = F)
  
  setwd('D:\\Rscript')
  
  #myData2 = read.table('input.txt',header = T,sep='\t',check.names=F)
  #facet_grid(City~., switch = 'x',scales = "free_y")  #加了个 ,scales = "free_y" y轴根据分页组拆分(y轴确实可拆分的情况下)
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
  
  P<-ggplot(myData, aes( Time,Name, fill = Transport)) +
    geom_col(position = "dodge") + #  stack # dodge
    scale_fill_viridis_d() +
    scale_x_continuous(limits = c(0, 10), breaks = 0:10) +
    labs(x = "Time", y = "City", title = "People") + 
    facet_grid(City ~ ., switch = 'y',scales = "free_y") +
    #theme_minimal(base_size = 14) +
    theme_classic(base_size = 14)+
    theme(strip.placement = 'outside',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = 'bold'),
          
          axis.ticks = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5, face = 'bold'))
  
  ggsave('gene_fpkm_bar_grouped.png',plot=P, height=12, width=12, units="in",dpi=600)
  ggsave('gene_fpkm_bar_stacked.pdf',plot=P, height=12, width=12, units="in",dpi=600)
  
}
# 补充：分页标签背景，angle 文字倾斜
if(FALSE){
  theme(
    strip.background = element_blank(), #element_rect( fill="#F8F8FF"), # 
    strip.text.y.left = element_text(size=10, angle = 0), #facet_grid 文字标签角度 # strip.text.x
  #strip.text.x.bottom,
  #strip.text.x.top,
  #strip.text.y,
  #strip.text.y.left,
  #strip.text.y.right,
    axis.text.y=element_text(angle=70,vjust=1,hjust=1,face="bold",size=8),
  )+guides(fill=guide_legend(title=""))
}
