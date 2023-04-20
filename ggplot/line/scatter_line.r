rm(list = ls())
setwd("D:\\Rscript\\line")

library(ggplot2)

#https://zhuanlan.zhihu.com/p/479361842

df <- read.csv("input.csv",header=T,check.names=F) #,sep = '\t'
#species   sample type  value
#A  a 1M  10
#B  b 2M  20
#C  c 3M  30
#D  d 4M  40


p0<-ggplot(data = df, 
           aes(x = type, 
               y = value, group = sample,
               color = sample,
               shape = species))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) #legend.box.background = element_rect(color = NA))   #为图例加边框线

p0
ggsave('all.v0.png',plot=p0, height=12, width=15, units="in",dpi=300)
ggsave('all.v0.pdf',plot=p0, height=12, width=15, units="in",dpi=600)

p1<-ggplot(data = df, 
               aes(x = type, 
                   y = value, group = sample,
                   color = sample,
                   shape = species))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) +#legend.box.background = element_rect(color = NA))   #为图例加边框线
  facet_grid(. ~ species, scales = "free_x")

p1

ggsave('all.png',plot=p1, height=12, width=15, units="in",dpi=300)
ggsave('all.pdf',plot=p1, height=12, width=15, units="in",dpi=600)



library(tidyverse)

df2 <- df %>% as_tibble() %>% filter(species=="A")

p2<-ggplot(data = df2, 
           aes(x = type, 
               y = value, group = sample,
               color = sample,
               shape = sample))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) #legend.box.background = element_rect(color = NA))   #为图例加边框线

p2
ggsave('A.png',plot=p2, height=12, width=15, units="in",dpi=300)
ggsave('A.pdf',plot=p2, height=12, width=15, units="in",dpi=600)




df3 <- df %>% as_tibble() %>% filter(species=="B")

p3<-ggplot(data = df3, 
           aes(x = type, 
               y = value, group = sample,
               color = sample,
               shape = sample))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) #legend.box.background = element_rect(color = NA))   #为图例加边框线

p3
ggsave('B.png',plot=p3, height=12, width=15, units="in",dpi=300)
ggsave('B.pdf',plot=p3, height=12, width=15, units="in",dpi=600)




df4 <- df %>% as_tibble() %>% filter(species=="C")

p4<-ggplot(data = df4, 
           aes(x = type, 
               y = value, group = sample,
               color = sample,
               shape = sample))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) #legend.box.background = element_rect(color = NA))   #为图例加边框线

p4
ggsave('C.png',plot=p4, height=12, width=15, units="in",dpi=300)
ggsave('C.pdf',plot=p4, height=12, width=15, units="in",dpi=600)


df5 <- df %>% as_tibble() %>% filter(species=="D")

p5<-ggplot(data = df5, 
           aes(x = type, 
               y = value, group = sample,
               color = sample,
               shape = sample))+
  geom_point(size=3)+
  geom_line(size=1)+
  xlab("")+
  ylab("Number of Gene")+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.minor =  element_blank(),    #以上代码用于去除网格线且保留坐标轴边框
        legend.position = c(.075, .655),    #更改图例位置，放至图内部的左上角
        legend.box.background = element_rect(fill=NA,color = NA)) #legend.box.background = element_rect(color = NA))   #为图例加边框线

p5
ggsave('D.png',plot=p4, height=12, width=15, units="in",dpi=300)
ggsave('D.pdf',plot=p4, height=12, width=15, units="in",dpi=600)

