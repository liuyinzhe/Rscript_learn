rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
setwd('D:\\Rscripts\\bar')



# 参考网址
# ggplot2:正负区分条形图及美化
# https://cloud.tencent.com/developer/article/1089081
# R语言绘制柱形图
# https://mp.weixin.qq.com/s?src=11&timestamp=1700917903&ver=4918&signature=LD*-ZJIZVQgaNRPQu8VebC*1rCACFkcTSDoMRMf5ck3ecM*HrxF9fYsfP7Yp39PX9eGS*HchaTme8NnroceGnUngYwbwBlKU6HIAzURdKw0lhCcdvvnk6jtFPT2prhZu&new=1
# Excel 实现轴标签正负数自动换位图（条形图&柱状图优化）
# https://zhuanlan.zhihu.com/p/337644712

dat <- read.table('lefse.res.plot.tsv',header = TRUE)


p <- ggplot(dat,aes(x = reorder(biomarker,LDA_score),y = LDA_score,fill = groups)) +
    geom_bar(stat = "identity",colour = "black", width = 0.78, position = position_dodge(0.7)) #position_dodge(width=0.9) 默认柱子宽度
# 使用reorder()函数，使得按照LDA_score取值大小，条形图进行排序

p <- p + 
    xlab("") +
    ylab(expression(LDA~score~(log~{10})))+ #log[10]^{x}
    #ylab("LDA score") + #设置 x,y 轴标签
    #scale_y_continuous(breaks=seq(-9,9,1)) +
    coord_flip() # 翻转坐标

p
# 美化图形
# 去除灰色背景
# 修改字体

p <- p +theme_bw()+
    theme(panel.grid.major.y = element_blank(),panel.grid.minor.y=element_blank())+ # 删除纵向网络(主网格线y轴，次网格线y轴)
    theme(axis.text.y = element_blank()) + # 删除刻度线(y轴)
    theme(axis.ticks = element_blank()) + # 删除刻度线
    theme(panel.border = element_blank()) + # 去除边框
    theme(legend.title = element_blank()) + # 去掉图例 标题
    #theme(axis.text.x = element_text(face = "bold")) + # x轴刻度文字加粗
    theme(axis.title.x = element_text(face = "bold")) + # x坐标轴标题文字加粗
    theme(legend.position="top")  #图例 位置放在上面
    # labs(fill = "Dose (mg)")  + #指定图例 标题内容
    # scale_fill_discrete(limits = c("SN", "QM")) # 图例顺序指定

#ggplot2图例修改详细介绍
#https://zhuanlan.zhihu.com/p/548703892?utm_id=0

# 修改图例的颜色、字体，并在条形图两侧加入标签文字
p

p <- p+geom_text(data=dat,
    aes(y=ifelse(LDA_score>0,-0.5,0.5),label = biomarker), # 标签放的坐标y轴固定位置
    fontface = 4, size=3.5, hjust=ifelse(dat$LDA_score>0,1,0)) + # hjust和vjust的值仅在0和1之间定义,0意味着左alignment,1表示右alignment
    # ifelse
    # > a <- c(5, 7, 2, 9)
    # > ifelse(a %% 2 == 0, "even", "odd")
    # [1] "odd"  "odd"  "even" "odd"
    
    # 设置图片中标签字体,位置等
    # ggplot绘图 003 美化
    # https://zhuanlan.zhihu.com/p/531364627
    # scale_fill_gradient2(low,mid,high,midpoint) # 2极渐变色
    #scale_fill_gradient2(low = "#36648B", high = "red",mid = "white", midpoint = 0, name="Difference\n") +
    # 设置图例颜色
    theme(legend.text = element_text(face="bold")) + # #,legend.title = element_text(face="bold"))
    theme(plot.margin = margin(l = 20, r = 20, t= 20, b= 20,unit = "pt")) # 绘图边距,输出过大画图过小时用
    # plot.margin	
    # margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
p


# R语言ggplot2 | 给图片添加上下标
# https://blog.csdn.net/qq_42830713/article/details/130298640

#   labs(x=expression(hwy~m^2~CO[2]^{-1}), #~表示空格
#       y=expression(displ~g~kg^{-1}))

#   annotate(geom ="text",x =4,y =7,
#            parse = TRUE,size=5,
#            label = "This^2~is~box[plot]",size =5)


ggsave(p,file="lefse_PN_bar.png",width=16,height=8,units='in',dpi=600) #
ggsave(p,file="lefse_PN_bar.pdf",width=16,height=8,units='in',dpi=300)
