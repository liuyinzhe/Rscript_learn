
library(data.table)
library(TwoSampleMR)

setwd("D:\\Rscript\\TwoSampleMR\\example")

#https://github.com/MRCIEU/TwoSampleMR/issues/495
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

pFilter=1e-05       #pvalue过滤条件  # 1*10-5  0.00001
# https://mibiogen.gcc.rug.nl/menu/main/home
inputFile="MBG.allHits.p1e4.txt"     #输入文件

Token="Token_xxxx"
# Token 注册获取
# https://api.opengwas.io
#ao<-available_outcomes(Token)
# https://cloufield.github.io/GWASTutorial/TwoSampleMR/
#write.table(ao,file = "./available_outcomes.txt", sep = "\t", quote = F, row.names = F)


#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)


#保留属水平的肠道菌
#grep()和grepl()是用来查找目标字符串或字符串向量中是否包含目标子串。grep()返回的是整数索引，grepl()返回的是布尔值索引。
rt=subset(rt, grepl("genus", bac))
rt=subset(rt, !grepl("unknowngenus", bac))
rt$bac=gsub("genus\\.+(.*?)\\.id.*", "\\1", rt$bac)

#根据pvalue<1e-05对SNP进行过滤
exposure_data=subset(rt, P.weightedSumZ<pFilter)
write.csv(exposure_data, file="exposure_data.csv", row.names=F)
# exposure_data$exposure 未暴露名

#关联性分析(p<1e-05)
#去除连锁不平衡的SNP
# clump_data(exposure_dat,clump_kb=10000, clump_r2=0.001)

#读取暴露数据
exposure_data=read_exposure_data(filename="exposure_data.csv",
                                 sep = ",",
                                 snp_col = "rsID",
                                 beta_col = "beta",
                                 se_col = "SE",
                                 pval_col = "P.weightedSumZ",
                                 effect_allele_col="eff.allele",
                                 other_allele_col = "ref.allele",
                                 phenotype_col = "bac",
                                 samplesize_col = "N",
                                 chr_col="chr", pos_col = "bp",
                                 clump=FALSE)

#去除连锁不平衡的SNP
exposure_dat_clumped=clump_data(exposure_data, clump_kb=10000, clump_r2=0.001,pop="EAS")
#exposure_dat_clumped=clump_data(exposure_data, clump_kb=500, clump_r2=0.01)
write.csv(exposure_dat_clumped, file="exposure_clumped.csv", row.names=F)



# bmi_exp <- extract_instruments(
#   outcomes='finn-b-K11_CONSTIPATION',
#   p1=0.00005,  # Significance threshold. The default is 5e-8
#   r2=0.001,
#   kb=5000,  # #10,000 kb window
#   opengwas_jwt=Token
# )
# write.table(bmi_exp,file = "./bmi_exp.txt", sep = "\t", quote = F, row.names = F)
# dim(bmi_exp)
# # [1] 80 15

# extract outcome data
t2d_out <- extract_outcome_data(
  snps=exposure_dat_clumped$SNP,
  outcomes='finn-b-K11_CONSTIPATION',
  proxies = FALSE,
  rsq = 0.001, #Minimum LD rsq(r square) value 
  opengwas_jwt=Token
)
write.table(t2d_out,file = "./t2d_out.txt", sep = "\t", quote = F, row.names = F)
dim(t2d_out)
# [1] 80 16


# SNPs with P value < 1 × 10− 5 and strong correlation to exposure were selected
# SNPs with linkage disequilibrium (r2 < 0.001, within a 10,000 kb window) were removed by using ‘extract outcome data’ function of ‘TwoSample MR’.
# r-square

mydata <- harmonise_data(
  exposure_dat=exposure_dat_clumped,
  outcome_dat=t2d_out
)

# MR 分析
result <- mr(mydata, method_list=c("mr_ivw"))#c("mr_egger_regression", "mr_ivw"))
write.table(result,file = "./res_table.txt", sep = "\t", quote = F, row.names = F)


# 异质性检验
mr_heterogeneity_table <- mr_heterogeneity(mydata)
write.table(mr_heterogeneity_table,file = "./mr_heterogeneity.txt", sep = "\t", quote = F, row.names = F)


# 水平多效性检验
mr_pleiotropy_table <- mr_pleiotropy_test(mydata)
write.table(mr_pleiotropy_table,file = "./mr_pleiotropy.txt", sep = "\t", quote = F, row.names = F)

# 逐个剔除检验
mr_leaveoneout_table<- mr_leaveoneout(mydata) 
write.table(mr_leaveoneout_table,file = "./mr_leaveoneout.txt", sep = "\t", quote = F, row.names = F)


# OR值
odds_ratios_table <- generate_odds_ratios(result)
write.table(odds_ratios_table,file = "./odds_ratios.txt", sep = "\t", quote = F, row.names = F)


# 散点图
p1 <- mr_scatter_plot(result, mydata)

png("mr_scatter.png",height =600,width =900,res=100,units = "px")
print(p1)
dev.off()


# 森林图
result_single <- mr_singlesnp(mydata)
p2 <- mr_forest_plot(result_single)

png("mr_forest.png",height =600,width =900,res=100,units = "px")
print(p2)
dev.off()

# 漏斗图
#result_single <- mr_singlesnp(mydata)
p3 <- mr_funnel_plot(result_single)

png("mr_funnel.png",height =600,width =900,res=100,units = "px")
print(p3)
dev.off()

# 留一图
result_loo <- mr_leaveoneout(mydata)
p4 <- mr_leaveoneout_plot(result_loo)

png("mr_leaveoneout.png",height =600,width =900,res=100,units = "px")
print(p4)
dev.off()




# TwoSampleMR：两样本孟德尔随机化
# https://www.jianshu.com/p/58c7d8541c84  # 分析过程
# 有相关性就有因果关系吗，教你玩转孟德尔随机化分析（mendelian randomization ）
# https://www.jianshu.com/p/253309a571aa   
# 孟德尔随机化 Mendelian randomization R语言实现
# https://www.bilibili.com/read/cv19152108/  # 绘图
# TwoSampleMR实战教程之计算并解读MR结果
# https://zhuanlan.zhihu.com/p/343214538

# TwoSampleMR：两样本孟德尔随机化
# https://www.jianshu.com/p/58c7d8541c84

# 孟德尔随机化（三）—— 再也不用担心网络或其他各种报错啦 | 从数据库下载数据到本地的数据处理方法
# https://blog.csdn.net/weixin_43843918/article/details/137961902

##################
# 孟德尔随机化分析生信套路
# https://www.cnblogs.com/lzryhc/category/2359045.html
# 02暴露数据过滤
# https://www.cnblogs.com/lzryhc/articles/17871031.html  # 数据读取
# https://www.cnblogs.com/lzryhc/category/2363576.html # 02暴露数据过滤 所在目录
##################
