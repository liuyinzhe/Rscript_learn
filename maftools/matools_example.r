
rm(list = ls())
options(stringsAsFactors = F)
setwd('D:\\Rscript\\maftools_example')


library(maftools)
#?annovarToMaf

# WES学习4：annovarToMaf以及maftools使用
#https://www.jianshu.com/p/1831f40df7df



var.annovar.maf = annovarToMaf(annovar = "all.annovar", 
                               Center = 'NA', 
                               refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', 
                               table = 'ensGene',
                               sep = "\t")

write.table(var.annovar.maf,file="all.annovar.maf",quote= F,sep="\t",row.names=F)



########  实际中使用 clinical.txt 对样品进行分组 var_annovar.maf 实际是多个样品合并数据

laml = read.maf(maf ='all.annovar.maf',clinicalData = 'clinical.txt')


vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Translation_Start_Site',
  'Nonstop_Mutation'
)

# vc_cols = RColorBrewer::brewer.pal(n = 11, name = 'Paired')

# names(vc_cols) = c(
#   "Silent",
#   'Multi_Hit',
#   "Frame_Shift_Del",
#   "In_Frame_Del",
#   "In_Frame_Ins",
#   "Frame_Shift_Ins",
#   "Splice_Site",
#   "Translation_Start_Site",
#   "Missense_Mutation",
#   "Nonstop_Mutation",
#   "Nonsense_Mutation"
#   # "IGR",
#   # "Splice_Region",
#   # "RNA",
#   # "5\'Flank",
#   # "5\'UTR",
#   # "Intron",
#   # "3\'UTR",
#   # "3\'Flank",
#   # "Unknown",
# )





# oncoplot
png("oncoplot.png", height=1000,width=2000, res=100)
oncoplot(maf=laml, top=20, colors = vc_cols,clinicalFeatures = 'Sample_organization',sortByAnnotation = TRUE,
         fontSize = 0.8,
         showTumorSampleBarcodes=FALSE) #borderCol=NULL #,draw_titv = TRUE
dev.off()

pdf("oncoplot.pdf")
oncoplot(maf=laml, top=20, colors = vc_cols,clinicalFeatures = 'Sample_organization',sortByAnnotation = TRUE) #borderCol=NULL #,draw_titv = TRUE
dev.off()


# somaticInteractions 绘制失败 No meaningful interactions found.
# png("somaticInteractions.png", height=1000,width=1000, res=100)
# somaticInteractions(maf = laml, top = 20, pvalue = c(0.05, 0.1))
# dev.off()

# pdf("somaticInteractions.pdf")
# somaticInteractions(maf = laml, top = 20, pvalue = c(0.05, 0.1))
# dev.off()

# plotmafSummary
png("plotmafSummary.png", height=1000,width=1000, res=100)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = "median", dashboard = TRUE,
               titvRaw = FALSE,showBarcodes = TRUE,textSize =0.6,color = vc_cols)
dev.off()

pdf("plotmafSummary.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = "median", dashboard = TRUE,
               titvRaw = FALSE,showBarcodes = TRUE, textSize =0.6,color = vc_cols)
dev.off()


if (TRUE){
# ti Tv
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary

png("plotTiTv.png", height=1000,width=1000, res=300)
plotTiTv(res = laml.titv)
dev.off()

pdf("plotTiTv.pdf", height=1000,width=1000)
plotTiTv(res = laml.titv)
dev.off()
}

## pathway

png("OncogenicPathways.png", height=1000,width=1000, res=100)
OncogenicPathways(maf = laml)
dev.off()

pdf("OncogenicPathways.pdf")
OncogenicPathways(maf = laml)
dev.off()

## PlotOncogenicPathways
png("RTK-RAS.Pathways.png", height=1000,width=1000, res=100)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS",showTumorSampleBarcodes=TRUE)
dev.off()


pdf("RTK-RAS.Pathways.pdf")
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS",showTumorSampleBarcodes=TRUE)
dev.off()

png("NOTCH.Pathways.png", height=1000,width=1000, res=100)
PlotOncogenicPathways(maf = laml, pathways = "NOTCH",showTumorSampleBarcodes=TRUE)
dev.off()

pdf("NOTCH.Pathways.pdf")
PlotOncogenicPathways(maf = laml, pathways = "NOTCH",showTumorSampleBarcodes=TRUE)
dev.off()


png("PI3K.Pathways.png", height=1000,width=1000, res=100)
PlotOncogenicPathways(maf = laml, pathways = "PI3K",showTumorSampleBarcodes=TRUE)
dev.off()

pdf("PI3K.Pathways.pdf")
PlotOncogenicPathways(maf = laml, pathways = "PI3K",showTumorSampleBarcodes=TRUE)
dev.off()


# TMB

png("TMB.png", height=1000,width=1000, res=100)
tmb(maf = laml,captureSize=60.45,logScale = TRUE)
dev.off()

pdf("TMB.pdf")
tmb(maf = laml,captureSize=60.45)
dev.off()


tmb_info = tmb(maf = laml,captureSize=60.45)
table(tmb_info)
tmb_info$Tumor_Sample_Barcode
tmb_info$total_perMB

write.table(tmb_info, file ="TMB.csv", sep=",",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# tcgaCompare
png("tcgaCompare.TMB.png", height=1000,width=1000, res=100)
tcgaCompare(maf=laml,cohortName="Tumors",capture_size=2.1,cohortFontSize = 0.7) #
dev.off()

pdf("tcgaCompare.TMB.pdf")
tcgaCompare(maf=laml,cohortName="Tumors",capture_size=2.1,cohortFontSize = 0.7) #
dev.off()

tcgaCompare_info <- tcgaCompare(maf=laml,cohortName="Tumors",capture_size=2.1,cohortFontSize = 0.7) #
print(tcgaCompare_info)
tcgaCompare_info$median_mutation_burden

write.table(tcgaCompare_info$median_mutation_burden, file ="median_mutation_burden.csv", sep=",",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


tcgaCompare_info$mutation_burden_perSample

write.table(tcgaCompare_info$mutation_burden_perSample, file ="mutation_burden_perSample.csv", sep=",",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

tcgaCompare_info$pairwise_t_test

write.table(tcgaCompare_info$pairwise_t_test, file ="pairwise_t_test.csv", sep=",",
            row.names = FALSE,col.names =TRUE, quote =TRUE)


# drugInteractions

png("drugInteractions.png", height=1000,width=1000, res=100)
drugInteractions(maf = laml, fontSize = 0.75)
dev.off()

pdf("drugInteractions.pdf")
drugInteractions(maf = laml, fontSize = 0.75)
dev.off()


# oncodrive.sig
laml.sig = oncodrive(maf = laml, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)

write.table(laml.sig, file ="oncodrive.sig.csv", sep=",",
            row.names = FALSE,col.names =TRUE, quote =TRUE)

# plotOncodrive

png("drugInteractions.png", height=1000,width=1000, res=100)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()

pdf("drugInteractions.pdf")
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()
