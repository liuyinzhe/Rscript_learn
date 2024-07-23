library(randomcoloR)
palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
palette <- distinctColorPalette(60) #差异明显的60种


colors=c('#CC92A5','#B28AAE','#B4A7E9','#5CBCE5','#4D869E','#4ED088','#57C663','#CDE86B','#93D767','#809F03')



library(RColorBrewer)
#433种
color433 = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(1,50), col=sample(color433, 50, replace = F))

#74种
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,60), col=sample(color74, 60, replace = F))

#37种
color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
pie(rep(1,37), col=sample(color37, 37))

#20种
color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
pie(rep(1,20), col=sample(color20, 20))

