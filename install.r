# 另外可参考https://www.bilibili.com/read/cv21650670/
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
BiocManager::install('ggplot2')
BiocManager::install('karyoploteR')
BiocManager::install("AnnotationForge",force = TRUE)
BiocManager::install(c("Hmisc","psych","igraph"))

# suppressWarnings
local({r <- getOption("repos")  
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 

# 依赖包列表安装与加载：
package_list <- c("ggplot2","RColorBrewer","randomForest","caret", "pROC","dplyr","ggrepel")

# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 指定网址
packageurl = 'https://cran.r-project.org/src/contrib/vctrs_0.5.2.tar.gz'
install.packages(packageurl, repos = NULL, type = 'source')

# 指定安装路径

# 本地安装
BiocManager::install("devtools")
library(devtools)
devtools::install_local("C:\\Downloads\\Compressed\\GenomeInfoDb_1.34.9.zip")
#install.packages('GenomeInfoDbData')

# 二进制安装，针对源码安装非0返回值问题
install.packages("pkgload",type="binary")
