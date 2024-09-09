options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("BiocManager")


############
library('BiocManager',quietly=TRUE)
BiocManager::install(c('ggplot2','devtools','ggpubr','argparse'))
BiocManager::install(c('clustree', 'harmony', 'R.utils', 'tidyverse'))
# https://satijalab.org/seurat/articles/install.html

remotes::install_github("satijalab/seurat", "seurat5", quiet = T)




remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)



remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("bnprks/BPCells/r")

# token xxxxxxxxxxxxxxxxxxx
# To increase your GitHub API rate limit
#- Use `usethis::create_github_token()` to create a Personal Access Token.
#- Use `gitcreds::gitcreds_set()` to add the token.
# usethis::create_github_token() 访问PAT 获取网页

# devtool Failed to install  unknown package from GitHub: HTTP error 401.
# https://blog.csdn.net/weixin_57360787/article/details/139073650
# usethis::create_github_token()
Sys.setenv(GITHUB_PAT = "token")
Sys.getenv("GITHUB_PAT")

#remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
update.packages(oldPkgs = c("withr", "rlang"))
remotes::install_github('satijalab/azimuth', ref = 'master')

#remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
