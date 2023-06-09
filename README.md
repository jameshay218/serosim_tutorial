# serosim_tutorial
 Scripts and additional files to run the serosim and serosolver tutorials

You will need to follow the following installation instructions. We're aiming to have `serosim` and `serosolver` installed. More information on these can be found on the [serosim](https://github.com/AMenezes97/serosim/) and [serosolver](https://github.com/seroanalytics/serosolver) webpages.

```r
## Note that you will need a Cpp tool chain through RTools. For example, my mac uses the clang compiler.

used_packages <- c("plyr","reshape2","data.table","foreach","doParallel","parallel",
                   "viridis","bayesplot","tidyverse","ggpubr","patchwork","Rcpp",
                   "coda","Matrix","GGally","MASS","deSolve","RcppArmadillo")

needed_packages <- used_packages[!(used_packages %in% installed.packages()[,"Package"])]
if(length(needed_packages)) install.packages(needed_packages)

devtools::install_github("AMenezes97/serosim")
devtools::install_github("seroanalytics/serosolver")
```