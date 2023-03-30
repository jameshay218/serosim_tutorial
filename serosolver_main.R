######################################################
## serosolver tutorial -- testing exposure history inference on serosim data
## Author: James Hay
## Last updated: 30 March 2023
######################################################
library(ggplot2)
library(coda)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(plyr)
library(dplyr)
library(readr)
library(viridis)
library(bayesplot)
library(tidyverse)
library(ggpubr)

#devtools::install_github("seroanalytics/serosolver")
library(serosolver)

run_name <- "serosim_recovery"
main_wd <- "~/Documents/GitHub/serosim_tutorial/"
source(paste0(main_wd,"extra_funcs.R"))

chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets <- 1 ## Ignore
prior_version <- 2 ## Which prior on the infection histories to use? Prior version 2 is generally preferred
n_chains <- 5 ## Number of MCMC chains to run

rerun <- TRUE ## Set to FALSE if you just want to load in previously run chains

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can be tweaked.
## It's fairly involved tweaking these to optimize the MCMC chains, but happy to chat through
## Main one to change is iterations
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=100,
               "iterations"=100000,
               "adaptive_period"=100000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

sero_data <- read.csv("data/simulated_serosurvey.csv")
sero_data <- convert_serodata_to_serosolver(sero_data)

true_inf_hist <- read.csv("data/simulated_serosurvey_exp_histories.csv")
## Total number of exposures
sum(true_inf_hist$value,na.rm=TRUE)
true_inf_hist <- convert_inf_hist_to_serosolver(true_inf_hist)

## Set up parameter table
par_tab <- read.csv("par_tab_base.csv",stringsAsFactors=FALSE)
## These are the alpha and beta priors on the Beta distribution prior on the per-time attack rate. 
## Set to 1/1 for uniform
par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- c(2,10)

## Setup some inputs for serosolver -- this just tells the function the vector of possible
## exposure times. The antigenic map is uninformative here, but is used in other examples
## to capture cross-reactivity
strain_isolation_times <- seq(1,max(sero_data$samples),by=1)
antigenic_map <- data.frame(x_coord=0,y_coord=0,inf_times=strain_isolation_times)

## Set up posterior function for later. Data_type=2 is for continuous data
f <- create_posterior_func(par_tab,sero_data,antigenic_map=antigenic_map,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=NULL,data_type=2
                           )
## Time runs and use dopar to run multiple chains in parallel
t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
if(rerun){
res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","serosolver")) %dopar% {
    index <- 1
    ## Try generating random starting values until we find a set that returns a finite likelihood
    lik <- -Inf
    inf_hist_correct <- 1
    while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
        start_tab <- generate_start_tab(par_tab)
        start_inf <- setup_infection_histories_total(sero_data,strain_isolation_times,2,100)
        
        inf_hist_correct <- sum(check_inf_hist(sero_data, strain_isolation_times, start_inf))
        y <- f(start_tab$values, start_inf)
        lik <- sum(y[[1]])
        index <- index + 1
    }
    ## Write starting conditions
    write.csv(start_tab, paste0(x, "_start_tab.csv"))
    write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
    write.csv(sero_data, paste0(x, "_titre_dat.csv"))
    
    ## Run serosolver!
    res <- serosolver::run_MCMC(start_tab, sero_data, antigenic_map=antigenic_map, 
                                strain_isolation_times = strain_isolation_times,
                                start_inf_hist=start_inf,filename=x,
                                CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                CREATE_PRIOR_FUNC = NULL,
                                version=prior_version,
                                mcmc_pars=mcmc_pars,
                                measurement_indices=NULL, ## NULL
                                measurement_random_effects = FALSE, ## FALSE
                                solve_likelihood=TRUE,data_type=2)
}
}
run_time_fast <- Sys.time() - t1
run_time_fast

## Read in chains for trace plot
chains <- load_mcmc_chains(chain_wd,par_tab=par_tab,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
## Save traceplots from coda package
pdf(paste0(save_wd,"/",run_name,"_chain.pdf"))
plot(as.mcmc.list(chains$theta_list_chains))
dev.off()

## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 

## Plot total number of infections
p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
ggsave(paste0(save_wd,"/",run_name,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)

## Plot individual infection history estimates
p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:25,
                                             real_inf_hist = as.matrix(true_inf_hist[,strain_isolation_times]),
                                             strain_isolation_times = strain_isolation_times,nsamp=100,
                                             number_col = 5)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists_dens.pdf"),p_cumu_infs[[2]],height=8,width=7,units="in",dpi=300)


## Plot attack rates
## Get number alive in each time point and true number of infections per time point
n_alive <- get_n_alive_group(sero_data,strain_isolation_times)
n_inf <- true_inf_hist[,strain_isolation_times] %>% colSums()
true_ar <- data.frame(j=strain_isolation_times,group=1,AR=n_inf/n_alive[1,])
p_ar <- plot_attack_rates(inf_chain, sero_data, strain_isolation_times,
                          pad_chain=FALSE,plot_den=TRUE,n_alive = n_alive,
                          true_ar=true_ar,
                          prior_pars = c("prior_version"=2,"alpha"=1,"beta"=10))
p_ar
ggsave(paste0(save_wd,"/",run_name,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)


## Plot model fits to titre data. We expand the sero_data object so that the function plots
## for all possible times, not just those with observation times.
sero_data_tmp <- expand_sero_data(sero_data)

titre_pred_p <- plot_infection_histories(chain = chain[chain$chain_no == 1,], 
                         infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                         titre_dat = sero_data_tmp, 
                         individuals = 1:25,
                         strain_isolation_times = strain_isolation_times,
                         nsamp = 100, # Needs to be smaller than length of sampled chain 
                         par_tab = par_tab,p_ncol=5, data_type=2)


ggsave(paste0(save_wd,"/",run_name,"_titre_fits.pdf"),titre_pred_p,height=7,width=8,units="in",dpi=300)

