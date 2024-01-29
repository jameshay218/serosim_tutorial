######################################################
## serosolver tutorial -- testing exposure history inference on serosim data
## Author: James Hay
## Last updated: 29/01/2024
######################################################
#devtools::install_github("seroanalytics/serosolver")
library(serosolver)
library(serosim)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(reshape2)
library(ggpubr)

set.seed(1234)

## Directly/file management
run_name <- "serosim_recovery"
main_wd <- "~/Documents/GitHub/serosim_tutorial/"
source(paste0(main_wd,"extra_funcs.R"))

chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))


## Read in simulated data from serosim
sero_data <- read.csv("data/simulated_serosurvey.csv")

## NOTE -- serosolver needs continuous individual IDs. This function re-indexes the individual IDs to be continuous from 1:N
## after removing individuals with no observed data
sero_data <- convert_serodata_to_serosolver(sero_data%>% drop_na()) 

## Simulated infection histories
true_inf_hist <- read.csv("data/simulated_serosurvey_exp_histories.csv")
## Total number of exposures
sum(true_inf_hist$value,na.rm=TRUE)
true_inf_hist <- convert_inf_hist_to_serosolver(true_inf_hist,sero_data=sero_data)

## Setup some inputs for serosolver -- this just tells the function the vector of possible
## exposure times. The antigenic map is uninformative here, but is used in other examples
## to capture cross-reactivity
possible_exposure_times <- seq(1,max(sero_data$sample_time),by=1)
biomarker_map <- data.frame(x_coord=0,y_coord=0,inf_times=possible_exposure_times)

## Look at serosim data before giving to serosolver
plot_antibody_data(sero_data, possible_exposure_times=seq(1,12,by=1),
                   infection_histories = true_inf_hist,
                   n_indivs=1:25,
                   study_design = "single-antigen") +
  facet_wrap(~individual)

## Set up parameter table
par_tab <- read.csv("par_tab.csv",stringsAsFactors=FALSE)

## These are the shape1 and shape2 priors on the Beta distribution prior on the per-time attack rate. 
## Set to 1/1 for uniform
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),"values"] <- c(2,10)

## Prior on waning rate parameter
prior_func <- function(par_tab){
  f <- function(pars){
    pr <- dlnorm(pars["wane_short"],log(0.003),1,log=TRUE)
    pr
  }
}

## Run serosolver and look at MCMC diagnostics
output <- serosolver::serosolver(par_tab, sero_data, biomarker_map,prior_func=prior_func,
                                 filename=paste0(chain_wd,"/",run_name), prior_version=2,n_chains=5,parallel=TRUE,
                                 mcmc_pars=c(adaptive_iterations=100000, iterations=100000),
                                 verbose=TRUE,data_type=2)
output

# Plot model predicted titres for a subset of individuals
chains <- load_mcmc_chains(location=chain_wd,par_tab=par_tab,burnin = 100000,unfixed=FALSE)
p_fits <- plot_model_fits(chain = chains$theta_chain,
                infection_histories = chains$inf_chain,
                known_infection_history = true_inf_hist,
                antibody_data = sero_data,individuals=1:25,
                antigenic_map=biomarker_map,
                par_tab=par_tab,expand_to_all_times=TRUE,
                p_ncol=5,
                orientation="longitudinal",data_type=2)
ggsave(paste0(save_wd,"/",run_name,"_fits.pdf"),p_fits,height=7,width=8,units="in",dpi=300)

## Plot estimated attack rates against truth
n_alive <- get_n_alive_group(sero_data,possible_exposure_times)
n_inf <- true_inf_hist[,possible_exposure_times] %>% colSums()
true_ar <- data.frame(j=possible_exposure_times,group=1,AR=n_inf/n_alive[1,])
p_ar <- plot_attack_rates_pointrange(infection_histories = chains$inf_chain,antibody_data = sero_data,
                             possible_exposure_times=possible_exposure_times,true_ar=true_ar,plot_den=TRUE,
                             prior_pars = c("prior_version"=2,"infection_model_prior_shape1"=2,"infection_model_prior_shape2"=10))
ggsave(paste0(save_wd,"/",run_name,"_attack_rates.pdf"),p_ar,height=4,width=7,units="in",dpi=300)

## Read in chains for other plots
chains <- load_mcmc_chains(chain_wd,par_tab=par_tab,convert_mcmc=FALSE,burnin = 100000,unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 

## Plot individual infection history estimates
p_cumu_infs <- plot_cumulative_infection_histories(inf_chain,indivs=1:25,
                                             real_inf_hist = as.matrix(true_inf_hist[,possible_exposure_times]),
                                             possible_exposure_times = possible_exposure_times,nsamp=100,
                                             number_col = 5)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists_dens.pdf"),p_cumu_infs[[2]],height=8,width=7,units="in",dpi=300)
