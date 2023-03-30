#######################
## serosim tutorial -- Imperial 31/03/2023
## Author: James Hay
## Last updated: 29/03/2023
#######################

## Install and load serosim 
## devtools::install_github("AMenezes97/serosim")
setwd("~/Documents/GitHub/serosim_tutorial/")
source("extra_funcs.R")

## Load additional packages required 
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(serosim)

set.seed(1234)

#####################################################################
## 1. Simulation settings
## Specify the number of time steps in the simulation
## Assume 10 years at monthly time steps
times <- seq(1,120,by=1) 

## Set simulation settings argument needed for runserosim 
simulation_settings <- list("t_start"=1,"t_end"=max(times))
#####################################################################

#####################################################################
## 2. Demography
## Generate the population demography tibble
## Specify the number of individuals in the simulation; N=1000
## Assume that 20% of individuals leave the population at some time; prob_removal=0.2
N <- 250

## We'll add in sex and location here. We'll make the FOI different 
## for the two locations, but won't use sex in the model
aux_demography <- list("Sex"=list("name"="sex","options"=c("male","female","other/unspecified"),"proportion"=c(0.475,0.475,0.05)),
                       "Group"=list("name"="location","options"=c("Urban","Rural"),"proportion"=c(0.7,0.3)))
demography <- generate_pop_demography(N=N, times=times, prob_removal=0.2, removal_min=70,aux=aux_demography,age_min=24)

## Note that an entry is provided for each time point of the simulation. 
## This allows characteristics to be changed over time.

## Examine the generated demography tibble
demography %>% select(-times) %>% distinct()
#####################################################################

#####################################################################
## 3. Biomarker map
# Create biomarker map
biomarker_map <- tibble(exposure_id=c("ifxn","vacc"),biomarker_id=c("IgG","IgG"))
biomarker_map

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map)
biomarker_map
reformat_biomarker_map(biomarker_map,exposure_key=c("infection","vaccination"),biomarker_key=c("IgG"))

## Alternative example
biomarker_map_alt <- tibble(exposure_id=c("Pathogen A","Pathogen B",
                                          "Vaccination","Vaccination"),
                            biomarker_id=c("Biomarker 1","Biomarker 2",
                                           "Biomarker 2"," Biomarker 2"))
#####################################################################

#####################################################################
## 4. Exposure model
## Create an empty array to store the force of exposure for all exposure types
# Dimension 1: location
# Dimension 2: time
# Dimension 3: exposure ID
foe_pars <- array(0, dim=c(2,max(times),n_distinct(biomarker_map$exposure_id)))

## Specify the force of exposure for Location 1, Exposure ID 1 which represents natural infection
foe_pars[1,,1] <- 0.01
## Specify the force of exposure for Location 2, Exposure ID 1 which represents natural infection
foe_pars[2,,1] <- 0.02

## Specify the force of exposure for exposure ID 2 which represents vaccination
foe_pars[1,,2] <- 0.25
foe_pars[2,,2] <- 0.01

## Specify a simple exposure model which calculates the probability of exposure 
## directly from the force of exposure at that time step. In this selected model, 
## the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.
exposure_model<-exposure_model_simple_FOE

## Examine the probability of exposure over time for the specified exposure model and foe_pars array
plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times, 
                    n_groups = 2, n_exposures = 2, foe_pars=foe_pars)
#####################################################################

#####################################################################
## 5. Immunity model
## Assume that immunity from infection is dependent on latent biomarker level at time of exposure
plot_biomarker_mediated_protection(seq(0,8,by=0.1), 2, 2)

## Tibble of model parameter controls related to infection immunity
model_pars_immunity <- tibble("exposure_id"="ifxn","biomarker_id"="IgG",
                              "name"=c("biomarker_prot_midpoint","biomarker_prot_width"),
                              "mean"=c(2,2),"sd"=NA,"distribution"="")

## Vaccination will be determined by age and doses, not dependent on biomarker quantity
max_events <- c(10,2)
vacc_exposures <- 2
vacc_age <- c(NA,12)

## See example in help file
immunity_model <- immunity_model_vacc_ifxn_biomarker_prot
#####################################################################


#####################################################################
## 6. Antibody model
## Specify the antibody model 
antibody_model<-antibody_model_monophasic

## Bring in the antibody parameters needed for the antibody model
model_pars_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)
model_pars_original 

model_pars <- reformat_biomarker_map(bind_rows(model_pars_original, model_pars_immunity))

## Specify the draw_parameters function
draw_parameters<-draw_parameters_random_fx

## Plot example biomarker trajectories given the specified antibody kinetics model, 
## model parameters and draw parameters function 
plot_antibody_model(antibody_model_monophasic, N=100, model_pars=model_pars %>% drop_na(),
                    draw_parameters_fn = draw_parameters_random_fx, 
                    biomarker_map=biomarker_map)
#####################################################################

#####################################################################
## 7. Observation model
## Specify the observation model 
observation_model<-observation_model_continuous_bounded_noise
#observation_model<-observation_model_discrete_noise

## Specify assay sensitivity and specificity needed for the observation model
model_pars_original[model_pars_original$name =="obs_sd","sd"] <- 0.5
sensitivity<-0.85
specificity<-0.95

bounds <- dplyr::tibble(biomarker_id=1,name=c("lower_bound","upper_bound"),value=c(1,10))
observation_model_continuous_bounded_noise(example_biomarker_states, model_pars=model_pars,bounds,
                                           sensitivity=sensitivity,specificity = specificity) %>% 
  pull(observed) %>% hist

#observation_model_discrete_noise(example_biomarker_states, model_pars=model_pars,cutoffs=matrix(seq(0,10,by=1),nrow=1),
# sensitivity=sensitivity,specificity = specificity) %>% pull(observed) %>% as.numeric() %>% hist


## Specify observation_times (serological survey sampling design) to observe
## biomarker 1 at two time points around t 80 and 110
observation_times<- tibble(i=rep(1:max(demography$i),2),
                           t=floor(c(rnorm(N,80,3), rnorm(N,110,3))),
                           b=1) %>% mutate(t = ifelse(t > max(times), max(times), t)) %>% 
  arrange(i,t)
#####################################################################

#####################################################################
## 8. Full model
## Run the core simulation and save outputs in "res"
res<- runserosim(
  simulation_settings,
  demography,
  observation_times,
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model,
  immunity_model,
  antibody_model,
  observation_model,
  draw_parameters,
  
  ## Other arguments needed
  bounds=bounds,
  #cutoffs=matrix(seq(0,10,by=1),nrow=1),
  max_events=max_events,
  vacc_exposures=vacc_exposures,
  vacc_age=vacc_age,
  sensitivity=sensitivity,
  specificity=specificity,

  VERBOSE=100
)
#####################################################################

#####################################################################
## 9. Some post processing and plots
## Plot biomarker kinetics and exposure histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states, res$exposure_histories_long, 
                                subset=10, demography)

## Plot individual exposure histories for all exposure types
plot_exposure_histories(res$exposure_histories_long)

## Plot true biomarker quantities for all individuals across the entire simulation period
plot_biomarker_quantity(res$biomarker_states)

## Plot the serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states)

## Save data for serosolver later
write_csv(res$observed_biomarker_states %>% left_join(demography %>% select(-times) %>% distinct()),file="data/simulated_serosurvey.csv")
write_csv(res$exposure_histories_long,file="data/simulated_serosurvey_exp_histories.csv")

output <- assess_seroconversion_threshold(threshold=2, observations=res$observed_biomarker_states, true_exposure_histories=res$exposure_histories_long, demography=demography)
calculate_sens_spec(output[[2]])

output[[1]]
output[[2]]
#####################################################################
