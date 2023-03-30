assess_seroconversion_threshold <- function(threshold, observations, true_exposure_histories, demography){
  
  ## How many people seroconverted between t=100 and t=120?
  seroconverted_data <- res$observed_biomarker_states %>% 
    arrange(i, t) %>%
    group_by(i) %>% mutate(t_point=1:n()) %>%
    pivot_wider(id_cols=c(i,b),names_from=t_point,values_from=observed) %>% 
    mutate(change=`2`-`1`) %>%  ## Get change in titer
    rename(t1=`1`,t2=`2`) %>%
    mutate(seroconv=change >= threshold) %>% ## Binary seroconversion
    left_join(demography_less) 
  
  ## Get true exposure state
  observation_times_wide <- observations %>% select(c(i,t)) %>% group_by(i) %>% 
    mutate(t_point=1:n()) %>% 
    mutate(t_point=paste0("t",t_point)) %>% 
    pivot_wider(id_cols=i, values_from=t, names_from=t_point)
  
  true_seroconversions <- res$exposure_histories_long %>% 
    left_join(observation_times_wide) %>%
    filter(t >= t1, t<=t2) %>% 
    group_by(i) %>% 
    summarize(true_exposures=sum(value,na.rm=TRUE))
  
  ## Combine estimated and true exposure events
  seroconverted_data <- seroconverted_data %>% left_join(true_seroconversions)
  
  ## Get number of individuals in the population during this time window, 
  ## and number of individuals who were actually sampled
  demography_summary <- demography %>% 
    filter(times <= removal, times >= birth) %>%
    #left_join(observation_times_wide)  %>% 
    filter(times >= min(observation_times_wide$t1), times <= max(observation_times_wide$t2)) %>% 
    select(i, location) %>% 
    distinct() %>% 
    group_by(location) %>% 
    tally()%>% 
    rename(N_true=n)
  
  ## Get number estimated to have seroconverted
  seroconverted_summary_obs <- seroconverted_data %>% 
    group_by(location) %>% 
    filter(!is.na(change)) %>%
    dplyr::summarize(seroconv=sum(seroconv,na.rm=TRUE), 
                     N=n()) %>%
    mutate(category="Estimate")
  
  ## Get number proportion exposed (call it seroconverted)
  seroconverted_summary_true <- seroconverted_data %>% 
    group_by(location) %>% 
    dplyr::summarize(seroconv=sum(true_exposures >= 1, na.rm=TRUE), N=n()) %>%
    mutate(category="Truth")
  
  ## Get proportions and binomial confidence intervals
  seroconverted_summary <- bind_rows(seroconverted_summary_obs,seroconverted_summary_true) %>%
    rowwise() %>%
    mutate(prop=seroconv/N, 
           lower=prop.test(x=seroconv,n=N)$conf.int[1],
           upper=prop.test(x=seroconv,n=N)$conf.int[2]) %>%
    mutate(group=paste0(location, " (", category,")"))
  
  ## Bar plot showing number sampled and number seroconverted
  p1 <- ggplot(seroconverted_summary %>% 
                 select(location, group,seroconv,N) %>% 
                 pivot_longer(-c(location,group))) + 
    geom_bar(aes(y=value,x=group,fill=name),stat="identity") +
    ylab("Number seroconverted") +
    xlab("Location") +
    facet_wrap(~location,nrow=1, scales="free_x") +
    ggtitle(paste0("Comparison of number and proportion seroconverted\n assuming seroconversion threshold=",threshold))
  
  ## Pointrange plot showing proportions and binomial confidence intervals
  p2 <- ggplot(seroconverted_summary) + 
    geom_pointrange(aes(y=prop,ymin=lower,ymax=upper,x=group,col=category),stat="identity") +
    ylab("Proportion seroconverted") +
    xlab("Location") +
    facet_wrap(~location,nrow=1, scales="free_x") +
    scale_y_continuous(limits=c(0,0.75))
  
  p_main <- p1/p2
  
  accuracy_table1 <- tibble(categorization = c("False negative","True negative","True positive","False positive"),n=0)
  
  accuracy_table <- seroconverted_data %>% 
    mutate(true_exposures = true_exposures > 0) %>%
    mutate(categorization=ifelse(true_exposures==TRUE & seroconv==TRUE,"True positive",
                                 ifelse(true_exposures==TRUE & seroconv==FALSE,"False negative",
                                        ifelse(true_exposures==FALSE & seroconv==TRUE,"False positive",
                                               "True negative")))) %>%
    group_by(categorization) %>% tally()
  accuracy_table <- bind_rows(accuracy_table, accuracy_table1)
  accuracy_table <- accuracy_table %>% drop_na() %>% group_by(categorization) %>% filter(n == max(n)) %>% ungroup()
  
  return(list(p_main, accuracy_table, seroconverted_data))
}

calculate_sens_spec <- function(output){
  output <- output %>% drop_na()
  TP <- output[output$categorization == "True positive","n"]
  FP <- output[output$categorization == "False positive","n"]
  TN <- output[output$categorization == "True negative","n"]
  FN <- output[output$categorization == "False negative","n"]
  sens <- unname(TP/(TP + FN))
  spec <- unname(TN/(TN + FP))
  return(c("sensitivity"=sens,"specificity"=spec))
}
