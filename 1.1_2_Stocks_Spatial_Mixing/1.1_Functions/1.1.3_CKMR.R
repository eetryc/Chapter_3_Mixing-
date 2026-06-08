#==================================##
#========= CKMR function ========####
#==================================##


library(tidyverse)
library(plyr)
library(dplyr)
library(lubridate)




## 1. Define function =====================================================####

### JAGS model  ===========================================================####

cat("model{


  mu ~ dunif(1,100000)

  

  sd ~ dunif(1,100000)

  

  

  # Survival 

  

  surv ~ dbeta(1, 1)



  # Adult numbers

  

  for(y in 1:years) {

      Nadult[y] ~ dnorm(mu, 1/(sd^2)) T(100, 1000000)

  }

    

  for(i in 1:nobs) {

    TruePairs[i] ~ dbinom((RObase[i] * (surv^AgeDif[i])) / Nadult[year[i]], Pair_viable_count[i])

  }





}", file = "HSPPOP.simplest.sim.jags")







CKMR_est <- function(samples){
  
  ### POPs filter 
  
  POpairs <- samples %>%
    
    rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns fro each indiv
    
    cross_join(samples %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
    
    # mutate to add the POP possibility (y/n) with the year and reproductive age filters
    
    mutate(HSPPOP_candidate = ifelse(cohort_1 < cohort_2 & 
                                       
                                       (cohort_2 - cohort_1 >= 5) & 
                                       
                                       (year_sampled_1 >= cohort_2),"1","0")) %>% 
    
    #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
    
    filter(id_1 != id_2) %>% 
    
    mutate(PairID = paste(id_1, id_2, sep = "_")) %>%
    
    mutate(older_first=ifelse(cohort_1 < cohort_2 |(cohort_1 == cohort_2 & id_1 < id_2), T,F) ) %>% 
    
    filter(older_first) %>%
    
    select(-older_first) %>% 
    
    mutate(AgeDif = cohort_2 - cohort_1)
  
  
  
  
  
  
  
  # count numnber of yes or no to get a better idea of how they are changing
  
  POpairs %>% 
    
    dplyr::count(HSPPOP_candidate)
  
  
  
  
  
  POPonly <- POpairs %>% 
    
    filter(HSPPOP_candidate == 1) %>% 
    
    mutate(RObase=ifelse(sex_1 == "U",2,1)) %>% 
    
    mutate(Age_dif_exp = 0)
  
  
  
  
  
  
  
  
  
  
  
  ### HSPs 
  
  HSpairs <- samples %>%
    
    rename_with(~ paste0(.x, "_1"), everything()) %>% # sort into col for each indiv
    
    cross_join(samples %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
    
    mutate(HSPPOP_candidate = ifelse(cohort_1 < cohort_2 & (cohort_1 != cohort_2),"1","0")) %>% 
    
    #get rid of the ones that are the same comp, eg 1 and 2 vs 2 and 1 
    
    filter(id_1 != id_2) %>% 
    
    mutate(PairID = paste(id_1, id_2, sep = "_")) %>%
    
    mutate(older_first=ifelse(cohort_1 < cohort_2 |(cohort_1 == cohort_2 & id_1 < id_2), T,F) ) %>% 
    
    filter(older_first) %>%
    
    select(-older_first) %>% 
    
    mutate(AgeDif = cohort_2 - cohort_1)
  
  
  
  
  
  
  
  
  
  # count number of yes or no to get a better idea of how they are changing
  
  HSpairs %>% 
    
    dplyr::count(HSPPOP_candidate)
  
  
  
  
  
  HSPonly  <- HSpairs %>% 
    
    filter(HSPPOP_candidate == 1) %>% 
    
    mutate(RObase= 4) %>% 
    
    mutate(Age_dif_exp = AgeDif)
  
  
  
  
  
  
  
  
  
  
  
  
  
  # join the two tables
  
  possible_pairs <- POPonly  %>% 
    
    rbind(data.frame=HSPonly) %>% 
    
    mutate(
      
      sex_num = case_when(
        
        sex_1 == "U" ~ NA_real_, # unknown to na
        
        sex_1 == "0" ~ 1, # male is 1
        
        sex_1 == "1" ~ 0 # female is 0
        
      )
      
    )
  
  
  
  possible_pairs %>% 
    
    dplyr::count(HSPPOP_candidate)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### Data to define for JAGS model ####
  
  
  
  
  
  # years being estimated (years that actually have potential)
  
  Cohort_years <- possible_pairs %>% 
    
    group_by(cohort_2) %>% 
    
    dplyr::summarise(
      
      Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
      
    ) %>% 
    
    filter(Pair_viable_count > 0) %>%     # keep only cohorts with >0
    
    arrange(cohort_2)
  
  
  
  # extract years (count to give to be i in JAGS)
  
  years <- nrow(Cohort_years)
  
  
  
  
  
  # group together by cohort year and age dif (differing probabilities)
  
  Cohort_years_cohort_age <- possible_pairs %>%
    
    group_by(cohort_2, AgeDif) %>%
    
    dplyr::summarise(
      
      Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE))
  
  
  
  
  
  # add the counts to each observation
  
  observations <- possible_pairs %>%
    
    left_join(Cohort_years_cohort_age, by = c("cohort_2", "AgeDif"))
  
  
  
  # year index to link kinship to years (just the number to match it, not the actual year)
  
  year_index <- factor(observations$cohort_2, levels = Cohort_years$cohort_2)
  
  observations <- observations %>%
    
    mutate(
      
      year_index = as.integer(year_index)
      
    )
  
  
  
  
  
  
  
  
  
  # Count the ones with known sex (U = 0)
  
  # Count females (F = 1)
  
  observations <- observations %>%
    
    mutate(
      
      knownSex = as.numeric(!is.na(sex_num)),
      
      IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0))  # if sex_num encoded 0 = female
      
    )
  
  
  
  
  
  
  
  
  
  # Match true pairs
  
  observations <- observations %>% 
    
    mutate(true_HSP = if_else(
      
      mother_1 == mother_2 | father_1 == father_2, 1, 0),
      
      true_POP = if_else(
        
        id_1 == mother_2 | id_1 == father_2, 1, 0),
      
      TruePairs = if_else(
        
        true_POP == 1 | true_HSP == 1, 1, 0))
  
  
  
  
  
  
  
  #### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
  
  
  
  group_vars <- c("cohort_2", "knownSex", "IsFemale", "RObase", "Age_dif_exp")
  
  
  
  
  
  collapsed <- observations %>% 
    
    drop_na(mother_1, father_1) %>% 
    
    group_by(year_index,across(all_of(group_vars))) %>%
    
    dplyr::summarise(
      
      Pair_viable_count_per_agedif = n(),                   
      
      TruePairs = sum(true_POP, na.rm=TRUE),
      
      .groups = "drop"
      
    ) %>%
    
    arrange(cohort_2) %>% 
    
    filter(TruePairs != 0)
  
  
  
  
  
  
  
  collapsed %>% 
    
    group_by(year_index) %>% 
    
    summarise(TotalPairs = sum(TruePairs, na.rm = TRUE)) %>% 
    
    ggplot(aes(x = year_index, y = TotalPairs)) +
    
    geom_line() 
  
  
  
  #### Just POP model
  
  data = list( years = years,  # number of cohorts,
               
               TruePairs = collapsed$TruePairs,
               
               RObase = collapsed$RObase,
               
               AgeDif = collapsed$Age_dif_exp,
               
               nobs = nrow(collapsed),
               
               Pair_viable_count = collapsed$Pair_viable_count_per_agedif,
               
               year=collapsed$year_index)
  
  
  
  # Initial values
  
  inits = function() {
    
    init_Nadult <- rnorm(years, mean = 1000, sd = 700)
    
    init_Nadult <- pmax(pmin(init_Nadult, 900000), 500)
    
    list(
      
      surv = runif(1,0.01, 0.99),
      
      Nadult = init_Nadult)
    
    # Nadult = pmin(pmax(rnorm(years, mean = 5000, sd = 2000), 0), 10000))
    
  }
  
  
  
  
  
  # Parameters to follow
  
  params = c("Nadult","surv")
  
  
  
  nburn <- 10000
  
  nchains <- 3
  
  niter <- 20000
  
  n.cores = 3
  
  
  
  
  
  Out_CKMR_sim1 = jagsUI::jags(data, inits, params, "HSPPOP.simplest.sim.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)
  
  
  
  
  
  
  
  #### Plot against true populations and compare
  
  # Extract posterior summary as a data frame
  
  nhat <- as.data.frame(Out_CKMR_sim1$summary)
  
  
  
  # Keep only the Nadult rows
  nhat <- nhat[grep("^Nadult", rownames(nhat)), ]
  
  nhat$Actual_Year <- Cohort_years$cohort_2
  
  # Select the useful columns
  nhat_df <- nhat %>%
    dplyr::select(Actual_Year, mean, `2.5%`, `50%`, `97.5%`,Rhat)
  
  return(nhat_df) 
}
  

