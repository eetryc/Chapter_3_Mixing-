#=========================================##
#========= CKMR median function ========####
#=========================================##


library(tidyverse)
library(plyr)
library(dplyr)
library(lubridate)



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

  for(y in 1:years) {
    w_Nadult[y] <- Nadult[y] * weights[y]
  }
  weighted_Nadult <- sum(w_Nadult[1:years])

}", file = "HSPPOP.simplest.sim.jags")


















CKMR_est.med <- function(samples){
  
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
    rbind(HSPonly) %>%        # Fixed: data.frame= was removed here as well to avoid empty list coercion
    mutate(
      sex_num = case_when(
        sex_1 == "U" ~ NA_real_, 
        sex_1 == "0" ~ 1, 
        sex_1 == "1" ~ 0 
      )
    )
  
  #### Data to define for JAGS model ####
  
  # years being estimated (years that actually have potential)
  Cohort_years <- possible_pairs %>% 
    group_by(cohort_2) %>% 
    dplyr::summarise(
      Pair_viable_count = sum(HSPPOP_candidate > 0, na.rm = TRUE)
    ) %>% 
    filter(Pair_viable_count > 0) %>%     
    arrange(cohort_2)
  
  # extract years (count to give to be i in JAGS)
  years <- nrow(Cohort_years)
  
  # group together by cohort year and age dif (differing probabilities)
  Cohort_years_cohort_age <- possible_pairs %>%
    group_by(cohort_2, AgeDif) %>%
    dplyr::summarise(
      Pair_viable_count_per_agedif = sum(HSPPOP_candidate > 0, na.rm = TRUE),
      .groups = "drop" # Prevent grouping warnings
    )
  
  # add the counts to each observation
  observations <- possible_pairs %>%
    left_join(Cohort_years_cohort_age, by = c("cohort_2", "AgeDif"))
  
  # year index to link kinship to years
  year_index_factor <- factor(observations$cohort_2, levels = Cohort_years$cohort_2)
  
  observations <- observations %>%
    mutate(
      year_index = as.integer(year_index_factor),
      knownSex = as.numeric(!is.na(sex_num)),
      IsFemale = as.numeric(!is.na(sex_num) & (sex_num == 0)) 
    )
  
  # Match true pairs
  observations <- observations %>% 
    mutate(
      true_HSP = if_else(mother_1 == mother_2 | father_1 == father_2, 1, 0),
      true_POP = if_else(id_1 == mother_2 | id_1 == father_2, 1, 0),
      TruePairs = if_else(true_POP == 1 | true_HSP == 1, 1, 0)
    )
  
  #### Now collapse based on cohort 2, age difference, RObase, if sex is known, and female or male 
  group_vars <- c("cohort_2", "knownSex", "IsFemale", "RObase", "Age_dif_exp")
  
  collapsed <- observations %>% 
    drop_na(mother_1, father_1) %>% 
    group_by(year_index, across(all_of(group_vars))) %>%
    dplyr::summarise(
      Pair_viable_count_per_agedif = n(),                   
      TruePairs = sum(TruePairs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cohort_2) #%>% 
    # filter(TruePairs != 0)
  
  
  
  
  # =========================================================================
  # CALCULATE WEIGHTS TOWARDS A TARGET YEAR
  # =========================================================================
  target_year <- Cohort_years$cohort_2[which.max(Cohort_years$Pair_viable_count)]     
  decay_rate  <- 7       
  
  actual_years <- Cohort_years$cohort_2
  
  # Calculate Gaussian weights: e^(-(year - target)^2 / (2 * decay_rate^2))
  raw_weights <- exp(-((actual_years - target_year)^2) / (2 * decay_rate^2))
  
  # Standardize weights so they sum to 1
  standardized_weights <- raw_weights / sum(raw_weights)
  # =========================================================================
  
  
  #### Filter out the years where there are not enough pairs to make real estimates
  
  
  true_kin_ratios <- collapsed %>% 
    mutate(ratio = TruePairs/Pair_viable_count_per_agedif)
  
  true_kin_ratios_sum <- true_kin_ratios %>%
    group_by(cohort_2) %>%
    dplyr::summarise(
      TotalTruePairs = sum(TruePairs),
      PercentTrue = sum(ratio),
      .groups = "drop"
    )
  
  drop_years <- true_kin_ratios_sum %>% 
    filter(TotalTruePairs != 0)
  
  collapsed <- collapsed %>% 
    semi_join(drop_years, by = "cohort_2")
  
  collapsed <- collapsed %>%
    mutate(year_index = as.integer(as.factor(cohort_2)))
  
  years <- nrow(drop_years)
  
  
  
  
  
  
  
  Cohort_years_filtered <- Cohort_years %>% 
    semi_join(drop_years, by = "cohort_2")
  
  actual_years <- Cohort_years_filtered$cohort_2
  
  
  
  
  raw_weights <- exp(-((actual_years - target_year)^2) / (2 * decay_rate^2))
  
  
  standardized_weights <- raw_weights / sum(raw_weights)
  
  #### Just POP model
  data = list( years = years,  # number of cohorts,
               TruePairs = collapsed$TruePairs,
               RObase = collapsed$RObase,
               AgeDif = collapsed$Age_dif_exp,
               nobs = nrow(collapsed),
               Pair_viable_count = collapsed$Pair_viable_count_per_agedif,
               year = collapsed$year_index,
               weights = standardized_weights) # <-- NEW: Pass weights to JAGS
  
  # Initial values
  inits = function() {
    
    list(
      
      surv = runif(1,0.01, 0.99),
      
      Nadult = runif(data$years,100, 100000))
    

  }
  
  # Parameters to follow
  params = c("Nadult","surv", "weighted_Nadult") # <-- NEW: Track the weighted summary
  
  nburn <- 10000
  nchains <- 3
  niter <- 20000
  n.cores = 3
  
  Out_CKMR_sim1 = jagsUI::jags(data, inits, params, "HSPPOP.simplest.sim.jags", 
                               n.burnin = nburn, n.chains = nchains, n.iter = niter, 
                               parallel = T, verbose = T)
  
  #### Plot against true populations and compare
  # Extract posterior summary as a data frame
  nhat_all <- as.data.frame(Out_CKMR_sim1$summary)  
  
  # Keep both individual Nadult rows AND the new weighted summary row
  nhat_years <- nhat_all[grep("^Nadult\\[", rownames(nhat_all)), ]
  nhat_years$Actual_Year <- as.character(Cohort_years_filtered$cohort_2)
  
  weighted_row_idx <- grep("^weighted_Nadult", rownames(nhat_all))
  
  if(length(weighted_row_idx) == 0) {
    stop("Could not locate 'weighted_Nadult' in the JAGS output summary. Check parameter names.")
  }
  
  nhat_weighted <- nhat_all[weighted_row_idx, ]
  nhat_weighted$Actual_Year <- paste0("Weighted_Focal_", target_year)
  
  nhat_combined <- rbind(nhat_years, nhat_weighted)
  
  nhat_df <- nhat_combined %>%
    dplyr::select(Actual_Year, mean, `2.5%`, `50%`, `97.5%`, Rhat)
  
  return(list(nhat = nhat_df, collapsed = collapsed))
}



test <- CKMR_est.med(samples1)





test_plot <- test$nhat %>%
  filter(Actual_Year != "Weighted_Focal_100") %>%
  mutate(Actual_Year = as.numeric(Actual_Year)) |> 
  mutate(dev_50 = 0.5 * (mean),
         upper=mean + dev_50,
         lower=mean - dev_50) |> 
  filter(`2.5%` > lower,
         `97.5%` < upper)
mean_nadult <- test$nhat$mean[21]

lower_nadult <- test$nhat$`2.5%`[21]

upper_nadult <- test$nhat$`97.5%`[21]


ggplot(test_plot, aes(x = Actual_Year)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = mean), colour = "red", linewidth = 1) +
  geom_line(data = true_abundance, aes(x = Actual_Year, y = True_Nadult), colour = "darkblue", linewidth = 1.1) +
  # coord_cartesian(xlim = c(93,105)) +
  geom_hline(yintercept = mean_nadult, linetype = 2, color = "darkgreen") +
  geom_ribbon(aes(ymin = lower_nadult, ymax = upper_nadult), fill = "darkgreen", alpha = 0.2)
  
