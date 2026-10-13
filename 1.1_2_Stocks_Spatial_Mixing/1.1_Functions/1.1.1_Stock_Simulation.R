#==================================##
#======= Stock Sim. Function ====####
#==================================##





# ^ = variable to change


## 0. Load libraries, init. storage =======================================####

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

# rm(pop_list, all_fish, current_pop,samples,empty_offspring) # or
# rm(list=ls())




## 1. Define function and constants =======================================####

### Defining constants for all stocks =====================================####

init.prop.female  <- 0.5  # assume equal amount of sexes

rep.age <- 5 # minimum age to reproduce

max.age <- 40 # max age they'll live to 

mean_offspring <- 150 # mean number of offspring per female ^



### Simulation timing =====================================================####
years <- 10 # number of years to sample and estimate ^

burn.in <- 100 # number of years to run the simulation before counting 

run.years <- burn.in + years # total years to run

sampling_window <- (run.years - years + 1):run.years




### Mortality related constants ===========================================####
surv_vec <- c(0.35, 0.6, 0.7, 0.75, 0.8, rep(0.9, max.age - 4))

p_catch <- 0.10 # percent mortality of catch ^ 

K <- 10000  # carrying capacity to stabilize population ^ 

p_larval_survival <- 0.15


### Stock specific information ============================================####

sample_rate <- 0.05 # change percent as needed


## 2. Function ============================================####

## Faster function? =======================================================####

stock.sim_fast <- function(init.pop.size) {
  
  # --- 0. Pre-allocate Vectors / Initial Population ---
  age    <- rep(4, times = init.pop.size)
  sex    <- sample(c('1', '0'), size = init.pop.size, prob = c(init.prop.female, 1 - init.prop.female), replace = TRUE)
  mother <- rep(NA_character_, init.pop.size)
  father <- rep(NA_character_, init.pop.size)
  cohort <- rep(0, times = init.pop.size)
  number <- seq_len(init.pop.size)
  
  id     <- sprintf("%d_%d_%s", cohort, number, sex)
  
  current_pop <- data.frame(id, number, cohort, age, sex, mother, father, stringsAsFactors = FALSE)
  
  pop_list    <- vector("list", run.years + 1)
  caught_list <- vector("list", run.years)
  pop_list[[1]] <- current_pop
  
  max_number <- init.pop.size
  # --- Loop ---
  for (yr in 1:run.years) {
    
    n_current <- nrow(current_pop)
    if (n_current == 0) {
      message(paste("Population extinct at year", yr))
      break
    }
    
    # 1. Age everyone 
    current_pop$age <- current_pop$age + 1
    
    # 2. Select mature individuals
    is_female_mature <- current_pop$sex == "1" & current_pop$age >= rep.age
    is_male_mature   <- current_pop$sex == "0" & current_pop$age >= rep.age
    
    n_females <- sum(is_female_mature)
    n_males   <- sum(is_male_mature)
    
    offspring_df <- NULL
    
    # 3. Generate offspring
    if (n_females > 0 && n_males > 0) {
      females <- current_pop[is_female_mature, ]
      
      # Determine mating capability and clutch sizes
      can_mate <- rbinom(n_females, 1, 0.8)
      num_offspring <- rep(0, n_females)
      if (sum(can_mate) > 0) {
        num_offspring[can_mate == 1] <- rpois(sum(can_mate), lambda = mean_offspring)
      }
      
      total_offspring <- sum(num_offspring)
      
      if (total_offspring > 0) {
        realized_larval_survival <- p_larval_survival / (1 + (n_current / K))
        
        surviving_idx <- which(rbinom(total_offspring, 1, realized_larval_survival) == 1)
        n_survivors <- length(surviving_idx)
        
        if (n_survivors > 0) {
          # Only generate IDs and records for surviving offspring
          mom_ids_all <- rep(females$id, times = num_offspring)
          mom_ids     <- mom_ids_all[surviving_idx]
          
          dad_ids_all <- sample(current_pop$id[is_male_mature], size = total_offspring, replace = TRUE)
          dad_ids     <- dad_ids_all[surviving_idx]
          
          sex_off     <- sample(c('1', '0'), size = n_survivors, replace = TRUE, 
                                prob = c(init.prop.female, 1 - init.prop.female))
          
          off_numbers <- seq_len(n_survivors) + max_number
          max_number  <- max_number + n_survivors
          
          offspring_df <- data.frame(
            id     = sprintf("%d_%d_%s", yr, off_numbers, sex_off),
            number = off_numbers,
            cohort = yr,
            age    = 0,
            sex    = sex_off,
            mother = mom_ids,
            father = dad_ids,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    # 4. Combine populations 
    if (!is.null(offspring_df) && nrow(offspring_df) > 0) {
      # rbind on base data.frames is much faster here than tracking dplyr overhead
      current_pop <- rbind(current_pop, offspring_df)
    }
    
    # 5 & 6. Natural survival & removing the old
    # Combine filtering steps to minimize data frame copying
    surv_prob <- surv_vec[pmin(current_pop$age, length(surv_vec))]
    survive_nat <- rbinom(nrow(current_pop), 1, surv_prob)
    
    alive_idx <- which(survive_nat == 1 & current_pop$age <= max.age)
    current_pop <- current_pop[alive_idx, , drop = FALSE]
    
    # 7. Fishing mortality
    if (nrow(current_pop) > 0) {
      catch_probabilities <- ifelse(current_pop$age > 4, p_catch, 0)
      caught <- rbinom(nrow(current_pop), 1, catch_probabilities)
      
      caught_idx <- which(caught == 1)
      surv_idx   <- which(caught == 0)
      
      caught_list[[yr]]  <- current_pop[caught_idx, , drop = FALSE]
      current_pop        <- current_pop[surv_idx, , drop = FALSE]
    } else {
      caught_list[[yr]]  <- current_pop[0, ]
    }
    
    # 9. Store population for this year
    pop_list[[yr + 1]] <- current_pop
  }
  
  return(list(pop_list = pop_list, caught_list = caught_list))
}




































stock.sim <- function(init.pop.size){
  
  ## 0. Init. objects
  ##------## Inits for the starting population ##------##
  age <- rep(4, times = init.pop.size) # set at 4 so these guys can mate when the loop starts post aging
  
  sex <- sample(c('1','0'), size = init.pop.size,
                prob = c(init.prop.female, 1 - init.prop.female),
                replace = TRUE) 
  
  mother <- rep(NA, init.pop.size)
  
  father <- rep(NA, init.pop.size)
  
  cohort <- rep(0, times = init.pop.size)
  
  number <- seq_len(init.pop.size)
  
  id <- paste(cohort, number, sex, sep = "_")
  
  pop.year0 <- data.frame(id, number, cohort, age, sex, mother, father, stringsAsFactors = FALSE) # bind for starting pop table
  
  
  current_pop <- pop.year0
  
  samples     <- data.frame()
  
  pop_list    <- vector("list", run.years + 1)
  
  caught_list <- vector("list", run.years)
  
  pop_list[[1]] <- current_pop
  
  empty_offspring <- pop.year0[0, ] %>%
    mutate(
      id      = as.character(id),
      number  = as.integer(number),
      cohort  = as.integer(cohort),
      age     = as.numeric(age),
      sex     = as.character(sex),
      mother  = as.character(mother),
      father  = as.character(father)
    )
  
  
  
  
  
  for (yr in 1:run.years) {
    
    offspring_df <- empty_offspring
    
    ## 1. Age everyone
    current_pop <- current_pop %>% mutate(age = age + 1)
    
    
    
    
    ## 2. Select mature females and males
    mothers <- current_pop %>%
      filter(sex == "1", age >= rep.age)
    
    if (nrow(mothers) > 0) {
      mothers$can_mate <- rbinom(nrow(mothers), 1, 0.8)
      mothers <- mothers[mothers$can_mate == 1, ]
      mothers$can_mate <- NULL
      if (nrow(mothers) > 0) {
        mothers$num_offspring <- rpois(nrow(mothers), lambda = mean_offspring)
      }
    }
    
    fathers <- current_pop %>%
      filter(sex == "0", age >= rep.age)
    
    
    
    
    ## 3. Generate offspring (Fully Vectorized)
    if (nrow(mothers) > 0 && nrow(fathers) > 0) {
      
      # Total number of babies to make this year
      total_offspring <- sum(mothers$num_offspring)
      
      if (total_offspring > 0) {
        
        # Vectorized: Repeat each mother's ID based on how many offspring she has
        mom_ids <- rep(mothers$id, times = mothers$num_offspring)
        
        # Vectorized: Grab a father for every single baby all at once
        dad_ids <- sample(fathers$id, size = total_offspring, replace = TRUE)
        
        # Vectorized: Determine sex for all babies at once
        sex_off <- sample(c('1','0'), size = total_offspring, replace = TRUE,
                          prob = c(init.prop.female, 1 - init.prop.female))
        
        # Create the complete offspring dataframe in one shot (No rbind loop!)
        offspring_df <- data.frame(
          age    = 0,
          sex    = sex_off,
          mother = mom_ids,
          father = dad_ids,
          cohort = yr,
          number = NA_integer_,
          stringsAsFactors = FALSE
        )
        
        # Assign unique IDs
        max_number <- if (nrow(current_pop) > 0) max(current_pop$number) else 0
        offspring_df$number <- seq_len(total_offspring) + max_number
        offspring_df$id <- paste(offspring_df$cohort, offspring_df$number, offspring_df$sex, sep = "_")
        
        # Reorder columns to match
        offspring_df <- offspring_df[, names(pop.year0)]
        
      } else {
        offspring_df <- empty_offspring
      }
    }
    
    
    
    
    ## 4. Combine current population + offspring
    current_pop <- bind_rows(current_pop, offspring_df)
    
    
    
    
    
    ## 5. Age-structured natural survival
    
    current_pop$survive_nat <- rbinom(
      n = nrow(current_pop),
      size = 1,
      prob = surv_vec[pmin(current_pop$age + 1, length(surv_vec))])
    
    # Take out the olds
    current_pop <- current_pop %>% filter(age <= max.age)
    
    
    
    
    
    ## 6. Natural survivors
    alive_after_nat <- current_pop %>% filter(survive_nat == 1)
    
    
    
    
    ## 7. Fishing mortality
    catch_probabilities <- ifelse(alive_after_nat$age > 4, p_catch, 0)
    
    alive_after_nat$caught <- rbinom(
      n = nrow(alive_after_nat), 
      size = 1, 
      prob = catch_probabilities)
    
    harvested_this_year   <- alive_after_nat %>% filter(caught == 1)
    
    caught_list[[yr]] <- harvested_this_year
    
    survivors_to_next_year <- alive_after_nat %>% filter(caught == 0)
    
    
    
    
    
    ## 8. Density dependence applied ONLY to new recruits
    total_before <- nrow(current_pop)
    total_off    <- nrow(offspring_df)
    
    if (total_before + total_off > K && total_off > 0) {
      excess  <- (total_before + total_off) - K
      kill_n  <- min(excess, total_off)
      kill_idx <- sample.int(total_off, size = kill_n, replace = FALSE)
      offspring_df <- offspring_df[-kill_idx, , drop = FALSE]
    }
    
    
    
    
    
    ## 9. STORE POPULATION FOR THIS YEAR
    pop_list[[yr + 1]] <- survivors_to_next_year
    
    
    
    
    ## 10. Update population for next year
    current_pop <- survivors_to_next_year
    current_pop$survive_nat <- NULL
    current_pop$caught <- NULL
    
    
    
    
    ## 11. Extinction check
    if (nrow(current_pop) == 0) {
      message(paste("Population extinct at year", yr))
      break
    }
  }
  
  
  
  
  # --- RETURN results for this stock ---
  return(list(
    pop_list = pop_list,
    caught_list = caught_list
  ))
  
}



















