# load packages


library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)



################################################################################
# Set up the simulation function ###############################################

#### 1. Define function and constants ####

##------## Defining constants for all stocks ##------##

init.prop.female  <- 0.5  # assume equal amount of sexes

rep.age <- 5 # minimum age to reproduce

max.age <- 40 # max age they'll live to 

mean_offspring <- 15 # mean number of offspring per female 



##------## Simulation timing ##------##
years <- 20 # number of years to sample and estiamte ^

burn.in <- 100 # number of years to run the simulation before counting 

run.years <- burn.in + years # total years to run

sampling_window <- (run.years - years + 1):run.years




##------## Mortality related constants ##------##
surv_vec <- c(0.35, 0.6, 0.7, 0.75, 0.8, rep(0.9, max.age - 4)) # age structured survival

p_catch <- 0.15 # percent mortality of catch ^ 

K <- 10000  # carrying capacity to stabilize population ^ 




##------## Stock specific information ##------##

sample_rate <- 0.03 # change percent as needed




## Function

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

  
  # put all that in a data frame
  pop.year0 <- data.frame(id, number, cohort, age, sex, mother, father, stringsAsFactors = FALSE) # bind for starting pop table

  next_id <- init.pop.size + 1L
  
  current_pop <- pop.year0
  samples     <- data.frame()
  pop_list    <- vector("list", run.years + 1)
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
    
    
## 3. Generate offspring (Fully Vectorized Rewrite)
    if (nrow(mothers) > 0 && nrow(fathers) > 0) {
      
      # Total number of new babies to generate this year
      total_offspring <- sum(mothers$num_offspring)
      
      if (total_offspring > 0) {
        
        # Repeat each mother's ID by her specific offspring count
        mom_ids <- rep(mothers$id, times = mothers$num_offspring)
        
        # Draw a father for every baby all at once
        dad_ids <- sample(fathers$id, size = total_offspring, replace = TRUE)
        
        # Assign sex to all offspring simultaneously 
        sex_off <- sample(c('1','0'), size = total_offspring, replace = TRUE,
                          prob = c(init.prop.female, 1 - init.prop.female))
        
        # Create the entire cohort's dataframe
        offspring_df <- data.frame(
          age    = 0,
          sex    = sex_off,
          mother = mom_ids,
          father = dad_ids,
          cohort = yr,
          number = NA_integer_,
          stringsAsFactors = FALSE
        )
        
        # Assign unique sequential numbers and IDs
    offspring_df$number <- seq.int(
      from = next_id,
      length.out = total_offspring
    )

    next_id <- next_id + total_offspring
        offspring_df$id <- paste(offspring_df$cohort, offspring_df$number, offspring_df$sex, sep = "_")
        
        # Reorder columns to match the baseline schema
        offspring_df <- offspring_df[, names(pop.year0)]
        
      } else {
        offspring_df <- empty_offspring
      }
    }
## 4. Density dependence applied ONLY to new recruits

total_before <- nrow(current_pop)
total_off <- nrow(offspring_df)

if (total_before + total_off > K && total_off > 0) {

  excess <- (total_before + total_off) - K

  kill_n <- min(excess, total_off)

  if (kill_n > 0) {

    kill_idx <- sample.int(
      total_off,
      size = kill_n,
      replace = FALSE
    )

    offspring_df <- offspring_df[-kill_idx, , drop = FALSE]
  }
}

## 5. Combine current population + surviving offspring
current_pop <- bind_rows(current_pop, offspring_df)
    
    ## 6. Age-structured natural survival
    
    current_pop$survive_nat <- rbinom(
      n = nrow(current_pop),
      size = 1,
      prob = surv_vec[pmin(current_pop$age + 1, length(surv_vec))])
    
    # Take out the olds
    current_pop <- current_pop %>% filter(age <= max.age)
    
    
    
    ## 7. Natural survivors
    alive_after_nat <- current_pop %>% filter(survive_nat == 1)
    
    ## 8. Fishing mortality
    alive_after_nat$caught <- rbinom(
    n = nrow(alive_after_nat),
    size = 1,
    prob = ifelse(
      alive_after_nat$age >= rep.age,
      p_catch,
      0
    )
  )
    
    harvested_this_year   <- alive_after_nat %>% filter(caught == 1)
    survivors_to_next_year <- alive_after_nat %>% filter(caught == 0)
    

  ## 9. Sampling (last years) – lethal sampling
  if (yr %in% sampling_window) {

    adults_alive <- survivors_to_next_year %>%
      filter(age >= rep.age)

    n_sample <- floor(sample_rate * nrow(adults_alive))

    sampled_fish <- adults_alive %>%
      slice_sample(
        n = min(n_sample, nrow(adults_alive)),
        replace = FALSE
      ) %>%
      mutate(year_sampled = yr)

    ## Store sampled individuals
    samples <- bind_rows(samples, sampled_fish)

    ## Remove sampled individuals from the population
    survivors_to_next_year <- survivors_to_next_year %>%
      filter(!id %in% sampled_fish$id)
  }
    
    ## 10. STORE POPULATION FOR THIS YEAR
    pop_list[[yr + 1]] <- survivors_to_next_year
    
    ## 11. Update population for next year
    current_pop <- survivors_to_next_year
    current_pop$survive_nat <- NULL
    current_pop$caught <- NULL
    
    ## 12. Extinction check
    if (nrow(current_pop) == 0) {
      message(paste("Population extinct at year", yr))
      break
    }
  }
  # --- RETURN results for this stock ---
  return(list(
    pop_list = pop_list,
    samples  = samples
  ))
  
}

################################################################################










################################################################################
# Simualte the stocks seperately ###############################################

#### Simulate GMZ1 ####


GMZ1 <- stock.sim(10000)

names(GMZ1$pop_list) <- 0:(length(GMZ1$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ1_years_to_estimate <- sort(unique(GMZ1$samples$cohort))


## Population size per year
GMZ1_true_abundance <- lapply(GMZ1$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ1_true_abundance <- as.data.frame(GMZ1_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ1_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ1_true_abundance$name <- as.integer(sub("X", "",GMZ1_true_abundance$name))
GMZ1_true_abundance$name <- GMZ1_true_abundance$name %>% 
  as.integer()
GMZ1_true_abundance <- GMZ1_true_abundance %>% 
  filter(name %in% GMZ1_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ1_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()



#### Simulate GMZ2 ####

GMZ2 <- stock.sim(10000)



names(GMZ2$pop_list) <- 0:(length(GMZ2$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ2_years_to_estimate <- sort(unique(GMZ2$samples$cohort))


## Population size per year
GMZ2_true_abundance <- lapply(GMZ2$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ2_true_abundance <- as.data.frame(GMZ2_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ2_true_abundance # check to make sure there isn't letters in the year names. Re run if true


GMZ2_true_abundance$name <- as.integer(sub("X", "",GMZ2_true_abundance$name))
GMZ2_true_abundance$name <- GMZ2_true_abundance$name %>% 
  as.integer()
GMZ2_true_abundance <- GMZ2_true_abundance %>% 
  filter(name %in% GMZ2_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ2_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()




#### Join samples based on sampling ####
# Add column to make it known where each fish came from
GMZ1.samp <- GMZ1$samples %>% 
  mutate(origin = "GMZ1" )

GMZ2.samp <- GMZ2$samples %>% 
  mutate(origin = "GMZ2" )


# Full mixing (equal proportions of both)
all_samples_full <- GMZ1.samp %>% 
  rbind(GMZ2.samp)
















################################################################################
# Do larger sized stocks ##########################################
# target stock is larger, secondary is smaller 
#### Simulate GMZ1 ####


GMZ1.large <- stock.sim(20000)

names(GMZ1.large$pop_list) <- 0:(length(GMZ1.large$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ1.large_years_to_estimate <- sort(unique(GMZ1.large$samples$cohort))


## Population size per year
GMZ1.large_true_abundance <- lapply(GMZ1.large$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ1.large_true_abundance <- as.data.frame(GMZ1.large_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ1.large_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ1.large_true_abundance$name <- as.integer(sub("X", "",GMZ1.large_true_abundance$name))
GMZ1.large_true_abundance$name <- GMZ1.large_true_abundance$name %>% 
  as.integer()
GMZ1.large_true_abundance <- GMZ1.large_true_abundance %>% 
  filter(name %in% GMZ1.large_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ1.large_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()



#### Simulate GMZ2 ####

GMZ2.large <- GMZ2



names(GMZ2.large$pop_list) <- 0:(length(GMZ2.large$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ2.large_years_to_estimate <- sort(unique(GMZ2.large$samples$cohort))


## Population size per year
GMZ2.large_true_abundance <- lapply(GMZ2.large$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ2.large_true_abundance <- as.data.frame(GMZ2.large_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ2.large_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ2.large_true_abundance$name <- as.integer(sub("X", "",GMZ2.large_true_abundance$name))
GMZ2.large_true_abundance$name <- GMZ2.large_true_abundance$name %>% 
  as.integer()
GMZ2.large_true_abundance <- GMZ2.large_true_abundance %>% 
  filter(name %in% GMZ2.large_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ2.large_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()




#### Join samples based on sampling ####
# Add column to make it known where each fish came from
GMZ1.large.samp <- GMZ1.large$samples %>% 
  mutate(origin = "GMZ1" )

GMZ2.large.samp <- GMZ2.large$samples %>% 
  mutate(origin = "GMZ2" )


# Full mixing (equal proportions of both)
all_samples_full_large <- GMZ1.large.samp %>% 
  rbind(GMZ2.large.samp)














# Do small stocks ##################################################################


#### Simulate GMZ1 ####


GMZ1.small <- GMZ1

names(GMZ1.small$pop_list) <- 0:(length(GMZ1.small$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ1.small_years_to_estimate <- sort(unique(GMZ1.small$samples$cohort))


## Population size per year
GMZ1.small_true_abundance <- lapply(GMZ1.small$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ1.small_true_abundance <- as.data.frame(GMZ1.small_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ1.small_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ1.small_true_abundance$name <- as.integer(sub("X", "",GMZ1.small_true_abundance$name))
GMZ1.small_true_abundance$name <- GMZ1.small_true_abundance$name %>% 
  as.integer()
GMZ1.small_true_abundance <- GMZ1.small_true_abundance %>% 
  filter(name %in% GMZ1.small_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ1.small_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()



#### Simulate GMZ2 ####

GMZ2.small <- GMZ1.large



names(GMZ2.small$pop_list) <- 0:(length(GMZ2.small$pop_list)-1) # rename the year for population so it starts at 0 (cohort will match year). this line should get rid of any weird letter

GMZ2.small_years_to_estimate <- sort(unique(GMZ2.small$samples$cohort))


## Population size per year
GMZ2.small_true_abundance <- lapply(GMZ2.small$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ2.small_true_abundance <- as.data.frame(GMZ2.small_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ2.small_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ2.small_true_abundance$name <- as.integer(sub("X", "",GMZ2.small_true_abundance$name))
GMZ2.small_true_abundance$name <- GMZ2.small_true_abundance$name %>% 
  as.integer()
GMZ2.small_true_abundance <- GMZ2.small_true_abundance %>% 
  filter(name %in% GMZ2.small_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


# Graph
ggplot(GMZ2.small_true_abundance, aes(x = year_index, y = value, group = 1)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  # scale_y_continuous(limits = c(0, 200)) + 
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  theme_bw()




#### Join samples based on sampling ####
# Add column to make it known where each fish came from
GMZ1.small.samp <- GMZ1.small$samples %>% 
  mutate(origin = "GMZ1" )

GMZ2.small.samp <- GMZ2.small$samples %>% 
  mutate(origin = "GMZ2" )


# Full mixing (equal proportions of both)
all_samples_full_small <- GMZ1.small.samp %>% 
  rbind(GMZ2.small.samp)










# Now, mix the stocks ######################################################################
# Equal ####################################################################################


## Sampled mixed proportionse
# 0 %
all_samples_mix0 <- GMZ1.samp

# 10%
mix10_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

mix10_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

all_samples_mix10 <- mix10_GMZ1 %>% 
  rbind(mix10_GMZ2)

# 20%
mix20_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

mix20_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

all_samples_mix20 <- mix20_GMZ1 %>% 
  rbind(mix20_GMZ2)

# 30% 
mix30_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.7, replace = F)

mix30_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

all_samples_mix30 <- mix30_GMZ1 %>% 
  rbind(mix30_GMZ2)


# 40% 
mix40_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

mix40_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

all_samples_mix40 <- mix40_GMZ1 %>% 
  rbind(mix40_GMZ2)

# 50% 
mix50_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

mix50_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

all_samples_mix50 <- mix50_GMZ1 %>% 
  rbind(mix50_GMZ2)

# 60% 
mix60_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

mix60_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

all_samples_mix60 <- mix60_GMZ1 %>% 
  rbind(mix60_GMZ2)



# 70% 
mix70_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

mix70_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.7, replace = F)

all_samples_mix70 <- mix70_GMZ1 %>% 
  rbind(mix70_GMZ2)

# 80% 
mix80_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

mix80_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

all_samples_mix80 <- mix80_GMZ1 %>% 
  rbind(mix80_GMZ2)


# 90% 
mix90_GMZ1 <- GMZ1.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

mix90_GMZ2 <- GMZ2.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

all_samples_mix90 <- mix90_GMZ1 %>% 
  rbind(mix90_GMZ2)

# 100% 


all_samples_mix100  <- GMZ2.samp  







# Larger stock ########################################################################
## Sampled mixed proportionse
# 0 %
all_samples_mix0.large <- GMZ1.large.samp

# 10%
mix10_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

mix10_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

all_samples_mix10.large <- mix10_GMZ1.large %>% 
  rbind(mix10_GMZ2.large)

# 20%
mix20_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

mix20_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

all_samples_mix20.large <- mix20_GMZ1.large %>% 
  rbind(mix20_GMZ2.large)

# 30% 
mix30_GMZ1.large <- GMZ1.large.samp %>% 
  slice_sample(prop = 0.7, replace = F)

mix30_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

all_samples_mix30.large <- mix30_GMZ1.large %>% 
  rbind(mix30_GMZ2.large)


# 40% 
mix40_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

mix40_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

all_samples_mix40.large <- mix40_GMZ1.large %>% 
  rbind(mix40_GMZ2.large)

# 50% 
mix50_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

mix50_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

all_samples_mix50.large <- mix50_GMZ1.large %>% 
  rbind(mix50_GMZ2.large)

# 60% 
mix60_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

mix60_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

all_samples_mix60.large <- mix60_GMZ1.large %>% 
  rbind(mix60_GMZ2.large)



# 70% 
mix70_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

mix70_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.7, replace = F)

all_samples_mix70.large <- mix70_GMZ1.large %>% 
  rbind(mix70_GMZ2.large)

# 80% 
mix80_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

mix80_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

all_samples_mix80.large <- mix80_GMZ1.large %>% 
  rbind(mix80_GMZ2.large)


# 90% 
mix90_GMZ1.large <- GMZ1.large.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

mix90_GMZ2.large <- GMZ2.large.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

all_samples_mix90.large <- mix90_GMZ1.large %>% 
  rbind(mix90_GMZ2.large)

# 100% 


all_samples_mix100.large  <- GMZ2.large.samp  




# Smaller #############################################################


## Sampled mixed proportionse
# 0 %
all_samples_mix0.small <- GMZ1.small.samp

# 10%
mix10_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

mix10_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

all_samples_mix10.small <- mix10_GMZ1.small %>% 
  rbind(mix10_GMZ2.small)

# 20%
mix20_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

mix20_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

all_samples_mix20.small <- mix20_GMZ1.small %>% 
  rbind(mix20_GMZ2.small)

# 30% 
mix30_GMZ1.small <- GMZ1.small.samp %>% 
  slice_sample(prop = 0.7, replace = F)

mix30_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

all_samples_mix30.small <- mix30_GMZ1.small %>% 
  rbind(mix30_GMZ2.small)


# 40% 
mix40_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

mix40_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

all_samples_mix40.small <- mix40_GMZ1.small %>% 
  rbind(mix40_GMZ2.small)

# 50% 
mix50_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

mix50_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.5, replace = F)

all_samples_mix50.small <- mix50_GMZ1.small %>% 
  rbind(mix50_GMZ2.small)

# 60% 
mix60_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.4, replace = F)

mix60_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.6, replace = F)

all_samples_mix60.small <- mix60_GMZ1.small %>% 
  rbind(mix60_GMZ2.small)



# 70% 
mix70_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.3, replace = F)

mix70_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.7, replace = F)

all_samples_mix70.small <- mix70_GMZ1.small %>% 
  rbind(mix70_GMZ2.small)

# 80% 
mix80_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.2, replace = F)

mix80_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.8, replace = F)

all_samples_mix80.small <- mix80_GMZ1.small %>% 
  rbind(mix80_GMZ2.small)


# 90% 
mix90_GMZ1.small <- GMZ1.small.samp  %>% 
  slice_sample(prop = 0.1, replace = F)

mix90_GMZ2.small <- GMZ2.small.samp  %>% 
  slice_sample(prop = 0.9, replace = F)

all_samples_mix90.small <- mix90_GMZ1.small %>% 
  rbind(mix90_GMZ2.small)

# 100% 


all_samples_mix100.small  <- GMZ2.small.samp  











# Set up CKMR #############################################################################

### Define a function to input the sample information into ####

CKMR_est <- function(samples){
### POPs filter 
POpairs <- samples %>%
  rename_with(~ paste0(.x, "_1"), everything()) %>% # add the columns for each indiv
  cross_join(samples %>% rename_with(~ paste0(.x, "_2"), everything())) %>% # remove self-pairings
  # mutate to add the POP possibility (y/n) with the year and reproductive age filters
  mutate(HSPPOP_candidate = ifelse(cohort_1 < cohort_2 & 
                                     (cohort_2 - cohort_1 >= rep.age), "1","0")) %>% 
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




# join the two tables
library(plyr)
possible_pairs <- POPonly  %>% 
  mutate(
    sex_num = case_when(
      sex_1 == "U" ~ NA_real_, # unknown to na
      sex_1 == "0" ~ 1, # male is 1
      sex_1 == "1" ~ 0 # female is 0
    )
  )






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

group_vars_POP <- c("cohort_2", "knownSex", "IsFemale", "RObase", "Age_dif_exp")


collapsed_POP <- observations %>% 
  filter(RObase == 1) %>% 
  drop_na(mother_1, father_1) %>% 
  group_by(year_index,across(all_of(group_vars_POP))) %>%
  dplyr::summarise(
    Pair_viable_count_per_agedif = n(),                   
    TruePairs = sum(true_POP, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(cohort_2) %>% 
  mutate(AgeDif = 0)


collapsed_POP %>%  
  ggplot(aes(x = year_index, y = TruePairs)) +
  geom_line()

#### Just POP model
data_POP = list( years = years,  # number of cohorts,
                 TruePairs = collapsed_POP$TruePairs,
                 RObase = collapsed_POP$RObase,
                 nobs = nrow(collapsed_POP),
                 Pair_viable_count_per_agedif = collapsed_POP$Pair_viable_count_per_agedif,
                 year=collapsed_POP$year_index)

# Initial values
inits_POP = function() {
  list(
    # surv = runif(1,0.01, 0.99),
    Nadult = pmin(pmax(rnorm(years, mean = 7500, sd = 1000), 0), 10000))
}


# Parameters to follow
params = c("Nadult")

nburn <- 10000
nchains <- 3
niter <- 20000
n.cores = 3


Out_CKMR_sim1_POP = jagsUI::jags(data_POP, inits_POP, params, "HSPPOP.simplest.sim.jags", n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = T, verbose = T)



#### Plot against true populations and compare
# Extract posterior summary as a data frame
nhat_POP <- as.data.frame(Out_CKMR_sim1_POP$summary)

# Keep only the Nadult rows
nhat_POP <- nhat_POP[grep("^Nadult", rownames(nhat_POP)), ]

# Add a Year index (1, 2, 3, ...)
nhat_POP$Year <- 1:nrow(nhat_POP)

# Select the useful columns
nhat_df_POP <- nhat_POP %>%
  dplyr::select(Year, mean, `2.5%`, `50%`, `97.5%`)
}



# JAGS

cat("model{

  # Adult numbers


  for(i in 1:years) {

  Nadult[i] ~ dunif(0, 10000) }

	
  for(j in 1:nobs) {
    
    TruePairs[j] ~ dbinom(
      (RObase[j]) /
      (Nadult[year[j]]),
      Pair_viable_count_per_agedif[j]
    )
  }}", file = "HSPPOP.simplest.sim.jags")



# Make CKMR estimates #######################################################################

### Run CKMR
mix.0 <- CKMR_est(all_samples_mix0) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 0)
mix.10 <- CKMR_est(all_samples_mix10)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 10)
mix.20 <- CKMR_est(all_samples_mix20) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 20)
mix.30 <- CKMR_est(all_samples_mix30) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 30)
mix.40 <- CKMR_est(all_samples_mix40)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 40)
mix.50 <- CKMR_est(all_samples_mix50) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 50)
mix.60 <- CKMR_est(all_samples_mix60)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 60)
mix.70 <- CKMR_est(all_samples_mix70) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 70)
mix.80 <- CKMR_est(all_samples_mix80)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 80)
mix.90 <- CKMR_est(all_samples_mix90)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 90)
mix.100 <- CKMR_est(all_samples_mix100)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 100)



# Large CKMR estimates ###########################################################


## Population size per year
GMZ1_true_abundance <- lapply(GMZ1.large$pop_list, function(df){
  df %>% 
    filter(age >=5) %>%
    nrow()})

GMZ1_true_abundance <- as.data.frame(GMZ1_true_abundance) %>% 
  pivot_longer(cols = everything())

GMZ1_true_abundance # check to make sure there isn't letters in the year names. Re run if true, idk why this is happening


GMZ1_true_abundance$name <- as.integer(sub("X", "",GMZ1_true_abundance$name))
GMZ1_true_abundance$name <- GMZ1_true_abundance$name %>% 
  as.integer()
GMZ1_true_abundance <- GMZ1_true_abundance %>% 
  filter(name %in% GMZ1_years_to_estimate) %>% 
  mutate(year_index = row_number(.))


### Run CKMR
mix.0.large <- CKMR_est(all_samples_mix0.large) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 0)
mix.10.large <- CKMR_est(all_samples_mix10.large)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 10)
mix.20.large <- CKMR_est(all_samples_mix20.large) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 20)
mix.30.large <- CKMR_est(all_samples_mix30.large) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 30)
mix.40.large <- CKMR_est(all_samples_mix40.large)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 40)
mix.50.large <- CKMR_est(all_samples_mix50.large) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 50)
mix.60.large <- CKMR_est(all_samples_mix60.large)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 60)
mix.70.large <- CKMR_est(all_samples_mix70.large) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 70)
mix.80.large <- CKMR_est(all_samples_mix80.large)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 80)
mix.90.large <- CKMR_est(all_samples_mix90.large)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 90)
mix.100.large <- CKMR_est(all_samples_mix100.large)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 100)




# Small #####################################################################
### Run CKMR
mix.0.small <- CKMR_est(all_samples_mix0.small) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 0)
mix.10.small <- CKMR_est(all_samples_mix10.small)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 10)
mix.20.small <- CKMR_est(all_samples_mix20.small) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 20)
mix.30.small <- CKMR_est(all_samples_mix30.small) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 30)
mix.40.small <- CKMR_est(all_samples_mix40.small)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 40)
mix.50.small <- CKMR_est(all_samples_mix50.small) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 50)
mix.60.small <- CKMR_est(all_samples_mix60.small)  %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 60)
mix.70.small <- CKMR_est(all_samples_mix70.small) %>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 70)
mix.80.small <- CKMR_est(all_samples_mix80.small)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 80)
mix.90.small <- CKMR_est(all_samples_mix90.small)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 90)
mix.100.small <- CKMR_est(all_samples_mix100.small)%>% 
  left_join(GMZ1_true_abundance, by = c("Year" = "year_index")) %>% 
  dplyr::rename(True_abundance = value) %>% 
  mutate(Difference = True_abundance - mean,
         Mixed = 100)

combined_estimates <- bind_rows(
  list(
    base  = bind_rows(mix.0, mix.10, mix.20, mix.30, mix.40, mix.50,
                      mix.60, mix.70, mix.80, mix.90, mix.100),
    
    large = bind_rows(mix.0.large, mix.10.large, mix.20.large, mix.30.large,
                      mix.40.large, mix.50.large, mix.60.large, mix.70.large,
                      mix.80.large, mix.90.large, mix.100.large),
    
    small = bind_rows(mix.0.small, mix.10.small, mix.20.small, mix.30.small,
                      mix.40.small, mix.50.small, mix.60.small, mix.70.small,
                      mix.80.small, mix.90.small, mix.100.small)
  ),
  .id = "size"
) %>% 
  mutate(conf.range = `97.5%`-`2.5%`)

combined_estimates <- combined_estimates %>% 
  filter(conf.range <= (mean)) 
  # filter(name >= 100)

ggplot(mix.50, aes(x = Year, y = mean)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(
    x = "Year (Cohort index)",
    y = expression(hat(N)[adult])) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = "lightblue") +
  theme_bw() +
  geom_point(data = GMZ1_true_abundance, aes(x = year_index, y = value, group = 1),color = "red") +
  geom_line(data = GMZ1_true_abundance, aes(x = year_index, y = value, group = 1),color = "red")

  
#### Graph ####
combined_estimates %>% 
  ggplot(aes(x = Mixed, y = Difference)) +
  geom_point() +
  theme_bw() +
  stat_smooth(method = "lm") +
  labs(x = "Percent fish sampled from non-studied stock", y = "Difference between true abundance and estimate")










