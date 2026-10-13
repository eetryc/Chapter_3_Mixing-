# Load plyr BEFORE tidyverse to prevent function masking conflicts
library(plyr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lubridate)

### 1. Define global constants first so they exist when needed =======####
init.prop.female <- 0.5   # assume equal amount of sexes
rep.age <- 5              # minimum age to reproduce
max.age <- 40             # max age they'll live to 
mean_offspring <- 100     # mean number of offspring per female 

### Simulation timing =====================================================####
years <- 20 
burn.in <- 100 
run.years <- burn.in + years 
sampling_window <- (run.years - years + 1):run.years

### Mortality related constants ===========================================####
surv_vec <- c(0.45, 0.65, 0.75, 0.85, 0.85, rep(0.94, max.age - 4))
p_catch <- 0.1

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
      current_pop <- rbind(current_pop, offspring_df)
    }
    
    # 5 & 6. Natural survival & removing the old
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

GMZ1.init <- 1000
GMZ2.init <- 1000
GMZ3.init <- 1000
GMZ4.init <- 1000
GMZ5.init <- 1000

sampler <- function(harvest_df, yr, sampling_window, sample_rate) {
  if (!(yr %in% sampling_window)) return(NULL)
  if (nrow(harvest_df) == 0) return(NULL)
  n_sample <- floor(sample_rate * nrow(harvest_df))
  if (n_sample == 0) return(NULL)
  
  sampled_fish <- harvest_df %>%
    slice_sample(n = min(n_sample, nrow(harvest_df)), replace = FALSE) %>%
    mutate(year_sampled = yr)
  
  return(sampled_fish)
}

sample_rate <- 0.1 

# --- Stock 1 ---
K <- 100000 
p_larval_survival <- 0.7
t1 <- stock.sim_fast(GMZ1.init)

# --- Stock 2 ---
mean_offspring <- 70 
p_larval_survival <- 0.4
t2 <- stock.sim_fast(GMZ2.init)

# --- Stock 3 ---
K <- 50000 
mean_offspring <- 30 
p_larval_survival <- 0.3
sample_rate <- 0.05 
t3 <- stock.sim_fast(GMZ3.init)

# --- Stock 4 ---
K <- 15000 
mean_offspring <- 25 
p_larval_survival <- 0.3
t4 <- stock.sim_fast(GMZ4.init)

# --- Stock 5 ---
K <- 15000 
mean_offspring <- 25 
p_larval_survival <- 0.2
t5 <- stock.sim_fast(GMZ5.init)

samp.props.real <- c(0.09, 0.44, 0.19, 0.17, 0.11)

set.seed(7777777)

samples1 <- bind_rows(lapply(seq_along(t1$caught_list), function(i) {
  sampler(harvest_df = t1$caught_list[[i]], yr = i, sampling_window = sampling_window, sample_rate = sample_rate)
}))
samples2 <- bind_rows(lapply(seq_along(t2$caught_list), function(i) {
  sampler(harvest_df = t2$caught_list[[i]], yr = i, sampling_window = sampling_window, sample_rate = sample_rate)
}))
samples3 <- bind_rows(lapply(seq_along(t3$caught_list), function(i) {
  sampler(harvest_df = t3$caught_list[[i]], yr = i, sampling_window = sampling_window, sample_rate = sample_rate)
}))
samples4 <- bind_rows(lapply(seq_along(t4$caught_list), function(i) {
  sampler(harvest_df = t4$caught_list[[i]], yr = i, sampling_window = sampling_window, sample_rate = sample_rate)
}))
samples5 <- bind_rows(lapply(seq_along(t5$caught_list), function(i) {
  sampler(harvest_df = t5$caught_list[[i]], yr = i, sampling_window = sampling_window, sample_rate = sample_rate)
}))

samples1 <- samples1 %>% mutate(origin = "GMZ1") 
samples2 <- samples2 %>% mutate(origin = "GMZ2")
samples3 <- samples3 %>% mutate(origin = "GMZ3")
samples4 <- samples4 %>% mutate(origin = "GMZ4")
samples5 <- samples5 %>% mutate(origin = "GMZ5")

total_target_n <- 5000
target_ns <- round(total_target_n * samp.props.real)

# Helper function to safely sample even if target_n exceeds available rows
safe_slice <- function(df, target_n) {
  if (nrow(df) == 0) return(df)
  use_replace <- nrow(df) < target_n
  df %>% slice_sample(n = target_n, replace = use_replace)
}

# Construct mixed_samp safely
mixed_samp <- bind_rows(
  safe_slice(samples1, target_ns[1]),
  safe_slice(samples2, target_ns[2]),
  safe_slice(samples3, target_ns[3]),
  safe_slice(samples4, target_ns[4]),
  safe_slice(samples5, target_ns[5])
)

sampled <- nrow(mixed_samp)
print(paste("Total mixed sample size:", sampled))

true_abundance.1 <- lapply(seq_along(t1$pop_list), function(i) {
  adults_stock1 <- t1$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
  data.frame(Actual_Year = i, True_Nadult = adults_stock1)
}) %>% bind_rows() |> mutate(Stock = 1)

true_abundance.2 <- lapply(seq_along(t2$pop_list), function(i) {
  adults_stock1 <- t2$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
  data.frame(Actual_Year = i, True_Nadult = adults_stock1)
}) %>% bind_rows() |> mutate(Stock = 2)

true_abundance.3 <- lapply(seq_along(t3$pop_list), function(i) {
  adults_stock1 <- t3$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
  data.frame(Actual_Year = i, True_Nadult = adults_stock1)
}) %>% bind_rows()|> mutate(Stock = 3)

true_abundance.4 <- lapply(seq_along(t4$pop_list), function(i) {
  adults_stock1 <- t4$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
  data.frame(Actual_Year = i, True_Nadult = adults_stock1)
}) %>% bind_rows()|> mutate(Stock = 4)

true_abundance.5 <- lapply(seq_along(t5$pop_list), function(i) {
  adults_stock1 <- t5$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
  data.frame(Actual_Year = i, True_Nadult = adults_stock1)
}) %>% bind_rows()|> mutate(Stock = 5)

full_lake_abundance <- true_abundance.1 %>% 
  bind_rows(true_abundance.2, true_abundance.3, true_abundance.4, true_abundance.5) %>% 
  group_by(Actual_Year) %>% 
  dplyr::summarize(Total_True_Nadult = sum(True_Nadult), .groups = "drop") 

full_lake_abundance_target <- full_lake_abundance
trunc <- true_abundance.1

# Extract Stock 3 samples
stock3_extracted_samples <- samples3

# Save simulation outputs
saveRDS(list(t1 = t1, t2 = t2, t3 = t3, t4 = t4, t5 = t5,
             true_abundance.1 = true_abundance.1, true_abundance.2 = true_abundance.2, 
             true_abundance.3 = true_abundance.3, true_abundance.4 = true_abundance.4, 
             true_abundance.5 = true_abundance.5, mixed_samp = mixed_samp, 
             trunc = trunc, full_lake_abundance_target = full_lake_abundance_target, 
             sampled = sampled, stock3_extracted_samples = stock3_extracted_samples), 
        file = "stock_sim_output_2.rds")
