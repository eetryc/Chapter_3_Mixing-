#==================================##
#===== Sampling function ========####
#==================================##


## 1. Sample Function =====================================================####
sampler <- function(harvest_df, yr, sampling_window, sample_rate) {
  
  # 1. Check if the current year falls within the sampling window
  if (!(yr %in% sampling_window)) return(NULL)
  
  # 2. Check if any fish were actually caught this year
  if (nrow(harvest_df) == 0) return(NULL)
  
  # 3. Calculate target sample size based on the catch
  n_sample <- floor(sample_rate * nrow(harvest_df))
  if (n_sample == 0) return(NULL)
  
  # 4. Draw samples without replacement from the harvested fish
  sampled_fish <- harvest_df %>%
    slice_sample(n = min(n_sample, nrow(harvest_df)), replace = FALSE) %>%
    mutate(year_sampled = yr)
  
  return(sampled_fish)
}




## 2. Apply to the pop_lists from the simulated objects list ==============####

# Define sampling rate per yearly population #

sample_rate <- 0.1 # change percent as needed


# Apply function #
samples <- bind_rows(
  lapply(seq_along(p$caught_list), function(i) {
    sampler(
      harvest_df      = p$caught_list[[i]],  
      yr              = i,                   
      sampling_window = sampling_window,
      sample_rate     = sample_rate          
    )
  })
)
