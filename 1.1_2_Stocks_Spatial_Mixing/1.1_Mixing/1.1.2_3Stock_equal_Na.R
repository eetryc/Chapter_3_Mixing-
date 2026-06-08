#====================================##
#======= Three Stock Mixing =======####
#====================================##


## Define init stock size for 3 zones =====================================####
  GMZ1.init <- 1000
  GMZ2.init <- 1000
  GMZ3.init <- 1000
  
  niter.stock <- 2


# Generate valid combinations of sampling proportions =====================####
  prop_grid <- expand.grid(
    p1 = seq(0, 1, by = 0.2),
    p2 = seq(0, 1, by = 0.2)
  ) %>%
    mutate(p3 = 1 - (p1 + p2)) %>%
    filter(p3 >= 0) # Keep only combinations that sum to <= 1
  
  head(prop_grid)



## Function for mixing stocks together  ===================================####
  
spatial_3stocks <- function(GMZ1.init, GMZ2.init, GMZ3.init, prop_grid) { 
  
  ### 1. Simulate 3 stocks ================================================####
  stock1 <- stock.sim_fast(GMZ1.init)
  stock2 <- stock.sim_fast(GMZ2.init)
  stock3 <- stock.sim_fast(GMZ3.init)
  
  ### 2. Sample stocks ====================================================####
  set.seed(7777777)
  
  # Stock 1 
  samples1 <- bind_rows(lapply(seq_along(stock1$caught_list), function(i) {
    sampler(harvest_df = stock1$caught_list[[i]], yr = i, 
            sampling_window = sampling_window, sample_rate = sample_rate)
  }))
  
  # Stock 2 
  samples2 <- bind_rows(lapply(seq_along(stock2$caught_list), function(i) {
    sampler(harvest_df = stock2$caught_list[[i]], yr = i,
            sampling_window = sampling_window, sample_rate = sample_rate)
  }))
  
  # Stock 3
  samples3 <- bind_rows(lapply(seq_along(stock3$caught_list), function(i) {
    sampler(harvest_df = stock3$caught_list[[i]], yr = i, 
            sampling_window = sampling_window, sample_rate = sample_rate)
  }))
  
  #### Assign stock origins ===============================================####
  samples1 <- samples1 %>% mutate(origin = "GMZ1") 
  samples2 <- samples2 %>% mutate(origin = "GMZ2")
  samples3 <- samples3 %>% mutate(origin = "GMZ3")
  
  ### 3. Mix samples and Run Iterations ===================================####
  results_list <- list()
  counter <- 1
  
  # Loop through rows of proportion grid
  for (i in 1:nrow(prop_grid)) {
    p1 <- prop_grid$p1[i]
    p2 <- prop_grid$p2[i]
    p3 <- prop_grid$p3[i]
    
    for (j in 1:niter.stock) {
      mixed_samp <- samples1 %>% slice_sample(prop = p1, replace = FALSE) %>% 
        rbind(slice_sample(samples2, prop = p2, replace = FALSE)) %>% 
        rbind(slice_sample(samples3, prop = p3, replace = FALSE))
      
      est_df <- CKMR_est(mixed_samp)
      
      est_df$p1        <- p1 
      est_df$p2        <- p2
      est_df$p3        <- p3
      est_df$iteration <- j
      
      results_list[[counter]] <- est_df
      counter <- counter + 1 
    }
  }
  
  ### 4. Process and Group Outputs ========================================####
  final_estimates <- bind_rows(results_list)
  
  # Group by all three proportions
  averaged_estimates <- final_estimates %>%
    group_by(p1, p2, p3, Actual_Year) %>%
    dplyr::summarise(
      avg_mean_Nadult = mean(mean, na.rm = TRUE),
      avg_lower_2.5   = mean(`2.5%`, na.rm = TRUE),
      avg_median_50   = mean(`50%`, na.rm = TRUE),
      avg_upper_97.5  = mean(`97.5%`, na.rm = TRUE),
      iters_counted   = dplyr::n(),
      rhat            = mean(Rhat, na.rm = TRUE),
      .groups         = "drop"
    )
  
  ### 5. Calculate true population abundance =============================####
  years_to_estimate <- sort(unique(samples1$cohort))
  
  true_abundance <- lapply(seq_along(stock1$pop_list), function(i) {
    adults_stock1 <- stock1$pop_list[[i]] %>% filter(age >= 5) %>% nrow()
    data.frame(Actual_Year = i, True_Nadult = adults_stock1)
  }) %>% bind_rows()
  
  true_abundance <- true_abundance %>% filter(Actual_Year %in% years_to_estimate)
  
  ### 6. Combine Truth and Estimates ======================================####
  comparison_df <- averaged_estimates %>%
    left_join(true_abundance, by = "Actual_Year")
  
  return(comparison_df)
}






## Run function ===========================================================####

stock.est.3 <- spatial_3stocks(GMZ1.init, GMZ2.init, GMZ3.init, prop_grid)

stock.est.3.trimmed <- stock.est.3 %>% 
  filter(rhat <= 1.05)



## Graph ==================================================================####

ggplot(stock.est.3.trimmed, aes(x = Actual_Year)) +
  geom_ribbon(aes(ymin = avg_lower_2.5, ymax = avg_upper_97.5), fill = "blue", alpha = 0.15) +
  geom_line(aes(y = avg_mean_Nadult, color = "CKMR Estimate"), size = 0.8) +
  geom_line(aes(y = True_Nadult, color = "True Population"), linetype = "dashed", size = 0.8) +
  # Rows represent Stock 1 prop, columns represent Stock 2 prop
  facet_grid(p1 ~ p2, labeller = label_both) + 
  labs(
    title = "CKMR Abundance Estimates Across 3-Stock Mixing Scenarios",
    subtitle = "Facets show proportions of Stock 1 (rows) and Stock 2 (columns)",
    x = "Model Year",
    y = "Adult Abundance (N)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 8))



### Violin plot ===========================================================####
stock.est.3_violin <- stock.est.3.trimmed %>% 
  mutate(
    rel.bias = ((avg_mean_Nadult - True_Nadult) / True_Nadult),
    mix_combo = paste0("P1:", p1, "|P2:", p2)
  ) 

ggplot(stock.est.3_violin, aes(x = mix_combo, y = rel.bias, fill = as.factor(p1))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  geom_violin(alpha = 0.6, show.legend = TRUE) +
  scale_fill_viridis_d(name = "Prop Stock 1") +
  labs(
    title = "Relative Bias by 3-Stock Sampling Proportions",
    x = "Sampling Combination Matrix",
    y = "Relative Bias"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))