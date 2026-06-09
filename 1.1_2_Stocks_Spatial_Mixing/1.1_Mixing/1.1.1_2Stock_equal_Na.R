#==================================##
#======= Two Stock Mixing =======####
#==================================##

library(ggplot2)


# This example is simply two equal sized stocks that are moving between 
# spatial management zones without engaging in natal straying (thus, no genetic
# mixture)


# You MUST define functions using scripts 1.1.2-1.1.4

## Define values for stock mixing function =================================####
# Define vector of stock mixing proportions
  # only for stock one as stock 2 will be 1-n

# samp.props <- c(seq(from= 0.1, to = 1, by = 0.1))
samp.props <- c(0.9,1)


# Define init stock size

GMZ1.init <- 1000

GMZ2.init <- 1000

niter.stock <- 2



## Function for mixing stocks together  ===================================####
spatial_2stocks <- function(GMZ1.init, GMZ2.init, samp.props) { 
  
  ### 1. Simulate stocks ==================================================####
  stock1 <- stock.sim_fast(GMZ1.init)
  stock2 <- stock.sim_fast(GMZ2.init)
  
  ### 2. Sample stocks  ===================================================####
  set.seed(7777777)
  
  #### Stock 1  ===========================================================####
  samples1 <- bind_rows(
    lapply(seq_along(stock1$caught_list), function(i) {
      sampler(
        harvest_df      = stock1$caught_list[[i]],  
        yr              = i,                                        
        sampling_window = sampling_window,
        sample_rate     = sample_rate          
      )
    })
  )
  
  #### Stock 2  ==========================================================####
  samples2 <- bind_rows(
    lapply(seq_along(stock2$caught_list), function(i) {
      sampler(
        harvest_df      = stock2$caught_list[[i]],  
        yr              = i,                                        
        sampling_window = sampling_window,
        sample_rate     = sample_rate          
      )
    })
  )
  
  #### Assign stock origins ===============================================####
  samples1 <- samples1 %>% mutate(origin = "GMZ1") 
  samples2 <- samples2 %>% mutate(origin = "GMZ2")
  
  ### 3. Mix samples and Run Iterations ===================================####
  results_list <- list()
  counter <- 1
  
  for (i in seq_along(samp.props)) {
    p <- samp.props[i]
    
    for (j in 1:niter.stock) {
      mixed_samp <- samples1 %>% 
        slice_sample(prop = p, replace = FALSE) %>% 
        rbind(slice_sample(samples2, prop = (1 - p), replace = FALSE))
      
      est_df <- CKMR_est(mixed_samp)
      
      est_df$proportion <- p 
      est_df$iteration  <- j
      
      results_list[[counter]] <- est_df
      counter <- counter + 1 
    }
  }
  
  ### 4. Process and Group Outputs ========================================####
  final_estimates <- bind_rows(results_list)
  
  averaged_estimates <- final_estimates %>%
    group_by(proportion, Actual_Year) %>%
    dplyr::summarise(
      avg_mean_Nadult = mean(mean, na.rm = TRUE),
      avg_lower_2.5   = mean(`2.5%`, na.rm = TRUE),
      avg_median_50   = mean(`50%`, na.rm = TRUE),
      avg_upper_97.5  = mean(`97.5%`, na.rm = TRUE),
      iters_counted   = dplyr::n(),
      rhat            = mean(Rhat, na.rm = TRUE),
      .groups         = "drop"
    )
  
  ### 5. Calculate true population abundance ==============================####
  years_to_estimate <- sort(unique(samples1$cohort))
  
  true_abundance <- lapply(seq_along(stock1$pop_list), function(i) {
    adults_stock1 <- stock1$pop_list[[i]] %>% 
      filter(age >= 5) %>% 
      nrow()
    
    data.frame(
      Actual_Year = i, 
      True_Nadult = adults_stock1
    )
  }) %>% 
    bind_rows()
  
  true_abundance <- true_abundance %>% 
    filter(Actual_Year %in% years_to_estimate)
  
  ### 6. Combine Truth and Estimates into One Data Frame ==================####
  comparison_df <- averaged_estimates %>%
    left_join(true_abundance, by = "Actual_Year")
  
  return(comparison_df)
}








## Run function ===========================================================####


stock.est.2 <- spatial_2stocks(GMZ1.init, GMZ2.init,samp.props)

stock.est.2.trimmed  <- stock.est.2 %>% 
  filter(rhat <=1.5)






# Plot individually =======================================================####
ggplot(stock.est.2.trimmed, aes(x = Actual_Year)) +
  geom_ribbon(aes(ymin = avg_lower_2.5, ymax = avg_upper_97.5), fill = "blue", alpha = 0.15) +
  geom_line(aes(y = avg_mean_Nadult, color = "CKMR Estimate"), size = 1) +
  geom_line(aes(y = True_Nadult, color = "True Population"), linetype = "dashed", size = 1) +
  facet_wrap(~ proportion, labeller = label_both) +
  labs(
    title = "CKMR Abundance Estimates vs True Baseline",
    x = "Model Year",
    y = "Adult Abundance (N)",
    color = "Legend"
  ) +
  theme_minimal()


# Plot as lines ===========================================================####
ggplot(stock.est.2.trimmed, aes(x = Actual_Year)) +
  geom_line(aes(y = avg_mean_Nadult, group = proportion, color = as.factor(proportion)), alpha = 0.7) +
  geom_line(aes(y = True_Nadult), color = "black", size = 1.2, linetype = "solid") +
  scale_color_viridis_d(name = "Sampling Prop (GMZ1)") +
  labs(
    title = "Influence of Mixing Bias on Population Recovery Estimates",
    subtitle = "Solid black line represents True Stock 1 Abundance",
    x = "Model Year",
    y = "Adult Abundance (N)"
  ) +
  theme_minimal()


# Violin plot of relative error ===========================================####
stock.est.2_violin <- stock.est.2.trimmed %>% 
  mutate(rel.bias = ((avg_mean_Nadult-True_Nadult)/True_Nadult)) 

ggplot(stock.est.2_violin, aes(x = as.factor(proportion), y = rel.bias, fill = as.factor(proportion))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  scale_fill_viridis_d() +
  labs(
    title = "Relative Bias by Sampling Proportion",
    x = "Sampling Proportion (p)",
    y = "Relative Bias",
    fill = "Proportion of GMZ1 Sampled"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# True abundance using non-mixed ==========================================####
no.mix.samp <- stock.est.2.trimmed %>% 
  filter(proportion == "1")

ggplot(no.mix.samp,aes(x = Actual_Year)) +
  geom_line(aes(y = avg_mean_Nadult ), colour = "red") +
  geom_ribbon(aes(ymin = avg_lower_2.5 , ymax = avg_upper_97.5), fill = "pink", alpha = 0.5) +
  geom_line(aes(y = True_Nadult), colour = "darkblue", linewidth = 1.1)





















