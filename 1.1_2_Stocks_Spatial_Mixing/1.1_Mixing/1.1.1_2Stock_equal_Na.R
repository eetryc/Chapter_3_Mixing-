#==================================##
#======= Two Stock Mixing =======####
#==================================##

# This example is simply two equal sized stocks that are moving between 
# spatial management zones without engaging in natal straying (thus, no genetic
# mixture)


## 1. Simulate stocks =====================================================####
 
  stock1 <- stock.sim_fast(1000)
  
  stock1 <- stock.sim_fast(1000)

  
  
  


## 2. Sample stocks  ======================================================####

  ### Stock 1  ============================================================####
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
  
  ### Stock 2  ============================================================####
  
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

  
  ### Assign stock values =================================================####
  
  # Stick with Genetic Management Zone (GMZ) nomenclature 
  
  samples1 <- samples1 %>% 
    mutate(origin = "GMZ1" ) 
  
  samples2 <- samples2 %>% 
    mutate(origin = "GMZ2" )
  
  
  
  
## 3. Mix samples =========================================================####
  
  
  
  


