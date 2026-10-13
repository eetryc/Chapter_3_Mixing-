# Load in libraries #==============================================#
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

stock_sim_RDS <- readRDS("C:/Users/eetry/Downloads/stock_sim_output_2 (3).rds")


nhat <- stock_sim_RDS$test

filter <- nhat |> 
  filter(Actual_Year != "Weighted_Focal_100") |> 
  mutate(Actual_Year = as.integer(Actual_Year))

Weighted_Focal_100 <- nhat |> 
  filter(Actual_Year == "Weighted_Focal_100") 

Weighted_Focal_100 <- Weighted_Focal_100$mean

stock3 <- stock_sim_RDS$true_abundance.3|> 
  filter(Actual_Year >=80)

full_lake_abundance_target <- stock_sim_RDS$full_lake_abundance_target|> 
  filter(Actual_Year >=80)





desired_order <- c("True Adult Abundance", "Annual Model Estimate", "Weighted Mean Estimate")



ggplot(filter, aes(x = Actual_Year)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#50C878", alpha = 0.2) +
  
  # Wrap the legend labels in factor() with explicit levels
  geom_line(aes(y = mean, 
                color = factor("Annual Model Estimate", levels = desired_order), 
                linetype = factor("Annual Model Estimate", levels = desired_order)), 
            linewidth = 1) +
  
  geom_line(data = stock3, 
            aes(x = Actual_Year, y = True_Nadult, 
                color = factor("True Adult Abundance", levels = desired_order), 
                linetype = factor("True Adult Abundance", levels = desired_order)), 
            linewidth = 1) +
  
  geom_hline(aes(yintercept = Weighted_Focal_100, 
                 color = factor("Weighted Mean Estimate", levels = desired_order), 
                 linetype = factor("Weighted Mean Estimate", levels = desired_order)), 
             linewidth = 1) +
  
  # Manual scales using the same named vector approach
  scale_color_manual(
    name = "Legend", 
    values = c(
      "Annual Model Estimate" = "#50C878", 
      "True Adult Abundance" = "darkblue", 
      "Weighted Mean Estimate" = "#008000"
    )
  ) +
  
  scale_linetype_manual(
    name = "Legend",
    values = c(
      "Annual Model Estimate" = 2, 
      "True Adult Abundance" = 1, 
      "Weighted Mean Estimate" = 3
    )
  ) +
  
  coord_cartesian(xlim = c(82, 113)) +
  theme_bw() +
  labs(x = "Year", y = "Adult Abundance") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) 





true_values <- stock3 |> 
  filter(Actual_Year <= 115) 
true_values <- true_values$True_Nadult


est_values <- filter |> 
  filter(Actual_Year >=80) 

est_values <- est_values$mean


relative_bias <-  100 * ((est_values - true_values)/ true_values)
  





upper_est_values <- filter |> 
  filter(Actual_Year >=80) 

upper_est_values <- upper_est_values$`97.5%`
  
  
upper_relative_bias <-  100 * ((upper_est_values - true_values)/ true_values)



lower_est_values <- filter |> 
  filter(Actual_Year >=80) 

lower_est_values <- lower_est_values$`2.5%`  

lower_relative_bias <-  100 * ((lower_est_values - true_values)/ true_values)


years <- filter |> 
  filter(Actual_Year >=80) 
  
years <-  years$Actual_Year


rel_bias_table <- tibble(years,relative_bias,upper_relative_bias,lower_relative_bias)

ggplot(rel_bias_table, aes(x = factor(years))) +
  geom_crossbar(
    aes(y = relative_bias, ymin = lower_relative_bias, ymax = upper_relative_bias),
    fill = "steelblue",
    alpha = 0.6,
    width = 0.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue", linewidth = 0.8) +
  theme_bw() +
  labs(x = "Year", y = "Relative Bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
       
