# ==============================================================================
# POWER ANALYSIS: Sensitivity, N-10%, N-15%, and N-20%
# Published in a paper: Digital Biomass, NDVI, NPCI
# ==============================================================================

# --- 1. SETUP & LIBRARIES ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(quantreg) 

# --- 2. GENERAL INPUTS ---
# Just place the CSV in the same folder as this script
input_file <- "Merged_dataset.csv"  
traits     <- list(
  list("Digital Biomass", "digital_biomass_mm3", "biomass_g"),
  list("NDVI",            "ndvi_average",        "avg_Chla"),
  list("NPCI",            "npci_average",        "avg_Chla")
)

# Load and basic clean
df0 <- read.csv(input_file)
df_clean <- df0 %>% filter(unique_id != "NN-4")

# --- 3. POWER CALCULATION FUNCTION ---
# Using Quantile Regression to model noise based on biological variance
calculate_pair_stats <- function(label, col_x, col_y, data) {
  clean_data <- na.omit(data.frame(x = data[[col_x]], y = data[[col_y]]))
  
  # A. Noise Model
  fit_med <- rq(x ~ y, tau = 0.5,  data = clean_data)
  fit_low <- rq(x ~ y, tau = 0.16, data = clean_data)
  fit_upp <- rq(x ~ y, tau = 0.84, data = clean_data)
  
  avg_y    <- median(clean_data$y)
  pred_med <- predict(fit_med, newdata = data.frame(y = avg_y))
  noise_sd <- abs(predict(fit_upp, newdata = data.frame(y = avg_y)) - 
                    predict(fit_low, newdata = data.frame(y = avg_y))) / 2
  
  # B. Simulation Helper
  get_power <- function(n_samp, pct_diff) {
    n_sims <- 400 
    y_shift <- avg_y * (1 + pct_diff/100)
    pred_shift <- predict(fit_med, newdata = data.frame(y = y_shift))
    
    ctrl <- matrix(rnorm(n_sims * n_samp, mean = pred_med,   sd = noise_sd), ncol = n_samp)
    trtm <- matrix(rnorm(n_sims * n_samp, mean = pred_shift, sd = noise_sd), ncol = n_samp)
    
    p_vals <- sapply(1:n_sims, function(i) t.test(ctrl[i,], trtm[i,])$p.value)
    return(mean(p_vals < 0.05))
  }
  
  # C. Running Repeats for Range/Error Bars
  n_repeats <- 30
  res_mde <- numeric(n_repeats); res_n10 <- numeric(n_repeats)
  res_n15 <- numeric(n_repeats); res_n20 <- numeric(n_repeats)
  
  message("Processing: ", label)
  for (i in 1:n_repeats) {
    for (eff in 1:100) { if (get_power(4, eff) >= 0.80) { res_mde[i] <- eff; break } }
    for (n in 3:60) { if (get_power(n, 10) >= 0.80) { res_n10[i] <- n; break } }
    for (n in 3:60) { if (get_power(n, 15) >= 0.80) { res_n15[i] <- n; break } }
    for (n in 3:60) { if (get_power(n, 20) >= 0.80) { res_n20[i] <- n; break } }
  }
  
  return(data.frame(
    Trait = label,
    MDE_Mean = mean(res_mde), MDE_Min = min(res_mde), MDE_Max = max(res_mde),
    N10_Mean = mean(res_n10), N10_Min = min(res_n10), N10_Max = max(res_n10),
    N15_Mean = mean(res_n15), N15_Min = min(res_n15), N15_Max = max(res_n15),
    N20_Mean = mean(res_n20), N20_Min = min(res_n20), N20_Max = max(res_n20)
  ))
}

# --- 4. EXECUTION ---
results_list <- lapply(traits, function(p) calculate_pair_stats(p[[1]], p[[2]], p[[3]], df_clean))
summary_df <- bind_rows(results_list)

# --- 5. PLOTTING ---
my_theme <- theme_minimal(base_size = 14) + theme(panel.grid.minor = element_blank())

p_plot <- function(y_var, title, y_lab, h_line = NULL) {
  p <- ggplot(summary_df, aes(x = Trait, y = !!sym(y_var), group = 1)) +
    geom_ribbon(aes(ymin = !!sym(gsub("Mean", "Min", y_var)), 
                    ymax = !!sym(gsub("Mean", "Max", y_var))), fill = "grey80", alpha = 0.5) +
    geom_line(color = "black", size = 1) +
    geom_point(size = 3) +
    labs(title = title, x = "", y = y_lab) +
    my_theme
  if(!is.null(h_line)) p <- p + geom_hline(yintercept = h_line, linetype = "dashed", color = "red")
  return(p)
}

p1 <- p_plot("MDE_Mean", "Sensitivity (n=4)", "% Difference")
p2 <- p_plot("N10_Mean", "N for 10% Diff", "Plants", 4)
p3 <- p_plot("N15_Mean", "N for 15% Diff", "Plants", 4)
p4 <- p_plot("N20_Mean", "N for 20% Diff", "Plants", 4)

final_panel <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# --- 6. EXPORT ---
ggsave("Power_Analysis_Plot.png", final_panel, width = 11, height = 8.5, dpi = 300)
write.csv(summary_df, "Power_Analysis_Results.csv", row.names = FALSE)

message("Analysis Finished. Results saved to your current folder.")