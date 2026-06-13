library(dplyr)
library(tidyr)

click_csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values clicks csv"
sweep_csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values sweeps csv"

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/rayleigh wave simulation figures"

# Setup
eta_R <- 0.9553^2
f_all <- c(113, 175, 350, 2200, 4000, 8700)

f_labels <- c(
  "Peak 1 ~113 Hz",
  "Peak 2 ~175 Hz",
  "Peak 3 ~350 Hz",
  "Peak 1 ~2200 Hz",
  "Peak 2 ~4000 Hz",
  "Peak 3 ~8700 Hz"
)

f_factor_levels <- f_labels

dist_order <- c("0um","50um","100um","200um","300um","400um",
                "600um","800um","1000um")

r_exp_um <- c(0, 50, 100, 200, 300, 400, 600, 800, 1000)

dist_to_num <- c(
  "0um"    =    0,
  "50um"   =   50,
  "100um"  =  100,
  "200um"  =  200,
  "300um"  =  300,
  "400um"  =  400,
  "600um"  =  600,
  "800um"  =  800,
  "1000um" = 1000
)

hankel0_abs <- function(x) {
  sqrt(besselJ(x, 0)^2 + besselY(x, 0)^2)
}

rho <- 1000
nu  <- 0.5

r_fine_um <- seq(5, 1000, by = 2)
r_fine    <- r_fine_um * 1e-6

vel_at_50um <- c(
  0.000338,
  0.000168,
  0.000407,
  0.000857,
  0.0000909,
  0.0000968
)

# Build model predictions at experimental distances only (E = 400 Pa)
df_model_at_exp <- data.frame()
for (fi in seq_along(f_all)) {
  E_i   <- 400
  G_i   <- E_i / (2 * (1 + nu))
  c_s_i <- sqrt(G_i / rho)
  c_R_i <- sqrt(eta_R) * c_s_i
  omega_i <- 2 * pi * f_all[fi]
  k_R_i   <- omega_i / c_R_i
  
  r_ref  <- 50e-6
  H_ref  <- hankel0_abs(k_R_i * r_ref)
  delta_0_scaled <- vel_at_50um[fi] / (2 * pi * f_all[fi] * 1000)
  A_norm <- delta_0_scaled / H_ref
  
  for (r_um in r_exp_um) {
    r_m <- r_um * 1e-6
    kRr <- k_R_i * r_m
    vel_exact <- A_norm * hankel0_abs(kRr) * 2 * pi * f_all[fi] * 1000
    
    df_model_at_exp <- rbind(df_model_at_exp, data.frame(
      r_um       = r_um,
      freq_label = factor(f_labels[fi], levels = f_factor_levels),
      vel_exact  = vel_exact
    ))
  }
}

# Read measured data
read_eqp_csv <- function(filepath, freq_label_val, exclude_400um = FALSE) {
  df <- read.csv(filepath, header = TRUE)
  df %>%
    filter(distance != "stinger") %>%
    { if (exclude_400um) filter(., distance != "400um") else . } %>%
    mutate(distance = factor(distance, levels = dist_order)) %>%
    select(distance, Gel1, Gel3, Gel4, Gel5) %>%
    pivot_longer(cols      = c(Gel1, Gel3, Gel4, Gel5),
                 names_to  = "gel",
                 values_to = "amplitude") %>%
    filter(!is.na(amplitude)) %>%
    group_by(distance) %>%
    summarise(
      mean_vel = mean(amplitude),
      sem_vel  = sd(amplitude) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    mutate(
      r_um       = dist_to_num[as.character(distance)],
      freq_label = factor(freq_label_val, levels = f_factor_levels)
    ) %>%
    filter(!is.na(r_um))
}

df_measured <- bind_rows(
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak1_values.csv"),
               "Peak 1 ~113 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak2_values.csv"),
               "Peak 2 ~175 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak3_values.csv"),
               "Peak 3 ~350 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak1_values.csv"),
               "Peak 1 ~2200 Hz", exclude_400um = TRUE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak2_values.csv"),
               "Peak 2 ~4000 Hz", exclude_400um = TRUE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak3_values.csv"),
               "Peak 3 ~8700 Hz", exclude_400um = TRUE)
)

# Calculate percent deviation
df_deviation <- df_measured %>%
  inner_join(df_model_at_exp, by = c("r_um", "freq_label")) %>%
  filter(r_um != 50) %>%
  mutate(
    deviation_pct = ((vel_exact - mean_vel) / vel_exact) * 100,
    stimulus = case_when(
      grepl("113|175|350",    as.character(freq_label)) ~ "Click",
      grepl("2200|4000|8700", as.character(freq_label)) ~ "Sweep",
      TRUE ~ NA_character_
    )
  )

# Print results
cat("\n--- Percent deviation of measured velocity from elastic Hankel model ---\n")
cat("Positive = measured decays faster than elastic model\n\n")
print(
  df_deviation %>%
    select(freq_label, r_um, mean_vel, vel_exact, deviation_pct) %>%
    arrange(freq_label, r_um),
  digits = 3
)

cat("\n--- Mean percent deviation per equivalent peak (100-1000 um) ---\n")
df_dev_summary <- df_deviation %>%
  group_by(freq_label) %>%
  summarise(
    mean_deviation_pct = mean(deviation_pct, na.rm = TRUE),
    sd_deviation_pct   = sd(deviation_pct,   na.rm = TRUE),
    n                  = n(),
    .groups = "drop"
  ) %>%
  arrange(freq_label)
print(df_dev_summary, digits = 3)

cat("\n--- Spearman correlation: mean deviation vs log10(frequency) ---\n")
df_spearman <- df_dev_summary %>%
  mutate(
    freq_hz  = c(113, 175, 350, 2200, 4000, 8700),
    log_freq = log10(freq_hz)
  )
sp_test <- cor.test(df_spearman$log_freq,
                    df_spearman$mean_deviation_pct,
                    method = "spearman")
print(sp_test)

# Save
write.csv(
  df_deviation %>% select(freq_label, r_um, mean_vel, vel_exact, deviation_pct),
  file.path(save_dir, "percent_deviation_from_elastic_model.csv"),
  row.names = FALSE
)
cat("\nSaved to:", file.path(save_dir, "percent_deviation_from_elastic_model.csv"), "\n")