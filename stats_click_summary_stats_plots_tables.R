library(tidyverse)

# set folder path 
click_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all click results"

# distance order for x-axis
dist_order <- c("0um", "50um", "100um", "200um", "300um", "400um", "600um", "800um", "1000um")

# pulls out distance and gel number from each filename
parse_filename <- function(fname) {
  distance <- str_extract(fname, "(?<=CLICK_)[^_]+")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  list(distance = distance, gel = gel)
}


# reads all click csv files from the folder
# the laser records in volts so we multiply by 5 to convert to mm/s
# only the time domain amplitude and noise need converting — fft values do not
load_click_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  map_dfr(files, function(f) {
    info <- parse_filename(basename(f))
    raw  <- read.csv(f, header = FALSE)
    colnames(raw) <- c("td_peak_amp", "fft_peak_amp", "fft_peak_freq", "noise")
    
    # convert time domain and noise from volts to mm/s
    laser_sens      <- 5
    raw$td_peak_amp <- raw$td_peak_amp * laser_sens
    raw$noise[1]    <- raw$noise[1] * laser_sens
    
    # time domain peak amplitude is in row 1, take absolute value
    td_peak_amp <- abs(raw$td_peak_amp[1])
    
    # for fft, remove zero-padded rows and find the dominant peak
    fft_data <- raw %>% filter(fft_peak_freq > 0)
    dom_row  <- fft_data[which.max(fft_data$fft_peak_amp), ]
    
    # noise is in row 1 only
    noise <- raw$noise[1]
    
    tibble(
      distance      = info$distance,
      gel           = info$gel,
      td_peak_amp   = td_peak_amp,
      fft_peak_freq = dom_row$fft_peak_freq,
      noise         = noise
    )
  })
}

# calculates summary statistics per distance
# mean, median, sd, sem, min, max and range for each variable
summarise_click <- function(df) {
  df$distance <- factor(df$distance, levels = dist_order)
  df %>%
    group_by(distance) %>%
    summarise(
      n = n(),
      
      td_pk_amp_mean   = mean(td_peak_amp,   na.rm = TRUE),
      td_pk_amp_median = median(td_peak_amp, na.rm = TRUE),
      td_pk_amp_sd     = sd(td_peak_amp,     na.rm = TRUE),
      td_pk_amp_sem    = sd(td_peak_amp, na.rm = TRUE) / sqrt(sum(!is.na(td_peak_amp))),
      td_pk_amp_min    = min(td_peak_amp,    na.rm = TRUE),
      td_pk_amp_max    = max(td_peak_amp,    na.rm = TRUE),
      td_pk_amp_range  = td_pk_amp_max - td_pk_amp_min,
      
      fft_pk_freq_mean   = mean(fft_peak_freq,   na.rm = TRUE),
      fft_pk_freq_median = median(fft_peak_freq, na.rm = TRUE),
      fft_pk_freq_sd     = sd(fft_peak_freq,     na.rm = TRUE),
      fft_pk_freq_min    = min(fft_peak_freq,    na.rm = TRUE),
      fft_pk_freq_max    = max(fft_peak_freq,    na.rm = TRUE),
      fft_pk_freq_range  = fft_pk_freq_max - fft_pk_freq_min,
      
      noise_mean   = mean(noise,   na.rm = TRUE),
      noise_median = median(noise, na.rm = TRUE),
      noise_sd     = sd(noise,     na.rm = TRUE),
      noise_sem    = sd(noise, na.rm = TRUE) / sqrt(sum(!is.na(noise))),
      noise_min    = min(noise,    na.rm = TRUE),
      noise_max    = max(noise,    na.rm = TRUE),
      noise_range  = noise_max - noise_min,
      
      .groups = "drop"
    ) %>%
    arrange(distance)
}

# plot theme — clean white background with angled x axis labels
theme_clean <- theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

# line plot function — shows mean with sem error bars across distances
make_lineplot <- function(df, y_var, y_label, title) {
  df$distance <- factor(df$distance, levels = dist_order)
  sum_df <- df %>%
    group_by(distance) %>%
    summarise(
      mean_val = mean(.data[[y_var]], na.rm = TRUE),
      sd_val   = sd(.data[[y_var]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[y_var]]))),
      .groups  = "drop"
    )
  ggplot(sum_df, aes(x = distance, y = mean_val, group = 1)) +
    geom_line(colour = "#2E6BAD", linewidth = 0.8) +
    geom_point(colour = "#2E6BAD", size = 3) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  width = 0.2, colour = "#2E6BAD") +
    labs(title = title, x = "Distance", y = y_label) +
    theme_clean
}

# load data, exclude stinger and Gel 2 for consistency with Friedman analysis
click_raw     <- load_click_files(click_folder)
click_raw     <- click_raw[click_raw$distance != "stinger", ]
click_raw     <- click_raw[click_raw$gel != 2, ]
click_summary <- summarise_click(click_raw)

cat("\nClick summary statistics:\n")
print(click_summary, width = Inf)

# ── Publishable supplementary table ──────────────────────────────────────────
pub_table <- click_summary %>%
  select(distance, n,
         td_pk_amp_mean, td_pk_amp_sd, td_pk_amp_sem,
         td_pk_amp_median, td_pk_amp_min, td_pk_amp_max,
         fft_pk_freq_mean, fft_pk_freq_sd,
         fft_pk_freq_median, fft_pk_freq_min, fft_pk_freq_max,
         noise_mean, noise_sd, noise_sem,
         noise_median, noise_min, noise_max) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  # transpose: distances become columns, statistics become rows
  pivot_longer(-distance, names_to = "Statistic", values_to = "value") %>%
  pivot_wider(names_from = distance, values_from = value) %>%
  # replace internal names with readable labels
  mutate(Statistic = recode(Statistic,
                            "n"                  = "n",
                            "td_pk_amp_mean"     = "TD Peak Amp – Mean (mm/s)",
                            "td_pk_amp_sd"       = "TD Peak Amp – SD",
                            "td_pk_amp_sem"      = "TD Peak Amp – SEM",
                            "td_pk_amp_median"   = "TD Peak Amp – Median",
                            "td_pk_amp_min"      = "TD Peak Amp – Min",
                            "td_pk_amp_max"      = "TD Peak Amp – Max",
                            "fft_pk_freq_mean"   = "FFT Peak Freq – Mean (Hz)",
                            "fft_pk_freq_sd"     = "FFT Peak Freq – SD",
                            "fft_pk_freq_median" = "FFT Peak Freq – Median",
                            "fft_pk_freq_min"    = "FFT Peak Freq – Min",
                            "fft_pk_freq_max"    = "FFT Peak Freq – Max",
                            "noise_mean"         = "Noise – Mean (mm/s)",
                            "noise_sd"           = "Noise – SD",
                            "noise_sem"          = "Noise – SEM",
                            "noise_median"       = "Noise – Median",
                            "noise_min"          = "Noise – Min",
                            "noise_max"          = "Noise – Max"
  ))

write.csv(
  pub_table,
  "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary/click_summary_stats.csv",
  row.names = FALSE,
  na        = ""
)
cat("Saved: click_summary_stats.csv\n")
# display interactive summary stats table
library(DT)
datatable(click_summary)

# noise bar chart with sem error bars
print(
  ggplot(click_raw %>%
           group_by(distance) %>%
           summarise(mean_noise = mean(noise, na.rm = TRUE),
                     sd_noise = sd(noise, na.rm = TRUE) / sqrt(sum(!is.na(noise))), .groups = "drop") %>%
           mutate(distance = factor(distance, levels = dist_order)),
         aes(x = distance, y = mean_noise)) +
    geom_bar(stat = "identity", fill = "#C2DEFF", colour = "#2E6BAD") +
    geom_errorbar(aes(ymin = mean_noise - sd_noise, ymax = mean_noise + sd_noise),
                  width = 0.2, colour = "#2E6BAD") +
    labs(title = "Click: noise vs distance", x = "Distance", y = "Noise (mm/s)") +
    theme_clean
)

# line plots
print(make_lineplot(click_raw, "td_peak_amp",
                    "Peak amplitude (mm/s)",
                    "Click: Peak amplitude (time domain) vs distance"))

print(make_lineplot(click_raw, "fft_peak_freq",
                    "Peak frequency (Hz)",
                    "Click: Peak frequency (FFT) vs distance"))
