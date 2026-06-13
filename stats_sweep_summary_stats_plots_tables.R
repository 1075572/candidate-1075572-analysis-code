# This code produces summary statistics and plots for sweep data
# it reads all sweep result csv files, converts voltage to velocity for the time domain,
# calculates summary statistics and produces a table and line plots
# it plots:
# average peak amplitude of all 8 sweeps vs distance in the time domain
# each of the 8 individual sweep peaks vs distance
# peak frequency vs distance (fft)
# and noise vs distance (in the time domain)
# excludes stinger data as it is not part of the distance continuum
# exclude stinger, Gel 2 and 400um for consistency with Friedman analysis on sweeps time amplitude peak data

library(tidyverse)
library(DT)
library(gridExtra)

# set folder path
sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"

# distance order for x-axis
dist_order <- c("0um", "50um", "100um", "200um", "300um", "600um", "800um", "1000um")

# pulls out distance and gel number from each filename
parse_filename <- function(fname) {
  distance <- str_extract(fname, "(?<=SWEEP_)[^_]+")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  list(distance = distance, gel = gel)
}

# reads all sweep csv files and extracts the averaged time domain amplitude,
# fft peak frequency and noise for each file
load_sweep_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  map_dfr(files, function(f) {
    info <- parse_filename(basename(f))
    raw  <- read.csv(f, header = FALSE)
    colnames(raw) <- c("td_peak_amp", "time", "fft_peak_amp", "fft_peak_freq", "noise")
    
    laser_sens      <- 5
    raw$td_peak_amp <- raw$td_peak_amp * laser_sens
    raw$noise[1]    <- raw$noise[1] * laser_sens
    
    td_mean  <- mean(abs(raw$td_peak_amp), na.rm = TRUE)
    fft_data <- raw %>% filter(fft_peak_freq > 0)
    dom_row  <- fft_data[which.max(fft_data$fft_peak_amp), ]
    noise    <- raw$noise[1]
    
    tibble(
      distance      = info$distance,
      gel           = info$gel,
      td_peak_amp   = td_mean,
      fft_peak_freq = dom_row$fft_peak_freq,
      noise         = noise
    )
  })
}

# reads all sweep csv files and extracts each of the 8 individual peaks separately
# excludes Gel 2 (noisy from 300um onwards) and 400um (incomplete data for Gel 3)
load_sweep_individual_peaks <- function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  # exclude Gel 2 and 400um entirely for individual peak analysis
  
  files <- files[!grepl("Gel2", files)]
  files <- files[!grepl("400um", files)]
  files <- files[!grepl("stinger", files)]   # add this line
  
  map_dfr(files, function(f) {
    info <- parse_filename(basename(f))
    raw  <- read.csv(f, header = FALSE)
    colnames(raw) <- c("td_peak_amp", "time", "fft_peak_amp", "fft_peak_freq", "noise")
    
    laser_sens      <- 5
    raw$td_peak_amp <- raw$td_peak_amp * laser_sens
    
    # remove zero padding, cap at 8 peaks and number each peak
    peaks <- abs(raw$td_peak_amp[raw$td_peak_amp != 0])
    peaks <- peaks[1:min(8, length(peaks))]
    
    map_dfr(seq_along(peaks), function(i) {
      tibble(
        distance  = info$distance,
        gel       = info$gel,
        peak_num  = i,
        amplitude = peaks[i]
      )
    })
  })
}

# calculates summary statistics per distance for averaged data
summarise_sweep <- function(df) {
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
      noise_min    = min(noise,    na.rm = TRUE),
      noise_max    = max(noise,    na.rm = TRUE),
      noise_range  = noise_max - noise_min,
      
      .groups = "drop"
    ) %>%
    arrange(distance)
}

# calculates summary statistics per distance for each individual peak
summarise_individual_peaks <- function(df) {
  df$distance <- factor(df$distance, levels = dist_order)
  df %>%
    group_by(peak_num, distance) %>%
    summarise(
      n        = n(),
      mean_amp = mean(amplitude,   na.rm = TRUE),
      median   = median(amplitude, na.rm = TRUE),
      sd       = sd(amplitude,     na.rm = TRUE),
      sem      = sd(amplitude, na.rm = TRUE) / sqrt(sum(!is.na(amplitude))),
      min      = min(amplitude,    na.rm = TRUE),
      max      = max(amplitude,    na.rm = TRUE),
      range    = max - min,
      .groups  = "drop"
    ) %>%
    arrange(peak_num, distance)
}

# plot theme
theme_clean <- theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

# line plot function — mean with sem error bars
make_lineplot <- function(df, y_var, y_label, title) {
  df$distance <- factor(df$distance, levels = dist_order)
  sum_df <- df %>%
    group_by(distance) %>%
    summarise(
      mean_val = mean(.data[[y_var]], na.rm = TRUE),
      sem_val  = sd(.data[[y_var]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[y_var]]))),
      .groups  = "drop"
    )
  ggplot(sum_df, aes(x = distance, y = mean_val, group = 1)) +
    geom_line(colour = "#AD6030", linewidth = 0.8) +
    geom_point(colour = "#AD6030", size = 3) +
    geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val),
                  width = 0.2, colour = "#AD6030") +
    labs(title = title, x = "Distance", y = y_label) +
    theme_clean
}


# load data and exclude stinger, Gel 2 and 400um for consistency with Friedman analysis
sweep_raw          <- load_sweep_files(sweep_folder)
sweep_raw          <- sweep_raw[sweep_raw$distance != "stinger", ]
sweep_raw          <- sweep_raw[sweep_raw$distance != "400um", ]
sweep_raw          <- sweep_raw[sweep_raw$gel != 2, ]

# check which files are failing to parse
files <- list.files(sweep_folder, pattern = "\\.csv$", full.names = TRUE)
bad <- files[is.na(str_extract(basename(files), "(?<=SWEEP_)[^_]+"))]
print(bad)


sweep_summary      <- summarise_sweep(sweep_raw)
sweep_peak_summary <- summarise_individual_peaks(sweep_peaks_long)


# summarise
sweep_peaks_long   <- load_sweep_individual_peaks(sweep_folder)
sweep_peaks_long   <- sweep_peaks_long[!is.na(sweep_peaks_long$distance), ]  # this first
dist_order_ind     <- dist_order[dist_order != "400um"]                       # then this
sweep_peaks_long$distance <- factor(sweep_peaks_long$distance, levels = dist_order_ind)  # then this

# print and display tables
cat("\nSweep summary statistics (averaged peaks):\n")
print(sweep_summary, width = Inf)
datatable(sweep_summary)

cat("\nSweep summary statistics (individual peaks):\n")
print(sweep_peak_summary, width = Inf)
datatable(sweep_peak_summary)

# add a label column to each table so you can tell them apart in the combined csv
sweep_summary$table       <- "averaged peaks"
sweep_peak_summary$table  <- "individual peaks"

# combine both tables into one csv — they have different columns so missing values
# will be filled with NA where columns dont overlap
combined_summary <- bind_rows(sweep_summary, sweep_peak_summary)

# save as two separate sheets 
write.csv(sweep_summary %>% select(-table),
          "/Users/kshamasaju/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/Summary stats/sweep_averaged_summary_stats.csv",
          row.names = FALSE)

write.csv(sweep_peak_summary %>% select(-table),
          "/Users/kshamasaju/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/Summary stats/sweep_individual_peak_summary_stats.csv",
          row.names = FALSE)
cat("Summary stats saved\n")

# line plots for averaged data
print(make_lineplot(sweep_raw, "td_peak_amp",
                    "Mean peak amplitude ± SEM (mm/s)",
                    "Sweep: Average peak amplitude (time domain) vs distance"))

print(make_lineplot(sweep_raw, "fft_peak_freq",
                    "Mean peak frequency ± SEM (Hz)",
                    "Sweep: Peak frequency (FFT) vs distance"))


# noise bar chart
print(
  ggplot(sweep_raw %>%
           group_by(distance) %>%
           summarise(mean_noise = mean(noise, na.rm = TRUE),
                     sem_noise  = sd(noise, na.rm = TRUE) / sqrt(sum(!is.na(noise))),
                     .groups = "drop") %>%
           mutate(distance = factor(distance, levels = dist_order)),
         aes(x = distance, y = mean_noise)) +
    geom_bar(stat = "identity", fill = "#AD6030", colour = "#AD6030", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_noise - sem_noise, ymax = mean_noise + sem_noise),
                  width = 0.2, colour = "black", linewidth = 0.5) +
    labs(title = "Sweep: noise vs distance", x = "Distance", y = "Mean noise ± SEM (mm/s)") +
    theme_clean
)

# 8 individual peak line plots arranged 4 on top 4 on bottom
# set distance order excluding 400um for individual peak plots
dist_order_ind <- dist_order[dist_order != "400um"]
sweep_peaks_long$distance <- factor(sweep_peaks_long$distance, levels = dist_order_ind)

library(cowplot)

individual_plots <- list()
for (pk in 1:8) {
  peak_data <- sweep_peaks_long %>% filter(peak_num == pk)
  if (nrow(peak_data) == 0) next
  
  sum_pk <- peak_data %>%
    group_by(distance) %>%
    summarise(
      mean_val = mean(amplitude, na.rm = TRUE),
      sem_val  = sd(amplitude, na.rm = TRUE) / sqrt(sum(!is.na(amplitude))),
      .groups  = "drop"
    )
  
  # only bottom row shows x tick labels, only left column shows y tick labels
  is_bottom <- pk %in% 5:8
  is_left   <- pk %in% c(1, 5)
  
  individual_plots[[pk]] <- ggplot(sum_pk, aes(x = distance, y = mean_val, group = 1)) +
    geom_line(colour = "#AD6030", linewidth = 0.7) +
    geom_point(colour = "#AD6030", size = 2) +
    geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val),
                  width = 0.2, colour = "#AD6030") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_clean +
    theme(
      axis.text.x  = if (is_bottom) element_text(angle = 45, hjust = 1, size = 7) else element_blank(),
      axis.ticks.x = if (is_bottom) element_line() else element_blank(),
      axis.text.y  = if (is_left) element_text(size = 7) else element_blank(),
      axis.ticks.y = if (is_left) element_line() else element_blank(),
      plot.margin  = margin(2, 2, 2, 2)
    )
}

valid_plots <- unname(individual_plots[!sapply(individual_plots, is.null)])

# build grid
pg <- plot_grid(plotlist = valid_plots, nrow = 2, ncol = 4, align = "hv", axis = "tblr")

# shared labels
y_label    <- ggdraw() + draw_label("Amplitude ± SEM (mm/s)", angle = 90, size = 15)
x_label <- ggdraw() + draw_label("Distance", size = 15, x = 0.53)
title_draw <- ggdraw() + draw_label("Sweep: individual peak amplitude vs distance",
                                    fontface = "plain", size = 12)

plot_with_y <- plot_grid(y_label, pg,        ncol = 2, rel_widths  = c(0.04, 1))
plot_with_x <- plot_grid(plot_with_y, x_label, nrow = 2, rel_heights = c(1, 0.04))
final_plot  <- plot_grid(title_draw, plot_with_x, nrow = 2, rel_heights = c(0.04, 1))

pg <- plot_grid(plotlist = valid_plots, nrow = 2, ncol = 4, align = "hv", axis = "tblr",
                rel_heights = c(1, 1))
options(repr.plot.width = 17, repr.plot.height = 7)
print(final_plot)

# plot all 8 individual peaks on one graph with the average overlaid
# each peak is a different colour and the average is thicker and black to stand out

# colour palette for the 8 individual peaks
peak_colours <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                  "#9467bd","#8c564b","#e377c2","#bcbd22")

# summarise each individual peak for plotting
all_peaks_summary <- sweep_peaks_long %>%
  group_by(peak_num, distance) %>%
  summarise(
    mean_val = mean(amplitude, na.rm = TRUE),
    sem_val  = sd(amplitude, na.rm = TRUE) / sqrt(sum(!is.na(amplitude))),
    .groups  = "drop"
  ) %>%
  mutate(peak_label = paste("Peak", peak_num))

# summarise the average across all 8 peaks for the overlay line
average_summary <- sweep_raw %>%
  group_by(distance) %>%
  summarise(
    mean_val = mean(td_peak_amp, na.rm = TRUE),
    sem_val  = sd(td_peak_amp, na.rm = TRUE) / sqrt(sum(!is.na(td_peak_amp))),
    .groups  = "drop"
  ) %>%
  mutate(peak_label = "Average")

# build the combined plot
print(
  ggplot() +
    # individual peak lines — one colour per peak, thinner lines
    geom_line(data = all_peaks_summary,
              aes(x = distance, y = mean_val,
                  group = peak_label, colour = peak_label),
              linewidth = 0.7, alpha = 0.8) +
    geom_point(data = all_peaks_summary,
               aes(x = distance, y = mean_val,
                   group = peak_label, colour = peak_label),
               size = 1.5, alpha = 0.8) +
 
    
    # average line — black, thicker, on top to stand out
    geom_line(data = average_summary,
              aes(x = distance, y = mean_val, group = 1),
              colour = "black", linewidth = 1) +
    geom_point(data = average_summary,
               aes(x = distance, y = mean_val, group = 1),
               colour = "black", size = 3) +
   
    
    # colours for the 8 individual peaks
    scale_colour_manual(values = setNames(peak_colours, paste("Peak", 1:8))) +
    
    labs(title  = "Sweep: all individual peaks and average amplitude vs distance",
         x      = "Distance",
         y      = "Mean amplitude ± SEM (mm/s)",
         colour = "Peak") +
    theme_clean +
    theme(legend.position = "right")
)
