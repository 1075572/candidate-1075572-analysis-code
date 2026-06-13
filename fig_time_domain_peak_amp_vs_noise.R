# This code creates 2 figures showing peak amplitude in the time domain 
# vs noise for clicks and sweeps
# converts the time-domain peak amplitude and noise floor from volts to mm/s
# removes the stinger position and Gel 2, then plots mean peak surface velocity ± SEM
# and mean background noise ± SEM against distance on a single graph with a legend

# creates 2 figures showing peak surface velocity vs noise for clicks and sweeps
# click figure includes 400um, sweep figure excludes it (incomplete Gel 3 data)

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(R.matlab)
library(ragg)
library(tidyverse)

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

laser_sens <- 5.0

# shared plot theme
stim_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    legend.title     = element_text(size = 11, family = "Helvetica"),
    legend.text      = element_text(size = 11, family = "Helvetica"),
    strip.text       = element_text(family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8),
    plot.tag         = element_text(family = "Helvetica", face = "bold", size = 20)
  )

# click uses blue/red, sweep uses orange tones
click_colours <- c("Peak surface velocity" = "#2E6BAD",
                   "Background noise (estimate)"      = "#636363")
sweep_colours <- c("Peak surface velocity" = "#AD6030",
                   "Background noise"      = "#636363")
series_lines  <- c("Peak surface velocity" = "solid",
                   "Background noise"      = "solid")

# clicks include 400um, sweeps exclude it
dist_order_click <- c("0um","50um","100um","200um","300um","400um","600um","800um","1000um")
dist_order_sweep <- c("0um","50um","100um","200um","300um","600um","800um","1000um")


# CLICK FUNCTIONS

# gets distance and gel from click filename
parse_click_filename <- function(fname) {
  distance <- str_extract(fname, "(?<=CLICK_)[^_]+")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  list(distance = distance, gel = gel)
}

# reads click csvs — 4 columns, takes row 1 amplitude and noise, converts to mm/s
load_click_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  map_dfr(files, function(f) {
    info <- parse_click_filename(basename(f))
    raw  <- read.csv(f, header = FALSE)
    colnames(raw) <- c("td_peak_amp", "fft_peak_amp", "fft_peak_freq", "noise")
    raw$td_peak_amp <- raw$td_peak_amp * laser_sens
    raw$noise[1]    <- raw$noise[1]    * laser_sens
    tibble(
      distance    = info$distance,
      gel         = info$gel,
      td_peak_amp = abs(raw$td_peak_amp[1]),
      noise       = raw$noise[1]
    )
  })
}


# SWEEP FUNCTIONS

# gets distance and gel from sweep filename
parse_sweep_filename <- function(fname) {
  distance <- str_extract(fname, "(?<=SWEEP_)[^_]+")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  list(distance = distance, gel = gel)
}

# reads sweep csvs — 5 columns, averages all 8 peaks, converts to mm/s
load_sweep_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  map_dfr(files, function(f) {
    info <- parse_sweep_filename(basename(f))
    raw  <- read.csv(f, header = FALSE)
    colnames(raw) <- c("td_peak_amp", "time", "fft_peak_amp", "fft_peak_freq", "noise")
    raw$td_peak_amp <- raw$td_peak_amp * laser_sens
    raw$noise[1]    <- raw$noise[1]    * laser_sens
    td_mean <- mean(abs(raw$td_peak_amp), na.rm = TRUE)
    tibble(
      distance    = info$distance,
      gel         = info$gel,
      td_peak_amp = td_mean,
      noise       = raw$noise[1]
    )
  })
}


# PLOTTING FUNCTION

# builds the combined amplitude + noise plot
# takes the data, distance order and colour palette as arguments
make_combined_plot <- function(data, dist_ord, colours) {
  
  amp_sum <- data %>%
    group_by(distance) %>%
    summarise(
      mean_val = mean(td_peak_amp, na.rm = TRUE),
      sem_val  = sd(td_peak_amp,   na.rm = TRUE) / sqrt(sum(!is.na(td_peak_amp))),
      .groups  = "drop"
    ) %>%
    mutate(dist_num = as.numeric(factor(distance, levels = dist_ord)),
           series   = "Peak surface velocity")
  
  noise_sum <- data %>%
    group_by(distance) %>%
    summarise(
      mean_val = mean(noise, na.rm = TRUE),
      sem_val  = sd(noise,   na.rm = TRUE) / sqrt(sum(!is.na(noise))),
      .groups  = "drop"
    ) %>%
    mutate(dist_num = as.numeric(factor(distance, levels = dist_ord)),
           series   = "Background noise")
  
  combined <- bind_rows(amp_sum, noise_sum)
  
  ggplot(combined, aes(x = dist_num, y = mean_val,
                       colour = series, linetype = series)) +
    geom_errorbar(aes(ymin = mean_val - sem_val,
                      ymax = mean_val + sem_val),
                  width = 0.2) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    scale_colour_manual(values = colours,      name = NULL) +
    scale_linetype_manual(values = series_lines, name = NULL) +
    scale_x_continuous(breaks = seq_along(dist_ord),
                       labels = dist_ord) +
    labs(x = "Distance", y = "Velocity (mm/s)") +
    stim_theme +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}


# fig 1: click - drop stinger and Gel 2

click_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all click results"

click_raw          <- load_click_files(click_folder)
click_raw          <- click_raw[click_raw$distance != "stinger", ]
click_raw          <- click_raw[click_raw$gel      != 2,         ]
click_raw$distance <- factor(click_raw$distance, levels = dist_order_click)

p_click <- make_combined_plot(click_raw, dist_order_click, click_colours)

agg_png(file.path(save_dir, "click_td_peak_amp_and_noise.png"),
        width = 8, height = 6, units = "in", res = 200)
print(p_click)
dev.off()
cat("Saved click_td_peak_amp_and_noise.png\n")


#fig 2: sweep, drop stinger, Gel 2 and 400um

sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"

sweep_raw          <- load_sweep_files(sweep_folder)
sweep_raw          <- sweep_raw[sweep_raw$distance != "stinger", ]
sweep_raw          <- sweep_raw[sweep_raw$distance != "400um",   ]
sweep_raw          <- sweep_raw[sweep_raw$gel      != 2,         ]
sweep_raw$distance <- factor(sweep_raw$distance, levels = dist_order_sweep)

p_sweep <- make_combined_plot(sweep_raw, dist_order_sweep, sweep_colours)

agg_png(file.path(save_dir, "sweep_td_peak_amp_and_noise.png"),
        width = 8, height = 6, units = "in", res = 200)
print(p_sweep)
dev.off()
cat("Saved sweep_td_peak_amp_and_noise.png\n")