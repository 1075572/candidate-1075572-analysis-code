# The code constructs the click stimulus waveform from first principles and reads both Gel 1 
# stinger .mat files to extract stimulus voltage and laser response columns, converts 
# responses to mm/s using the sensitivity factor, and saves two figures 

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(R.matlab)
library(ragg)

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters 
fs            <- 100000
Vclick        <- 0.0136
t_prestim     <- 0.2
t_poststim    <- 2.0
t_stimS_click <- 0.00010
laser_sens    <- 5.0
ds            <- 100



# Shared theme 
stim_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8),
    plot.tag         = element_text(family = "Helvetica", face = "bold", size = 20)
  )



# FIGURE 1: CLICK
n_pre         <- round(t_prestim     * fs)
n_click       <- round(t_stimS_click * fs)
n_post        <- round(t_poststim    * fs)
n_total_click <- n_pre + n_click + n_post

click_signal  <- c(rep(0, n_pre), rep(Vclick, n_click), rep(0, n_post))
click_time    <- seq(1/fs, n_total_click/fs, by = 1/fs)

df_click      <- tibble(time = click_time, voltage = click_signal)
df_click_plot <- df_click[seq(1, nrow(df_click), by = ds), ]

p_click_stim <- ggplot(df_click_plot, aes(x = time, y = voltage)) +
  geom_line(colour = "#2E6BAD", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 2.5, by = 0.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, Vclick * 1.4), expand = c(0, 0)) +
  labs(x = "Time (s)", y = "Voltage (V)", tag = "(A)") +
  stim_theme

click_mat  <- readMat("~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Gels 1-5/Gel 1/clicks dataout gel 1/dataout_26_04_08_12_57_47_CLICK_stinger_Gel1.mat")
click_raw  <- click_mat$dataout
click_data <- click_raw[3:nrow(click_raw), ]
valid_c    <- !is.nan(click_data[, 1])
click_data <- click_data[valid_c, ]

df_click_resp      <- tibble(
  time     = click_data[, 1],
  velocity = click_data[, 2] * laser_sens
)
df_click_resp_plot <- df_click_resp

p_click_resp <- ggplot(df_click_resp_plot, aes(x = time, y = velocity)) +
  geom_line(colour = "#2E6BAD", linewidth = 0.4, alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, ceiling(max(df_click_resp$time)), by = 0.5),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "Time (s)", y = "Velocity (mm/s)", tag = "(B)") +
  stim_theme

fig1 <- arrangeGrob(
  p_click_stim, p_click_resp,
  ncol = 1
)

agg_png(file.path(save_dir, "click_stimulus_and_response_signal.png"),
        width = 10, height = 10, units = "in", res = 200)
grid.draw(fig1)
dev.off()
cat("Saved click_stimulus_and_response_signal.png\n")

# FIGURE 2: SWEEP
stinger      <- readMat("~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Gels 1-5/Gel 1/sine sweep dataout gel 1/dataout_26_04_08_11_50_42_SWEEP_stinger_Gel1.mat")
stinger_data <- stinger$dataout[3:nrow(stinger$dataout), ]
valid_s      <- !is.nan(stinger_data[, 1])
stinger_data <- stinger_data[valid_s, ]

df_sweep_stim      <- tibble(
  time    = stinger_data[, 1],
  voltage = stinger_data[, 2]
)
df_sweep_stim_plot <- df_sweep_stim[seq(1, nrow(df_sweep_stim), by = 10), ]

p_sweep_stim <- ggplot(df_sweep_stim_plot, aes(x = time, y = voltage)) +
  geom_line(colour = "#AD6030", linewidth = 0.3, alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, 14, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "Time (s)", y = "Voltage (V)", tag = "(A)") +
  stim_theme

df_freq <- tibble(
  time      = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
  inst_freq = c(1, 9500, 1, 9500, 1, 9500, 1, 9500, 1, 9500, 1)
)

p_sweep_freq <- ggplot(df_freq, aes(x = time, y = inst_freq)) +
  geom_line(colour = "#AD6030", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2000), expand = c(0.02, 0.02)) +
  labs(x = "Time (s)", y = "Frequency (Hz)", tag = "(B)") +
  stim_theme

df_sweep_resp      <- tibble(
  time     = stinger_data[, 1],
  velocity = stinger_data[, 3] * laser_sens
)

df_sweep_resp_plot <- df_sweep_resp

p_sweep_resp <- ggplot(df_sweep_resp_plot, aes(x = time, y = velocity)) +
  geom_line(colour = "#AD6030", linewidth = 0.15) +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "Time (s)", y = "Velocity (mm/s)", tag = "(C)") +
  stim_theme

fig2 <- arrangeGrob(
  p_sweep_stim,
  p_sweep_freq,
  p_sweep_resp,
  ncol = 1
)

agg_png(file.path(save_dir, "sweep_stimulus_and_response_signal.png"),
        width = 14, height = 14, units = "in", res = 200)
grid.draw(fig2)
dev.off()
cat("Saved sweep_stimulus_and_response_signal.png\n")