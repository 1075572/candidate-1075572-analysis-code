# This code creates a figure representing how noise was quantified during MATLAB analysis
## Creates a two-panel figure showing the pre-stimulus noise estimation window
# (rows 20-200, first ~2 ms) for click (top) and sweep (bottom) recordings at
# 1000 um, Gel 1. Each panel shows the raw signal, the half-max threshold zone
# (grey band), dotted threshold lines, red dashed average noise level lines,
# a green triangle for the single highest peak, and bright yellow triangles
# for all other peaks above the half-max threshold.

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(R.matlab)
library(ragg)
library(cowplot)   # <-- added for panel labelling

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

laser_sens <- 5.0  # mm/s per volt

noise_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8)
  )

find_noise_peaks <- function(amp_window) {
  abs_window <- abs(amp_window)
  max_val    <- max(abs_window)
  thresh     <- max_val / 2
  noise_av   <- mean(abs_window[abs_window > thresh])
  max_idx    <- which.max(abs_window)
  list(noise_av = noise_av, thresh = thresh, max_idx = max_idx)
}

# CLICK 1000 um 

click_mat  <- readMat("~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Clicks (Gel 1-5)/dataout_26_04_08_01_31_49_CLICK_1000um_Gel1.mat")
click_data <- click_mat$dataout[3:nrow(click_mat$dataout), ]
click_data <- click_data[!is.nan(click_data[, 1]), ]

df_click <- tibble(
  time     = click_data[, 1],
  velocity = click_data[, 2] * laser_sens
)

noise_window_c <- df_click$velocity[20:200]
result_c       <- find_noise_peaks(noise_window_c)
noise_av_c     <- result_c$noise_av
thresh_c       <- result_c$thresh

df_click_zoom <- tibble(
  time_ms  = df_click$time[20:200] * 1000,
  velocity = df_click$velocity[20:200]
)

abs_c          <- abs(df_click_zoom$velocity)
max_idx_c      <- which.max(abs_c)
above_thresh_c <- which(abs_c > thresh_c)
other_idx_c    <- above_thresh_c[above_thresh_c != max_idx_c]

df_click_max    <- df_click_zoom[max_idx_c, ]
df_click_others <- df_click_zoom[other_idx_c, ]

p_click_zoom <- ggplot(df_click_zoom, aes(x = time_ms, y = velocity)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -thresh_c, ymax = thresh_c,
           fill = "grey90", alpha = 0.8) +
  geom_hline(yintercept =  thresh_c, linetype = "dotted",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -thresh_c, linetype = "dotted",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept =  noise_av_c, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  geom_hline(yintercept = -noise_av_c, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  geom_line(colour = "#2E6BAD", linewidth = 0.5) +
  geom_point(data = df_click_others, aes(x = time_ms, y = velocity),
             colour = "#FFE800", size = 2.5, shape = 17) +
  geom_point(data = df_click_max, aes(x = time_ms, y = velocity),
             colour = "#00C000", size = 3.5, shape = 17) +
  annotate("text", x = df_click_max$time_ms + 0.01, y = df_click_max$velocity + 0.02,
           label = "Maximum point", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "#00C000") +
  annotate("text", x = 0.77, y = noise_av_c + 0.012,
           label = "Average noise level", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "red") +
  annotate("text", x = 0.77, y = thresh_c + 0.011,
           label = "Half-max threshold", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "grey40") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  labs(x = "Time (ms)", y = "Velocity (mm/s)") +
  noise_theme

# SWEEP 1000 um 

sweep_mat  <- readMat("~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Sine sweeps (Gel 1-5)/dataout_26_04_08_12_38_47_SWEEP_1000um_Gel1.mat")
sweep_data <- sweep_mat$dataoutSweep[3:nrow(sweep_mat$dataoutSweep), ]
sweep_data <- sweep_data[!is.nan(sweep_data[, 1]), ]

df_sweep <- tibble(
  time     = sweep_data[, 1],
  velocity = sweep_data[, 2] * laser_sens
)

noise_window_s <- df_sweep$velocity[20:200]
result_s       <- find_noise_peaks(noise_window_s)
noise_av_s     <- result_s$noise_av
thresh_s       <- result_s$thresh

df_sweep_zoom <- tibble(
  time_ms  = df_sweep$time[20:200] * 1000,
  velocity = df_sweep$velocity[20:200]
)

abs_s          <- abs(df_sweep_zoom$velocity)
max_idx_s      <- which.max(abs_s)
above_thresh_s <- which(abs_s > thresh_s)
other_idx_s    <- above_thresh_s[above_thresh_s != max_idx_s]

df_sweep_max    <- df_sweep_zoom[max_idx_s, ]
df_sweep_others <- df_sweep_zoom[other_idx_s, ]

p_sweep_zoom <- ggplot(df_sweep_zoom, aes(x = time_ms, y = velocity)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -thresh_s, ymax = thresh_s,
           fill = "grey90", alpha = 0.8) +
  geom_hline(yintercept =  thresh_s, linetype = "dotted",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -thresh_s, linetype = "dotted",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept =  noise_av_s, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  geom_hline(yintercept = -noise_av_s, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  geom_line(colour = "#AD6030", linewidth = 0.5) +
  geom_point(data = df_sweep_others, aes(x = time_ms, y = velocity),
             colour = "#FFE800", size = 2.5, shape = 17) +
  geom_point(data = df_sweep_max, aes(x = time_ms, y = velocity),
             colour = "#00C000", size = 3.5, shape = 17) +
  annotate("text", x = df_sweep_max$time_ms + 0.01, y = df_sweep_max$velocity - 0.02,
           label = "Maximum point", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "#00C000") +
  annotate("text", x = 0.77, y = noise_av_s + 0.011,
           label = "Average noise level", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "red") +
  annotate("text", x = 0.77, y = thresh_s + 0.009,
           label = "Half-max threshold", hjust = 0,
           size = 3.2, family = "Helvetica", colour = "grey40") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  labs(x = "Time (ms)", y = "Velocity (mm/s)") +
  noise_theme

# Stack vertically with (A) and (B) panel labels outside the plot area


fig_noise <- plot_grid(
  p_click_zoom,
  p_sweep_zoom,
  ncol        = 1,
  labels      = c("(A)", "(B)"),
  label_size  = 19,
  label_fontfamily = "Helvetica",
  label_fontface   = "bold",   # change to "bold" if preferred
  hjust       = -0.5,           # negative pulls label left of the panel edge
  vjust       = 1.2             # >1 pushes label above the panel edge
)

agg_png(file.path(save_dir, "noise_quantification_click_sweep_1000um.png"),
        width = 10, height = 10, units = "in", res = 200)
print(fig_noise)
dev.off()
cat("Saved noise_quantification_click_sweep_1000um.png\n")