# Builds the sweep delivered-velocity envelope figure for the calibration section.
# Shows delivered velocity across frequency against the 0.5 mm/s target.
# Time maps to frequency because the sweep is a linear chirp.

library(ggplot2)
library(dplyr)
library(grid)
library(R.matlab)
library(ragg)

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# settings
fs         <- 100000
fmin       <- 1
fmax       <- 9500
t_sweep    <- 2            # length of one up or down sweep
laser_sens <- 5.0          # mm/s per volt
target     <- 0.5          # the velocity we were aiming for
rms_win    <- 500          # window for the envelope (~5 ms)
nbins      <- 200          # bins for smoothing the envelope
up_starts  <- c(2, 6, 10)  # the three up-sweeps start at these times

# shared look (light grid)
stim_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8)
  )

# smooth amplitude envelope using a sliding RMS, then scale to a peak value
sliding_rms <- function(x, win) {
  half <- floor(win / 2)
  n    <- length(x)
  out  <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- sqrt(mean(x[lo:hi]^2))
  }
  out * sqrt(2)
}

# turns a sweep response into delivered velocity against frequency
envelope_vs_freq <- function(time, vel) {
  rows <- list()
  for (t0 in up_starts) {
    idx <- which(time >= t0 & time < (t0 + t_sweep))
    if (length(idx) == 0) next
    v   <- vel[idx]
    tl  <- time[idx] - t0
    env <- sliding_rms(v, rms_win)
    f   <- fmin + (fmax - fmin) * (tl / t_sweep)
    keep <- tl > 0.005 * t_sweep & tl < 0.995 * t_sweep   # trim only the exact turnaround samples
    rows[[length(rows) + 1]] <- tibble(freq = f[keep], env = env[keep])
  }
  all <- bind_rows(rows)
  brks <- seq(min(all$freq), max(all$freq), length.out = nbins + 1)
  grp  <- cut(all$freq, breaks = brks, labels = FALSE)
  tibble(freq = as.numeric(tapply(all$freq, grp, mean)),
         vel  = as.numeric(tapply(all$env,  grp, mean)))
}

# load the sweep recording (column 3 = laser response)
stinger      <- readMat("~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Sine sweeps (Gel 1-5)/dataout_26_04_08_11_50_42_SWEEP_stinger_Gel1.mat")
stinger_data <- stinger$dataout[3:nrow(stinger$dataout), ]
stinger_data <- stinger_data[!is.nan(stinger_data[, 1]), ]

df_resp <- tibble(time = stinger_data[, 1], velocity = stinger_data[, 3] * laser_sens)

# delivered velocity vs frequency
df_env <- envelope_vs_freq(df_resp$time, df_resp$velocity)

p_env <- ggplot(df_env, aes(x = freq, y = vel)) +
  geom_hline(yintercept = target, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 3900, y = target,
           label = "target 0.5 mm/s", vjust = -0.6, hjust = 0,
           size = 4, colour = "grey30") +
  geom_line(colour = "#AD6030", linewidth = 0.5) +
  scale_x_continuous(limits = c(1, fmax), breaks = seq(0, fmax, by = 1000),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "Frequency (Hz)", y = "Delivered velocity (mm/s)") +
  stim_theme

agg_png(file.path(save_dir, "sweep_delivered_velocity_envelope.png"),
        width = 10, height = 5, units = "in", res = 200)
grid.draw(p_env)
dev.off()
cat("Saved sweep_delivered_velocity_envelope.png\n")