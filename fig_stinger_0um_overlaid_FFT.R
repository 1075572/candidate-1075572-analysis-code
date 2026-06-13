# Delivery (stinger arm) vs gel response (0 um) FFT comparison
#
# ONE figure, two stacked panels:
#   (A) CLICK  - 4 gels across, x = 0-600 Hz, own y axis
#   (B) SWEEP  - 4 gels across, x = 0-9500 Hz, own y axis
# Gel titles kept with transparent strip backgrounds.
# Light display-binning makes the spectra read as clean curves;
# peak-finding for the table uses the RAW spectra (unbinned).
#
# Colours: CLICK blue family, SWEEP orange family.
#   darker = delivery (stinger), lighter = response (0 um)
#
# CLICK .mat: response in column 2 ; SWEEP .mat: response in column 3
# Base R only (no dplyr). Saves figure + table CSV to Figures folder.

library(ggplot2)
library(R.matlab)
library(gridExtra)
library(grid)

# Parameters
fs         <- 100000
laser_sens <- 5.0
gels       <- c(1, 3, 4, 5)
click_peaks <- c(113, 175, 350)
sweep_peaks <- c(2200, 4000, 8700)
click_xmax <- 600
sweep_xmax <- 9500
search_hw_click <- 40    # how far either side of a nominal peak to search (Hz)
search_hw_sweep <- 400

# how close delivery and response peaks need to be to count as a match
match_thresh_click <- 20
match_thresh_sweep <- 100

# bin width for display only - makes spectra look like smooth curves rather than noisy lines
click_bin <- 2
sweep_bin <- 15

# colours
click_delivery <- "#0A2E5C"; click_response <- "#4FC3F7"
sweep_delivery <- "#6E2C00"; sweep_response <- "#F5A623"

fig_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

base <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files"
click_dir <- file.path(base, "Clicks (Gel 1-5)")
sweep_dir <- file.path(base, "Sine sweeps (Gel 1-5)")

# Helpers
find_file <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) stop(paste("No file matching", pattern, "in", dir))
  hits[1]
}

load_velocity <- function(path, resp_col) {
  m <- readMat(path)
  # the data matrix may be stored under different variable names
  varname <- intersect(c("dataout", "dataoutSweep", "dataoutC"), names(m))[1]
  mat <- m[[varname]]
  raw <- mat[3:nrow(mat), ]
  ok <- !is.nan(raw[, 1]); raw <- raw[ok, ]
  list(time = raw[, 1], vel = raw[, resp_col] * laser_sens)
}

compute_fft <- function(time, vel, lo, hi) {
  idx <- which(time >= lo & time <= hi); v <- vel[idx]; v <- v - mean(v)
  n <- length(v); w <- 0.5 - 0.5 * cos(2*pi*(0:(n-1))/(n-1))
  f <- fft(v*w); mag <- Mod(f)[1:floor(n/2)]/n
  freq <- (0:(floor(n/2)-1))*fs/n
  data.frame(freq = freq, mag = mag)
}

nearest_peak <- function(spec, nominal, hw) {
  # find the frequency with the highest magnitude within hw Hz of the nominal peak
  win <- spec[spec$freq >= nominal-hw & spec$freq <= nominal+hw, ]
  if (nrow(win) == 0) return(NA_real_)
  win$freq[which.max(win$mag)]
}

# average magnitude within frequency bins - only used for display, not peak finding
bin_spectrum <- function(spec, binwidth, fmax) {
  spec <- spec[spec$freq <= fmax, ]
  brks <- seq(0, fmax + binwidth, by = binwidth)
  grp  <- cut(spec$freq, breaks = brks, labels = FALSE)
  fb   <- tapply(spec$freq, grp, mean)
  mb   <- tapply(spec$mag,  grp, mean)
  data.frame(freq = as.numeric(fb), mag = as.numeric(mb))
}

# Build
click_rows <- list(); sweep_rows <- list(); table_rows <- list()

for (g in gels) {
  # clicks use column 2 for the response signal
  cs <- load_velocity(find_file(click_dir, paste0("CLICK_stinger_Gel", g, "\\.mat$")), 2)
  c0 <- load_velocity(find_file(click_dir, paste0("CLICK_0um_Gel", g, "\\.mat$")), 2)
  cs_fft <- compute_fft(cs$time, cs$vel, 0.2, 0.6)
  c0_fft <- compute_fft(c0$time, c0$vel, 0.2, 0.6)
  
  for (pk in click_peaks) {
    dpk <- nearest_peak(cs_fft, pk, search_hw_click)
    rpk <- nearest_peak(c0_fft, pk, search_hw_click)
    diff_hz <- abs(dpk - rpk)
    matched <- ifelse(is.na(diff_hz), "NA",
                      ifelse(diff_hz <= match_thresh_click, "Match", "No match"))
    table_rows[[length(table_rows)+1]] <- data.frame(
      Stimulus="Click", Gel=g, Nominal_peak_Hz=pk,
      Delivery_peak_Hz=round(dpk,1), Response_peak_Hz=round(rpk,1),
      Delivery_offset_Hz=round(dpk-pk,1), Response_offset_Hz=round(rpk-pk,1),
      Delivery_minus_Response_Hz=round(dpk-rpk,1),
      Abs_difference_Hz=round(diff_hz,1),
      Match=matched, stringsAsFactors=FALSE)
  }
  
  csb <- bin_spectrum(cs_fft, click_bin, click_xmax); csb$trace <- "Delivery (stinger)"
  c0b <- bin_spectrum(c0_fft, click_bin, click_xmax); c0b$trace <- "Response (0 um)"
  cdf <- rbind(csb, c0b); cdf$gel <- paste("Gel", g)
  click_rows[[as.character(g)]] <- cdf
  
  # sweeps use column 3 for stinger and column 2 for 0 um response
  ss <- load_velocity(find_file(sweep_dir, paste0("SWEEP_stinger_Gel", g, "\\.mat$")), 3)
  s0 <- load_velocity(find_file(sweep_dir, paste0("SWEEP_0um_Gel", g, "\\.mat$")), 2)
  ss_fft <- compute_fft(ss$time, ss$vel, 2, 12)
  s0_fft <- compute_fft(s0$time, s0$vel, 2, 12)
  
  for (pk in sweep_peaks) {
    dpk <- nearest_peak(ss_fft, pk, search_hw_sweep)
    rpk <- nearest_peak(s0_fft, pk, search_hw_sweep)
    diff_hz <- abs(dpk - rpk)
    matched <- ifelse(is.na(diff_hz), "NA",
                      ifelse(diff_hz <= match_thresh_sweep, "Match", "No match"))
    table_rows[[length(table_rows)+1]] <- data.frame(
      Stimulus="Sweep", Gel=g, Nominal_peak_Hz=pk,
      Delivery_peak_Hz=round(dpk,1), Response_peak_Hz=round(rpk,1),
      Delivery_offset_Hz=round(dpk-pk,1), Response_offset_Hz=round(rpk-pk,1),
      Delivery_minus_Response_Hz=round(dpk-rpk,1),
      Abs_difference_Hz=round(diff_hz,1),
      Match=matched, stringsAsFactors=FALSE)
  }
  
  ssb <- bin_spectrum(ss_fft, sweep_bin, sweep_xmax); ssb$trace <- "Delivery (stinger)"
  s0b <- bin_spectrum(s0_fft, sweep_bin, sweep_xmax); s0b$trace <- "Response (0 um)"
  sdf <- rbind(ssb, s0b); sdf$gel <- paste("Gel", g)
  sweep_rows[[as.character(g)]] <- sdf
}

click_df <- do.call(rbind, click_rows)
sweep_df <- do.call(rbind, sweep_rows)

# scale each trace to its own maximum so all gels sit on the same 0-1 y axis
normalise <- function(d) {
  d$key <- paste(d$gel, d$trace); mx <- tapply(d$mag, d$key, max, na.rm=TRUE)
  d$mag_norm <- d$mag / mx[d$key]
  d$gel <- factor(d$gel, levels = paste("Gel", gels)); d
}
click_df <- normalise(click_df); sweep_df <- normalise(sweep_df)

# Table
match_table <- do.call(rbind, table_rows)
match_table <- match_table[order(match_table$Stimulus, match_table$Nominal_peak_Hz, match_table$Gel), ]
cat("\nPEAK FREQUENCY MATCHING TABLE\n")
print(match_table, row.names = FALSE)
write.csv(match_table, file.path(fig_dir, "peak_frequency_matching.csv"), row.names = FALSE)

# shared theme - transparent strip backgrounds so gel titles don't have a coloured box
strip_theme <- theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text   = element_text(size = 15, face = "bold"),
        axis.title   = element_text(size = 18),
        axis.text    = element_text(size = 12),
        legend.text  = element_text(size = 16),
        legend.key.width = unit(1.4, "cm"),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 20))

# Panel A: clicks
p_click <- ggplot(click_df, aes(freq, mag_norm, colour = trace)) +
  geom_vline(data = data.frame(peak = click_peaks), aes(xintercept = peak),
             linetype = "dashed", colour = "grey60", linewidth = 0.5, inherit.aes = FALSE) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  facet_wrap(~ gel, nrow = 1) +
  scale_colour_manual(values = c("Delivery (stinger)" = click_delivery,
                                 "Response (0 um)" = click_response)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 3))) +
  labs(x = "Frequency (Hz)", y = "Normalised FFT magnitude", colour = NULL, tag = "(A)") +
  strip_theme

# Panel B: sweeps
p_sweep <- ggplot(sweep_df, aes(freq, mag_norm, colour = trace)) +
  geom_vline(data = data.frame(peak = sweep_peaks), aes(xintercept = peak),
             linetype = "dashed", colour = "grey60", linewidth = 0.5, inherit.aes = FALSE) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  facet_wrap(~ gel, nrow = 1) +
  scale_colour_manual(values = c("Delivery (stinger)" = sweep_delivery,
                                 "Response (0 um)" = sweep_response)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 3))) +
  labs(x = "Frequency (Hz)", y = "Normalised FFT magnitude", colour = NULL, tag = "(B)") +
  strip_theme

# stack the two panels and save
png(file.path(fig_dir, "delivery_vs_response_comparison.png"),
    width = 18, height = 11, units = "in", res = 200)
grid.arrange(p_click, p_sweep, ncol = 1)
dev.off()
cat("\nSaved delivery_vs_response_comparison.png (A=click, B=sweep)\n")
cat("Saved peak_frequency_matching.csv\n")