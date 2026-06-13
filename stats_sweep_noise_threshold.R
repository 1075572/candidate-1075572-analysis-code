# SNR noise screening for sweep recordings
# 5 gels x 9 distances = 45 files minus stinger files and no result 400um Gel 3
# noise = row 1 col 5 x laser sensitivity. SNR = signal / noise.
# Shapiro-Wilk is run on the raw SNR first. Was normal (p >= 0.05), the 2SD
# threshold was applied directly on the raw scale.
# Gel 2 is excluded from statistical analysis and figures but retained in the
# exported CSV with NA for status.

library(tidyverse)
library(ragg)
library(gridExtra)
library(grid)

# Paths

sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"
results_path <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary/sweep_snr_results.csv"
save_dir     <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

laser_sens <- 5  # mm/s per volt

# File loading

files <- list.files(sweep_folder, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("stinger", ignore.case = TRUE, files)]

n_files    <- length(files)
n_expected <- 45

cat("File count check\n")
cat("Files found (after stinger exclusion):", n_files, "\n")
cat("Expected:", n_expected, "\n")

if (n_files == n_expected) {
  cat("File count looks correct.\n\n")
} else {
  cat("WARNING: expected", n_expected, "files but found", n_files,
      "- check for duplicates or missing files in the folder.\n")
  cat("Listing all files being read:\n")
  cat(paste(" ", basename(files)), sep = "\n")
  cat("\n")
}

# Read each file and parse gel/distance from filename
# Format: SWEEP_<distance>um_Gel<N>_results

records <- list()
for (f in files) {
  fname    <- tools::file_path_sans_ext(basename(f))
  raw      <- read.csv(f, header = FALSE)
  td_peaks <- raw[raw[, 1] != 0, 1]
  td_amp   <- mean(abs(td_peaks), na.rm = TRUE) * laser_sens
  noise    <- abs(raw[1, 5]) * laser_sens
  snr      <- if (noise > 0) td_amp / noise else NA
  gel      <- as.integer(str_match(fname, regex("gel(\\d+)", ignore_case = TRUE))[, 2])
  distance <- as.integer(str_match(fname, "(\\d+)um")[, 2])
  label    <- paste0("Gel ", gel, " ", distance, "um")
  records[[length(records) + 1]] <- data.frame(
    label = label, gel = gel, distance = distance,
    td_amp = td_amp, noise = noise, snr = snr
  )
}

df <- bind_rows(records) %>%
  arrange(gel, distance)

# Split into analysis set (Gel 2 excluded) and full set (all gels for CSV)

is_gel2     <- df$gel == 2
df_analysis <- df[!is_gel2, ]

cat("Gel 2 exclusion\n")
cat("Total files loaded:", nrow(df), "\n")
cat("Gel 2 files excluded from analysis:", sum(is_gel2), "(",
    paste(df$label[is_gel2], collapse = ", "), ")\n")
cat("Files entering statistical analysis:", nrow(df_analysis), "\n\n")

# Normality testing (analysis set only)

shapiro_raw <- shapiro.test(df_analysis$snr)
cat("Shapiro-Wilk: raw SNR (Gel 2 excluded)\n")
cat("W =", round(shapiro_raw$statistic, 4), " p =", round(shapiro_raw$p.value, 4), "\n")

if (shapiro_raw$p.value >= 0.05) {
  cat("Raw SNR is approximately normal - no transformation needed, 2SD applied on raw scale.\n\n")
  use_log <- FALSE
} else {
  cat("Raw SNR not normal - applying log transformation.\n\n")
  df_analysis$snr_log <- log(df_analysis$snr + 0.001)
  shapiro_log <- shapiro.test(df_analysis$snr_log)
  cat("Shapiro-Wilk: log SNR (Gel 2 excluded)\n")
  cat("W =", round(shapiro_log$statistic, 4), " p =", round(shapiro_log$p.value, 4), "\n")
  if (shapiro_log$p.value >= 0.05) {
    cat("Log-transformed values are approximately normal - 2SD threshold is justified.\n\n")
  } else {
    cat("Still not normal after log transform - treat the 2SD threshold with caution.\n\n")
  }
  use_log <- TRUE
}

# Threshold (derived from and applied to analysis set only)

if (!use_log) {
  raw_mean  <- mean(df_analysis$snr, na.rm = TRUE)
  raw_sd    <- sd(df_analysis$snr,   na.rm = TRUE)
  threshold <- raw_mean - 2 * raw_sd
} else {
  log_mean      <- mean(df_analysis$snr_log, na.rm = TRUE)
  log_sd        <- sd(df_analysis$snr_log,   na.rm = TRUE)
  log_threshold <- log_mean - 2 * log_sd
  threshold     <- exp(log_threshold)
}

df_analysis$status <- ifelse(df_analysis$snr < threshold, "noisy", "clean")

cat("Threshold (Gel 2 excluded)\n")
cat("Threshold (original SNR units):", round(threshold, 3), "\n")
cat("Clean:", sum(df_analysis$status == "clean"),
    " | Noisy:", sum(df_analysis$status == "noisy"), "\n\n")

if (any(df_analysis$status == "noisy")) {
  cat("Noisy files:\n")
  cat(paste(" ", df_analysis$label[df_analysis$status == "noisy"]), sep = "\n")
  cat("\n")
}

# Shared theme

snr_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 12, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8),
    plot.tag         = element_text(family = "Helvetica", face = "bold", size = 20)
  )

# Histogram (analysis set only — Gel 2 excluded)

p_hist <- ggplot(df_analysis, aes(x = snr, fill = status)) +
  geom_histogram(bins = 20, boundary = floor(min(df_analysis$snr, na.rm = TRUE)),
                 colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = threshold,
             colour = "#C41E3A", linewidth = 0.8, linetype = "dashed") +
  annotate("text",
           x = threshold + 0.05, y = Inf,
           label = paste0("threshold = ", round(threshold, 2)),
           hjust = 0, vjust = 1.5, size = 3.5,
           colour = "#C41E3A", family = "Helvetica") +
  scale_fill_manual(
    values = c(clean = "#AD6030", noisy = "#4B4B4B"),
    labels = c(clean = "Clean", noisy = "Noisy")
  ) +
  labs(x = "SNR (raw)", y = "Number of files", fill = NULL) +
  snr_theme +
  theme(legend.position = "none")

agg_png(file.path(save_dir, "sweep_snr_distribution.png"),
        width = 6, height = 5, units = "in", res = 200)
print(p_hist)
dev.off()
cat("Saved sweep_snr_distribution.png\n")

# Export results
# All files included. Gel 2 rows have NA for status as they were not analysed.

df$status  <- NA_character_
df$status[!is_gel2] <- df_analysis$status

if (use_log) {
  df$snr_log <- NA_real_
  df$snr_log[!is_gel2] <- df_analysis$snr_log
  df_out <- df[, c("label", "gel", "distance", "td_amp", "noise", "snr", "snr_log", "status")]
} else {
  df_out <- df[, c("label", "gel", "distance", "td_amp", "noise", "snr", "status")]
}

write.csv(df_out, results_path, row.names = FALSE)
cat("Saved sweep_snr_results.csv\n")

cat("\nSummary\n")
cat("Total files in CSV:               ", nrow(df), "(expected", n_expected, ")\n")
cat("Gel 2 files (status = NA in CSV): ", sum(is_gel2), "\n")
cat("Files analysed:                   ", nrow(df_analysis), "\n")
cat("  Clean:                          ", sum(df_analysis$status == "clean"), "\n")
cat("  Noisy:                          ", sum(df_analysis$status == "noisy"), "\n")