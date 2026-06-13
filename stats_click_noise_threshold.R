# SNR noise screening for click recordings
#
# 5 gels x 9 distances = 45 files minus stinger files and 3 no results in gel 2
# noise = row 1 col 5 x laser sensitivity. SNR = signal / noise.
#
# For each file: signal = row 1 col 1 × laser sensitivity, noise = row 1 col 4 × same.
# SNR = signal / noise. Raw SNR turned out right-skewed (Shapiro-Wilk p < 0.05),
# so log-transform before setting the threshold. Log SNR is approximately normal
# (W ≈ 0.95, p ≈ 0.077), which means a mean − 2SD cutoff is sensible. The threshold
# was back-transformed via exp() so its the same units as the raw data.
# Gel 2 is excluded from statistical analysis and figures but retained in the
# exported CSV with NA for status.

library(tidyverse)
library(ragg)
library(gridExtra)
library(grid)
library(cowplot)

# Paths

click_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all click results"
results_path <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary/click_snr_results.csv"
save_dir     <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

laser_sens <- 5  # mm/s per volt

# File loading

files <- list.files(click_folder, pattern = "\\.csv$", full.names = TRUE)
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
# Format: CLICK_<distance>um_Gel<N>_results

records <- list()
for (f in files) {
  fname    <- tools::file_path_sans_ext(basename(f))
  raw      <- read.csv(f, header = FALSE)
  td_amp   <- abs(raw[1, 1]) * laser_sens
  noise    <- abs(raw[1, 4]) * laser_sens
  snr      <- if (noise > 0) td_amp / noise else NA
  gel      <- as.integer(str_match(fname, regex("gel(\\d+)", ignore_case = TRUE))[, 2])
  distance <- as.integer(str_match(fname, "(\\d+)um")[, 2])
  label    <- paste0("Gel ", gel, " ", distance, "um")
  records[[length(records) + 1]] <- data.frame(
    label = label, td_amp = td_amp, noise = noise, snr = snr
  )
}

df <- bind_rows(records) %>%
  arrange(gel = as.integer(str_match(label, "Gel (\\d+)")[, 2]),
          distance = as.integer(str_match(label, "(\\d+)um")[, 2]))

# Split into analysis set (Gel 2 excluded) and full set (all gels for CSV)

is_gel2     <- grepl("gel\\s*2", df$label, ignore.case = TRUE)
df_analysis <- df[!is_gel2, ]

cat("Gel 2 exclusion\n")
cat("Total files loaded:", nrow(df), "\n")
cat("Gel 2 files excluded from analysis:", sum(is_gel2), "(",
    paste(df$label[is_gel2], collapse = ", "), ")\n")
cat("Files entering statistical analysis:", nrow(df_analysis), "\n\n")

# Normality testing (analysis set only)

shapiro_raw <- shapiro.test(df_analysis$snr)
cat("Shapiro-Wilk: raw SNR (Gel 2 excluded)\n")
cat("W =", round(shapiro_raw$statistic, 4), " p =", round(shapiro_raw$p.value, 4), "\n\n")

df_analysis$snr_log <- log(df_analysis$snr + 0.001)
shapiro_log <- shapiro.test(df_analysis$snr_log)
cat("Shapiro-Wilk: log SNR (Gel 2 excluded)\n")
cat("W =", round(shapiro_log$statistic, 4), " p =", round(shapiro_log$p.value, 4), "\n")

if (shapiro_log$p.value >= 0.05) {
  cat("Log-transformed values are approximately normal - 2SD threshold is justified.\n\n")
} else {
  cat("Still not normal after log transform - treat the 2SD threshold with caution.\n\n")
}

# Threshold (derived from and applied to analysis set only)

log_mean      <- mean(df_analysis$snr_log, na.rm = TRUE)
log_sd        <- sd(df_analysis$snr_log,   na.rm = TRUE)
log_threshold <- log_mean - 2 * log_sd
threshold     <- exp(log_threshold)

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

# Panel A: raw SNR distribution (analysis set only)

p_raw <- ggplot(df_analysis, aes(x = snr, fill = status)) +
  geom_histogram(bins = 20, boundary = 2,
                 colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = threshold,
             colour = "#C41E3A", linewidth = 0.8, linetype = "dashed") +
  annotate("text",
           x = threshold + 0.05, y = Inf,
           label = paste0("threshold = ", round(threshold, 2)),
           hjust = 0, vjust = 1.5, size = 3.5,
           colour = "#C41E3A", family = "Helvetica") +
  scale_fill_manual(
    values = c(clean = "#2E6BAD", noisy = "#4B4B4B")
  ) +
  labs(x = "SNR (raw)", y = "Number of files", tag = "(A)") +
  snr_theme +
  theme(legend.position = "none")

# Panel B: log SNR distribution (analysis set only)

p_log <- ggplot(df_analysis, aes(x = snr_log, fill = status)) +
  geom_histogram(bins = 20, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = log_threshold,
             colour = "#C41E3A", linewidth = 0.8, linetype = "dashed") +
  scale_fill_manual(
    values = c(clean = "#2E6BAD", noisy = "#4B4B4B")
  ) +
  labs(x = "Log(SNR)", y = "Number of files", tag = "(B)") +
  snr_theme +
  theme(legend.position = "none")

# Panel C: Q-Q plot (analysis set only)

p_qq <- ggplot(df_analysis, aes(sample = snr_log)) +
  stat_qq(shape = 21, size = 1.5,
          fill = "#2E6BAD", colour = "black", stroke = 0.4) +
  stat_qq_line(colour = "#C41E3A", linewidth = 0.8) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", tag = "(C)") +
  snr_theme

# Layout

top_row    <- arrangeGrob(p_raw, p_log, ncol = 2)
bottom_row <- arrangeGrob(
  nullGrob(), p_qq, nullGrob(),
  ncol   = 3,
  widths = c(0.25, 0.5, 0.25)
)

fig <- arrangeGrob(top_row, bottom_row, nrow = 2)

# Save

agg_png(file.path(save_dir, "click_snr_distribution_transformation.png"),
        width = 10, height = 10, units = "in", res = 200)
grid.draw(fig)
dev.off()
cat("Saved click_snr_distribution_transformation.png\n")

# Export results
# All files included. Gel 2 rows have NA for status and snr_log as they were
# not part of the analysis.

df$status  <- NA_character_
df$status[!is_gel2] <- df_analysis$status

df$snr_log <- NA_real_
df$snr_log[!is_gel2] <- df_analysis$snr_log

df_out <- df[, c("label", "td_amp", "noise", "snr", "snr_log", "status")]
write.csv(df_out, results_path, row.names = FALSE)
cat("Saved click_snr_results.csv\n")

cat("\nSummary\n")
cat("Total files in CSV:               ", nrow(df), "(expected", n_expected, ")\n")
cat("Gel 2 files (status = NA in CSV): ", sum(is_gel2), "\n")
cat("Files analysed:                   ", nrow(df_analysis), "\n")
cat("  Clean:                          ", sum(df_analysis$status == "clean"), "\n")
cat("  Noisy:                          ", sum(df_analysis$status == "noisy"), "\n")