# This code conducts the Friedman Test on the average Peak Amplitude of 8 sweeps in
# the time domain and produces Kendall's W to see the strentgh of the effect

# what does this code do?
# reads all sweep result csv files from a folder,
# excludes all 400um files because Gel 3 did not produce enough sweeps at this distance,
# this also excludes Gel 2 ( as it wasnoisy from 300um onwards),
# stigner also excluded
# pulls out the time domain peak amplitudes from each file (up to 8 peaks per file),
# averages those 8 peaks into one value per file,
# converts the values from volts to mm/s using the laser sensitivity,
# and then runs a Friedman test to ask whether distance has a significant effect on the amplitude.
# non-parametric equivalent of a repeated measures ANOVA.

library(tidyverse)
library(rcompanion)

# set folder paths
sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"
results_path <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/td friedman test results/td_friedman_average_sweeps_results.csv"

# the laser records in volts, multiply by laser sensitivity to get mm/s
laser_sens <- 5

# distances in the order we want them to appear (400um excluded)
dist_order <- c("0um","50um","100um","200um",
                "300um","600um","800um","1000um")

# find all csv files in the folder, then remove all 400um files
# 400um is excluded because Gel 3 did not produce enough sweeps at this distance
# Gel 2 is excluded because it has noisy data from 300um onwards
files <- list.files(sweep_folder, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("400um", files)]
files <- files[!grepl("Gel2", files)]

cat("Files loaded (400um and Gel 2 excluded):\n")
for (f in files) cat(" ", basename(f), "\n")
cat("\n")

# goes through each file and pulls out the time domain peak amplitudes
# sweep files have up to 8 peaks in column 1 — zeros are padding and are removed
# the 8 peaks are averaged into one value per file, then converted to mm/s
records <- list()
for (f in files) {
  fname    <- tools::file_path_sans_ext(basename(f))
  dist     <- str_extract(fname, "(?<=SWEEP_)(stinger|\\d+um)")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  raw      <- read.csv(f, header = FALSE)
  
  # remove zero rows (padding), take absolute value and average, then convert to mm/s
  td_peaks <- raw[raw[,1] != 0, 1]
  td_amp   <- mean(abs(td_peaks), na.rm = TRUE) * laser_sens
  
  records[[length(records) + 1]] <- data.frame(
    distance = dist,
    gel      = gel,
    td_amp   = td_amp
  )
}

# combine all records into one dataframe and set distance order
# combine all records, exclude stinger and set distance order
df <- bind_rows(records)
df <- df[df$distance != "stinger", ]
df$distance <- factor(df$distance, levels = dist_order)

cat("Parsed data:\n")
print(df)
cat("\n")

# build the matrix that the Friedman test needs
# rows are gels (the repeated blocks) and columns are distances (the treatments)
all_gels <- sort(unique(df$gel))
mat <- matrix(NA,
              nrow = length(all_gels),
              ncol = length(dist_order),
              dimnames = list(paste0("Gel", all_gels), dist_order))

for (i in 1:nrow(df)) {
  gel_idx  <- which(all_gels == df$gel[i])
  dist_idx <- which(dist_order == as.character(df$distance[i]))
  mat[gel_idx, dist_idx] <- df$td_amp[i]
}

cat("Data matrix (rows are gels, columns are distances):\n")
print(mat)
cat("\n")

# check every gel has a value at every distance before running the test
# Friedman will fail if any values are missing
incomplete <- which(apply(mat, 1, function(x) any(is.na(x))))
if (length(incomplete) > 0) {
  cat("Warning: these gels have missing distances and will cause the test to fail:\n")
  cat(rownames(mat)[incomplete], "\n\n")
} else {
  cat("All blocks are complete, Friedman test can proceed\n\n")
}

# run the Friedman test
cat("Friedman test: sweep time domain peak amplitude\n\n")

result <- tryCatch(
  friedman.test(mat),
  error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  }
)

# print the result, interpretation, Kendall's W and confidence intervals
if (!is.null(result)) {
  print(result)
  p <- result$p.value
  cat("\nInterpretation: ")
  if (!is.na(p) && p < 0.05) {
    cat("Distance has a significant effect on time domain peak amplitude (p =",
        round(p, 4), ")\n")
  } else {
    cat("No significant effect of distance on time domain peak amplitude (p =",
        round(p, 4), ")\n")
  }
  
  # Kendall's W measures the strength of the effect (effect size)
  # thresholds for k = 9 distances (from rcompanion handbook)
  n_gels     <- nrow(mat)
  k_dist     <- ncol(mat)
  kendalls_w <- result$statistic / (n_gels * (k_dist - 1))
  cat("\nEffect size (Kendall's W):", round(kendalls_w, 4))
  cat("\nInterpretation (k =", k_dist, "groups): ")
  if (kendalls_w < 0.1)  cat("negligible effect\n")
  if (kendalls_w >= 0.1 && kendalls_w < 0.2) cat("small effect\n")
  if (kendalls_w >= 0.2) cat("large effect\n")
  
  # confidence intervals for Kendall's W using bootstrap
  # with small sample sizes these will likely be wide — this is expected
  ci_result <- kendallW(mat, correct = TRUE, ci = TRUE)
  cat("Kendall's W 95% CI: [", round(ci_result$lower.ci, 4),
      ",", round(ci_result$upper.ci, 4), "]\n")
  
  # save the results to a csv file
  results_df <- data.frame(
    test                = "Friedman",
    stimulus            = "Sweep",
    variable            = "Time domain peak amplitude",
    gels_used           = paste(paste0("Gel", all_gels), collapse = ", "),
    statistic           = round(result$statistic, 4),
    df                  = result$parameter,
    p_value             = round(result$p.value, 6),
    significant         = ifelse(result$p.value < 0.05, "YES", "NO"),
    kendalls_w          = round(kendalls_w, 4),
    kendalls_w_lower_ci = round(ci_result$lower.ci, 4),
    kendalls_w_upper_ci = round(ci_result$upper.ci, 4)
  )
  write.csv(results_df, results_path, row.names = FALSE)
  cat("\nResults saved to:", results_path, "\n")
}