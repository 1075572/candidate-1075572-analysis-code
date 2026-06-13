# This code conducts the Friedman Test on each of the 8 individual sweep peaks
# separately in the time domain
#
# what does this code do?
# reads all sweep result csv files from a folder,
# excludes 400um (Gel 3 incomplete) and Gel 2 (noisy from 300um onwards),
# extracts each of the 8 individual peaks separately per file,
# converts values from volts to mm/s,
# runs a separate Friedman test for each of the 8 peak positions,
# asking whether distance has a significant effect on that peak's amplitude.
# applies Bonferroni correction (p < 0.00625) because 8 tests are run on the same data.
# non-parametric equivalent of a repeated measures ANOVA.

library(tidyverse)
library(rcompanion)

# set folder paths
sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"
results_path <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/td friedman test results/td_friedman_sweeps_individual_peaks_results.csv"

# the laser records in volts, multiply by laser sensitivity to get mm/s
laser_sens <- 5

# Bonferroni corrected significance threshold for 8 tests
bonferroni_threshold <- 0.05 / 8

# distances in the order we want them to appear (400um excluded)
dist_order <- c("0um","50um","100um","200um",
                "300um","600um","800um","1000um")

# find all csv files, exclude 400um and Gel 2
files <- list.files(sweep_folder, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("400um", files)]
files <- files[!grepl("Gel2",  files)]

cat("Files loaded (400um and Gel 2 excluded):\n")
for (f in files) cat(" ", basename(f), "\n")
cat("\n")

# goes through each file and extracts all 8 individual peaks separately
# zeros are padding and are removed, peaks are capped at 8
records <- list()
for (f in files) {
  fname <- tools::file_path_sans_ext(basename(f))
  dist  <- str_extract(fname, "(?<=SWEEP_)(stinger|\\d+um)")
  gel   <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  raw   <- read.csv(f, header = FALSE)
  
  td_peaks <- raw[raw[,1] != 0, 1]
  td_peaks <- abs(td_peaks[1:min(8, length(td_peaks))]) * laser_sens
  
  for (pk in seq_along(td_peaks)) {
    records[[length(records) + 1]] <- data.frame(
      distance = dist,
      gel      = gel,
      peak_num = pk,
      td_amp   = td_peaks[pk]
    )
  }
}

# combine into one dataframe, exclude stinger and set distance order
df <- bind_rows(records)
df <- df[df$distance != "stinger", ]
df$distance <- factor(df$distance, levels = dist_order)

cat("Total records:", nrow(df), "\n\n")

# store results for all 8 peaks
all_results <- list()

# run a separate Friedman test for each of the 8 peak positions
for (pk in 1:8) {
  
  cat("=====================================================\n")
  cat("Peak", pk, "\n\n")
  
  pk_df <- df %>% filter(peak_num == pk)
  all_gels <- sort(unique(pk_df$gel))
  
  # build the matrix for this peak
  mat <- matrix(NA,
                nrow = length(all_gels),
                ncol = length(dist_order),
                dimnames = list(paste0("Gel", all_gels), dist_order))
  
  for (i in 1:nrow(pk_df)) {
    gel_idx  <- which(all_gels == pk_df$gel[i])
    dist_idx <- which(dist_order == as.character(pk_df$distance[i]))
    mat[gel_idx, dist_idx] <- pk_df$td_amp[i]
  }
  
  cat("Data matrix:\n")
  print(mat)
  cat("\n")
  
  # check for complete blocks
  incomplete <- which(apply(mat, 1, function(x) any(is.na(x))))
  if (length(incomplete) > 0) {
    cat("Warning: incomplete blocks for gels:", rownames(mat)[incomplete], "\n\n")
    next
  }
  
  # run Friedman test
  result <- tryCatch(
    friedman.test(mat),
    error = function(e) { cat("Error:", e$message, "\n"); return(NULL) }
  )
  
  if (!is.null(result)) {
    print(result)
    p <- result$p.value
    
    cat("\nInterpretation (Bonferroni threshold p <", bonferroni_threshold, "):\n")
    if (!is.na(p) && p < bonferroni_threshold) {
      cat("SIGNIFICANT effect of distance on peak", pk, "amplitude (p =", round(p, 6), ")\n")
    } else {
      cat("No significant effect of distance on peak", pk, "amplitude (p =", round(p, 6), ")\n")
    }
    
    # Kendall's W effect size — thresholds for k = 9 distances
    n_gels     <- nrow(mat)
    k_dist     <- ncol(mat)
    kendalls_w <- result$statistic / (n_gels * (k_dist - 1))
    cat("\nEffect size (Kendall's W):", round(kendalls_w, 4))
    cat("\nInterpretation (k =", k_dist, "groups): ")
    if (kendalls_w < 0.1)  cat("negligible effect\n")
    if (kendalls_w >= 0.1 && kendalls_w < 0.2) cat("small effect\n")
    if (kendalls_w >= 0.2) cat("large effect\n")
    
    # confidence intervals for Kendall's W
    ci_result <- kendallW(mat, correct = TRUE, ci = TRUE)
    cat("Kendall's W 95% CI: [", round(ci_result$lower.ci, 4),
        ",", round(ci_result$upper.ci, 4), "]\n\n")
    
    # store results
    all_results[[pk]] <- data.frame(
      test                = "Friedman",
      stimulus            = "Sweep",
      variable            = paste("Individual peak", pk, "time domain amplitude"),
      gels_used           = paste(paste0("Gel", all_gels), collapse = ", "),
      peak_num            = pk,
      statistic           = round(result$statistic, 4),
      df                  = result$parameter,
      p_value             = round(p, 6),
      bonferroni_threshold = bonferroni_threshold,
      significant         = ifelse(p < bonferroni_threshold, "YES", "NO"),
      kendalls_w          = round(kendalls_w, 4),
      kendalls_w_lower_ci = round(ci_result$lower.ci, 4),
      kendalls_w_upper_ci = round(ci_result$upper.ci, 4)
    )
  }
}

# combine and save all results
results_df <- bind_rows(all_results)
write.csv(results_df, results_path, row.names = FALSE)
cat("All results saved to:", results_path, "\n")