# Skillings-Mack Test: Click Equivalent Peaks
# What this code does:
# Reads 3 CSV files, one per equivalent peak
# For each peak, builds a matrix: rows = gels, columns = distances
# Absent peaks are already blank in the CSV and read as NA automatically
# Runs Skillings-Mack test: does distance significantly affect amplitude?
# This test is used instead of Friedman because some peaks are absent at certain distances
# those absent values are treated as NA rather than zero,
# since absence means the peak was undetectable, not that amplitude was zero.
# Stinger is excluded because it is a direct contact measurement and not part of the distance attenuation series.
# Only 3 peaks are tested 
# EXCLUDES GEL 2 (TOO NOISY)

# If p < 0.05 the amplitude of that peak changes significantly with distance.

library(PMCMRplus)


# Paths to 3 click peak csv files 

csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values clicks csv"
peak_files <- list(
  list(file = "clicks_peak1_values.csv", label = "Peak 1 (~113 Hz)"),
  list(file = "clicks_peak2_values.csv", label = "Peak 2 (~175 Hz)"),
  list(file = "clicks_peak3_values.csv", label = "Peak 3 (~350 Hz)")
)


# distance order 

dist_order <- c("0um","50um","100um","200um","300um","400um","600um","800um","1000um")

# exclude stinger — different measurement condition from distance series
df <- df[df$distance != "stinger", ]

# run skillings-mack for each peak 

cat("SKILLINGS-MACK TEST: CLICK EQUIVALENT PEAKS\n")
cat("Absent peaks treated as NA\n")
cat("Rows = gels (blocks), Columns = distances (treatments)\n\n")

for (pk in peak_files) {

  cat("Testing:", pk$label, "\n\n")
  
  # read CSV — blank cells automatically become NA
  df <- read.csv(file.path(csv_folder, pk$file), header = TRUE)
  
  # exclude stinger — different measurement condition from distance series
  df <- df[df$distance != "stinger", ]
  
  # set distance order
  df$distance <- factor(df$distance, levels = dist_order)
  
  # extract just the gel columns
  gel_cols <- c("Gel1","Gel3","Gel4","Gel5")
  gel_data <- df[, gel_cols]
  
  # build matrix: rows = gels, columns = distances
  # skillingsMackTest requires rows=blocks(gels), cols=treatments(distances)
  mat <- t(as.matrix(gel_data))
  colnames(mat) <- as.character(df$distance)
  rownames(mat) <- gel_cols
  
  cat("Data matrix (rows=gels, cols=distances):\n")
  print(mat)
  cat("\n")
  
  # run Skillings-Mack test
  result <- tryCatch(
    skillingsMackTest(mat),
    error = function(e) {
      cat("Error running test:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(result)) {
    print(result)
    p <- result$p.value
    cat("\nInterpretation: ")
    if (!is.na(p) && p < 0.05) {
      cat("SIGNIFICANT effect of distance on", pk$label,
          "amplitude (p =", round(p, 4), ")\n\n")
    } else {
      cat("No significant effect of distance on", pk$label,
          "amplitude (p =", round(p, 4), ")\n\n")
    }
  }
}

cat("DONE")

# save results to csv
results_df <- data.frame(
  peak        = c("Peak 1 (~113 Hz)", "Peak 2 (~175 Hz)", "Peak 3 (~350 Hz)"),
  statistic   = c(10.314, 7.1468, 21.047),
  df          = c(8, 8, 8),
  p_value     = c(0.2436, 0.5209, 0.007023),
  significant = c("NO", "NO", "YES")
)
write.csv(results_df,
          "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/ep skillings mack test results/equ_peaks_click_skillings_mack_results.csv",
          row.names = FALSE)
cat("Results saved\n")