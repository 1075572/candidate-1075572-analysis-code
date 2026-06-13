
# Skillings-Mack Test: Sweep Equivalent Peaks
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

# If p < 0.05 the amplitude of that peak changes significantly with distance.

library(PMCMRplus)

# set the path to your sweep equivalent peak csv files
csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values sweeps csv"

# three peaks to test — csv files have been renamed to match
peak_files <- list(
  list(file = "sweep_peak1_values.csv", label = "Peak 1 (~2200 Hz)"),
  list(file = "sweep_peak2_values.csv", label = "Peak 2 (~4000 Hz)"),
  list(file = "sweep_peak3_values.csv", label = "Peak 3 (~8700 Hz)")
)

# distances to include — stinger excluded
dist_order <- c("0um","50um","100um","200um","300um","400um","600um","800um","1000um")

cat("Skillings-Mack Test: Sweep Equivalent Peaks\n")
cat("Absent peaks treated as NA\n")
cat("Rows = gels (blocks), Columns = distances (treatments)\n\n")

for (pk in peak_files) {
  
  cat("-----------------------------------------------------\n")
  cat("Testing:", pk$label, "\n\n")
  
  # read the csv, remove stinger and set distance order
  df <- read.csv(file.path(csv_folder, pk$file), header = TRUE)
  df <- df[df$distance != "stinger", ]
  df$distance <- factor(df$distance, levels = dist_order)
  df <- df[order(df$distance), ]
  
  # build the matrix — rows are gels (blocks), columns are distances (treatments)
  gel_cols <- c("Gel1","Gel2","Gel3","Gel4","Gel5")
  gel_data <- df[, gel_cols]
  mat <- t(as.matrix(gel_data))
  colnames(mat) <- as.character(df$distance)
  rownames(mat) <- gel_cols
  
  cat("Data matrix (rows=gels, cols=distances):\n")
  print(mat)
  cat("\n")
  
  # run the Skillings-Mack test
  result <- tryCatch(
    skillingsMackTest(mat),
    error = function(e) {
      cat("Error running test:", e$message, "\n")
      return(NULL)
    }
  )
  
  # print the result and a plain english interpretation
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

cat("Done\n\n")

# automatically extract results from the test and save to csv
results_list <- list()

for (i in seq_along(peak_files)) {
  pk <- peak_files[[i]]
  
  # read csv, remove stinger and build matrix same as above
  df <- read.csv(file.path(csv_folder, pk$file), header = TRUE)
  df <- df[df$distance != "stinger", ]
  df$distance <- factor(df$distance, levels = dist_order)
  df <- df[order(df$distance), ]
  
  gel_data <- df[, c("Gel1","Gel2","Gel3","Gel4","Gel5")]
  mat <- t(as.matrix(gel_data))
  colnames(mat) <- as.character(df$distance)
  rownames(mat) <- c("Gel1","Gel2","Gel3","Gel4","Gel5")
  
  result <- tryCatch(skillingsMackTest(mat), error = function(e) NULL)
  
  if (!is.null(result)) {
    results_list[[i]] <- data.frame(
      peak        = pk$label,
      statistic   = ifelse(is.na(result$statistic), NA, round(result$statistic, 4)),
      df          = result$parameter,
      p_value     = ifelse(is.na(result$p.value), NA, round(result$p.value, 6)),
      significant = ifelse(is.na(result$p.value), "UNTESTABLE",
                           ifelse(result$p.value < 0.05, "YES", "NO"))
    )
  }
}

results_df <- do.call(rbind, results_list)
write.csv(results_df,
          "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/ep skillings mack test results/equ_peaks_sweep_skillings_mack_results.csv",
          row.names = FALSE)
cat("Results saved\n")