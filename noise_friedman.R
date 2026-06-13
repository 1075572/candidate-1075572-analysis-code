# Friedman test on background noise across distance for click and sweep recordings
# Tests whether background noise floor varies significantly with distance from source
# Follows same structure as click/sweep time domain amplitude Friedman tests
# Gel 2 excluded throughout; stinger excluded throughout
# Click: 4 gels x 9 distances (complete block)
# Sweep: 4 gels x 8 distances (400um excluded as Gel 3 missing)

library(tidyverse)
library(rcompanion)

# Paths
click_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all click results"
sweep_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/all sweep results"
results_path <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary/noise_friedman_results.csv"

laser_sens <- 5

dist_order_click <- c("0um","50um","100um","200um","300um","400um","600um","800um","1000um")
dist_order_sweep <- c("0um","50um","100um","200um","300um","600um","800um","1000um")

# run Friedman + Kendall's W on a noise matrix 
run_friedman_noise <- function(mat, stimulus_label) {
  cat("\n--- Friedman test:", stimulus_label, "background noise ---\n")
  cat("Data matrix (rows = gels, columns = distances):\n")
  print(round(mat, 5))
  cat("\n")
  
  result <- tryCatch(
    friedman.test(mat),
    error = function(e) { cat("Error:", e$message, "\n"); return(NULL) }
  )
  
  if (is.null(result)) return(NULL)
  
  print(result)
  p          <- result$p.value
  n_gels     <- nrow(mat)
  k_dist     <- ncol(mat)
  kendalls_w <- result$statistic / (n_gels * (k_dist - 1))
  
  cat("\nInterpretation: ")
  if (!is.na(p) && p < 0.05) {
    cat("Distance has a significant effect on background noise (p =", round(p, 4), ")\n")
  } else {
    cat("No significant effect of distance on background noise (p =", round(p, 4), ")\n")
  }
  
  cat("Kendall's W:", round(kendalls_w, 4))
  cat("  (k =", k_dist, "groups)\n")
  
  ci_result <- kendallW(mat, correct = TRUE, ci = TRUE)
  cat("95% CI: [", round(ci_result$lower.ci, 4), ",", round(ci_result$upper.ci, 4), "]\n")
  
  data.frame(
    test                = "Friedman",
    stimulus            = stimulus_label,
    variable            = "Background noise",
    gels_used           = paste(rownames(mat), collapse = ", "),
    statistic           = round(result$statistic, 4),
    df                  = result$parameter,
    p_value             = round(result$p.value, 6),
    significant         = ifelse(p < 0.05, "YES", "NO"),
    kendalls_w          = round(kendalls_w, 4),
    kendalls_w_lower_ci = round(ci_result$lower.ci, 4),
    kendalls_w_upper_ci = round(ci_result$upper.ci, 4)
  )
}

# ── Click noise ───────────────────────────────────────────────────────────────
# noise = column 4, row 1 of each click results CSV (same as click SNR script)

files_click <- list.files(click_folder, pattern = "\\.csv$", full.names = TRUE)
files_click <- files_click[!grepl("Gel2|stinger", files_click, ignore.case = TRUE)]

records_click <- list()
for (f in files_click) {
  fname    <- tools::file_path_sans_ext(basename(f))
  dist     <- str_extract(fname, "(?<=CLICK_)(stinger|\\d+um)")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  raw      <- read.csv(f, header = FALSE)
  noise    <- abs(raw[1, 4]) * laser_sens
  records_click[[length(records_click) + 1]] <- data.frame(
    distance = dist, gel = gel, noise = noise
  )
}

df_click <- bind_rows(records_click) %>%
  filter(distance != "stinger") %>%
  mutate(distance = factor(distance, levels = dist_order_click))

all_gels_click <- sort(unique(df_click$gel))
mat_click <- matrix(NA,
                    nrow = length(all_gels_click),
                    ncol = length(dist_order_click),
                    dimnames = list(paste0("Gel", all_gels_click), dist_order_click)
)
for (i in 1:nrow(df_click)) {
  gi <- which(all_gels_click == df_click$gel[i])
  di <- which(dist_order_click == as.character(df_click$distance[i]))
  mat_click[gi, di] <- df_click$noise[i]
}

res_click <- run_friedman_noise(mat_click, "Click")

# ── Sweep noise ───────────────────────────────────────────────────────────────
# noise = column 5, row 1 of each sweep results CSV (same as sweep SNR script)
# 400um excluded (Gel 3 missing → incomplete block)

files_sweep <- list.files(sweep_folder, pattern = "\\.csv$", full.names = TRUE)
files_sweep <- files_sweep[!grepl("Gel2|stinger|400um", files_sweep, ignore.case = TRUE)]

records_sweep <- list()
for (f in files_sweep) {
  fname    <- tools::file_path_sans_ext(basename(f))
  dist     <- str_extract(fname, "(?<=SWEEP_)(stinger|\\d+um)")
  gel      <- as.integer(str_extract(fname, "(?<=Gel)\\d+"))
  raw      <- read.csv(f, header = FALSE)
  noise    <- abs(raw[1, 5]) * laser_sens
  records_sweep[[length(records_sweep) + 1]] <- data.frame(
    distance = dist, gel = gel, noise = noise
  )
}

df_sweep <- bind_rows(records_sweep) %>%
  filter(distance != "stinger") %>%
  mutate(distance = factor(distance, levels = dist_order_sweep))

all_gels_sweep <- sort(unique(df_sweep$gel))
mat_sweep <- matrix(NA,
                    nrow = length(all_gels_sweep),
                    ncol = length(dist_order_sweep),
                    dimnames = list(paste0("Gel", all_gels_sweep), dist_order_sweep)
)
for (i in 1:nrow(df_sweep)) {
  gi <- which(all_gels_sweep == df_sweep$gel[i])
  di <- which(dist_order_sweep == as.character(df_sweep$distance[i]))
  mat_sweep[gi, di] <- df_sweep$noise[i]
}

res_sweep <- run_friedman_noise(mat_sweep, "Sweep")

# Save combined results 
results_combined <- bind_rows(res_click, res_sweep)
write.csv(results_combined, results_path, row.names = FALSE)
cat("\nResults saved to:", results_path, "\n")
print(results_combined)