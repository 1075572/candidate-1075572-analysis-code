# reads three sweep equivalent peak csv files and plots mean FFT amplitude ± SEM
# vs distance for each peak as stacked panels, saved as a single PNG
# excludes stinger, gel 2 and 400um 

library(tidyverse)
library(ragg)
library(gridExtra)

csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values sweeps csv"
save_dir   <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# sweeps exclude 400um
dist_order <- c("0um","50um","100um","200um","300um","600um","800um","1000um")

peak_info <- list(
  list(file = "sweep_peak1_values.csv", tag = "(A)"),
  list(file = "sweep_peak2_values.csv", tag = "(B)"),
  list(file = "sweep_peak3_values.csv", tag = "(C)")
)

# colourblind-safe Okabe-Ito colours — warm family for sweeps, no repeats from click peaks
peak_colours <- c("#296339", "#FF007F", "#5D3999")

# shared theme — no title
stim_theme <- theme_bw(base_size = 11) +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_blank(),
    axis.title       = element_text(size = 13, family = "Helvetica"),
    axis.text        = element_text(size = 11, family = "Helvetica"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8),
    plot.tag         = element_text(family = "Helvetica", face = "bold", size = 20),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# build one plot per peak
plots <- lapply(seq_along(peak_info), function(i) {
  pk  <- peak_info[[i]]
  col <- peak_colours[i]
  
  df <- read.csv(file.path(csv_folder, pk$file), header = TRUE)
  colnames(df)[1:11] <- c("distance","Gel1","Gel2","Gel3","Gel4","Gel5",
                          "mean","median","SD","SEM","range")
  df <- df[df$distance != "stinger", ]
  df <- df[df$distance != "400um",   ]
  df$distance <- factor(df$distance, levels = dist_order)
  
  ggplot(df, aes(x = distance, y = mean, group = 1)) +
    geom_line(colour = col, linewidth = 0.8) +
    geom_point(colour = col, size = 3) +
    geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM),
                  width = 0.2, colour = col) +
    labs(x   = "Distance",
         y   = "Amplitude (mm/s)",
         tag = pk$tag) +
    stim_theme
})

# stack all 3 vertically and save
fig <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], ncol = 1)

agg_png(file.path(save_dir, "sweep_equivalent_peaks.png"),
        width = 8, height = 14, units = "in", res = 200)
grid.draw(fig)
dev.off()
cat("Saved sweep_equivalent_peaks.png\n")