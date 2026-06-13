# Reads raw click .mat files and saves one figure per gel (1, 3, 4, 5):
# a 3x3 grid of FFT spectra, one panel per distance (0-1000um),
# line colour graduating light to dark blue with increasing distance.
# Distance labels are black with a coloured line beside them matching the line.
# Peak legend kept below the grid (enlarged). Panel B (amplitude vs distance)
# has been REMOVED. Stinger and Gel 2 excluded throughout.
#
# Gel 1 saves to Main text folder; Gels 3, 4, 5 save to Supplementary folder.

library(R.matlab)
library(tidyverse)
library(gridExtra)
library(pracma)
library(cowplot)
library(grid)
library(ragg)

click_folder   <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/cleaned data files/Clicks (Gel 1-5)"
csv_folder     <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/EQP freq values csv/EQP freq values clicks csv"

plot_folder_main <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
plot_folder_supp <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary"
dir.create(plot_folder_main, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_folder_supp, recursive = TRUE, showWarnings = FALSE)

laser_sens <- 5

peak_info <- list(
  list(file = "clicks_peak1_freq_values.csv", label = "Peak 1 (~113 Hz)", colour = "#00E676"),
  list(file = "clicks_peak2_freq_values.csv", label = "Peak 2 (~175 Hz)", colour = "#FF9ECD"),
  list(file = "clicks_peak3_freq_values.csv", label = "Peak 3 (~350 Hz)", colour = "#AA00FF")
)

dist_order   <- c("0um","50um","100um","200um","300um","400um","600um","800um","1000um")
dist_colours <- colorRampPalette(c("#9DC4F0", "#1E4F85"))(length(dist_order))
names(dist_colours) <- dist_order

read_click_fft <- function(filepath) {
  mat        <- readMat(filepath)
  dataout    <- mat$dataout
  last_row   <- max(which(!is.na(dataout[, 1])))
  data_clean <- dataout[3:last_row, , drop = FALSE]
  time       <- data_clean[, 1]
  amp_tm     <- data_clean[, 2] * laser_sens
  m          <- length(amp_tm)
  NFFT       <- 2^ceiling(log2(m))
  delta_t    <- time[3] - time[2]
  sf         <- 1 / delta_t
  Y          <- fft(c(amp_tm, rep(0, NFFT - m)))
  f          <- sf * (0:(NFFT/2)) / NFFT
  P          <- Mod(Y / NFFT)
  Mag        <- 2 * P[1:(NFFT/2 + 1)]
  Mag_smooth <- movavg(Mag, 5, type = "s")
  n_bins     <- min(2500, length(f))
  list(f = f[1:n_bins], mag = Mag_smooth[1:n_bins])
}

files <- list.files(click_folder, pattern = "\\.mat$", full.names = TRUE)
file_info <- tibble(
  filepath = files,
  fname    = tools::file_path_sans_ext(basename(files)),
  distance = str_extract(fname, "(?<=CLICK_)(stinger|\\d+um)"),
  gel      = as.integer(str_extract(fname, "(?<=Gel)\\d+"))
)
file_info$distance[is.na(file_info$distance)] <- "800um"

get_peak_freqs <- function(csv_file, gel_col) {
  path <- file.path(csv_folder, csv_file)
  if (!file.exists(path)) stop("Cannot find freq CSV: ", path)
  df <- read.csv(path, header = TRUE)
  df %>%
    dplyr::select(distance, dplyr::all_of(gel_col)) %>%
    dplyr::filter(distance != "stinger") %>%
    dplyr::rename(freq_val = dplyr::all_of(gel_col)) %>%
    dplyr::filter(!is.na(freq_val)) %>%
    dplyr::mutate(freq_val = as.numeric(freq_val))
}

get_local_peak <- function(fft_data, target_freq, window_hz = 20) {
  freq_step  <- fft_data$f[2] - fft_data$f[1]
  idx_centre <- which.min(abs(fft_data$f - target_freq))
  half_win   <- max(1, round(window_hz / freq_step))
  idx_lo     <- max(1, idx_centre - half_win)
  idx_hi     <- min(length(fft_data$f), idx_centre + half_win)
  idx_peak   <- idx_lo + which.max(fft_data$mag[idx_lo:idx_hi]) - 1
  list(freq = fft_data$f[idx_peak], mag = fft_data$mag[idx_peak])
}

make_peak_legend_grob <- function(pk) {
  df <- tibble(x = 1, y = 1, label = pk$label, colour = pk$colour)
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(colour = label, fill = label), shape = 23, size = 9) +
    scale_colour_manual(values = setNames(pk$colour, pk$label), name = NULL) +
    scale_fill_manual(values   = setNames(pk$colour, pk$label), name = NULL) +
    guides(colour = guide_legend(title = NULL),
           fill   = guide_legend(title = NULL)) +
    theme_void() +
    theme(
      legend.position   = "bottom",
      legend.text       = element_text(size = 30, family = "Helvetica"),
      legend.key.size   = unit(2.2, "cm"),
      legend.margin     = margin(t = 1, b = 2, l = 0, r = 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  cowplot::get_legend(p)
}

make_title_grob <- function(dist_label, line_colour) {
  txt_gp   <- gpar(fontsize = 24, fontface = "bold",
                   fontfamily = "Helvetica", col = "black")
  gap      <- unit(0.025, "npc")
  line_len <- unit(0.04,  "npc")
  lwd_val  <- 10
  
  shift <- (line_len + gap) * 0.5 + unit(0.06, "npc")
  txt_centred <- textGrob(dist_label,
                          x  = unit(0.5, "npc") + shift,
                          y  = unit(0.5, "npc"),
                          gp = txt_gp)
  half_txt2 <- unit(0.5, "grobwidth", txt_centred)
  dash_x1 <- unit(0.5, "npc") + shift - half_txt2 - gap
  dash_x0 <- dash_x1 - line_len
  left_line <- segmentsGrob(
    x0 = dash_x0, x1 = dash_x1,
    y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
    gp = gpar(col = line_colour, lwd = lwd_val)
  )
  grobTree(left_line, txt_centred)
}

x_title_centre <- 0.525

fft_y_label <- textGrob("Magnitude (mm/s)",
                        rot = 90, hjust = 0.5, vjust = 0.275,
                        gp = gpar(fontsize = 26, fontfamily = "Helvetica"))
fft_x_label <- textGrob("Frequency (Hz)",
                        x = unit(x_title_centre, "npc"), hjust = 0.5, vjust = 1.4,
                        gp = gpar(fontsize = 26, fontfamily = "Helvetica"))

fft_bottom_row <- c(7, 8, 9)


# Main loop

for (gel_num in c(1, 3, 4, 5)) {
  
  gel_col <- paste0("Gel", gel_num)
  cat("\nProcessing Gel", gel_num, "...\n")
  
  gel_files <- file_info %>%
    dplyr::filter(gel == gel_num, distance != "stinger") %>%
    dplyr::arrange(factor(distance, levels = dist_order))
  
  peak_freqs <- lapply(peak_info, function(pk) {
    freqs <- get_peak_freqs(pk$file, gel_col)
    list(freqs = freqs, label = pk$label, colour = pk$colour)
  })
  
  fft_list <- list()
  for (i in 1:nrow(gel_files)) {
    fft_list[[i]] <- tryCatch(
      read_click_fft(gel_files$filepath[i]),
      error = function(e) { cat("Could not read:", gel_files$fname[i], "\n"); NULL }
    )
  }
  
  global_y_max <- max(sapply(fft_list, function(d) {
    if (is.null(d) || !is.list(d) || is.null(d$mag)) return(0)
    max(d$mag, na.rm = TRUE)
  }))
  global_x_max <- max(sapply(fft_list, function(d) {
    if (is.null(d) || !is.list(d) || is.null(d$f)) return(0)
    max(d$f, na.rm = TRUE)
  }))
  diamond_offset <- global_y_max * 0.03
  
  blank_panel <- ggplot() + theme_void()
  panels <- setNames(vector("list", length(dist_order)), dist_order)
  
  for (i in 1:nrow(gel_files)) {
    dist   <- gel_files$distance[i]
    result <- fft_list[[i]]
    if (is.null(result)) { panels[[dist]] <- blank_panel; next }
    
    slot_idx  <- which(dist_order == dist)
    is_bottom <- slot_idx %in% fft_bottom_row
    line_col  <- dist_colours[dist]
    df_plot   <- tibble(freq = result$f, mag = result$mag)
    
    lp_list <- lapply(peak_freqs, function(pk) {
      freq_row <- pk$freqs %>% dplyr::filter(distance == dist)
      if (nrow(freq_row) > 0 && !is.na(freq_row$freq_val[1])) {
        get_local_peak(result, freq_row$freq_val[1], window_hz = 20)
      } else NULL
    })
    
    title_grob <- make_title_grob(dist, line_col)
    
    p <- ggplot(df_plot, aes(x = freq, y = mag)) +
      geom_line(colour = line_col, linewidth = 0.7)
    
    for (j in seq_along(peak_freqs)) {
      lp <- lp_list[[j]]
      if (!is.null(lp)) {
        p <- p +
          annotate("segment", x = 0, xend = global_x_max,
                   y = lp$mag, yend = lp$mag,
                   colour = peak_freqs[[j]]$colour,
                   linewidth = 0.4, linetype = "dashed") +
          geom_point(data = tibble(x = lp$freq, y = lp$mag + diamond_offset),
                     aes(x = x, y = y), shape = 23, size = 5,
                     colour = peak_freqs[[j]]$colour,
                     fill = peak_freqs[[j]]$colour, inherit.aes = FALSE)
      }
    }
    
    p <- p +
      labs(x = NULL, y = NULL) +
      scale_x_continuous(limits = c(0, global_x_max),
                         expand = expansion(mult = c(0, 0.03))) +
      coord_cartesian(ylim = c(0, global_y_max * 1.08)) +
      theme_bw() +
      theme(
        text             = element_text(family = "Helvetica"),
        plot.title       = element_blank(),
        axis.text.x      = if (is_bottom) element_text(size = 18, family = "Helvetica") else element_blank(),
        axis.ticks.x     = if (is_bottom) element_line() else element_blank(),
        axis.text.y      = element_text(size = 18, family = "Helvetica"),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
      )
    
    panels[[dist]] <- arrangeGrob(title_grob, p, nrow = 2, heights = c(0.12, 1))
  }
  
  panels <- lapply(panels, function(p) if (is.null(p)) blank_panel else p)
  
  # enlarged peak legend (per-gel labels)
  gel_peak_info <- lapply(peak_info, function(pk) {
    lbl <- sub("^(Peak \\d+)", paste0("\\1 across Gel ", gel_num), pk$label)
    list(label = lbl, colour = pk$colour)
  })
  peak_legend_grobs <- lapply(gel_peak_info, make_peak_legend_grob)
  
  legend_row <- arrangeGrob(
    peak_legend_grobs[[1]], peak_legend_grobs[[2]], peak_legend_grobs[[3]],
    ncol = 3
  )
  
  out_folder   <- if (gel_num == 1) plot_folder_main else plot_folder_supp
  out_filename <- paste0("click_fft_Gel", gel_num, "_all_distances_PANELA.png")
  out_path     <- file.path(out_folder, out_filename)
  
  # taller/bigger FFT panels: bigger canvas, FFT takes most of the height
  ragg::agg_png(out_path, width = 3400, height = 3600, res = 150)
  
  gridExtra::grid.arrange(
    arrangeGrob(
      arrangeGrob(grobs = panels, ncol = 3),
      left   = fft_y_label,
      bottom = fft_x_label
    ),
    legend_row,
    nrow    = 2,
    heights = c(0.9, 0.1)
  )
  
  dev.off()
  cat("Saved:", out_path, "\n")
}

cat("\nAll gels done.\n")