#  Rayleigh wave simulation in gelatin hydrogel
#  Figures:
#    Fig 1 - predicted velocity amplitude decay (central E = 400 Pa)
#    Fig 2 - sensitivity to Young's modulus estimate
#    Fig 3 - predicted vs measured velocity amplitude (data read from CSVs)
#    Fig 4 - near/far field approximation error
#    Fig 5 - master curve: u_z*sqrt(r) vs r/lambda_R

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggh4x)

#  Resolve namespace conflicts: bind these names to dplyr
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
rename <- dplyr::rename
arrange <- dplyr::arrange
summarise <- dplyr::summarise
bind_rows <- dplyr::bind_rows

save_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/rayleigh wave simulation figures"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

click_csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values clicks csv"
sweep_csv_folder <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/statistics/equivalent peaks analysis/EQP values sweeps csv"

#Setup

# Rayleigh wave speed ratio for incompressible material (nu = 0.5)
# exact solution to the Rayleigh secular equation
eta_R <- 0.9553^2   # c_R = 0.9553 * c_s, so eta_R = (c_R/c_s)^2

# Young's modulus values to test (Pa)
E_vals   <- c(200, 400, 600)
E_labels <- c("E = 200 Pa", "E = 400 Pa", "E = 600 Pa")

# Frequencies (Hz) - click peaks then sweep peaks
f_all <- c(113, 175, 350, 2200, 4000, 9000)

f_labels <- c(
  "Peak 1 ~113 Hz",
  "Peak 2 ~175 Hz",
  "Peak 3 ~350 Hz",
  "Peak 1 ~2200 Hz",
  "Peak 2 ~4000 Hz",
  "Peak 3 ~9000 Hz"
)

# factor level order controls plot panel order (clicks left, sweeps right)
f_factor_levels <- f_labels

# experimental distances (excluding stinger, excluding 400um for sweeps)
# used for vertical guide lines and near-field error calculation
dist_order <- c("0um","50um","100um","200um","300um","400um",
                "600um","800um","1000um")

r_exp_um <- c(0, 50, 100, 200, 300, 400, 600, 800, 1000)
r_exp    <- r_exp_um * 1e-6

# named numeric lookup: distance label -> micrometres
dist_to_num <- c(
  "0um"    =    0,
  "50um"   =   50,
  "100um"  =  100,
  "200um"  =  200,
  "300um"  =  300,
  "400um"  =  400,
  "600um"  =  600,
  "800um"  =  800,
  "1000um" = 1000
)

# Hankel function magnitude: |H0^(1)(x)| = sqrt(J0(x)^2 + Y0(x)^2)
# This is the exact Rayleigh wave amplitude decay function
hankel0_abs <- function(x) {
  sqrt(besselJ(x, 0)^2 + besselY(x, 0)^2)
}

# Far-field approximation: sqrt(2 / (pi * x))
# Valid when x = k_R * r >> 1 (i.e. many wavelengths from source)
farfield_abs <- function(x) {
  sqrt(2 / (pi * x))
}

# Material parameters - 20 mg/ml gelatin
rho     <- 1000
nu      <- 0.5

# Build model decay data
r_fine_um <- seq(5, 1000, by = 2)   # capped at 1000um
r_fine    <- r_fine_um * 1e-6

# measured mean velocity at 50um for each frequency (mm/s)
# used to scale delta_0 so model passes through measured value at 50um
# this anchors the model to the data at the closest measurement point
vel_at_50um <- c(
  0.000338,   # Click Peak 1 ~113 Hz
  0.000168,   # Click Peak 2 ~175 Hz
  0.000407,   # Click Peak 3 ~350 Hz
  0.000857,   # Sweep Peak 1 ~2200 Hz
  0.0000909,  # Sweep Peak 2 ~4000 Hz
  0.0000968   # Sweep Peak 3 ~9000 Hz
)

df_decay <- data.frame()
for (ei in seq_along(E_vals)) {          # loop through each stiffness value (Young's modulus)
  E_i   <- E_vals[ei]                    # get the current stiffness (E)
  G_i   <- E_i / (2 * (1 + nu))         # calculates shear modulus (resistance to shape change)
  c_s_i <- sqrt(G_i / rho)              # calculates the speed of shear waves in the gel
  c_R_i <- sqrt(eta_R) * c_s_i          # calculates the speed of the Rayleigh surface ripples
  
  for (fi in seq_along(f_all)) {        # loop through each vibration frequency (clicks and sweeps)
    omega_i  <- 2 * pi * f_all[fi]      # convert standard frequency (Hz) to angular frequency
    k_R_i    <- omega_i / c_R_i         # calculate how many wave cycles fit in a meter
    lam_R_i  <- 2 * pi / k_R_i          # calculates wavelength (distance between two ripple crests)
    
    r_ref  <- 50e-6   # set a reference distance of 50 micrometers
    # The exact theoretical model uses a Bessel function (Y0) to calculate the wave's height.
    # The mathematical quirk of this function is that at a distance of exactly 0 (right under the needle),
    # the calculated wave height blows up to infinity (a logarithmic singularity).
    # Because you cannot calculate or plot "infinity," a reference point slightly away from the source
    # must be chosen to start the calculation
    H_ref  <- hankel0_abs(k_R_i * r_ref)   # find the exact mathematical wave height at this 50um mark
    
    # scale delta_0 so model velocity at 50um matches measured value
    # vel at 50um = delta_0 * 2 * pi * f * 1000
    # so delta_0 = vel_at_50um / (2 * pi * f * 1000)
    # this scaling only applies to E = 400 Pa (central estimate)
    # for other E values the shape is preserved but anchored at the 400 Pa measured amplitude
    delta_0_scaled <- vel_at_50um[fi] / (2 * pi * f_all[fi] * 1000)
    
    A_norm   <- delta_0_scaled / H_ref              # create a scaling factor so our model matches the source height
    kRr      <- k_R_i * r_fine                      # multiply wavenumber by our distances (creates a dimensionless scale)
    
    amp_exact <- A_norm * hankel0_abs(kRr)          # calculate exact wave heights for all distances
    amp_ff    <- A_norm * H_ref * farfield_abs(kRr) / farfield_abs(k_R_i * r_ref)  # calculate shortcut far-field heights
    
    df_decay <- rbind(df_decay, data.frame(         # save all these calculated numbers into data table
      r_um         = r_fine_um,                     # distance in micrometers
      amp_exact    = amp_exact,                     # exact wave height (meters)
      amp_ff       = amp_ff,                        # shortcut wave height (meters)
      vel_exact    = amp_exact * 2 * pi * f_all[fi] * 1000,  # exact velocity amplitude (mm/s) - displacement * 2*pi*f
      vel_ff       = amp_ff    * 2 * pi * f_all[fi] * 1000,  # shortcut velocity amplitude (mm/s)
      E_label      = E_labels[ei],                  # name of the stiffness used
      E_val        = E_i,                           # number value of stiffness
      freq         = f_all[fi],                     # number value of frequency
      freq_label   = factor(f_labels[fi], levels = f_factor_levels),  # name of frequency (keeps plot order nice)
      lambda_R_um  = lam_R_i * 1e6,                 # wavelength in micrometers
      nearfield_um = 0.6 * lam_R_i * 1e6            # boundary where the far-field shortcut becomes accurate (<1% error)
    ))
  }
}

df_nf <- df_decay %>%                 # create a mini-table just for the near-field boundary lines
  dplyr::filter(E_val == 400) %>%          # only grab the boundaries for the middle stiffness (400 Pa)
  select(freq_label, nearfield_um) %>% # keep only the frequency names and boundary distances
  distinct()                           # remove duplicate rows

# Build near-field error data for Figure 4
# calculated fresh at exactly the experimental distances rather than filtering
# df_decay (which uses r_fine_um and can cause floating point mismatches)
df_nf_err <- data.frame()
for (fi in seq_along(f_all)) {
  E_i     <- 400                            # central estimate only
  G_i     <- E_i / (2 * (1 + nu))
  c_s_i   <- sqrt(G_i / rho)
  c_R_i   <- sqrt(eta_R) * c_s_i
  omega_i <- 2 * pi * f_all[fi]
  k_R_i   <- omega_i / c_R_i
  
  r_ref          <- 50e-6
  H_ref          <- hankel0_abs(k_R_i * r_ref)
  delta_0_scaled <- vel_at_50um[fi] / (2 * pi * f_all[fi] * 1000)
  A_norm         <- delta_0_scaled / H_ref
  
  r_use <- r_exp_um[r_exp_um > 0] * 1e-6   # experimental distances excluding source point
  kRr   <- k_R_i * r_use
  
  vel_exact_i <- A_norm * hankel0_abs(kRr) * 2 * pi * f_all[fi] * 1000
  vel_ff_i    <- A_norm * H_ref * farfield_abs(kRr) / farfield_abs(k_R_i * r_ref) * 2 * pi * f_all[fi] * 1000
  
  df_nf_err <- rbind(df_nf_err, data.frame(
    r_um       = r_exp_um[r_exp_um > 0],    # experimental distances in micrometers
    ff_error   = abs(vel_ff_i - vel_exact_i) / vel_exact_i * 100,  # percentage error
    freq_label = case_when(
      f_all[fi] == 113  ~ "Click: Peak 1 ~113 Hz",
      f_all[fi] == 175  ~ "Click: Peak 2 ~175 Hz",
      f_all[fi] == 350  ~ "Click: Peak 3 ~350 Hz",
      f_all[fi] == 2200 ~ "Sweep: Peak 1 ~2200 Hz",
      f_all[fi] == 4000 ~ "Sweep: Peak 2 ~4000 Hz",
      f_all[fi] == 9000 ~ "Sweep: Peak 3 ~9000 Hz"
    )
  ))
}

# Read equivalent peaks CSVs and compute mean and SEM velocity across replicates
# Each CSV is expected to have a 'distance' column (character, e.g. "50um") and
# one numeric column per replicate containing velocity amplitude values (mm/s).
# exclude_400um: set TRUE for sweeps where 400um was not a valid measurement point
read_eqp_csv <- function(path, label, exclude_400um = FALSE) {
  df_raw <- read.csv(path, stringsAsFactors = FALSE)
  
  # rename first column to 'distance' regardless of what it was called in the CSV
  colnames(df_raw)[1] <- "distance"
  
  # optionally drop the 400um row (sweeps only)
  if (exclude_400um) {
    df_raw <- df_raw %>% dplyr::filter(distance != "400um")
  }
  
  # convert distance label to numeric micrometres using the lookup table
  df_raw <- df_raw %>%
    dplyr::filter(distance %in% names(dist_to_num)) %>%
    mutate(r_um = as.numeric(dist_to_num[distance]))
  
  # pivot all replicate columns to long format, then summarise
  rep_cols <- setdiff(colnames(df_raw), c("distance", "r_um"))
  
  df_long <- df_raw %>%
    pivot_longer(cols = all_of(rep_cols),
                 names_to  = "replicate",
                 values_to = "vel") %>%
    dplyr::filter(!is.na(vel))
  
  df_out <- df_long %>%
    group_by(r_um) %>%
    summarise(
      mean_vel = mean(vel,  na.rm = TRUE),
      sem_vel  = sd(vel, na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    mutate(freq_label = factor(label, levels = f_factor_levels))
  
  return(df_out)
}
# in df_measured, change "Peak 3 ~8700 Hz" → "Peak 3 ~9000 Hz" to match f_labels
df_measured <- bind_rows(
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak1_values.csv"),
               "Peak 1 ~113 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak2_values.csv"),
               "Peak 2 ~175 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(click_csv_folder, "clicks_peak3_values.csv"),
               "Peak 3 ~350 Hz",  exclude_400um = FALSE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak1_values.csv"),
               "Peak 1 ~2200 Hz", exclude_400um = TRUE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak2_values.csv"),
               "Peak 2 ~4000 Hz", exclude_400um = TRUE),
  read_eqp_csv(file.path(sweep_csv_folder, "sweep_peak3_values.csv"),
               "Peak 3 ~9000 Hz", exclude_400um = TRUE)   # was "Peak 3 ~8700 Hz"
)

# Shared settings
x_lim <- c(0, 1000)

custom_strips <- strip_themed(
  background_x = elem_list_rect(fill = c(
    "#C2DEFF", "#C2DEFF", "#C2DEFF",
    "#FFD9BC", "#FFD9BC", "#FFD9BC"
  ))
)

sim_theme <- theme_bw(base_size = 14) +
  theme(
    text            = element_text(family = "Helvetica"),
    plot.title      = element_text(family = "Helvetica", face = "bold",
                                   hjust = 0.5, size = 18),
    plot.subtitle   = element_text(family = "Helvetica", hjust = 0.5, size = 14),
    axis.title.x    = element_text(family = "Helvetica", size = 15,
                                   margin = margin(t = 8)),
    axis.title.y    = element_text(family = "Helvetica", size = 15),
    axis.text       = element_text(family = "Helvetica", size = 13),
    legend.text     = element_text(family = "Helvetica", size = 13),
    legend.title    = element_text(family = "Helvetica", face = "bold", size = 15),
    strip.text      = element_text(family = "Helvetica", face = "bold", size = 13),
    legend.position = "bottom"
  )

model_colour_scale <- scale_colour_manual(
  name   = NULL,
  values = c(
    "Exact (Hankel)"                 = "black",
    "Far-field shortcut (1/sqrt(r))" = "#CC0000",
    "Near-field boundary"            = "#0000FF"
  )
)

model_linetype_scale <- scale_linetype_manual(
  name   = NULL,
  values = c(
    "Exact (Hankel)"                 = "solid",
    "Far-field shortcut (1/sqrt(r))" = "longdash",
    "Near-field boundary"            = "dashed"
  )
)

# Figure 1
fig1 <- ggplot(df_decay %>% filter(E_val == 400), aes(x = r_um)) +
  geom_vline(data = data.frame(r_exp_um_val = r_exp_um),
             aes(xintercept = r_exp_um_val),
             colour = "grey92", linetype = "solid", linewidth = 0.3) +
  geom_line(aes(y = vel_exact,
                colour   = "Exact (Hankel)",
                linetype = "Exact (Hankel)"),
            linewidth = 0.9) +
  geom_line(aes(y = vel_ff,
                colour   = "Far-field shortcut (1/sqrt(r))",
                linetype = "Far-field shortcut (1/sqrt(r))"),
            linewidth = 0.9) +
  geom_vline(data = df_nf,
             aes(xintercept = nearfield_um,
                 colour     = "Near-field boundary",
                 linetype   = "Near-field boundary"),
             linewidth = 0.7) +
  model_colour_scale +
  model_linetype_scale +
  scale_x_continuous(limits = x_lim, breaks = c(0, 200, 400, 600, 800, 1000)) +
  facet_wrap2(~ freq_label, scales = "free_y", ncol = 3, strip = custom_strips) +
  labs(
    x = "Distance from source (um)",
    y = "Velocity amplitude (mm/s)"
  ) +
  sim_theme +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

# Figure 2 
fig2 <- ggplot(df_decay,
               aes(x = r_um, y = vel_exact,
                   colour = E_label, linetype = E_label)) +
  geom_vline(data = data.frame(r_exp_um_val = r_exp_um),
             aes(xintercept = r_exp_um_val),
             colour = "grey92", linetype = "solid", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#0055AA", "black", "#CC0000")) +
  scale_x_continuous(limits = x_lim, breaks = c(0, 200, 400, 600, 800, 1000)) +
  facet_wrap2(~ freq_label, scales = "free_y", ncol = 3, strip = custom_strips) +
  labs(
    x        = "Distance from source (um)",
    y        = "Velocity amplitude (mm/s)",
    colour   = "Young's modulus",
    linetype = "Young's modulus"
  ) +
  sim_theme

# Figure 3 - predicted vs measured velocity amplitude
# freq_label stays as-is (matches model) for correct faceting.
# Each measured series gets a fixed colour string based on stimulus type,
# so the legend shows e.g. "Measured (Click)" and "Measured (Sweep)" once each,
# not six separate entries.
measured_levels <- c(
  "Measured: Peak 1 (Click)",
  "Measured: Peak 2 (Click)",
  "Measured: Peak 3 (Click)",
  "Measured: Peak 1 (Sweep)",
  "Measured: Peak 2 (Sweep)",
  "Measured: Peak 3 (Sweep)"
)

df_measured_fig3 <- df_measured %>%
  mutate(
    freq_label   = factor(as.character(freq_label), levels = f_factor_levels),
    series_label = factor(
      case_when(
        grepl("113",  as.character(freq_label)) ~ "Measured: Peak 1 (Click)",
        grepl("175",  as.character(freq_label)) ~ "Measured: Peak 2 (Click)",
        grepl("350",  as.character(freq_label)) ~ "Measured: Peak 3 (Click)",
        grepl("2200", as.character(freq_label)) ~ "Measured: Peak 1 (Sweep)",
        grepl("4000", as.character(freq_label)) ~ "Measured: Peak 2 (Sweep)",
        grepl("9000", as.character(freq_label)) ~ "Measured: Peak 3 (Sweep)"
      ),
      levels = measured_levels
    )
  )

fig3 <- ggplot(df_decay %>% filter(E_val == 400), aes(x = r_um)) +
  geom_vline(data = data.frame(r_exp_um_val = r_exp_um),
             aes(xintercept = r_exp_um_val),
             colour = "grey92", linetype = "solid", linewidth = 0.3) +
  geom_line(aes(y        = vel_exact,
                colour   = "Exact (Hankel)",
                linetype = "Exact (Hankel)"),
            linewidth = 0.9) +
  geom_line(aes(y        = vel_ff,
                colour   = "Far-field shortcut (1/sqrt(r))",
                linetype = "Far-field shortcut (1/sqrt(r))"),
            linewidth = 0.9) +
  geom_vline(data = df_nf,
             aes(xintercept = nearfield_um,
                 colour     = "Near-field boundary",
                 linetype   = "Near-field boundary"),
             linewidth = 0.7) +
  geom_line(data = df_measured_fig3,
            aes(x = r_um, y = mean_vel,
                colour   = series_label,
                linetype = series_label),
            linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(data = df_measured_fig3,
             aes(x = r_um, y = mean_vel,
                 colour = series_label),
             shape = 23, size = 3, fill = "white", inherit.aes = FALSE) +
  geom_errorbar(data = df_measured_fig3,
                aes(x      = r_um,
                    ymin   = mean_vel - sem_vel,
                    ymax   = mean_vel + sem_vel,
                    colour = series_label),
                width = 20, linewidth = 0.5, inherit.aes = FALSE) +
  scale_colour_manual(
    name   = NULL,
    values = c(
      "Exact (Hankel)"                 = "black",
      "Far-field shortcut (1/sqrt(r))" = "#CC0000",
      "Near-field boundary"            = "#0000FF",
      "Measured: Peak 1 (Click)"       = "#00E676",
      "Measured: Peak 2 (Click)"       = "#FF9ECD",
      "Measured: Peak 3 (Click)"       = "#AA00FF",
      "Measured: Peak 1 (Sweep)"       = "#296339",
      "Measured: Peak 2 (Sweep)"       = "#FF007F",
      "Measured: Peak 3 (Sweep)"       = "#5D3999"
    )
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(
      "Exact (Hankel)"                 = "solid",
      "Far-field shortcut (1/sqrt(r))" = "longdash",
      "Near-field boundary"            = "dashed",
      "Measured: Peak 1 (Click)"       = "solid",
      "Measured: Peak 2 (Click)"       = "solid",
      "Measured: Peak 3 (Click)"       = "solid",
      "Measured: Peak 1 (Sweep)"       = "solid",
      "Measured: Peak 2 (Sweep)"       = "solid",
      "Measured: Peak 3 (Sweep)"       = "solid"
    )
  ) +
  facet_wrap2(~ freq_label, scales = "free_y", ncol = 3, strip = custom_strips) +
  scale_x_continuous(limits = x_lim, breaks = c(0, 200, 400, 600, 800, 1000)) +
  labs(
    x = "Distance from source (um)",
    y = "Velocity amplitude (mm/s)"
  ) +
  sim_theme +
  guides(colour   = guide_legend(ncol = 3, byrow = FALSE),
         linetype = guide_legend(ncol = 3, byrow = FALSE))

# Figure 4
fig4 <- ggplot(df_nf_err,
               aes(x = r_um, y = ff_error, colour = freq_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "#E69F00") +
  annotate("text", x = 950, y = 6.5,
           label  = "5% threshold", size = 5,
           family = "Helvetica", colour = "grey40") +
  annotate("text", x = 950, y = 2.2,
           label  = "1% threshold", size = 5,
           family = "Helvetica", colour = "#E69F00") +
  scale_colour_manual(
    name = "Frequency",
    values = c(
      "Click: Peak 1 ~113 Hz"  = "#1B9E77",
      "Click: Peak 2 ~175 Hz"  = "#D95F02",
      "Click: Peak 3 ~350 Hz"  = "#7570B3",
      "Sweep: Peak 1 ~2200 Hz" = "#E7298A",
      "Sweep: Peak 2 ~4000 Hz" = "#66A61E",
      "Sweep: Peak 3 ~9000 Hz" = "#E6AB02"
    ),
    breaks = c(
      "Click: Peak 1 ~113 Hz",  "Click: Peak 2 ~175 Hz",  "Click: Peak 3 ~350 Hz",
      "Sweep: Peak 1 ~2200 Hz", "Sweep: Peak 2 ~4000 Hz", "Sweep: Peak 3 ~9000 Hz"
    )
  ) +
  scale_x_continuous(limits = x_lim, breaks = c(0, 200, 400, 600, 800, 1000)) +
  labs(
    x      = "Distance from source (um)",
    y      = "Far-field approximation error (%)",
    colour = "Frequency"
  ) +
  sim_theme +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

fig4

# save block
main_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Main text"
supp_dir <- "~/Library/CloudStorage/OneDrive-Nexus365/Year 4/Data analysis/Figures/Supplementary"

if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

ggsave(file.path(main_dir, "rayleigh_sim_amplitude_decay.png"),
       fig1, width = 14, height = 8, dpi = 200)
cat("Saved Fig 1 → Main text\n")

ggsave(file.path(main_dir, "rayleigh_sim_measured_vs_predicted.png"),
       fig3, width = 14, height = 8, dpi = 200)
cat("Saved Fig 3 → Main text\n")

ggsave(file.path(supp_dir, "rayleigh_sim_E_sensitivity.png"),
       fig2, width = 14, height = 8, dpi = 200)
cat("Saved Fig 2 → Supplementary\n")

ggsave(file.path(supp_dir, "rayleigh_sim_nearfield_error.png"),
       fig4, width = 12, height = 6, dpi = 200)
cat("Saved Fig 4 → Supplementary\n")