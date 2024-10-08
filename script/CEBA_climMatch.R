# R script that uses different methods to do climate matching for C. basicorne
# Author: Brittany Barker, Oregon IPM Center, Oregon State University
# Associated publication: "Climate matching models for the yellow starthistle 
# rosette weevil (Coleoptera: Apionidae) for the western U.S."
# Method 1: CLIMEX approach (regional)
# Method 2: Climatch algorithm using the "climatchR" package

# Load packages
library(here)
library(tidyverse)
library(ggspatial) # Spatial options for ggplot2
library(rnaturalearth) # Polygons for countries/states
library(sf) # Working with spatial features
library(terra) # Working with rasters
library(tidyterra) # Convert rasters to data frames
library(Euclimatch) # Climatch
library(GA) # Color palettes
library(cowplot) # Combining plots

#### Plot theme ----

mytheme <- theme(legend.text = element_text(size = 14), 
                 legend.justification = "center",
                 legend.position = "bottom",
                 plot.title = element_text(size = 20, face = "bold", vjust = 0.9),
                 panel.border = element_rect(colour = "black", linewidth = 0.75, fill = NA),
                 plot.background = element_rect(colour = "white", fill = "white"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(fill = "white", colour = "white"), 
                 axis.title.x = element_blank(), 
                 axis.title.y = element_blank(),
                 axis.text = element_text(size = 13, face = "bold", colour = "black"),
                 plot.margin = margin(t = 0.5, r = 1, b = 0, l = 1, unit = "cm"),
                 legend.key.width  =  unit(1.75, "cm"),
                 legend.title = element_text(size = 18, hjust = 0, vjust = 0.95, face = "bold"))


#### Functions ----

# Calculate Match Indices for temperature and soil moisture.
# Equations for indices are provided on pp. 99-100 of CLIMEX manual (vs. 2, 2015).
# Itmin = Tmin index (weight = W1) = exp(-kT * Tdmin)
# Itmax = Tmax index (weight = W2) = exp(-kT * Tdmax)
# Ismst = Soil moisture index (weight = W7) = exp((-kSM * Sd))
Index_calc <- function(v) {
  
  # Select columns
  if (v == "ppt") {
    home <- select(home_clim, contains(v), "ppt_ann")
    away <- select(away_clim, contains(v), "ppt_ann")  
  } else {
    home <- select(home_clim, contains(v))
    away <- select(away_clim, contains(v))  
  }
  
  # Calculate absolute differences between home and away locations for each week
  diff_weeks <- map(1:length(weeks), function(i) {
    
    # Only 1 data point for ppt (annual for year)
    # weekly (12 data points) for other variables
    # It's pointless to run ppt 12 times, but results in less code
    if (v == "ppt_ann") {
      home_m <- as.numeric(home)
      away_m <- away
    } else {
      home_m <- as.numeric(select(home, paste0(v, weeks[i])))
      away_m <- select(away, paste0(v, weeks[i]))
    }
    
    if (v == "ppt_ann") {
      # Annual total precipitation
      # Difference in annual precipitation, adjusted so that differences in precipitation 
      # are more significant for locations with lower precipitation
      away_m$diff <- abs(home_m - away_m[,1]) /
        (1 + (0.001*(home_m + away_m[,1])))
      # Temperature
    } else {
      away_m$diff <- abs(home_m - away_m[,1])
    }
    
    away_m <- select(away_m, diff)
    names(away_m) <- paste0(v, weeks[i])
    
    return(away_m)
    
  })
  
  # Combine results and calculate average of differences for each week
  all_diffs <- Reduce(cbind, diff_weeks)
  all_diffs$avg <- rowMeans(all_diffs)
  
  # Calculate the index
  # Constants are the same as CLIMEX except for soil moisture
  # (constant in CLIMEX assumes a relative scale - i.e. 0 to 1)
  # Precipitation index
  if (v == "ppt_ann") {
    all_diffs$ind <- exp(-0.004*all_diffs$avg)
    # Soil moisture index
  } else if (v == "sm") {
    all_diffs$ind <- exp(-0.0008*all_diffs$avg)
    # Temperature index
  } else {
    all_diffs$ind <- exp(-0.1*all_diffs$avg)
  }
  
  return(all_diffs)
}

# Plot function to create maps for CLIMEX-based match indices
# "ind" is the index for which map is being created
CMI_base_map <- function(ind) {
  ggplot() +
    geom_sf(data = na_states, fill = "gray90") +
    geom_raster(data = filter(all_ind_l, name == ind), aes(x = x, y = y, fill = value)) +
    geom_sf(data = na_states_l, color="black") +
    annotation_spatial_hline(intercept = kilkis[2], crs = 4326, color = "white", size = 1.25, alpha=0.85) +
    annotation_spatial_hline(intercept = kilkis[2], crs = 4326, color = "black", size = 0.5, linetype="dashed") +
    mytheme 
}

# Plot function to create maps for Climatch scores
# Input is variable(s) to be mapped
Climatch_base_map <- function(v) {
  ggplot() +
    geom_sf(data = na_states, fill = "gray90", linewidth = 0.5, color = "gray20") +
    geom_raster(data = filter(match_outs_df, variable == v), aes(x = x, y = y, fill = value)) +
    geom_sf(data = na_states_l, color = "gray15") +
    annotation_spatial_hline(intercept = kilkis[2], crs = 4326, color = "white", size = 1.25, alpha=0.85) +
    annotation_spatial_hline(intercept = kilkis[2], crs = 4326, color = "black", size = 0.5, linetype="dashed") +
    mytheme
}

#### Home location and spatial features ----

# Home climate (Kilkis)
kilkis <- c(22.844, 40.994)

# North America features
sf_use_s2(FALSE) # Avoids issues with spherical geometries

# North America features
na_states <- ne_states(returnclass = "sf") %>%
  filter(geonunit %in% c("United States of America", "Canada", "Mexico")) 
na_states <- st_crop(na_states, ext(-125.0, -102, 31.1892, 49.4))
na_states_l <- st_cast(na_states, "MULTILINESTRING")

#### Climate data ----

# Western US (away) climate
# Daymet data are 30-year weekly averages (1991-2020) at 1-km2 resolution
# Data were exported to CSV format from daily Daymet gridded data at 1km2 resolution
# Processing Daymet rasters each time was too computationally intensive
# Data are too large to store for free at GitHub
# Column names are variable plus week number
# (e.g., "tmin01" = minimum temperature of week 1)
# However, the first, column, "ppt_ann" is total annual precipitation
away_clim <- read.csv(here("data", "weekly", "away", "all_away_dat.csv"))

# X and Y coords for away
away_xy <- select(away_clim, x, y)

# Remove XY coords for analysis
away_clim <- select(away_clim, -x, -y)

# Home climate
# E-OBS data that have been downscaled to a 1-km2 resolution
home_clim <- read.csv(here("data", "weekly", "home", "kilkis.csv"))
home_clim <- home_clim[names(away_clim)]

#### (1) CLIMEX  approach  ----

### (1a) Calculate indices ----
# Uses equations from the CLIMEX manual for the regional climate matching model. 
# See Appendix 1 of the JEE publication for details.

# Weeks in a year
weeks <- sprintf('%02d', c(1:52))

# Calculate indices for all variables
rtot_ind <- Index_calc("ppt_ann")
sm_ind <- Index_calc("sm")
tmax_ind <- Index_calc("tmax")
tmin_ind <- Index_calc("tmin")

# Combined indices (see pg. 100 of manual)
w1 <- 1 # tmin
w2 <- 1 # tmax
w3 <- 1 # rtot
w4 <- 1 # sm
wt <- 1 # max of w1 and w2
all_ind <- data.frame(
  tempi = ((tmax_ind$ind * w2) + (tmin_ind$ind * w1))/ (w1 + w2),
  rtoti = rtot_ind$ind * w3,
  smi = sm_ind$ind * w4,
  x = away_xy$x, 
  y = away_xy$y) 
all_ind <-  all_ind %>% 
  mutate(cmi = (tempi^wt * rtoti^w3 * smi^w4)^(1/(wt + w3 + w4))) %>%
  select(x, y, everything())

# CMI plot

# Create maps for all indices
indices <- c("tempi", "smi", "rtoti", "cmi")

# Get max values for indices to help define cutoffs
maxvals <- data.frame(
  rbind(
    lapply(1:length(indices) + 2, function(i) {
      m <- plyr::round_any(max(all_ind[,i]), 0.1, ceiling)
    }))
)

names(maxvals) <- indices
maxvals <- tidyr::unnest(maxvals, cols = c(tempi, smi, rtoti, cmi))

# Cut-off threshold (values below this aren't shown)
cutoff <- 0.5
all_ind_l <- all_ind %>%
  pivot_longer(cols = indices) %>%
  filter(value >= cutoff)

# (1b) Indices plots ----

# Apply plotting function to climate matching model outputs
# Composite match index
cmi_plot <- CMI_base_map(ind = "cmi") +
  labs(title = "Composite Match Index (CMI)") +
  scale_fill_gradientn(colours = GA::jet.colors(10), name = "CMI",
                       breaks = seq(cutoff, maxvals$cmi, 0.1), 
                       limits = c(cutoff, maxvals$cmi)) 

# Temperature index
tempi_plot <- CMI_base_map(ind = "tempi") +
  labs(title = "(A) Temperature Index") +
  scale_fill_gradientn(colours = rev(pals::kovesi.linear_kry_5_95_c72(10)), 
                       name = expression(bolditalic(I[t])),
                       breaks = seq(cutoff, maxvals$tempi, 0.1), 
                       limits = c(cutoff, 0.9)) 
# Moisture index
rtoti_plot <- CMI_base_map(ind = "rtoti") +
  labs(title = "(A) Moisture Index") +
  scale_fill_gradientn(colours = rev(pals::kovesi.linear_blue_5_95_c73(10)), 
                       name = expression(bolditalic(I[m])),
                       breaks = seq(cutoff, maxvals$rtoti,0.1), 
                       limits = c(cutoff,1.0))
# Soil moisture index
smi_plot <- CMI_base_map(ind = "smi") +
  labs(title = "(A) Soil Moisture Index") +
  scale_fill_gradientn(colours = rev(pals::ocean.solar(10)),  
                       name = expression(bolditalic(I[sm])),
                       breaks = seq(cutoff, maxvals$smi, 0.1), 
                       limits = c(cutoff,1.0)) 

# Save all plots individually
# Add date to file name
dat <- Sys.Date()
dat <- format(dat, "%m-%d-%y")

# List of files to save
indmaps <- list(tempi_plot, smi_plot, rtoti_plot, cmi_plot)
map(1:length(indices), function(i) {
  ggsave(indmaps[[i]], 
         filename = here("plots", "CLIMEX_custom", paste0(indices[i], "_", dat, ".png")), 
         width = 7, height = 7, units = c('in'), dpi = 300)
})

#### (2) Euclimatch ----
# An interface for performing climate matching using the Euclidean "Climatch" algorithm. 
# Functions provide a vector of climatch scores (0-10) for each location (i.e., grid cell) 
# within the recipient region, the percent of climatch scores >= a threshold value, 
# and mean climatch score. See: https://cran.r-project.org/web/packages/Euclimatch/index.html.

# Bioclim variables used in model
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO12 = Annual Precipitation
# BIO29 = Highest weekly moisture index
# BIO30 = Lowest weekly moisture index

# Get rasters for all bioclim variables
home_bio <- rast(
  list.files(here("data", "bioclim", "home"), pattern = ".tif$", full.names = TRUE))

away_bio <- rast(
  list.files(here("data", "bioclim", "away"),pattern = ".tif$", full.names = TRUE))

# Remove exponential notation and increase decimal lengths in R
options("scipen"=100, "digits"=10)

# Convert raster to a data frame and select desired columns
away_bio_df <- as.data.frame(away_bio_wus, xy=TRUE)
away_xy2 <- select(away_bio_df, x, y) # Need coordinates but then remove below
away_bio_df <- away_bio_df %>%
  select(bio5, bio6, bio12, bio29, bio30, -x, -y)

# Extract home (Kilkis)
home_extr <- extract(home_bio, matrix(kilkis, ncol=2)) %>%
  select(bio5, bio6, bio12, bio29, bio30)

# Calculate global variance of recipient climate data (all columns)
all_var <- away_bio_df %>% 
  summarise_all(list(var), na.rm = TRUE)
globvar <- unname(unlist(c(all_var)))

### (2a) Run Climatch ----
matchout <- climatch_vec(
  recipient = away_clim, 
  source = home_clim,
  globvar = globvar)

# Conduct analysis for individual bioclim variables and all variables
vars <- list("all", c("bio5","bio6"), "bio12", c("bio29","bio30"))

match_outs <- map(1:length(vars), function(i) {
  
  if (any(vars[[i]] == "all")) {
    # Calculate global variance of recipient climate data 
    globvar <- away_bio_df %>% 
      summarise_all(list(var), na.rm = TRUE)
    # Climatch algorithm
    matchout <- climatch_vec(
      recipient = away_bio_df, 
      source = home_extr,
      globvar = unlist(globvar))
  } else {
    # Calculate global variance of recipient climate data 
    globvar <- away_bio_df %>% 
      select(contains(vars[[i]])) %>%
      summarise_all(list(var), na.rm = TRUE)
    # Climatch algorithm
    matchout <- climatch_vec(
      recipient = select(away_bio_df, contains(vars[[i]])), 
      source = select(home_extr, contains(vars[[i]])),
      globvar = unname(unlist(c(globvar))))
  }
  
  # Convert values to categories
  match_df <- away_xy2 %>%
    mutate(value = matchout) %>%
    # Why are some match values crazy negative values??? -2147483648
    filter(value >= 0) %>%
    mutate(Climatch_fact = case_when(
      value == 0 ~ "0",
      value > 0 & value <= 1 ~ "0.1-1.0",
      value >= 1 & value <= 2 ~ "1.1-2.0",
      value >= 2 & value <= 3 ~ "2.1-3.0",
      value >= 3 & value <= 4 ~ "3.1-4.0",
      value > 4 & value <= 5 ~ "4.1-5.0",
      value > 5 & value <= 6 ~ "5.1-6.0",
      value > 6 & value <= 7 ~ "6.1-7.0",
      value > 7 & value <= 8 ~ "7.1-8.0",
      value > 8 & value <= 9 ~ "8.1-9.0",
      value > 9 & value <= 10 ~ "9.1-10"
    )) 
  # Factorize Climatch score
  match_df$Climatch_fact <- factor(
    match_df$Climatch_fact, levels = unique(match_df$Climatch_fact[order(match_df$value)]))
  
  return(match_df)
  
})

names(match_outs) <- c("all", "temp", "ppt",  "sm")

# Combine into a single data frame
cutoff2 <- 0.1
match_outs_df <- bind_rows(match_outs, .id = "variable") %>%
  filter(complete.cases(.)) %>%
  filter(value >= cutoff2) %>%
  mutate(Climatch = round(value, 2))

# (2b) Climatch plots ----

# All variables
allMatch_plot <- Climatch_base_map("all") +
  labs(title = "(B) Climatch - All Variables") +
  scale_fill_gradientn(
    colours = GA::jet.colors(10), 
    name = "Score", 
    breaks = c(cutoff2, 2, 4, 6, 8), 
    limits = c(cutoff2, 9)) 

# Temperature variables (bio5, bio6)
tempMatch_plot <- Climatch_base_map("temp") +
  labs(title = "(A) Temperature") +
  scale_fill_gradientn(
    colours = rev(pals::kovesi.linear_kry_5_95_c72(10)), 
    name = "Score", 
    breaks = c(cutoff2, seq(2,10,2)), 
    limits = c(cutoff2,10)) 

# Total precipitation (bio12)
rtotMatch_plot <- Climatch_base_map("ppt") +
  labs(title = "(B) Precipitation") +
  scale_fill_gradientn(
    colours = rev(pals::kovesi.linear_blue_5_95_c73(10)), 
    name = "Score", 
    breaks = c(cutoff2, seq(2,10,2)), 
    limits = c(cutoff2,10))

# Soil moisture variables (bio29, bio30)
smMatch_plot <- Climatch_base_map("sm") +
  labs(title = "(C) Soil Moisture") +
  scale_fill_gradientn(
    colours = rev(pals::ocean.solar(10)), 
    name = "Score", 
    breaks = c(cutoff2, seq(2,10,2)), 
    limits = c(cutoff2,10))

# Save all plots individually
matchmaps <- list(allMatch_plot, tempMatch_plot, rtotMatch_plot, smMatch_plot)
flnams <- c("allvars", "temp", "ppt", "sm")
map(1:length(flnams), function(i) {
  ggsave(matchmaps[[i]], 
         filename = here("plots", "Climatch", 
                         paste0("Climatch_", flnams[i], "_", dat, ".png")), 
         width = 7, height = 7, units = c('in'), dpi = 300)
})

### (3) Combine all plots ----

# Combine plots into one plot
match_1 <- plot_grid(tempMatch_plot, rtotMatch_plot)
match_2 <- plot_grid(smMatch_plot, allMatch_plot)
match_all <- plot_grid(match_1, match_2, nrow = 2) 
ggsave(match_all, 
       filename = here("plots", "Climatch", paste0("Climatch_4grid_", dat, ".png")), 
       width = 12, height = 8, units = c('in'), dpi = 300, bg = "white")

# Combine CMI plot and Climatch plot w/ all variables
final_all <- plot_grid(
  cmi_plot + 
    ggtitle("(A)") +
    theme(plot.margin = margin(t = 0.1, r = 0, b = 0, l = 0, unit = "cm"),
          plot.title = element_text(size = rel(2), vjust = 1.25, hjust = -0.15, face = "bold")), 
  allMatch_plot + 
    ggtitle("(B)") +
    theme(plot.margin = margin(t = 0.1, r = 0, b = 0, l = 0, unit = "cm"),
          plot.title = element_text(size = rel(2), vjust = 1.25, hjust = -0.15, face = "bold")), 
  nrow = 1)
ggsave(final_all, 
       filename = here("plots", 
                       paste0("CMI_Climatch_final_", dat, ".png")), 
       width = 11, height = 6, units = c('in'), dpi = 300, bg = "white")
# Plot for publication (Figure 1)
# ggsave(final_all, 
#        filename = here("plots", paste0("Figure2.tif")), 
#        width = 11, height = 6, units = c('in'), dpi = 300, bg = "white")

# Combine temperature index plot and Climatch plot for temperature variables
# Temperature
final_temp <- plot_grid(tempi_plot + ggtitle("(A)"), 
                        tempMatch_plot + ggtitle("(B)"), nrow = 1)
# ggsave(final_temp, 
#        filename = here("plots", paste0("TempIt_Climatch_final_", dat, ".png")), 
#        width = 12, height = 6, units = c('in'), dpi = 300, bg = "white")

# Combine moisture index plot and Climatch plot for annual precipitation 
# precipitation/moisture
final_rtot <- plot_grid(rtoti_plot + ggtitle("(A)"), 
                        rtotMatch_plot + ggtitle("(B)"), nrow = 1)
# ggsave(final_rtot, 
#        filename = here("plots", paste0("RainIm_Climatch_final_", dat, ".png")), 
#        width = 12, height = 6, units = c('in'), dpi = 300, bg = "white")

# Combine soil moisture index plot and Climatch plot for soil moisture variables
# Soil moisture
final_sm <- plot_grid(smi_plot + ggtitle("(A)"), 
                      smMatch_plot + ggtitle("(B)"), nrow = 1)
# ggsave(final_sm, 
#        filename = here("plots", paste0("SMIsm_Climatch_final_", dat, ".png")), 
#        width = 12, height = 6, units = c('in'), dpi = 300, bg = "white")

# Combine CLIMEX-based plots and Climatch plots
mytheme2 <- theme(legend.position = "right",
                  legend.key.height  =  unit(1.25, "cm"),
                  legend.key.width = unit(0.75, "cm"),
                  plot.title = element_text(size = 16, face = "bold", hjust = -0.25, vjust = -2),
                  axis.text = element_text(size = 10, face = "bold", colour = "black"),
                  plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0.5, unit = "cm"),
                  legend.title = element_text(size = 16, hjust = 0, vjust = 0.95, face = "bold"))

# Combine all 6 plots for publication
# All 3 climate factors (panels A to F)
final_temp2 <- plot_grid(tempi_plot + ggtitle("(A)") + mytheme2,
                         tempMatch_plot + ggtitle("(B)") + mytheme2, nrow = 1)
final_rtot2 <- plot_grid(rtoti_plot + ggtitle("(C)") + mytheme2,
                         rtotMatch_plot + ggtitle("(D)") + mytheme2, nrow = 1)
final_sm2 <- plot_grid(smi_plot + ggtitle("(E)") + mytheme2,
                       smMatch_plot + ggtitle("(F)") + mytheme2, nrow = 1)
final_3fact <- plot_grid(final_temp2, final_rtot2, final_sm2, nrow = 3)
ggsave(final_3fact, 
       filename = here("plots", paste0("ThreeFactors_final_", dat, ".png")), 
       width = 11, height = 12, units = c('in'), dpi = 300, bg = "white")
# TIF file for publication
# ggsave(final_3fact, 
#        filename = here("plots", paste0("Figure3.tif")), 
#        width = 11, height = 12, units = c('in'), dpi = 300, bg = "white")

# (4) Convert results to rasters ----

# Template
template <- rast(extent = ext(away_bio_wus), 
                 ncol=ncol(away_bio_wus), 
                 nrow = nrow(away_bio_wus))

# Loop through output list, rasterize, and save
inds <- names(select(all_ind, -x, -y))

# Custom CMI
for (ind in inds) {
  # Rasterize
  df <- select(all_ind, x, y, ind)
  r <- rast(df, extent = ext(template))
  # Save raster
  writeRaster(
    r, filename = here("raster_outputs", "CLIMEX_custom",
                       paste0(ind, "_CLIMEX_custom_", dat, ".tif")), overwrite=TRUE)
}

# Climatch
for (i in 1:length(match_outs)) {
  nam <- names(match_outs[i])
  df <- match_outs[[i]] %>% 
    select(-Climatch_fact)
  r <- rast(df, extent = ext(template))
  # Save raster
  writeRaster(
    r, filename = here("raster_outputs", "Climatch",
                       paste0(nam, "_Climatch_", dat, ".tif")), overwrite=TRUE)
}
