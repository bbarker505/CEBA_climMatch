# R script to create a map of locality records for C. basicorne in the native range
# Author: Brittany Barker, Oregon IPM Center, Oregon State University
# Associated publication: Barker (2025) Climate matching models for Ceratapion 
# basicorne (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle. 
# Journal of Economic Entomology.
#
# Load packages
library(here)
library(tidyverse)
library(ggspatial) # project hline in ggplot
library(rnaturalearth)
library(sf)

# Records
pts_sf <- readxl::read_excel(
  here("data", "records", "CEBA_all_records.xlsx")) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
kilkis <- filter(pts_sf, Site == "Kilkis")

# Europe features
sf_use_s2(FALSE)
europe_c <- st_crop(ne_countries(returnclass = "sf", scale = 10), 
                xmin = -12, xmax = 52.5, ymin = 32.5, ymax = 65)

# Map
# Color vs. black and white
p <- ggplot() +
  geom_sf(data = europe_c, color = "gray10", fill = "gray93") +  
  #geom_sf(data = pts_sf,  color = "white", fill = "gray30", size = 2, shape = 21) +
  #geom_sf(data = kilkis, color = "white", fill = "black", size = 3, shape = 23) +
  geom_sf(data = pts_sf,  color = "black", fill = "#007FFF", size = 1.5, shape = 21) +
  geom_sf(data = kilkis, color = "black", fill = "red", size = 2, shape = 24) +
  coord_sf(crs = "EPSG:3857") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", 
                                    linewidth = 0.75, fill = NA),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, face = "bold", colour = "black"))

# Save
ggsave(p, 
       filename = here("plots", "Fig1_CEBA_range_map.png"), 
       width = 7, height = 5, units = c('in'), dpi = 300, bg = "white")
knitr::plot_crop(here("Records",  "plots", "All_records_9-13-24.png"))
# Figure for publication
# ggsave(p, 
#        filename = here("plots",  "Figure1.tiff"), 
#        width = 7, height = 5, units = c('in'), dpi = 300, bg = "white")
# knitr::plot_crop(here("Records",  "plots", "Figure1.tiff"))
