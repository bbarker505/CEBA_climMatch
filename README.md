# Climate matching models for *Ceratapion basicorne* (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle

👥 Brittany Barker 

Questions? 📧 bbarker505@gmail.com or brittany.barker@oregonstate.edu  

# Purpose

The code and files in this repository reproduce the outputs presented in a climate matching modeling study of *Ceratapion basicorne* (manuscript in review). 

# Abstract

*Ceratapion basicorne* (Illiger) (Coleoptera: Apionidae), a weevil native to Europe and western Asia, shows promise for enhancing the control of yellow starthistle (*Centaurea solstitialis* L.), an invasive annual forb in the western U.S. However, a paucity of data on this biocontrol agent’s environmental constraints has made it difficult to assess the suitability of potential release locations. Climate matching models were developed for *C. basicorne* to help identify areas of the western U.S. with similar climates to the source area of breeding colonies being used for releases (home location). The models used climate variables derived from daily estimates of minimum temperature, maximum temperature, precipitation, and soil moisture for a 30-year period spanning 1991−2020 at 1-km2 resolution. Models indicated that the Central California Foothills, Eastern Cascades Foothills, Columbia Plateau, and mountainous parts of northcentral Utah had the most similar climates to the home location. Of these areas, the Eastern Cascades foothills in northeastern California and Wasatch Range in Utah occurred at a similar latitude as the home location, which may be important to consider if *C. basicorne* has photoperiodic diapause. The least similar climates occurred in very wet coastal areas, high-elevation (cold) mountains, and hot deserts. The development of process-based models for predicting the establishment of this agent will require a more detailed understanding of the agent’s requirements for development and survival.

# Inputs and outputs

Model inputs are climate data in the "data" subfolder. However, a very large data file must be downloaded from Zenodo because it's too big for GitHub (see instructions below). The CLIMEX-based model used weekly climate data whereas the Climatch model used six bioclimatic variables. For each model type, the data are split according to the 'home' (Kilkis, Greece) and 'away' (western U.S.) locations. Weekly data were converted from raster to CSV format prior to modeling to reduce computational load. Conversely, bioclim data are still in raster (GeoTIFF) format. Details about the climatic data and variables used in models can be found in the publication for this study.

Model outputs are static image files (.PNG) and rasters (GeoTIFF) that are saved to the "plots" and "raster_outputs" subfolders, respectively. Outputs are further organized by model type ("CLIMEX_custom" vs. "Climatch"). 

# Required R packages

The follwing packages must be intalled:
- here
- tidyverse
- ggspatial
- rnaturalearth
- sf
- terra
- tidyterra
- Euclimatch
- GA
- cowplot

# Instructions

(1) Clone the repository. Don't move or delete any subfolders or data. 
(2) Download the climate data for the western U.S. ("all_away_dat.csv") from Zenodo: https://zenodo.org/records/13770391.  
(3) Place "all_away_dat.csv" in the following directory of the repository: /data/data/weekly/away/  
(4) Open the R project ("CEBA_climMatch.Rproj") in RStudio/Posit.   
(5) Open R script named "CEBA_climMatch.R" in the "script" subfolder.  
(6) Install any necessary packages as listed above.    
(7) Run the R script.

# Acknowledgements

This work was funded by the U.S. Department of Defense Strategic Environmental Research and Development Program (U.S. Army Corps of Engineers, contract no. RC23-3611) and the USDA National Institute of Food and Agriculture (NIFA) Agriculture and Food Research Initiative (AFRI) grant no. 2022-68013-37138.

# References

Barker, B. S. 2024. Climate matching models for *Ceratapion basicorne* (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle. In review.

