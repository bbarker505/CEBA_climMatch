# Climate matching models for *Ceratapion basicorne* (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle

👥 Brittany Barker 

Questions? 📧 bbarker505@gmail.com or brittany.barker@oregonstate.edu  

# Purpose

The code and files in this repository produce two climate matching models for *Ceratapion basicorne* (Barker et al. 2025 In Press). 

# Abstract

*Ceratapion basicorne* (Illiger) (Coleoptera: Apionidae), a weevil native to Europe and western Asia, shows promise for enhancing the control of yellow starthistle (*Centaurea solstitialis* L.), an invasive annual forb in the western U.S. However, a paucity of data on this biocontrol agent’s environmental constraints has made it difficult to assess the suitability of potential release locations. Climate matching models were developed for *C. basicorne* to help identify areas of the western U.S. with similar climates to the source area of breeding colonies being used for releases (home location). The models used climate variables derived from daily estimates of minimum temperature, maximum temperature, precipitation, and soil moisture for a 30-year period spanning 1991−2020 at 1-km2 resolution. Models indicated that the Central California Foothills, Eastern Cascades Foothills, Columbia Plateau, and mountainous parts of northcentral Utah had the most similar climates to the home location. Of these areas, the Eastern Cascades foothills in northeastern California and Wasatch Range in Utah occurred at a similar latitude as the home location, which may be important to consider if *C. basicorne* has photoperiodic diapause. The least similar climates occurred in very wet coastal areas, high-elevation (cold) mountains, and hot deserts. The development of process-based models for predicting the establishment of this agent will require a more detailed understanding of the agent’s requirements for development and survival.

# Inputs and outputs

Data inputs for the study ("data" subfolder) include:  
(1) Climate data for modeling: the CLIMEX-based model uses weekly climate data ("/data/weekly/*") whereas the Climatch model uses six bioclimatic variables ("/data/bioclim/*"). For each model type, the data are split according to the 'home' (Kilkis, Greece) and 'away' (western U.S.) locations. Data for the 'home' location were extracted from rasters for Europe whereas data for the 'away' locations (western U.S.) are in raster format. Details about the climatic data and variables used in models can be found in the publication for this study.  
(2) Locality records: records for *C. basicorne* in Europe and western Asia ("/data/records/CEBA_all_records.xlsx") were derived from the literature, GBIF, and field studies. The "Estimated" column indicates whether coordinate information for the record were precise (Estimated = 0) or approximated from geographic information, such as a city name (Estimated = 1).  
(3) County-level detections for yellow starthistle ("/data/YST_counties/YST_counties_9-19-24.shp") were estimated from GBIF and EDDMaPS data.  

Model outputs are static image files (.PNG) and rasters (GeoTIFF) that are saved to the "plots" and "raster_outputs" subfolders, respectively. Outputs are further organized by model type ("CLIMEX_custom" vs. "Climatch"). 

# Required R packages

The following packages must be intalled:
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
(2) Open the R project ("CEBA_climMatch.Rproj") in RStudio/Posit.   
(3) Open R script named "CEBA_climMatch.R" in the "script" subfolder.  
(4) Install any necessary packages as listed above.    
(5) Run the R script, which produces the two climate matching models and associated figures.  
(6) (Optional) The "CEBA_localityMap_nativeRange.R" script produces a map of locality records for *C. basicorne*.  

# Acknowledgements

This work was funded by the U.S. Department of Defense Strategic Environmental Research and Development Program (U.S. Army Corps of Engineers, contract no. RC23-3611).

# References

Barker, B. S. 2025. Climate matching models for *Ceratapion basicorne* (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle. Journal of Economic Entomology. In Press.
