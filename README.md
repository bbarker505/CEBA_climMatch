# **Climate matching models for *Ceratapion basicorne* (Coleoptera: Apionidae), a biocontrol agent of yellow starthistle**

Access this dataset on Dryad: <https://doi.org/10.5061/dryad.vmcvdnd2w>

👥 Brittany Barker

Questions? 📧
[bbarker505\@gmail.com](mailto:bbarker505@gmail.com){.email} or
[brittany.barker\@oregonstate.edu](mailto:brittany.barker@oregonstate.edu){.email}

## **Purpose**

The code and files in this repository produce climate matching models
for *Ceratapion basicorne* for the western U.S. (Barker et al. 2025).

## **Abstract**

*Ceratapion basicorne* (Illiger) (Coleoptera: Apionidae), a weevil
native to Europe and western Asia, shows promise for enhancing the
control of yellow starthistle (*Centaurea solstitialis* L.), an invasive
annual forb in the western U.S. However, a paucity of data on this
biocontrol agent’s environmental constraints has made it difficult to
assess the suitability of potential release locations. Climate matching
models were developed for *C. basicorne* to help identify areas of the
western U.S. with similar climates to the source area of breeding
colonies being used for releases (home location). The models used
climate variables derived from daily estimates of minimum temperature
(Tmin), maximum temperature (Tmax), precipitation, and soil moisture for
a 30-year period spanning 1991−2020 at 1-km^2^ resolution. Models
indicated that the Central California Foothills, Eastern Cascades
Foothills, Columbia Plateau, and mountainous parts of northcentral Utah
had the most similar climates to the home location. Of these areas, the
Eastern Cascades foothills in northeastern California and Wasatch Range
in Utah occurred at a similar latitude as the home location, which may
be important to consider if *C. basicorne* has photoperiodic diapause.
The least similar climates occurred in very wet coastal areas,
high-elevation (cold) mountains, and hot deserts. The development of
process-based models for predicting the establishment of this agent will
require a more detailed understanding of the agent’s requirements for
development and survival.

## **Description of the data and file structure**

### **Subfolders**

The main folder (directory) contains the following subfolders:

-   `data`: all data needed for modeling (see descriptions below).
-   `plots`: image files (`.png` files) produced by models are saved
    here. Images produced by the CLIMEX-based model and Climatch models
    are saved to the `CLIMEX_custom` and `Climatch` folder,
    respectively.
-   `raster_outputs`: rasters (`.tif` files) produced by models are
    saved here. Rasters produced by the CLIMEX-based model and Climatch
    models are saved to the `CLIMEX_custom` and `Climatch` subfolder,
    respectively.
-   `script`: running the R script `CEBA_climMatch.R` re-produces all
    models, figures, and rasters presented in the manuscript (see
    below).

### **Data**

Climate data for modeling: the CLIMEX-based model uses weekly climate
data (`/data/weekly/`) whereas the Climatch model uses six bioclimatic
variables (`/data/bioclim/`). For each model type, the data are split
according to the 'home' location (Kilkis, Greece) and 'away' locations
in the western U.S. Averages of climate for a 30-year period spanning
`1991-2020` were derived from the Daymet and E-OBS datasets. In raster
datasets, all locations in the western U.S. (extent:
`xmin = -125.0024, xmax = -101.9945,ymin = 31.18652, ymax = 49.40551`)
are at `1-km2` resolution and the coordinate reference system is WGS84
(EPSG:4326). Units are degrees Celsius for temperature (Tmin, Tmax),
millimeters for precipitation, and `m3 m-3` for soil moisture. Further
details about the climatic data and variables used in models can be
found in the publication for this study.

A summary of these datasets is below:

-   `/data/weekly/home/kilkis.csv`: Data for the 'home' location used
    for the CLIMEX-based model. The columns correspond to averages of
    weekly climate for Tmin, Tmax, and soil moisture (sm) (`52` weeks
    `x` `3` variables `=` `156` columns) (e.g., `ppt01` = precipitation
    of the first week of the year), as well as average annual total
    precipitation (mm) (`ppt_ann`).
-   `/data/weekly/away/*tif:` Data for the 'away' locations used for
    CLIMEX-based models. There are `156` rasters corresponding to weekly
    estimates of Tmin, Tmax, and soil moisture (sm) (`52` weeks `x` `3`
    variables `=` `156` rasters), as well as an annual estimate of
    precipitation (`ppt_ann`).
-   `/data/bioclim/home/kilkis.csv`: Data for the 'home' location used
    for the Climatch models. The columns correspond to bioclimatic
    variables derived from 30-year averages of climate. These included
    Tmax of the warmest month (`bio5`), Tmin of the coldest month
    (`bio6`), annual precipitation (`bio12`), highest monthly soil
    moisture (`bio29`), and lowest monthly soil moisture (`bio30`). The
    row corresponds to the source location for *C. basicorne* in Kilkis,
    Greece (`latitude = 22.844, longitude = 40.994`).
-   `/data/bioclim/away/*tif:` Data for 'away' locations used for
    Climatch models (same bioclimatic variables as the home location).
-   `/data/YST_counties/YST_counties_9-19-24.shp:` Geospatial vector
    data of counties in the western U.S. where yellow starthistle has
    been observed. Observations were obtained from the Early Detection
    and Distribution Mapping System (EDDMapS, 2024) and the Global
    Biodiversity Information Facility (GBIF.org, 2024). A single record
    was retained for each county, for a total of `182` records from `11`
    states. The geometry type is a multipolygon (extent:
    `xmin = -125.0024, xmax = -101.9945,ymin = 31.18652, ymax = 49.40551`)
    and the coordinate reference system is NAD83. The columns include:
    `STATEFP` (state FIPS code), `COUNTYFP` (U.S. county FIPS code),
    `NAME` (county name), `geometry` (longitude, latitude), `n_number`
    (number of records), and `present` [if starthisle has been detected
    (present) or not (absent)].

## **Sharing/Access information**

Data were derived from the following sources:

-   [Daymet dataset for North America](https://daymet.ornl.gov/getdata)
-   [E-OBS dataset for Europe](https://surfobs.climate.copernicus.eu)
-   [SiTHv2 global soil moisture
    dataset](https://data.tpdc.ac.cn/en/data/bc51e1b0-494c-4cd5-ae4d-eba6b9d2322c)
-   [EDDMapS database](https://www.eddmaps.org/)
-   [GBIF \| Global Biodiversity Information
    Facility](https://www.gbif.org/)

Citations for datasets can be found in the publication (Barker et al.
2025)

## **Code/Software**

The R statistical software (version 4.3.2) was used to produce models.

### **Required R packages**

The following packages must be installed:
`here, tidyverse, ggspatial, rnaturalearth, sf, terra, tidyterra, Euclimatch, GA, cowplot.`

### **Instructions**

(1) Clone the repository (or download the directory from Dryad). Don't
    move or delete any subfolders or datasets.\
(2) Open the R project (`CEBA_climMatch.Rproj`) in RStudio/Posit.\
(3) Open R script named `CEBA_climMatch.R` in the `script` subfolder.\
(4) Install any necessary packages as listed above.\
(5) Run the R script, which produces the climate matching models and
    associated outputs.

## **Acknowledgements**

This work was funded by the U.S. Department of Defense Strategic
Environmental Research and Development Program (U.S. Army Corps of
Engineers, contract no.
[RC23-3611](https://demo.serdp-estcp.mil/projects/details/48787c60-7e33-4cda-bd02-5bc31b402265/managing-yellow-starthistle-using-a-new-biocontrol-agent-an-integrative-experimental-and-geo-climatic-modeling-approach)).

## **References**

Barker, B. S. 2025. Climate matching models for *Ceratapion basicorne*
(Coleoptera: Apionidae), a biocontrol agent of yellow starthistle.
Journal of Economic Entomology. <http:/doi.org/10.1093/jee/toae299>.

EDDMapS. 2024. Early Detection & Distribution Mapping System. The
University of Georgia - Center for Invasive Species and Ecosystem
Health. Available online at: <http://www.eddmaps.org> (accessed 19
September 2024).

GBIF.org. 2024. (19 September 2024) GBIF Occurrence Download
<https://doi.org/10.15468/dl.sqmyns>.
