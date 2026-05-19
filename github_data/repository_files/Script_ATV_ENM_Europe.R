# ------------------------------------------------------------------------------
# Script name: Script_ATV_ENM_Europe
#
# Purpose: Ecological Niche modelling analyses for "Escalating human exposure to 
#          Aedes mosquito-borne viruses in Europe" Study 
#
# Author(s): Kyla Serres
# ------------------------------------------------------------------------------

# Load or install necessary libraries ------------------------------------------

# Geospatial libraries
if (!require(ncdf4)) install.packages("ncdf4")
if (!require(exactextractr)) install.packages("exactextractr")
if (!require(raster)) install.packages("raster")
if (!require(sf)) install.packages("sf")


# Data manipulation libraries
if (!require(dplyr)) install.packages("dplyr")
if (!require(purrr)) install.packages("purrr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(stringr)) install.packages("stringr")

# Statistical modeling libraries
if (!require(gmodels)) install.packages("gmodels")
if (!require(devtools)) install.packages("devtools")
if (!require(dismo)) install.packages("dismo")
if (!require(gbm)) install.packages("gbm")
if (!require(blockCV)) install.packages("blockCV")

# Plotting libraries
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(scales)) install.packages("scales")

# Parallel processing libraries
if (!require(parallel)) install.packages("parallel")
if (!require(foreach)) install.packages("foreach")
if (!require(future.apply)) install.packages("future.apply")

# ------------------------------------------------------------------------------
#                                   SUMMARY
# ------------------------------------------------------------------------------

    ## I. DATA PREPARATION
        # 1. Load optimized European administrative area shapefile (VectorNet)  
        # 2. Aedes-borne (DENV, CHIKV, ZIKV) viral occurrence data curation

    ## II. BRT MODEL TRAINING AND EVALUATION
        # 1. Pseudo-absence sampling
        # 2. Environmental data extraction for BRT training
        # 3. Investigating the environmental variables to include in the model
          # 3.1 Building correlation matrix among predictors (Pearson correlations)
          # 3.2 Investagating spatial autocorrelation and defining block size for spatial cross-validation
          # 3.3 Comparing the relative importance (RI) of each variable with the RI of its shadow version
        # 4. BRT model training
          # 4.1 hyperparameter tuning on selected set of variables
          # 4.2 Final BRT training, with selected variables and hyperparameters, and spatial cross-validation (SCV)
          # 4.3 Predictive performance table

    ## III. VARIABLE RELATIVE IMPORTANCE AND RESPONSE CURVES
        # 1. Relative importances 
        # 2. Response Curves

    ## IV. CURRENT RISK MAP (2020-2024)
        # 1. Data extraction
        # 2. Current predictions (t0)

    ## V. PAST PROJECTIONS
        # 1. Data extraction 
        # 2. Historical predictions
        # 3. Saving predictions
        # 4. Historical predictions figures

    ## VI. FUTURE PREDICTIONS
        # 1. subset and aggregate future landuse
        # 2. future data extraction
        # 3. future projections
          # 3.1 Saving future projections
          # 3.2 extracting population counts and saving future population counts per NUTS3 area
        # 4. future predictions figures
        # 5. Quantifying the population at risk of exposure
        # 6. Estimating the ecological suitability variation at the continent scale
        # 7. Estimating the population at risk variation at the continent scale


# I. DATA PREPARATION ----------------------------------------------------------

  # 1. Load optimized European administrative area shapefile (VectorNet)  
  # --------------------------------------------------------------------

nutsM = crop(shapefile("Shapefiles/MOOD_EUROPE_NUTS3_WillyW/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # modified NUT3 shapefile
# Only keep EU areas
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
# remove the small islands where land-use data is not available
nutsM = subset(nutsM, !LocationCo%in%c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533","ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65")) 

nutsM_sf = st_as_sf(nutsM)
nutsM_sf$area = nutsM_sf %>% st_area %>% as.numeric/(1000*1000) ; #st_write(nutsM_sf, "Shapefiles/nutsM.gpkg")

# Correspondences from original NUTS3 to optimized admin. areas (NUTSM)
correspondences = shapefile("Shapefiles/MOOD_EUROPE_NUTS3_WillyW/vnMOODdatamapcodejoinpoly.shp")@data[,c("DATLOCODE","MAPLOCODE")]

# Subset polygon IDs
nutsM_IDs = as.data.frame(matrix(nrow = dim(nutsM@data)[1], ncol = 3))
colnames(nutsM_IDs) = c("country","region_name", "nutM_code"); nutsM_IDs$nutM_code = nutsM@data[,"LocationCo"]
nutsM_IDs$country = nutsM@data[,"CountryNam"]; nutsM_IDs$region_name = nutsM@data[,"LocationNa"]



  # 2. Aedes-borne (DENV, CHIKV, ZIKV) viral occurrence data curation
  # ------------------------------------------------------------------

curatingData = FALSE
if(curatingData== TRUE) {
  ### ECDC sources: the data comes from 2 sources, a 2024 ECDC data request ###
  #       and a 2025 ECDC ATLAS surveillance platform data download

  # ECDC 2024 request data: dengue, chikungunya, west nile
  denv = read.csv("Occurrence_data/ECDC/Request_jul_24/DENGUE.csv", header = T, sep = ","); denv = denv[which(denv[,"Imported"] == "N"), ] # remove imported cases
  germany = c("DE216", "DE40C"); denv = denv[!apply(denv, 1, function(row) any(row %in% germany)), ] # these are not locally acquired vector borne cases. 
  denv$Place.of.Infection[denv$Place.of.Infection == "ES"] = denv$PlaceOfNotification[denv$Place.of.Infection == "ES"] # Replace "with the corresponding NUTS3 region from PlaceOfNotification, (since "ES" is a general label for Spain and the notification has the specific region)
  denv_subset = data.frame(year = denv$DateUsedForStatisticsYear, disease = "dengue", region_name = NA, region_code = denv$Place.of.Infection, case_numbers = NA, stringsAsFactors = FALSE)

  chikv = read.csv("Occurrence_data/ECDC/Request_jul_24/CHIK.csv", header = T, sep = ","); chikv = chikv[which(chikv[,"Imported"] == "N"), ] # remove imported cases
  chikv_subset = data.frame(year = chikv$DateUsedForStatisticsYear, disease = "chikungunya", region_name = NA, region_code = chikv$Place.of.Infection, case_numbers = NA, stringsAsFactors = FALSE)

  ### ECDC ATLAS platform data: dengue, chikungunya, west nile 
  Denv = read.csv("Occurrence_data/ECDC/Atlas_feb_25/ECDC_surveillance_data_Dengue.csv", header = T, sep = ",")
  Denv = Denv[!apply(Denv, 1, function(row) any(row %in% germany)), ] # these are not locally acquired vector borne cases. 
  Denv_subset = data.frame(year = Denv$Time, disease = "dengue", region_name = Denv$RegionName, region_code = Denv$RegionCode, case_numbers = Denv$NumValue, stringsAsFactors = FALSE)
  
  Chikv = read.csv("Occurrence_data/ECDC/Atlas_feb_25/ECDC_surveillance_data_Chikungunya_virus_disease.csv", header = T, sep = ",")
  Chikv_subset = data.frame(year = Chikv$Time, disease = "chikungunya", region_name = Chikv$RegionName, region_code = Chikv$RegionCode, case_numbers = Chikv$NumValue, stringsAsFactors = FALSE)

  ecdc_pres = rbind(denv_subset, Denv_subset, chikv_subset, Chikv_subset); ecdc_pres$source = "ECDC"
  
  #### Literature review source: dengue, chikungunya, zika
    # Data from the literature will also be added to the occurrences. Specifically for DENV, CHIKV and ZIKV
    # Lit. data: Canttaneo et al. (2025) 
  
  literature = read.csv("Occurrence_data/Cattaneo_2025/Nuts3_cattaneo.csv", header = T, sep = ",") 
  literature[nrow(literature) + 1, ] = NA
  # add Croatian local case, article published in march 2025 by Medić et al. 2025 microorganisims 
  literature[nrow(literature), "country"] = "Croatia";literature[nrow(literature), "disease"] = "Dengue"
  literature[nrow(literature), "year"] = 2024; literature[nrow(literature), "location"] = NA 
  literature[nrow(literature), "NUTS3"] = "HR033"; literature[nrow(literature), "case_number"] = NA; literature[nrow(literature), "authors"] = "Medić et al. 2025"
  
  literature2 = literature[,c("year", "disease", "location","NUTS3","case_number")]
  colnames(literature2) = c("year","disease", "region_name", "region_code", "case_numbers")
  literature2$source = "Cattaneo et al. 2025"

  #### Santé Publique France (SPF) aedes viruses data: dengue, chikungunya, zika
  atv.spf = read.csv("Occurrence_data/SPF/SPF_ATV_nuts3_cases.csv", header = T) 
  atv.spf = atv.spf[, c("year","disease", "nuts3_name", "nuts3_codes", "number_of_cases")]
  colnames(atv.spf) = c("year","disease", "region_name", "region_code", "case_numbers"); atv.spf$source = "Santé Publique France"
  
  #### PADI-web
    # Padi-web curation

  # load included articles
articles = read.csv("Included_articles.csv", header = T, sep=",")
ids = articles$article_id

  # load location extracted info
denv_loc = read.csv("location_extracted_information_den_europe.csv", header = T, sep=",")
denv_loc$disease = "Dengue"
chik_loc = read.csv("location_extracted_information_chik_europe.csv", header = T, sep = ",")
chik_loc$disease = "Chikungunya"

  # subset locations to match only included articles
denv = subset(denv_loc, article %in% ids)
chik = subset(chik_loc, article %in% ids)

padi_web_loc = rbind(denv, chik)

  # curate locations to include only those in departments
all_locations = unique(padi_web_loc$location)
toremove = c("France", "the West Indies", "Santé", "Guadeloupe", "Actu", "Yogyakarta", 
             "UK", "Kraemer", "the United Kingdom", "Netherlands", "the Netherlands", 
             "Zika", "this?Europe", "Europe", "Berlin", "Rome", "Italy", "Zurich", "Le", 
             "Dengue", "Pacific", "Australia", "America", "Latin", "the Middle East", 
             "India", "China", "Southeast Asia", "Den-1",  "Sicily", "The West Indies", 
             "North", "E.L.", "E.N.", "F.M.", "F.V.", "A.C.", "A.D.", "G.S.", "E.S.", 
             "D.L.", "C.M.", "S.V.", "Appendix", "Remove", "US", "All Saints", 
             "Saint-Pierre", "the United States", "Martinique", "Loire-Atlantique", 
             "West Indies", "Val-de-Marne", "Nandy", "Outremer", "Normandy", "Brittany", 
             "Asia", "Casinalbo", "ed.", "Roman", "West Nile", "Mozambique", "Tanzania", 
             "knuckles", "United States", "Africa",  "Ville", "the District Areas", 
             "Mediterranean",  "Indonesia", "Latin America", "Arezzo", "Veneto", "Cesena", 
             "Trieste", "Garda", "Lake Garda", "dengue", "Thaïlande", "Americas", 
             "Caribbean", "Baho", "Florac", "Aimargues", "Montpellier-Pérols", 
             "Central and South America", "Nature", "Eur", "Virginia", "Baden-Württemberg",
             "Bavaria", "Germany", "German", "South America", "south-east", "Thailand", 
             "the Region\n\nCases", "Brazil", "Giuseppe", "Serenissima", "Via", "Malatesta",
             "Erbe", "Benvenuti", "San Carlo", "Pertini", "Pays-de-la-Loire", "Nouvelle-Aquitaine", 
             "Valneva", "Kimakonde", "Petra", "Makonde", "Largo", "Raggi", "Lillà", "Rome Virginia")

all_locations = subset(padi_web_loc, !(location %in% toremove))
coordinates_df = data.frame(year = all_locations$publication_date, location = all_locations$location, x= all_locations$longitude, y = all_locations$latitude, disease = all_locations$disease)
coords = coordinates_df
coordinates(coords) = ~x + y
points_sf = st_as_sf(coordinates_df, coords = c("x", "y"), crs = st_crs(nutsM_sf))
nutsM_joined = st_join(points_sf, nutsM_sf, join = st_intersects)

# force these regions to the correct Nuts3 area
polygon_ITI43 = nutsM_sf %>% filter(LocationCo == "ITI43") %>% st_drop_geometry()
polygon_FRL03 = nutsM_sf %>% filter(LocationCo == "FRL03") %>% st_drop_geometry()
cols = names(polygon_ITI43)

nutsM_joined = nutsM_joined %>%
  mutate(across(all_of(cols), 
                ~ if_else(location %in% c("Anzio", "the Municipality of Rome"), polygon_ITI43[[cur_column()]][1],
                          if_else(location == "Paca", polygon_FRL03[[cur_column()]][1], .))))

padi_web_nuts = data.frame(year = nutsM_joined$year, disease= nutsM_joined$disease, location =nutsM_joined$LocationNa, nuts_code = nutsM_joined$LocationCo)
padi_web_nuts = padi_web_nuts %>% mutate(year = paste0("20", substr(year, 1, 2)))
write.csv(padi_web_nuts, "padi_wed_nutsM.csv", row.names = FALSE)

padi_web_nuts$case_numbers = NA
padi_web_nuts = padi_web_nuts[,c("year", "disease", "location","nuts_code","case_numbers")]
colnames(padi_web_nuts) = c("year","disease", "region_name", "region_code", "case_numbers"); padi_web_nuts$source = "PADI-Web"


#### Compile aedes transmitted viruses together and combine all wnv presences
Aedes_presences = rbind(ecdc_pres, literature2, atv.spf, padi_web_nuts)
Aedes_presences$disease = tools::toTitleCase(tolower(Aedes_presences$disease))

# Some region codes are old versions and need to be updated to the new ones
Aedes_presences$region_code[Aedes_presences$region_code == "FR812"] = "FRJ12"
Aedes_presences$region_code[Aedes_presences$region_code == "FR813"] = "FRJ13"
Aedes_presences$region_code[Aedes_presences$region_code == "FR815"] = "FRJ15"
Aedes_presences$region_code[Aedes_presences$region_code == "FR823"] = "FRL03"
Aedes_presences$region_code[Aedes_presences$region_code == "FR824"] = "FRL04"
Aedes_presences$region_code[Aedes_presences$region_code == "FR825"] = "FRL05"
Aedes_presences$region_code[Aedes_presences$region_code == "FR826"] = "FRL06"
Aedes_presences$region_code[Aedes_presences$region_code == "FR831"] = "FRM01"
Aedes_presences$region_code[Aedes_presences$region_code == "HU101"] = "HU110"
Aedes_presences$region_code[Aedes_presences$region_code == "HU102"] = "HU120"
Aedes_presences = Aedes_presences[Aedes_presences$region_code != "NULL", ] # remove where there is no info
Aedes_presences = Aedes_presences[Aedes_presences$region_code != "ES531", ] # remove because it's an island, too small no environmental data for this area


#### transform presences to modified nuts 3 regions (vectornet map)
Aedes_presences$nutM_code = NA  
if (!"nutM_code" %in% colnames(Aedes_presences)) {Aedes_presences$nutM_code = NA}
  valid_rows = logical(nrow(Aedes_presences))  
  
  for (j in seq_len(nrow(Aedes_presences))) {
    nut3_code = Aedes_presences[j, "region_code"]
    
# Check if nut3_code is not NA and exists in nutsM_IDs
    if (!is.na(nut3_code) && nut3_code %in% nutsM_IDs) {
      index = which(nutsM_IDs == nut3_code)
      if (length(index) == 1) {
        Aedes_presences[j, "nutM_code"] = nutsM@data[index, "LocationCo"]
        valid_rows[j] = TRUE
      }
    } else {
      # Check against the correspondences table
      nutsM_ID = correspondences[correspondences$DATLOCODE == nut3_code, "MAPLOCODE"]
      
      if (length(nutsM_ID) > 0 && !is.na(nutsM_ID)) {
        Aedes_presences[j, "nutM_code"] = nutsM_ID
        valid_rows[j] = TRUE
      } else {
        # If no match is found
        print(paste("No match found for NUT3:", nut3_code))
      }
    }
  }
  
  # Keep only rows with valid nutM_code and remove rows with NA values
  Aedes_presences = Aedes_presences[valid_rows, ]
  # Remove the region_code column
  Aedes_presences = Aedes_presences[, !(colnames(Aedes_presences) == "region_code")]
  Aedes_presences$year = as.numeric(Aedes_presences$year)
  write.csv(Aedes_presences, "github_data/Occurrence_data/Aedes_viruses_European_presences.csv", row.names = FALSE)
  } else {
  Aedes_presences = read.csv("github_data/Occurrence_data/Aedes_viruses_European_presences.csv", header = T)
}

# clean environment
keep = c("Aedes_presences", "nutsM", "contour","nutsM_sf","nutsM_IDs") 
rm(list = setdiff(ls(), keep))


# II. BRT MODEL TRAINING AND EVALUATION ----------------------------------------

  # 1. Pseudo-absence sampling 
  # ---------------------------

# Group presence data by year ranges
# To avoid over weighting admin areas with repeated presences across years, 
# we keep only the year range per location. This allows averaging 
# environmental variables over that period.


# create the grouping function
years_by_nuts = function(df) {
  df %>%
    group_by(nutM_code) %>%
    summarize(years = list(sort(unique(year))), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      groups = list({
        y = years
        print(paste("Processing ID:", nutM_code))
        print(paste("All years:", paste(y, collapse = ", ")))
        
        if (length(y) == 1) {
          print("Single year")
          list(y)
        } else {
          breaks = c(0, which(diff(y) >= 10), length(y))
          print(paste("Year gap positions:", paste(breaks, collapse = ", ")))
          lapply(seq_along(breaks[-1]), function(i) y[(breaks[i] + 1):breaks[i + 1]])
        }
      })
    ) %>%
    unnest(groups) %>%
    mutate(
      n_years = lengths(groups),
      year = ifelse(n_years == 1, groups, NA_integer_),
      year_start = ifelse(n_years > 1, map_int(groups, min), NA_integer_),
      year_end = ifelse(n_years > 1, map_int(groups, max), NA_integer_)
    ) %>%
    select(nutM_code, year, year_start, year_end) %>%
    {
      .
    }
}

# apply function to the Aedes presences
ae.presences_agg = years_by_nuts(Aedes_presences)
ae.presences_agg$response = 1 # add response column
# write.csv(ae.presences_agg,"github_data/Occurrence_data/European_brt_training_presences.csv")


# Creating the pseudo-absence shapefile: for pseudo-absence (PA) sampling, remove polygons 
# adjacent to presence sites in order to avoid sampling pseudo-absences in areas too close to known presences.

# Identify adjacent polygons for each presence location
ae.pres_geoms = which(nutsM_sf$LocationCo %in% ae.presences_agg$nutM_code)
ae.pres_geoms = nutsM_sf[ae.pres_geoms, ]
ae.adjacent = st_intersects(ae.pres_geoms$geometry, nutsM_sf$geometry)
ae.adj = unlist(ae.adjacent)
ae.adj_geoms = nutsM_sf[ae.adj, ]


# Remove adjacent polygons from PA sampling shapefile
absences = setdiff(nutsM_IDs$nutM_code, ae.presences_agg) # select all NUTM areas that are not in presences
Aedes_abs = data.frame(
  nutM_code = absences, 
  year = NA
)

remove_rows = which(Aedes_abs$nutM_code %in% ae.adj_geoms$LocationCo) 
Aedes_abs = Aedes_abs[(-remove_rows), ]


# Pseudo-absences must be assigned to a specific year,
# as the model is trained temporally. Sampled years are
# weighted according to the frequency of presence years.

# get presence outbreak dates (year of infection)
Aedes_presences.unique = Aedes_presences[!duplicated(Aedes_presences[c("year", "nutM_code")]), ] # remove duplicated presences
presence_years = Aedes_presences.unique$year

year_counts = table(presence_years)
year_probs = year_counts / sum(year_counts)
years = as.integer(names(year_probs))


## Surveillance capability pseudo-absence sampling:
# Pseudo-absences are sampled from areas with high surveillance suitability 
# and no recorded arbovirus presence.
# Load GeoTIFF surveillance raster Lim et al 2025
rast_surv = raster("Rasters/Lim_et_al_2025_data/outputs/Rasters/Surveillance_map_wmean.tif")
rast_surv_cropped = crop(rast_surv, nutsM) # crop to Europe

# extract surveillance capability values per Nuts area
surv_nutsM = exact_extract(rast_surv_cropped, nutsM, fun="mean")
surv.nutsM_data = nutsM_IDs
surv.nutsM_data$surv = surv_nutsM 
surv.nutsM_data$surv[is.na(surv.nutsM_data$surv)] = 0

# Compute weighted probabilities
Aedes_abs$surv = NA
Aedes_abs$Probabilities = NA
for (k in seq_len(nrow(Aedes_abs))) {
  index = which(surv.nutsM_data$nutM_code == Aedes_abs[k, "nutM_code"])
  Aedes_abs[k, "surv"] = surv.nutsM_data[index, "surv"]
}

Aedes_abs$Probabilities = Aedes_abs$surv / sum(Aedes_abs$surv)


  # 2. Environmental data extraction for BRT training
  # -------------------------------------------------

temperature = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"); temperature2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_2022_2024_monmean.nc")
temperature = addLayer(temperature,temperature2)

precipitation = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"); precipitation2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_2022_2024_monmean.nc")
precipitation = addLayer(precipitation, precipitation2)

relativehumidity = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"); relativehumidity2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_2022_2024_monmean.nc")
relativehumidity= addLayer(relativehumidity,relativehumidity2) 

croplands = brick("Rasters/Environmental_rasters/Land-use/croplands_LUH2-GCB2024_1901_2024.nc")
pastures = brick("Rasters/Environmental_rasters/Land-use/pastures_LUH2-GCB2024_1901_2024.nc")
rangelands = brick("Rasters/Environmental_rasters/Land-use/rangeland_LUH2-GCB2024_1901_2024.nc")
urbanAreas = brick("Rasters/Environmental_rasters/Land-use/urbanAreas_LUH2-GCB2024_1901_2024.nc")
primaryForest = brick("Rasters/Environmental_rasters/Land-use/primaryForest_LUH2-GCB2024_1901_2024.nc")
primaryNonForest = brick("Rasters/Environmental_rasters/Land-use/primaryNonF_LUH2-GCB2024_1901_2024.nc")
secondaryForest = brick("Rasters/Environmental_rasters/Land-use/secondaryForest_LUH2-GCB2024_1901_2024.nc")
secondaryNonForest = brick("Rasters/Environmental_rasters/Land-use/secondaryNonF_LUH2-GCB2024_1901_2024.nc")
secondaryLand = brick("Rasters/Environmental_rasters/Land-use/secondaryLand_LUH2-GCB2024_1901_2024.nc") # merged secondary forested and non-forested areas

population_density = brick("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_1900_2025_historical.nc")
population_density = population_density[[2:125]] # 1901 to 2024

# cell grid areas 
areas = brick("Rasters/Environmental_rasters/clm45_area.nc4")
areas = areas/(1000*1000) # conversion to km2

# convert population density to total population count 
total_population = population_density * areas

# add all variables
covariates = list(temperature, precipitation, relativehumidity, croplands, 
                  pastures, rangelands, urbanAreas, primaryForest, primaryNonForest, 
                  secondaryForest, secondaryNonForest, secondaryLand, total_population)

names(covariates) = c("temperature", "precipitation", "relative_humidity", "croplands", 
                    "pastures", "rangelands", "urbanAreas", "primaryForest", "primaryNonForest", 
                    "secondaryForest", "secondaryNonForest", "secondaryLand", "total_population")

DataExtraction = FALSE
if (DataExtraction == TRUE) {
  # Data extraction #
  plan(multisession, workers = parallel::detectCores() - 1) # set parallele process
  nruns = 100 # this will be used as the number of repetitions for the brt training
  data_train = list()
  
  for(n in 1:nruns) {
    # pseudo-absence sampling
    indices = seq_len(nrow(Aedes_abs))
    selected_absences = sample(indices, 1*length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$Probabilities)
    absences = Aedes_abs[selected_absences, ]

    absences$year = NA_integer_
    absences$year_start = NA_integer_
    absences$year_end = NA_integer_

    # Assign pseudo-years
    absences$year = sample(years, size = length(absences$year), replace = TRUE, prob = year_probs)

    # Add response variable
    absences$response = 0
    absences = absences[, !(colnames(absences) %in% c("surv", "Probabilities"))]
    Aedes_abs = Aedes_abs[, !(colnames(Aedes_abs) %in% c("surv", "Probabilities"))]
    data = rbind(ae.presences_agg, absences) # bind presences and pseudo-absences


    ## data extraction ##
    data_for_brt = as.data.frame(matrix(NA,nrow=nrow(data), ncol = 32))
    data_for_brt[,c(1:5)] = data[,c("nutM_code","year","year_start","year_end", "response")]
    
    colnames(data_for_brt) = c("nutM_code","year","year_start","year_end", "response","annual_temperature","annual_precipitation",
                               "annual_relative_humidity","temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                                "seasonal_variation","precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                                "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                                "croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonForest","secondaryForest","secondaryNonForest", "secondaryLand","total_population","population_density")
    
    
    geometry_matches = nutsM_sf$geometry[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
    area_matches = nutsM_sf$area[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
    nutsM_train = st_sf(data_for_brt, geometry = geometry_matches); nutsM_train$area = area_matches
    names_env = names(data_for_brt[,c(6:31)])
      
    cat("Starting data extraction", n, "of", nruns)
    
    # parallel data extraction  
    results = future_lapply(1:nrow(data_for_brt), function(k) {
      if (!is.na(data_for_brt$year[k])) { 
        # For temporal training, split climate variables into seasonal subsets
        message("temporal training on single year")
        year = as.numeric(data_for_brt[k,"year"])
        year_index = year - 1900
        range1  = ((year - 1901)*12+1)-1; range2 = range1 + 11 # annual indices december of previous year to november of current year
        ranges = list(range1:range2)
        
        month_start = ((year-1) - 1901)*12+1
        winter_indices = sort(c(month_start + 11, month_start + 12, month_start + 13)) # December (previous year), January, February
        spring_indices = sort(c(month_start + 14, month_start + 15, month_start + 16)) # March, April, May
        summer_indices = sort(c(month_start + 17, month_start + 18, month_start + 19)) # June, July, August
        fall_indices   = sort(c(month_start + 20, month_start + 21, month_start + 22)) # September, October, November
        
        } else {
          message("temporal training on average year range")
          year_start = as.numeric(data_for_brt[k,"year_start"]); year_end = as.numeric(data_for_brt[k,"year_end"])
          year_index = (year_start:year_end) - 1900
          
          winter_indices=c()
          spring_indices=c()
          summer_indices=c()
          fall_indices  =c()
          ranges = list()
          for (yr in year_start:year_end) {
            range1  = ((yr - 1901)*12+1)-1; range2 = range1 + 11 # annual indices december of previous year to november of current year
            ranges = c(ranges, list(range1:range2))
            month_start=((yr - 1) - 1901) * 12 + 1
            winter_indices=c(winter_indices, month_start + 11, month_start + 12, month_start + 13) # December (previous year), January, February
            spring_indices=c(spring_indices, month_start + 14, month_start + 15, month_start + 16) # March, April, May
            summer_indices=c(summer_indices, month_start + 17, month_start + 18, month_start + 19) # June, July, August
            fall_indices  =c(fall_indices,   month_start + 20, month_start + 21, month_start + 22) # September, October, November
          }
        }
        
      # Extract environmental variables for the given year or year range
      covariates_temp = list()
      covariates_temp[["temp_annual"]] = mean(stack(lapply(ranges, function(x) mean(covariates[[1]][[x]])))) - 273.15
      covariates_temp[["prec_annual"]] = mean(stack(lapply(ranges, function(x) mean(covariates[[2]][[x]]))))* 60 * 60 * 24 # conversion to kg/m2/day
      covariates_temp[["relh_annual"]] = mean(stack(lapply(ranges, function(x) mean(covariates[[3]][[x]]))))
      
      covariates_temp[["temp_winter"]] = mean(covariates[[1]][[winter_indices]]) - 273.15 # conversion to Celsius
      covariates_temp[["temp_spring"]] = mean(covariates[[1]][[spring_indices]]) - 273.15
      covariates_temp[["temp_summer"]] = mean(covariates[[1]][[summer_indices]]) - 273.15
      covariates_temp[["temp_fall"]]   = mean(covariates[[1]][[fall_indices]]) - 273.15
      
      covariates_temp[["seasonal_var"]] = covariates_temp[["temp_summer"]] - covariates_temp[["temp_winter"]]
      covariates_temp[["prec_winter"]] = mean(covariates[[2]][[winter_indices]]) * 60 * 60 * 24 # conversion to kg/m2/day
      covariates_temp[["prec_spring"]] = mean(covariates[[2]][[spring_indices]]) * 60 * 60 * 24
      covariates_temp[["prec_summer"]] = mean(covariates[[2]][[summer_indices]]) * 60 * 60 * 24
      covariates_temp[["prec_fall"]]   = mean(covariates[[2]][[fall_indices]]) * 60 * 60 * 24
      
      covariates_temp[["relh_winter"]] = mean(covariates[[3]][[winter_indices]])
      covariates_temp[["relh_spring"]] = mean(covariates[[3]][[spring_indices]])
      covariates_temp[["relh_summer"]] = mean(covariates[[3]][[summer_indices]])
      covariates_temp[["relh_fall"]]   = mean(covariates[[3]][[fall_indices]])
      
      covariates_temp[["crops"]] = mean(covariates[[4]][[year_index]])
      covariates_temp[["pastures"]] = mean(covariates[[5]][[year_index]])
      covariates_temp[["rangelands"]] = mean(covariates[[6]][[year_index]])
      covariates_temp[["urbanAreas"]] = mean(covariates[[7]][[year_index]])
      covariates_temp[["primForest"]] = mean(covariates[[8]][[year_index]])
      covariates_temp[["primNonForest"]] = mean(covariates[[9]][[year_index]])
      covariates_temp[["secForest"]] = mean(covariates[[10]][[year_index]])
      covariates_temp[["secNonForest"]] = mean(covariates[[11]][[year_index]])
      covariates_temp[["secLand"]] = mean(covariates[[12]][[year_index]])
      covariates_temp[["totalPop"]] = mean(covariates[[13]][[year_index]])
      
        vals = numeric(length(covariates_temp))
        for (i in seq_along(covariates_temp)) {
          if (i == length(covariates_temp)) {
            vals[i] = exact_extract(covariates_temp[[i]], nutsM_train[k,], fun = "sum") 
          } else {
            vals[i] = exact_extract(covariates_temp[[i]], nutsM_train[k,], fun = "mean")
          }
        }
        return(vals)
      })
      
    data_for_brt[, names_env] = do.call(rbind, results)
    data_for_brt[,32] = data_for_brt$total_population/nutsM_train$area # population density by NUTS3
    data_for_brt$year = as.integer(unlist(data_for_brt$year))
    data_train[[n]] = data_for_brt
    
    write.csv(data_for_brt,paste0("github_data/Brt_training_data/European_BRT_data/20crv3-era5_env_data_replicate",n,".csv"), row.names = F)  
    
    cat("Data extraction run", n, "of", nruns, "completed.\n")
  } 
  
  plan(sequential)  
  saveRDS(data_train, file = "github_data/Brt_training_data/European_BRT_data/brt_training_data_170426.rds")
} else {
  data_train = readRDS("github_data/Brt_training_data/European_BRT_data/brt_training_data_170426.rds")
}

  # 3. Investigating the environmental variables to include in the model (shadow variables approach)
  # ------------------------------------------------------------------------------------------------
nutsM_centroids = coordinates(nutsM); dataframes = list()
all_var_names  = c("temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                      "seasonal_variation", "precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                      "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                      "croplands","pastures","rangelands","urbanAreas", "primaryForest","primaryNonForest","secondaryForest",
                      "secondaryNonForest", "total_population", "population_density" )

# remove annual variables  
for (i in 1:length(data_train)){
  df = data_train[[i]]
  df = df[, c("nutM_code","year","year_start","year_end", "response",all_var_names), drop = FALSE]
  data_train[[i]] = df
}

# get centroids
for (i in 1:length(data_train))
{
  dataframe = data_train[[i]]
  dataframe$x = rep(NA, dim(dataframe)[1]); dataframe$y = rep(NA, dim(dataframe)[1])
  for (j in 1:dim(dataframe)[1])
  {
    centroid = nutsM_centroids[which(nutsM@data[,"LocationCo"]==dataframe[j,"nutM_code"]),]
    dataframe[j,"x"] = centroid[1]; dataframe[j,"y"] = centroid[2]
  }
  dataframes[[i]] = dataframe
}

    # 3.1 Building correlation matrix among predictors (Pearson correlations)

rP_mean = matrix(nrow=length(all_var_names), ncol=length(all_var_names)); colnames(rP_mean) = all_var_names
rP_lowerCI = matrix(nrow=length(all_var_names), ncol=length(all_var_names)); colnames(rP_lowerCI) = all_var_names
rP_upperCI = matrix(nrow=length(all_var_names), ncol=length(all_var_names)); colnames(rP_upperCI) = all_var_names
for (i in 1:length(all_var_names))
{
  for (j in 1:length(all_var_names))
  {
    if (i != j)
    {
      vS = rep(NA, length(dataframes))
      for (k in 1:length(dataframes))
      {
        vS[k] = cor(dataframes[[k]][,all_var_names[i]], dataframes[[k]][,all_var_names[j]], method="pearson")
      }
      rP_mean[i,j] = mean(vS); rP_lowerCI[i,j] = t.test(vS)$conf.int[1]; rP_upperCI[i,j] = t.test(vS)$conf.int[2]
    }
  }
}
for (i in 2:length(all_var_names))
{
  for (j in 1:(i-1))
  {
    to_report = FALSE
    if ((rP_mean[i,j] < 0)&(rP_upperCI[i,j] <= -0.7)) to_report = TRUE
    if ((rP_mean[i,j] > 0)&(rP_lowerCI[i,j] >= 0.7)) to_report = TRUE
    if (to_report)
    {
      cat(all_var_names[i]," - ",all_var_names[j],": rP = ",round(rP_mean[i,j],2),", 95% CI = [",round(rP_lowerCI[i,j],2),", ",round(rP_upperCI[i,j],2),"]\n",sep="")
    }
  }
}
# temperature_spring - temperature_winter: rP = 0.88, 95% CI = [0.88, 0.89]
# temperature_summer - temperature_spring: rP = 0.88, 95% CI = [0.87, 0.88]
# temperature_fall - temperature_winter: rP = 0.91, 95% CI = [0.91, 0.91]
# temperature_fall - temperature_spring: rP = 0.92, 95% CI = [0.92, 0.92]
# temperature_fall - temperature_summer: rP = 0.85, 95% CI = [0.85, 0.86]
# relative_humidity_fall - temperature_summer: rP = -0.72, 95% CI = [-0.73, -0.72]
# secondaryNonForest - secondaryForest: rP = -0.71, 95% CI = [-0.72, -0.71]
# population_density - urbanAreas: rP = 0.97, 95% CI = [0.97, 0.98] 


    # 3.2 Investagating spatial autocorrelation and defining block size for spatial cross-validation

# nutsM polygon centroids
centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  

plottingCorrelogram = FALSE
if (plottingCorrelogram == TRUE)
{
  par(mar=c(2,2,1,1), oma=c(0.5,0.5,0,0), lwd=0.3)
  correlograms = list(); values = rep(NA, 100)
  for (i in 1:100)
  {
    data = training_df[[i]]
    matching_centroids = as.data.frame(matrix(nrow = nrow(data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
    matching_centroids$NUTM = data$nutM_code 
    matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
    
    correlograms[[i]] = ncf::correlog(matching_centroids[,2], matching_centroids[,3], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)
  }
  plot(correlograms[[1]]$mean.of.class, correlograms[[1]]$correlation, ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-1,1))
  for (i in 1:100)
  {
    values[i] = correlograms[[i]]$mean.of.class[which(correlograms[[i]]$correlation <= 0)[1]]
    lines(correlograms[[i]]$mean.of.class, correlograms[[i]]$correlation, lwd=0.1, col="gray60")
  }
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=0.5, tck=-0.05, col.axis="gray30", mgp=c(0,-0.15,0), seq(0,4000,500))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=0.5, tck=-0.05, col.axis="gray30", mgp=c(0,0.18,0))
  title(xlab="Distance (m)", cex.lab=0.65, mgp=c(0.9,0,0), col.lab="gray30")
  title(ylab="Correlation", cex.lab=0.65, mgp=c(1.2,0,0), col.lab="gray30")
  print(mean(unlist(lapply(correlograms, `[[`, "x.intercept")))) # 641.3303, use block size of 650
}  # range size 650 km2

    # 3.3 Comparing the relative importance (RI) of each variable with the RI of its shadow version
brt_analysis = function(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, 
                        bag.fraction, n.trees, theRange, crs, 
                        fold.vector = NULL) { # 1. Added argument here
  
  theRanges = c(theRange, theRange) * 1000
  
  if(is.null(fold.vector)) {
    spdf = SpatialPointsDataFrame(dataframe[c("x","y")], 
                                  dataframe[, c("response", gbm.x)], 
                                  proj4string = crs)
    myblocks = cv_spatial(spdf, column = "response", k = 5, size = theRanges[1], 
                          selection = "random", plot = F, progress = F, report = F)
    fold.vector = myblocks$folds_ids
  }
  
  offset = NULL
  site.weights = rep(1, dim(dataframe)[1])
  var.monotone = rep(0, length(gbm.x))
  n.folds = 5
  prev.stratify = TRUE
  family = "bernoulli"
  step.size = 5
  max.trees = 10000
  tolerance.method = "auto"
  tolerance = 0.001
  plot.main = FALSE
  plot.folds = FALSE
  verbose = FALSE
  silent = TRUE
  keep.fold.models = FALSE
  keep.fold.vector = FALSE
  keep.fold.fit = FALSE
  
  brt_model = gbm.step(dataframe, gbm.x, gbm.y, offset, fold.vector, 
                       tree.complexity, learning.rate, bag.fraction, site.weights,
                       var.monotone, n.folds, prev.stratify, family, n.trees, 
                       step.size, max.trees, tolerance.method, tolerance,
                       plot.main, plot.folds, verbose, silent, 
                       keep.fold.models, keep.fold.vector, keep.fold.fit)
  
  return(brt_model)
}

gbm.y = "response"
gbm.x = all_var_names
tree.complexity = 5
learning.rate = 0.001
bag.fraction = 0.80
n.trees = 10
theRange = 650 # km
n.folds = 5
crs = crs(nutsM)

tab = matrix(0, nrow=length(all_var_names), ncol=1); row.names(tab) = all_var_names
for (h in 1:10)
{
  for (i in 1:length(dataframes))
  {
    # Create shadow variables (randomised versions of predictors)
    # These break any real ecological signal
    shadows = dataframes[[i]][,gbm.x]; colnames(shadows) = paste0(colnames(shadows), "_shadow")
    # shuffle each variable independently
    for (k in 1:dim(shadows)[2]) shadows[,k] = shadows[sample(1:dim(shadows)[1],dim(shadows)[1],replace=F),k]
    # combine real + shadow predictors
    dataframe_shadow = cbind(dataframes[[i]], shadows); gbm.x_shadow = c(gbm.x, colnames(shadows))
    # Fit BRT model including shadow variables
    brt_model = brt_analysis(dataframe_shadow, gbm.y, gbm.x_shadow, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs)
    # Extract relative influence (RI) of each variable
    relativeImportances = matrix(nrow=length(gbm.x_shadow), ncol=1); row.names(relativeImportances) = gbm.x_shadow
    for (j in 1:length(gbm.x_shadow)) relativeImportances[j,1] = summary(brt_model)[gsub("-","\\.",gbm.x_shadow)[j],"rel.inf"]
    # Compare real variables vs their shadow versions
    for (j in 1:length(gbm.x))
    {
      if (relativeImportances[gbm.x[j],1] < relativeImportances[paste0(gbm.x[j],"_shadow"),1]) tab[j,] = tab[j,]+1
    }
  }
}
# Convert counts into proportions (proxy for p-values)
tab = tab/(length(dataframes)*10); colnames(tab) = c("p-value") # to get p-values

# Final variables to include, those with p-value < 0.05 
# (i.e., with RI higher than their shadow version in at least 95% of the runs): 
included_vars = which(tab < 0.05) # variables with p-value < 0.05
included_vars = tab[included_vars,]  
print(included_vars)
# temperature_winter: 0.011         
# temperature_fall:  0.007    
# precipitation_spring: 0.027      
# precipitation_fall:  0.000
# relative_humidity_winter: 0.005 
# relative_humidity_fall: 0.004               
# pastures: 0.005              
# urbanAreas: 0.011         
# primaryNonForest: 0.023


  # 4. BRT model training
  # ----------------------

# final inlcuded variables
training_df = list()
for (i in 1:length(dataframes)) {
  df = dataframes[[i]]
  shadow_005 = df %>% dplyr::select(-c(temperature_spring, temperature_summer, precipitation_winter,
                                       precipitation_summer, relative_humidity_spring, relative_humidity_summer,seasonal_variation,rangelands,croplands,
                                       primaryForest, secondaryForest, secondaryNonForest, total_population, population_density))
  training_df[[i]] = shadow_005
}

var_names = names(training_df[[1]][6:14])

    # 4.1 hyperparameter tuning on selected set of variables
n_cores = parallel::detectCores() - 3
cl = makeCluster(n_cores)
registerDoParallel(cl)

# tuning Grid
learning_rates    = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
tree_complexities = c(2, 3, 5, 6, 7, 10)

# fixed parameters
gbm.y = "response"
gbm.x = var_names
bag.fraction = 0.75
theRange = 650 
n.trees = 10
crs_obj = crs(nutsM) 

results_list = foreach(j = 1:100, .combine = rbind, 
                       .packages = c("dismo", "gbm", "sp", "blockCV"),
                       .export = c("training_df", "learning_rates", "tree_complexities",
                                   "bag.fraction", "theRange", "gbm.y", "gbm.x", 
                                   "crs_obj", "brt_analysis")) %dopar% {
                                     
                                     current_data = as.data.frame(training_df[[j]])
                                     
                                     
                                     # Generate folds once here so all tuning combos use the SAME spatial split
                                     spdf = SpatialPointsDataFrame(current_data[c("x","y")], 
                                                                   current_data[, c(gbm.y, gbm.x)], 
                                                                   proj4string = crs_obj)
                                     
                                     myblocks = cv_spatial(spdf, column = gbm.y, k = 5, size = theRange * 1000, 
                                                           selection = "random", plot = F, progress = F, report = F)
                                     current_folds = myblocks$folds_ids
                                     
                                     local_results = data.frame()
                                     
                                     # Loop through Tuning Grid
                                     for (lr in learning_rates) {
                                       for (tc in tree_complexities) {
                                         
                                         # Use tryCatch correctly around the function call
                                         model <- tryCatch({
                                           brt_analysis(
                                             dataframe = current_data, 
                                             gbm.y = gbm.y, 
                                             gbm.x = gbm.x, 
                                             tree.complexity = tc, # Use tc directly
                                             learning.rate = lr,   # Use lr directly
                                             bag.fraction = bag.fraction, 
                                             n.trees = n.trees, 
                                             theRange = theRange, 
                                             crs = crs_obj
                                           )
                                         }, error = function(e) NULL)
                                         
                                         if (!is.null(model)) {
                                           local_results = rbind(local_results, data.frame(
                                             replicate       = j,
                                             learning_rate   = lr,
                                             tree_complexity = tc,
                                             bag_fraction    = bag.fraction,
                                             n_trees         = model$gbm.call$best.trees,
                                             cv_deviance     = model$cv.statistics$deviance.mean,
                                             AUC             = model$cv.statistics$discrimination.mean
                                           ))
                                         }
                                       }
                                     }
                                     local_results
                                   }

stopCluster(cl)


# Group results, and mean over replicates for each combination of learning rate and tree complexity
summary_results = results_list %>%
  group_by(learning_rate, tree_complexity) %>%
  summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    mean_dev = mean(cv_deviance, na.rm = TRUE),
    sd_dev   = sd(cv_deviance, na.rm = TRUE),
    ci_low   = t.test(cv_deviance)$conf.int[1],
    ci_high  =  t.test(cv_deviance)$conf.int[2],
    mean_nt  = mean(n_trees, na.rm = TRUE),
    .groups = "drop"
  )

# find best parameter combination based on mean deviance, while ensuring that the nt is > 1000 and 
# that tc is not too low to allow for complex interactions.

best = summary_results %>% arrange(mean_dev)
best = best %>% filter(mean_nt > 1000)
head(best, 10) # we chose lr = 0.001, tc = 3 and bag fraction = 0.75.



    # 4.2 Final BRT training, with selected variables and hyperparameters, and spatial cross-validation (SCV)

# fixed parameters
gbm.y = "response"
gbm.x = var_names
learning.rate =  0.001
tree.complexity = 3
bag.fraction = 0.75
theRange = 650 
n.trees_start = 10
crs = crs(nutsM) 
nruns = 100

brt_model_scvs = list() # brt with spatial cross-validations (SCVs)
AUCs = matrix(nrow = nruns, ncol = 1); colnames(AUCs) = "AUCs"
SIppcs = matrix(nrow = nruns, ncol = 1); colnames(SIppcs) = "SIpcc"
thresholds = matrix(nrow = nruns, ncol = 1); colnames(thresholds) = "thresholds"

for(i in 1:length(training_df)){
  current_data = training_df[[i]]
  brt_model_scvs[[i]] = brt_analysis(current_data, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees_start, theRange, crs)
  
  # Model predictive performance evaluation
  # AUC 
  AUCs[i,] = brt_model_scvs[[i]]$cv.statistics$discrimination.mean # Mean test AUC
  
  # Calculation of the sorensen index to establish predictive performance of our models
  # Sources:
  # - computation performed according to the formulas of Leroi et al. (2018, J. Biogeography)
  # - optimisation of the threshold with a 0.01 step increment according to Li & Guo (2013, Ecography)

  tmp = matrix(nrow = 101, ncol = 2); tmp[,1] = seq(0,1,0.01)
  df = brt_model_scvs[[i]]$gbm.call$dataframe
  responses = df$response;
  data = df[,6:dim(df)[2]]
  n.trees = brt_model_scvs[[i]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
  prediction = predict.gbm(brt_model_scvs[[i]], data, n.trees, type, single.tree)
  N = dim(data)[1]; P = sum(responses==1); A = sum(responses==0)
  prev = P/(P+A) # proportion of recorded sites where the species is present
  x = (P/A)*((1-prev)/prev);
  sorensen_ppc = 0
  
  for (threshold in seq(0,1,0.01)) {
    TP = length(which((responses==1)&(prediction>=threshold))) # true positives
    FN = length(which((responses==1)&(prediction<threshold))) # false negatives
    FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
    sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
    
    if (sorensen_ppc < sorensen_ppc_tmp) {
      sorensen_ppc = sorensen_ppc_tmp
      optimised_threshold = threshold
      }
    }
  
  SIppcs[i,] = sorensen_ppc
  thresholds[i,] = optimised_threshold
}

saveRDS(brt_model_scvs, file = "github_data/Brt_training_data/European_BRT_data/brt_trained_model.rds")

  # 4.3 Predictive performance table
  # --------------------------------

metrics = cbind(AUCs, SIppcs, thresholds)
colnames(metrics) = c("Area Under the Curve (AUC)", 
                      "Prevalence-pseudo-absence calibrated Sorensen index (SIppc)",
                      "Threshold value maximising the SIppc")

metrics[1:100, 1:2] = round(metrics[1:100, 1:2], 3); metrics = as.data.frame(metrics)
for (i in 1:3) {
  x = metrics[1:100, i]
  mean_val = round(mean(x, na.rm = TRUE), 3)
  se = sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  ci = round(1.96 * se, 3)
  CI_lower = round(mean_val - ci, 3)
  CI_upper = round(mean_val + ci, 3)
  
  metrics[101, i] = paste0(mean_val, " (", CI_lower, "-", CI_upper, ")")
}
# Mean AUC: 0.886 (0.878-0.894); Mean SIppc: 0.978 (0.975-0.981); Mean SIppc threshold: 0.48 (0.466-0.494)
		
write.csv(metrics, "github_data/Brt_training_data/European_BRT_data/brt_model_performances.csv", row.names = FALSE)


# III. VARIABLE RELATIVE IMPORTANCES AND RESPONSE CURVES ------------------------

  # 1. Relative importances 
  # -----------------------

brt_model_scvs = readRDS("github_data/Brt_training_data/European_BRT_data/brt_trained_model.rds")
brt_model_scv = brt_model_scvs

var_names = brt_model_scvs[[1]]$var.names

nruns = 100
relativeInfluences = matrix(0, nrow=length(var_names), ncol=3) #nrow depends on the number of covariates
relativeInfluences_all = matrix(0, nrow=length(var_names), ncol=nruns) #nrow depends on the number of covariates


# get relative influence of each covariate
for (j in 1:length(brt_model_scv)) {
  for (k in 1:length(var_names)) {
      relativeInfluences[k,1] = relativeInfluences[k] + summary(brt_model_scv[[j]])[var_names[k],"rel.inf"]
      relativeInfluences_all[k,j] = summary(brt_model_scv[[j]], plot = FALSE)[var_names[k],"rel.inf"]
      
    }
  }
  
row.names(relativeInfluences) = var_names
relativeInfluences[,1] = relativeInfluences[,1]/length(brt_model_scv)
relativeInfluences[,2] = as.data.frame(t(apply(as.matrix(relativeInfluences_all), 1, function(x) ci(x))))$`CI lower`
relativeInfluences[,3] = as.data.frame(t(apply(as.matrix(relativeInfluences_all), 1, function(x) ci(x))))$`CI upper`
write.csv(relativeInfluences, file = "github_data/Brt_training_data/European_BRT_data/Relative_influences.csv", row.names = TRUE)
write.csv(relativeInfluences_all, file = "github_data/Brt_training_data/European_BRT_data/All_relative_influences.csv", row.names = TRUE)

  # 2. Response curves  
  # -------------------
relativeInfluences = as.data.frame(relativeInfluences)
colnames(relativeInfluences) = c("mean", "CI_lower", "CI_upper")
relativeInfluences = arrange(relativeInfluences, desc(mean))
var_names = row.names(relativeInfluences)

data_curv = brt_model_scv[[1]]$gbm.call$dataframe
data_curv = data_curv[data_curv[, "response"] == 1, ] 
envVariableValues = matrix(nrow=3, ncol=length(var_names))
row.names(envVariableValues) = c("median","minV","maxV")
colnames(envVariableValues) = var_names
  
  for (j in 1:length(var_names))
  {
    minV = min(data_curv[,var_names[j]], na.rm=T)
    maxV = max(data_curv[,var_names[j]], na.rm=T)
    medianV = median(data_curv[,var_names[j]], na.rm=T)
    envVariableValues[,j] = rbind(medianV, minV, maxV)
  }
  
envVariableValues_list = envVariableValues
                           

# pretty names
var_names1 = c("Fall precipitation","Humidity fall", "Fall temperature", "Humidity winter", "Urban areas",
                      "Pastures", "Winter temperature", "Primary non-forested areas", "Spring precipitation")

# plotting response curve code
pdf("github_data/Brt_training_data/European_BRT_data/response_curves2.pdf", width = 4.5, height = 3.5)
par(mfrow=c(3,3), oma=c(0.5,1,0.5,1), mar=c(2,1.38,0.2,0.5), lwd=0.2, col="gray30", bg="white")

# Colors
col_fullrange  = adjustcolor("#e76e55", alpha.f = 0.20)  # light band: min-max
col_median     = "#de4427"                         # bold median line

for (i in 1:length(var_names))
{
  valuesInterval = (envVariableValues_list["maxV",i] - envVariableValues_list["minV",i]) / 100
  df = data.frame(matrix(nrow = length(seq(envVariableValues_list["minV",i],
                                           envVariableValues_list["maxV",i],
                                           valuesInterval)),
                         ncol = length(var_names)))
  colnames(df) = var_names
  
  for (k in 1:length(var_names))
  {
    valuesInterval_k = (envVariableValues_list["maxV",k] - envVariableValues_list["minV",k]) / 100
    if (i == k) df[,var_names[k]] = seq(envVariableValues_list["minV",k],
                                               envVariableValues_list["maxV",k],
                                               valuesInterval_k)
    if (i != k) df[,var_names[k]] = rep(envVariableValues_list["median",k], nrow(df))
  }
  
  pred_matrix = matrix(NA, nrow = nrow(df), ncol = length(brt_model_scv))
  
  for (j in 1:length(brt_model_scv))
  {
    n.trees = brt_model_scv[[j]]$gbm.call$best.trees
    pred_matrix[,j] = predict.gbm(brt_model_scv[[j]], newdata=df,
                                  n.trees, "response", FALSE)
  }
  
  pred_median = apply(pred_matrix, 1, median)
  pred_min    = apply(pred_matrix, 1, min)
  pred_max    = apply(pred_matrix, 1, max)
  
  x_vals = df[, var_names[i]]
  
  # Base plot (empty, sets axes)
  plot(x_vals, pred_median, type="n", ann=F, axes=F,
       xlim=range(x_vals), ylim=c(0,1))
  
  # full min-max range
  polygon(c(x_vals, rev(x_vals)), c(pred_min, rev(pred_max)),
          col=col_fullrange, border=NA)
  
  
  # Median line
  lines(x_vals, pred_median, col=col_median, lwd=0.7)
  
  # Axes
  ticks <- axTicks(1)
  axis(side=1, at = ticks,lwd.tick=0.2, cex.axis=0.4, lwd=0, tck=-0.06,labels = round(ticks, 2),
       col.axis="gray30", mgp=c(0,0.15,0))
  axis(side=2, lwd.tick=0.2, cex.axis=0.4, lwd=0, tck=-0.060,
       col.axis="gray30", mgp=c(0,0.15,0))
  
  # Labels
  env_label = var_names1[i]
  ri_label  = paste0("RI = ", round(relativeInfluences[i,1],1), "% [",
                     round(relativeInfluences[i,2],2), "-",
                     round(relativeInfluences[i,3],2), "]")
  
  n_col = 3; panel_id = i
  if ((panel_id-1) %% n_col == 0) {
    title(ylab="Predicted value", cex.lab=0.58, mgp=c(0.8,0,0), col.lab="gray30")
  }
  title(xlab=env_label, cex.lab=0.58, mgp=c(0.8,0,0), col.lab="gray30")
  #title(xlab=ri_label,  cex.lab=0.55, mgp=c(1.3,0,0), col.lab="gray50")
  text(x = mean(range(x_vals)), y = 0.25, 
       labels = ri_label, 
       cex = 0.55, 
       col = "gray30")
  
  
  box(lwd=0.2, col="gray30")
}

dev.off()


  
# IV. CURRENT RISK MAP (2020-2024) ---------------------------------------------

  # 1.data extraction
  # ------------------

years = 2020:2024
all_winter_indices <- c()
all_spring_indices <- c()
all_summer_indices <- c()
all_fall_indices   <- c()

for (yr in years) {
  month_start = (yr - 1901) * 12 + 1
  
  # 2. Calculate current year indices
  current_winter = c(month_start - 1, month_start, month_start + 1) # Dec(prev), Jan, Feb
  current_spring = c(month_start + 2, month_start + 3, month_start + 4)
  current_summer = c(month_start + 5, month_start + 6, month_start + 7)
  current_fall   = c(month_start + 8, month_start + 9, month_start + 10)
  
  all_winter_indices <- c(all_winter_indices, current_winter)
  all_spring_indices <- c(all_spring_indices, current_spring)
  all_summer_indices <- c(all_summer_indices, current_summer)
  all_fall_indices   <- c(all_fall_indices,   current_fall)
}

year_index = 120:124 # indices for land use and population variables (annual data)

# Extract environmental variables
covariates_2024 = list()
covariates_2024[["temp_winter"]] = mean(covariates[["temperature"]][[all_winter_indices]]) - 273.15 # conversion to Celsius
covariates_2024[["temp_spring"]] = mean(covariates[["temperature"]][[all_spring_indices]]) - 273.15
covariates_2024[["temp_summer"]] = mean(covariates[["temperature"]][[all_summer_indices]]) - 273.15
covariates_2024[["temp_fall"]]   = mean(covariates[["temperature"]][[all_fall_indices]]) - 273.15

covariates_2024[["seasonal_var"]] = covariates_2024[["temp_summer"]] - covariates_2024[["temp_winter"]]

covariates_2024[["prec_winter"]] = mean(covariates[["precipitation"]][[all_winter_indices]]) * 60 * 60 * 24 # conversion to kg/m2/day
covariates_2024[["prec_spring"]] = mean(covariates[["precipitation"]][[all_spring_indices]]) * 60 * 60 * 24
covariates_2024[["prec_summer"]] = mean(covariates[["precipitation"]][[all_summer_indices]]) * 60 * 60 * 24
covariates_2024[["prec_fall"]]   = mean(covariates[["precipitation"]][[all_fall_indices]]) * 60 * 60 * 24

covariates_2024[["relh_winter"]] = mean(covariates[["relative_humidity"]][[all_winter_indices]])
covariates_2024[["relh_spring"]] = mean(covariates[["relative_humidity"]][[all_spring_indices]])
covariates_2024[["relh_summer"]] = mean(covariates[["relative_humidity"]][[all_summer_indices]])
covariates_2024[["relh_fall"]]   = mean(covariates[["relative_humidity"]][[all_fall_indices]])

covariates_2024[["crops"]] = mean(covariates[["croplands"]][[year_index]])
covariates_2024[["pastures"]] = mean(covariates[["pastures"]][[year_index]])
covariates_2024[["rangelands"]] = mean(covariates[["rangelands"]][[year_index]])
covariates_2024[["urbanAreas"]] = mean(covariates[["urbanAreas"]][[year_index]])
covariates_2024[["primForest"]] = mean(covariates[["primaryForest"]][[year_index]])
covariates_2024[["primNonForest"]] = mean(covariates[["primaryNonForest"]][[year_index]])
covariates_2024[["secForest"]] = mean(covariates[["secondaryForest"]][[year_index]])
covariates_2024[["secNonForest"]] = mean(covariates[["secondaryNonForest"]][[year_index]])
covariates_2024[["secLand"]] = mean(covariates[["secondaryLand"]][[year_index]])
covariates_2024[["totalPop"]] = mean(covariates[["total_population"]][[year_index]])

data2024 = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(all_var_names) + 1))
colnames(data2024) = c("NUTM", all_var_names); data2024$NUTM = nutsM_IDs$nutM_code

for (i in seq_along(covariates_2024)) {
  if (i == length(covariates_2024)) {
    extracted_values = exactextractr::exact_extract(covariates_2024[[i]], nutsM, fun = "sum")
  } else {
    extracted_values = exact_extract(covariates_2024[[i]], nutsM, fun = "mean")
    }
  data2024[, i + 1] = extracted_values  # + 1 because first column is NUTM
  }
data2024$population_density = data2024$total_population/nutsM_sf$area
included_vars_2024 = data2024[,c("NUTM",var_names)]
write.csv(data2024,"github_data/Brt_projections_outputs/Current_projections/20crv3-era5_europe_data_2020_2024_130526.csv", row.names = F)
write.csv(included_vars_2024,"github_data/Brt_projections_outputs/Current_projections/20crv3-era5_europe_included_vars_2020_2024_130526.csv", row.names = F)

  # 2. current predictions (t0)
  # ------------------------

eu_2024 = as.data.frame(matrix(nrow=nrow(present_df_kyla), ncol=100))
for(j in 1:length(brt_model_scvs)){
  object = brt_model_scvs[[j]]
  n.trees = brt_model_scvs[[j]]$gbm.call$best.trees 
  type = "response"; single.tree = FALSE
  prediction = predict.gbm(object, present_df_kyla, n.trees, type, single.tree)
  eu_2024[,j] = prediction
}
colnames(eu_2024)[1:100] = paste0("rep_", 1:100)
write.csv(eu_2024,"github_data/Brt_projections_outputs/Current_projections/current_predictions_130526.csv", row.names = F)


# V. PAST PREDICTIONS ----------------------------------------------------------

  # 1. Data extraction 
  # ------------------
# Past time intervals
years_start = c(1901,1925,1950,1975,2000)
years_end = c(1924,1949,1974,1999,2021)

# start past data extraction and predictions
models = c("obsclim","counterclim")
var_names = brt_model_scvs[[1]]$var.names

options(future.globals.maxSize = 8 * 1024^3)
plan(multisession, workers = parallel::detectCores() - 1) # parallel across models

past_predictions_list = future_lapply(seq_along(models), function(m) {
  past_predictions = list()

  temperature = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_tas_1901_2021","_monmean.nc"))
  precipitation = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_pr_1901_2021","_monmean.nc"))
  relativehumidity = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_hurs_1901_2021","_monmean.nc"))
  
  for(y in 1:length((years_start))) {
    past_predictions[[y]] = list()
    start = years_start[[y]] ; end = years_end[[y]] 
    year_index = (start-1900):(end-1900)
    
    if (y==1) {month_start = (start - 1901)*12+1; month_end = (end - 1901)*12 + 12
    }else{month_start = ((start-1) - 1901)*12+1; month_end = (end - 1901)*12 + 12}
    
    temperature_temp = temperature[[month_start:month_end]]
    precipitation_temp = precipitation[[month_start:month_end]]
    relativehumidity_temp = relativehumidity[[month_start:month_end]]
    
    # remove 1901 december
    if (y==1) {
      winter_indices = sort(c(12, 1, 2)) + rep((start:end - start) * 12, each = 3); winter_indices = winter_indices[-c(1,2,length(winter_indices))]
      spring_indices = sort(c(3,4,5) + rep((start:end - start)*12, each= 3))
      summer_indices = sort(c(6,7,8) + rep((start:end - start)*12, each= 3))
      fall_indices = sort(c(9,10,11) + rep((start:end - start)*12, each= 3))
    }else{
      winter_indices = sort(c(12, 1, 2)) + rep((start:end - start) * 12, each = 3); winter_indices = winter_indices[-c(1,2,length(winter_indices))]
      spring_indices = sort(c(3,4,5) + rep((start:end - start)*12, each= 3)); spring_indices = spring_indices[-c(1,2,3)]
      summer_indices = sort(c(6,7,8) + rep((start:end - start)*12, each= 3)); summer_indices = summer_indices[-c(1,2,3)]
      fall_indices = sort(c(9,10,11) + rep((start:end - start)*12, each= 3)); fall_indices = fall_indices[-c(1,2,3)]
    }
    
    
    tempm_winter =  mean(temperature_temp[[winter_indices]]) - 273.15
    tempm_spring = mean(temperature_temp[[spring_indices]]) - 273.15
    tempm_summer = mean(temperature_temp[[summer_indices]]) - 273.15
    tempm_fall = mean(temperature_temp[[fall_indices]]) - 273.15
    
    precp_winter = mean(precipitation_temp[[winter_indices]]) * 60 * 60 * 24
    precp_spring = mean(precipitation_temp[[spring_indices]]) * 60 * 60 * 24
    precp_summer = mean(precipitation_temp[[summer_indices]]) * 60 * 60 * 24
    precp_fall = mean(precipitation_temp[[fall_indices]]) * 60 * 60 * 24
    
    relh_winter = mean(relativehumidity_temp[[winter_indices]])
    relh_spring = mean(relativehumidity_temp[[spring_indices]]) 
    relh_summer = mean(relativehumidity_temp[[summer_indices]])
    relh_fall = mean(relativehumidity_temp[[fall_indices]])
    
    envVariables = list()
    envVariables[[1]] = tempm_winter
    envVariables[[2]] = tempm_spring
    envVariables[[3]] = tempm_summer
    envVariables[[4]] = tempm_fall
    
    envVariables[[5]] = tempm_summer - tempm_winter # seasonal variation
    
    envVariables[[6]] = precp_winter
    envVariables[[7]] = precp_spring
    envVariables[[8]] = precp_summer
    envVariables[[9]] = precp_fall
    
    envVariables[[10]] =  relh_winter
    envVariables[[11]] =  relh_spring
    envVariables[[12]] =  relh_summer
    envVariables[[13]] =  relh_fall
    
    envVariables[[14]] =  mean(covariates[[4]][[year_index]]) # croplands
    envVariables[[15]] =  mean(covariates[[5]][[year_index]]) # pastures
    envVariables[[16]] =  mean(covariates[[6]][[year_index]]) # rangelands
    envVariables[[17]] =  mean(covariates[[7]][[year_index]]) # urban areas
    envVariables[[18]] =  mean(covariates[[8]][[year_index]]) # primary forest
    envVariables[[19]] =  mean(covariates[[9]][[year_index]]) # primary non forest
    envVariables[[20]] =  mean(covariates[[10]][[year_index]]) # secondary forest
    envVariables[[21]] =  mean(covariates[[11]][[year_index]]) # secondary non forest
    envVariables[[22]] =  mean(covariates[[12]][[year_index]]) # secondary land
    envVariables[[23]] =  mean(covariates[[13]][[year_index]]) # total population
    
    
    nutsM_data_past = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(names_env) + 1))
    colnames(nutsM_data_past) = c("NUTM", names_env)
    nutsM_data_past$NUTM = nutsM_IDs$nutM_code
    
    for (i in seq_along(envVariables)) { 
      if (i == length(envVariables)) {
        extracted_values = log(exactextractr::exact_extract(envVariables[[i]], nutsM, fun = "sum") + 1)
      } else {
        extracted_values = exact_extract(envVariables[[i]], nutsM, fun = "mean")
      }
      nutsM_data_past[, i + 1] = extracted_values  # + 1 because first column is NUTM
    }
    nutsM_data_past$population_density = nutsM_data_past$total_population/nutsM_sf$area
    
    included_data_past = nutsM_data_past[,c("NUTM",var_names)]
    write.csv(nutsM_data_past, paste0("github_data/Brt_projections_outputs/Past_projections/extracted_env_data/20crv3-era5_",models[[m]],"_",years_start[[y]],"-",years_end[[y]],".csv"), row.names = F)
    write.csv(included_data_past, paste0("github_data/Brt_projections_outputs/Past_projections/extracted_env_data/20crv3-era5_",models[[m]],"_",years_start[[y]],"-",years_end[[y]],"_included_vars.csv"), row.names = F)
    
    
      # 2. Historical predictions
      # --------------------------
  
     for(j in 1:length(brt_model_scvs)) {
       object = brt_model_scvs[[j]]
       n.trees = brt_model_scvs[[j]]$gbm.call$best.trees;
       type = "response"; single.tree = FALSE
       prediction = predict.gbm(object, included_data_past[, -1], n.trees, type, single.tree)
       past_predictions[[y]][[j]] = prediction
     }
  }
  saveRDS(past_predictions, paste0("github_data/Brt_projections_outputs/Past_projections/Past_predictions_", models[m], ".rds")) 
  return(past_predictions)
})


  # 3. Saving predictions
  # --------------------------
year_intervals = c("1901-1924","1925-1949","1950-1974","1975-1999","2000-2021")
path = "github_data/Brt_projections_outputs/Past_projections/"

for (m in 1:length(models)) {
  for (y in 1:length(year_intervals)) {
    df = as.data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = 100))
    for (j in 1:100) {
      df[,j] = past_predictions_list[[m]][[y]][[j]]
    }
    df = cbind(NUTSM = nutsM_IDs$nutM_code, df); colnames(df)[2:101] = paste0("rep_",1:100)
    write.csv(df, paste0(path,"Past_preds_",models[[m]],"_",year_intervals[[y]],"_20CRv3-ERA5_300426.csv"), row.names = F)
  }
}

  # 4. Historical predictions figures
  # ----------------------------------

coasts = st_union(nutsM_sf) # europe outline for plotting
thresholds_SI = read.csv("github_data/Brt_training_data/European_BRT_data/brt_model_performances.csv", head=T) # read SIppc thresholds 
colnames(thresholds_SI) = c("AUC", "SIppc", "threshold")
thresholds_SI = thresholds_SI[-101,]

# load predictions for all scenarios and periods 
brts = list(); scenarios = c("obsclim","counterclim"); periods = c("1901-1924","1925-1949","1950-1974","1975-1999","2000-2021")
for (i in 1:length(scenarios))
{
  buffer = list()
  for (j in 1:length(periods))
  {
    buffer[[j]] = read.csv(paste0("github_data/Brt_projections_outputs/Past_projections/Past_preds_",scenarios[i],"_",periods[j],"_20CRv3-ERA5_300426.csv"), head=F)
    buffer[[j]] = buffer[[j]][,2:dim(buffer[[j]])[2]]
  }
  brts[[i]] = buffer
}

# plotting parameters for the different maps (A to E)
suffix = "A"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE # baseline map
suffix = "B"; plottingStDevs = TRUE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE # standard deviations map
suffix = "C"; plottingStDevs = FALSE; plottingLower95CI = TRUE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE # lower 95% confidence interval map
suffix = "D"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = TRUE; plottingBinaryMap = FALSE # upper 95% confidence interval map
suffix = "E"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = TRUE # binary map (proportion of predictions above the SIppc threshold)

# plotting
pdf(paste0("github_data/Figures/Figure_1",suffix,"_010526.pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
scenario_names = c("Historical","Counterfactual")
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111]) # main figure color scale
# colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
# colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # Set 85% opacity ("D9") 
# colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") 
# colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5" 
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18]) # difference map color scale
colourScale2[9] = paste0(colourScale2[9],"00") # "00" = 0% midpoint (index 76), "white" for values near zero
colourScale3 = colourScale2; colourScale3[9] = "#E5E5E5" # Same as Scale 2, but instead of transparency, use a solid light gray to explicitly mark the "no change" or zero-value areas.
colourScale4 = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[11:111] # standard deviation color palette
colourScale = colourScale1
if (plottingStDevs) { colourScale = colourScale4 }
for (i in 1:length(scenarios))
{
  for (j in 1:length(periods))
  {
    brt = as.matrix(brts[[i]][[j]]);colnames(brt) = brt[1, ]; brt = brt[-1, ]; brt = apply(brt, 2, as.numeric)  
    means = rep(NA, dim(brt)[1]); stdvs = rep(NA, dim(brt)[1])
    lowerCI95 = rep(NA, dim(brt)[1]); upperCI95 = rep(NA, dim(brt)[1])
    binaries = rep(0, dim(brt)[1])
    for (k in 1:dim(brt)[1])
    {
      means[k] = mean(brt[k,]); stdvs[k] = sd(brt[k,])
      lowerCI95[k] = t.test(brt[k,])$conf.int[1]
      upperCI95[k] = t.test(brt[k,])$conf.int[2]
      temp = rep(0, dim(brt)[2])
      temp[which(brt[k,]>=thresholds_SI[,"threshold"])] = 1
      binaries[k] = round(mean(temp))
    }
    cols = colourScale[(((means-0)/(1-0))*100)+1]; # print(max(stdvs))
    if (plottingStDevs) cols = colourScale[(((stdvs-0)/(0.25-0))*100)+1]
    if (plottingLower95CI) cols = colourScale[(((lowerCI95-0)/(1-0))*100)+1]
    if (plottingUpper95CI) cols = colourScale[(((upperCI95-0)/(1-0))*100)+1]
    if (plottingBinaryMap) cols = c("gray90","#de4327")[binaries+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA)
    plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
    # if ((i == 1)&(j == 1)) mtext(expression(bold(A)), side=3, line=-1.45, at=-8.5, cex=0.70, col="gray30")
    # if ((i == 2)&(j == 1)) mtext(expression(bold(B)), side=3, line=-1.45, at=-8.5, cex=0.70, col="gray30")
    mtext(scenario_names[i], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
    if (plottingStDevs) rast = raster(as.matrix(c(0,0.25)))
    if ((i == length(scenarios))&(j == length(periods))&(!plottingBinaryMap))
    {
      plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
           axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
    }								
  }
}

# difference maps
vS = c()
for (j in 1:length(periods))
{
  brt1 = as.matrix(brts[[1]][[j]]); colnames(brt1) = brt1[1, ]; brt1 = brt1[-1, ]; brt1 = apply(brt1, 2, as.numeric)  

  brt2 = as.matrix(brts[[2]][[j]]); colnames(brt2) = brt2[1, ]; brt2 = brt2[-1, ]; brt2 = apply(brt2, 2, as.numeric)
  for (k in 1:dim(brt1)[1])
  {
    vS = c(vS, mean(brt1[k,]-brt2[k,]))
  }
}
vS = vS[which(!is.na(vS))]; minVS = min(vS); maxVS = max(vS)
if (abs(minVS) < abs(maxVS)) minVS = -maxVS
if (abs(maxVS) < abs(minVS)) maxVS = -minVS
maxVS = 0.75; minVS = -0.75
for (j in 1:length(periods))
{
  brt1 = as.matrix(brts[[1]][[j]]); colnames(brt1) = brt1[1, ]; brt1 = brt1[-1, ]; brt1 = apply(brt1, 2, as.numeric)
  brt2 = as.matrix(brts[[2]][[j]]); colnames(brt2) = brt2[1, ]; brt2 = brt2[-1, ]; brt2 = apply(brt2, 2, as.numeric)
  means1 = rep(NA, dim(brt1)[1]); means2 = rep(NA, dim(brt2)[1]); differences = rep(NA, dim(brt1)[1])
  for (k in 1:dim(brt1)[1])
  {
    means1[k] = mean(brt1[k,]); means2[k] = mean(brt2[k,]); differences[k] = mean(brt1[k,]-brt2[k,])
  }
  cols = colourScale2[(((differences-minVS)/(maxVS-minVS))*length(colourScale2))+1]		
  plot(coasts, lwd=0.4, border="gray30", col=NA)
  plot(coasts, add=T, border=NA, col="#E5E5E5")
  plot(nutsM, col=cols, border=NA, add=T)
  mtext("Difference", side=3, line=-1.7, at=1, cex=0.50, col="gray30")
  mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
  if ((i == length(scenarios))&(j == length(periods)))
  {
    rast = raster(as.matrix(c(minVS,maxVS))); par(lwd=0.1, fg="white")
    plot(rast, legend.only=T, add=T, col="#E5E5E5", legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96),
         legend.args=list(text="", cex=0.7, line=0.0, col=NA, col.axis=NA), col.tick=NA, adj=3, alpha=1, side=3,
         axis.args=list(cex.axis=0.65, lwd=0, col=NA, lwd.tick=0, col.tick=NA, tck=-1.2, col.axis=NA, line=0, mgp=c(0,0.45,0)))
    par(lwd=0.2, fg="gray30")
    plot(rast, legend.only=T, add=T, col=colourScale3, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
         axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, 
                        mgp=c(0,0.45,0), at=c(-0.6,-0.3,0,0.3,0.6)), alpha=1, side=3)
  }								
}
dev.off()


# VI. FUTURE PREDICTIONS -------------------------------------------------------

  # 1. subset and aggregate future landuse
  # --------------------------------------
scenarios = c("ssp126","ssp370","ssp585")
scenarios2 = c("SSP1","SSP3","SSP5")

# subset land-use data per land cover variable
LAND_USE = list()
years = c(2015:2100)
for (i in 1:length(scenarios)){
  land_cover = nc_open(paste0("Rasters/Environmental_rasters/Land-use/Future_landuse_annual/",scenarios[i],"/landcover_",scenarios[i],"_2015-2100.nc"))
  landCoverVariableIDs = names(land_cover$var); land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
  landCoverVariableNames = as.character(read.csv("Rasters/Environmental_rasters/Land-use/LC_vars.csv")[1:12,2])
  for (j in 1:12)
  {
    land_covers1[[j]] = brick(paste0("Rasters/Environmental_rasters/Land-use/Future_landuse_annual/",scenarios[i],"/landcover_",scenarios[i],"_2015-2100.nc"), varname=landCoverVariableIDs[j])
  }
  variable_codes = c("croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
  variable_names = c("crops","pasture","rangeland","urban land","forested primary land","non-forested primary land",
                     "potentially forested secondary land","potentially non-forested secondary land")
  for (j in 1:length(variable_names))
  {
    names = gsub("\\."," ",landCoverVariableNames)
    indices = which(landCoverVariableNames==variable_names[j])
    if (length(indices) == 0) indices = which(grepl(variable_names[j],names))
    land_cover = land_covers1[[indices[1]]]
    names(land_cover) = paste0(variable_codes[j], "_", years)
    if (length(indices) > 1)
    {
      for (k in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[k]]][]
    }
    land_covers2[[j]] = land_cover
  }
  land_covers2 = lapply(land_covers2, function(x) aggregate(x, fact = 2)) # aggregate to get same resolution as climate data
  LAND_USE[[i]] = land_covers2
}


# average subseted variables per year interval
year_intervals = list(2025:2049,2050:2074,2075:2100) # future year intervals

mean_land_cover = list()
for (i in seq_along(scenarios)) {
  lc = LAND_USE[[i]]       
  mean_land_cover[[i]] = list()
  
  for (y in seq_along(year_intervals)) {
    lc_mean = lapply(lc, function(x) {
      yrs = as.numeric(sub(".*_", "", names(x)))
      
      idx = which(yrs %in% year_intervals[[y]])
      mean(x[[idx]])   
    })
    mean_land_cover[[i]][[y]] = lc_mean # final list per scenario, year interval and variable
  }
}

  # 2. future data extraction
  # -------------------------
areas = brick("Rasters/Environmental_rasters/clm45_area.nc4"); areas = areas/(1000*1000) # conversion to km2 # total_population = population_density*areas
years_start = c(2025,2050,2075)
years_end = c(2049,2074,2100)
year_intervals = c("2025-2049","2050-2074","2075-2100")

MODELS_isimip3b = c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2","CanESM5-2",
                    "CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3", "MIROC6")
models_isimip3b = c("gfdl-esm4","ipsl-cm6a-lr","mpi-esm1-2-hr","mri-esm2-0","ukesm1-0-ll", "canesm5",
                    "cnrm-cm6-1","cnrm-esm2-1","ec-earth3","miroc6")
nmodels = length(models_isimip3b)


options(future.globals.maxSize = 8 * 1024^3)
plan(multisession, workers = parallel::detectCores() - 4) # parallel across scenarios

future_lapply(seq_along(scenarios), function(s) {
  
  # Load population raster per scenario
  Pop_D = brick(paste0("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_2026_2100_",scenarios2[[s]], ".nc"))
  totalPopulation = Pop_D * areas
  
  scenario_predictions = list()
  log_file = paste0("github_data/Brt_projections_outputs/log_file_", scenarios[s], ".txt")
  
  for (y in seq_along(years_start)) {
    year_start = years_start[y]; year_end = years_end[y]
    year_predictions = list()
    

    # Month indicies 
    month_start = ((year_start - 1) - 2015) * 12 + 1
    month_end = (year_end - 2015) * 12 + 12
    
    for (m in seq_along(models_isimip3b)) {
      
      output_file <- paste0("github_data/Brt_projections_outputs/Future_projections/data_",scenarios[[s]], "_", year_intervals[[y]], "_", MODELS_isimip3b[[m]], "_120526.csv")
      # if env. data exisits use extacted data
      if (file.exists(output_file)) {
        message("Using existing file: ", output_file)
        nutsM_future = read.csv(output_file)
      } else {
        # Load climate rasters
        temp_files = list.files(
          path = paste0("Rasters/Environmental_rasters/ISIMIP3b/", scenarios[s], "/", models_isimip3b[m]),
          pattern = "tas.*\\.nc$", full.names = TRUE)
        temperature = do.call(addLayer, lapply(temp_files, brick))
        
        pr_files = list.files(
          path = paste0("Rasters/Environmental_rasters/ISIMIP3b/", scenarios[s], "/", models_isimip3b[m]),
          pattern = "pr.*\\.nc$", full.names = TRUE)
        precipitation = do.call(addLayer, lapply(pr_files, brick))
        
        rh_files = list.files(
          path = paste0("Rasters/Environmental_rasters/ISIMIP3b/", scenarios[s], "/", models_isimip3b[m]),
          pattern = "hurs.*\\.nc$", full.names = TRUE)
        relative_humidity = do.call(addLayer, lapply(rh_files, brick))
        

        # Select relevant months
        temperature_temp = temperature[[month_start:month_end]]
        precipitation_temp = precipitation[[month_start:month_end]]
        relativehumidity_temp = relative_humidity[[month_start:month_end]]
        
        # Seasonal indices
        winter_indices = sort(c(12,1,2)) + rep((year_start:year_end - year_start) * 12, each = 3)
        winter_indices = winter_indices[-c(1,2,length(winter_indices))]
        spring_indices = sort(c(3,4,5) + rep(((year_start+1):year_end - year_start) * 12, each = 3))
        summer_indices = sort(c(6,7,8) + rep(((year_start+1):year_end - year_start) * 12, each = 3))
        fall_indices = sort(c(9,10,11) + rep(((year_start+1):year_end - year_start) * 12, each = 3))
        
        # Seasonal means
        tempm_winter = mean(temperature_temp[[winter_indices]]) - 273.15
        tempm_spring = mean(temperature_temp[[spring_indices]]) - 273.15
        tempm_summer = mean(temperature_temp[[summer_indices]]) - 273.15
        tempm_fall   = mean(temperature_temp[[fall_indices]]) - 273.15
        
        seasonal_var = tempm_summer - tempm_winter
        
        precp_winter = mean(precipitation_temp[[winter_indices]]) * 60*60*24
        precp_spring = mean(precipitation_temp[[spring_indices]]) * 60*60*24
        precp_summer = mean(precipitation_temp[[summer_indices]]) * 60*60*24
        precp_fall   = mean(precipitation_temp[[fall_indices]]) * 60*60*24
        
        relh_winter = mean(relativehumidity_temp[[winter_indices]])
        relh_spring = mean(relativehumidity_temp[[spring_indices]])
        relh_summer = mean(relativehumidity_temp[[summer_indices]])
        relh_fall   = mean(relativehumidity_temp[[fall_indices]])
        
        # Land-use and population
        croplands_temp = mean_land_cover[[s]][[y]][[1]]
        pastures_temp  = mean_land_cover[[s]][[y]][[2]]
        rangelands_temp = mean_land_cover[[s]][[y]][[3]]
        urbanAreas_temp = mean_land_cover[[s]][[y]][[4]]
        primaryForest_temp = mean_land_cover[[s]][[y]][[5]]
        primaryNonForest_temp = mean_land_cover[[s]][[y]][[6]]
        secondaryForest_temp = mean_land_cover[[s]][[y]][[7]]
        secondaryNonForest_temp = mean_land_cover[[s]][[y]][[8]]

        
        year_start_temp = if(y==1) 1 else year_start - 2025
        year_end_temp   = year_end - 2025
        population = mean(totalPopulation[[year_start_temp:year_end_temp]])
        
        # Combine environmental variables
        envVariables = list(
          tempm_winter, tempm_spring, tempm_summer, tempm_fall,
          seasonal_var,precp_winter, precp_spring, precp_summer, precp_fall,
          relh_winter, relh_spring, relh_summer, relh_fall,
          croplands_temp, pastures_temp, rangelands_temp, urbanAreas_temp,
          primaryForest_temp, primaryNonForest_temp, secondaryForest_temp, secondaryNonForest_temp,population
        )
        
        # Prepare data.frame
        nutsM_future = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(names_env)+1))
        colnames(nutsM_future) = c("NUTM", names_env)
        nutsM_future$NUTM = nutsM_IDs$nutM_code
        
        # Extract values
        for (v in seq_along(envVariables)) {
          fun_type = if (v == length(envVariables)) "sum" else "mean"
          nutsM_future[, v+1] <- exactextractr::exact_extract(envVariables[[v]], nutsM, fun = fun_type)
        }
        nutsM_future$population_density = nutsM_future$total_population / nutsM_sf$area
        included_vars_future = nutsM_future[,c("NUTM",var_names)]
        
        # Save CSV
        write.csv(nutsM_future, output_file, row.names = FALSE)
        write.csv(included_vars_future, paste0("github_data/Brt_projections_outputs/Future_projections/included_vars_",
                                               scenarios[[s]],"_",year_intervals[[y]], "_", MODELS_isimip3b[[m]],"_120526.csv"), row.names = FALSE)
        
      } 
      
      # 3. future projections
      # ----------------------
      model_predictions = list()
      for (j in seq_along(brt_model_scvs)) {
        object = brt_model_scvs[[j]]
        n.trees = object$gbm.call$best.trees
        type = "response"
        single.tree = FALSE

        msg = paste("Processing", scenarios[[s]], "year interval", year_intervals[[y]],
                    "model", m, "of", length(models_isimip3b), "BRT run", j)
        write(msg, file = log_file, append = TRUE)

        model_predictions[[j]] = predict.gbm(object, nutsM_future[, -1], n.trees, type, single.tree)
      }
      year_predictions[[m]] = model_predictions
    } # end models_isimip3b
    scenario_predictions[[y]] = year_predictions
  } # end years
  
  saveRDS(scenario_predictions, file = paste0("github_data/Brt_projections_outputs/Future_projections/future_predictions_", scenarios[s], "_120526.rds"))
  return(paste("Scenario", scenarios[s], "done"))
})


  # 3.1 Saving future projections
ssp1 = readRDS("github_data/Brt_projections_outputs/Future_projections/future_predictions_ssp126_120526.rds")
ssp3 = readRDS("github_data/Brt_projections_outputs/Future_projections/future_predictions_ssp370_120526.rds")
ssp5 = readRDS("github_data/Brt_projections_outputs/Future_projections/future_predictions_ssp585_120526.rds")
predictions_future = list (ssp1,ssp3,ssp5)

year_intervals = c("2025-2049","2050-2074","2075-2100")
path = "github_data/Brt_projections_outputs/Future_projections/"
for (s in 1:length(scenarios)) {
  for (y in 1:length(year_intervals)) {
    for (m in 1:length(models_isimip3b)){
      df = as.data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = 100))
      for (j in 1:100) {
      df[,j] = predictions_future[[s]][[y]][[m]][[j]]
    }
    df = cbind(NUTSM = nutsM_IDs$nutM_code, df); colnames(df)[2:101] = paste0("rep_",1:100)
    write.csv(df, paste0(path,"Future_preds_",scenarios[[s]],"_",year_intervals[[y]],"_",MODELS_isimip3b[[m]],"_120526.csv"), row.names = F)
    }
  }
}

  # 3.2 extracting population counts and saving future population counts per NUTS3 area
pop = brick("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_1900_2025_historical.nc")
pop_2024 = pop[[125]] ; pop_2024 = pop_2024 * areas
pop_2025 = pop[[126]] 

future_population = list()
for (s in 1:length(scenarios)) {
fp = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = 5 ))
colnames(fp) = c("NUTS_ID", "2024", "2025-2049", "2050-2074", "2075-2100")
fp$NUTS_ID = nutsM_IDs$nutM_code
fp$"2024" = exactextractr::exact_extract(pop_present, nutsM, fun = "sum")

Pop_D = brick(paste0("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_2026_2100_",scenarios2[[s]], ".nc"))
Pop_D = addLayer(pop_2025, Pop_D) # add 2025 population density to the brick of future population densities
totalPopulation = Pop_D * areas

  for (y in 1:length(year_intervals)) {
    year_start = years_start[y]; year_end = years_end[y]
    start_idx = year_start - 2025 +1
    end_idx  = year_end - 2025 + 1
    subset = mean(totalPopulation[[start_idx:end_idx]])
    fp[,y+2] = exactextractr::exact_extract(subset, nutsM, fun = "sum")
  
  }
future_population[[s]] = fp
write.csv(fp,paste0("github_data/Environmental_data/Population_density_data/NUTS3_population_",scenarios[[s]],"_010526.csv"), row.names = F)
}

  # 4. future predictions figures
  # -----------------------------

present = read.csv(paste0("github_data/Brt_projections_outputs/Current_projections/current_predictions_120526.csv"), head=T) # load 2024 predictions 
brts = list(); pops = list(); scenarios = c("ssp126","ssp370","ssp585"); periods = c("2025-2049","2050-2074","2075-2100")
models = c("CanESM5-2","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2")
models = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")

scenarios3 = c("ssp1-2.6","ssp3-7.0","ssp5-8.5")
for (i in 1:length(scenarios)) {
  brts[[scenarios[i]]] = list() 
  for (j in 1:length(periods)) {
    buffer_models = list()
    for (k in 1:length(models)) {
      path = paste0("~/Desktop/simon_enm/All_BRT_projection_outputs/Future_projections/Future_", scenarios3[i], "_", periods[j], "_", models[k], ".csv") #github_data/Brt_projections_outputs/Future_projections/
      dat = read.csv(path, head=T)
      buffer_models[[models[k]]] = dat[, 2:ncol(dat)]
    }
    brts[[scenarios[i]]][[periods[j]]] = buffer_models
  }
}

# load population counts per scenario
for (i in 1:length(scenarios)) {
  pops[[scenarios[i]]] = read.csv(paste0("github_data/Environmental_data/Population_density_data/NUTS3_population_",scenarios[i],"_010526.csv"), head=T)
}

# Plotting parameters
nberOfReplicates = 100
suffix = "A"; plottingDifference = TRUE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "B"; plottingDifference = FALSE; plottingStDevs = TRUE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "C"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "D"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = TRUE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "E"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = TRUE; plottingBinaryMap = FALSE
suffix = "F"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = TRUE

# plotting
pdf(paste0("github_data/Figures/Figure_2",suffix,"_010526.pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")		
scenario_names = c("SSP1-2.6","SSP3-7.0","SSP5-8.5")
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111]) # main figure color scale
# colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
# colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # "D9" = 85%
# colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") # "00" = 0%
# colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5"
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18]) # difference map color scale
colourScale2[9] = paste0(colourScale2[9],"00") # "00" = 0% transparency for the middle color (white) 
colourScale3 = colourScale2; colourScale3[9] = "#E5E5E5" 
colourScale4 = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[11:111]
I1 = 1; I2 = length(scenarios); c = 0; colourScale = colourScale1
if (plottingStDevs) { colourScale = colourScale4 }
if (plottingDifference) { I1 = 2; I2 = 2 }
for (i in I1:I2)
{
  c = c+1
  plot.new()
  if (c == 1)
  {
    brt = as.matrix(present[,2:dim(present)[2]]); brt1 = brt
    means = rep(NA, dim(brt)[1]); stdvs = rep(NA, dim(brt)[1])
    lowerCI95 = rep(NA, dim(brt)[1]); upperCI95 = rep(NA, dim(brt)[1])
    binaries = rep(0, dim(brt)[1])
    for (j in 1:dim(brt)[1])
    {
      means[j] = mean(brt[j,]); stdvs[j] = sd(brt[j,])
      lowerCI95[j] = t.test(brt[j,])$conf.int[1]
      upperCI95[j] = t.test(brt[j,])$conf.int[2]
      temp = rep(0, dim(brt)[2])
      temp[which(brt[j,]>=thresholds_SI[,"threshold"])] = 1
      binaries[j] = round(mean(temp))
    }
    cols = colourScale[(((means-0)/(1-0))*100)+1]; # print(max(stdvs))
    if (plottingStDevs) cols = colourScale[(((stdvs-0)/(0.25-0))*100)+1]
    if (plottingLower95CI) cols = colourScale[(((lowerCI95-0)/(1-0))*100)+1]
    if (plottingUpper95CI) cols = colourScale[(((upperCI95-0)/(1-0))*100)+1]
    if (plottingBinaryMap) cols = c("gray90","#de4327")[binaries+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
    mtext("Present time", side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext("(2020-2024)", side=3, line=-2.5, at=1, cex=0.50, col="gray30"); rast = raster(as.matrix(c(0,1)))
  }	else	{
    plot.new()
  }
  for (j in 1:length(periods))
  {
    means1 = rep(NA, length(means)); stdvs1 = rep(NA, length(means))
    lowerCI95 = rep(NA, dim(brt)[1]); upperCI95 = rep(NA, dim(brt)[1])
    binaries1 = rep(NA, dim(brt)[1])
    for (k in 1:length(means))
    {
      means2 = rep(NA, length(models)); stdvs2 = rep(NA, length(models))
      vS = matrix(nrow=nberOfReplicates, ncol=length(models))
      binaries2 = rep(0, length(models))
      for (l in 1:length(models))
      {
        brt = as.matrix(brts[[i]][[j]][[l]]); vS[,l] = brt[k,]
        means2[l] = mean(brt[k,]); stdvs2[l] = sd(brt[k,])
        temp = rep(0, dim(brt)[2])
        temp[which(brt[k,]>=thresholds_SI[,"threshold"])] = 1
        binaries2[l] = mean(temp)
      }
      means1[k] = mean(means2); stdvs1[k] = mean(stdvs2)
      lowerCI95[k] = t.test(vS)$conf.int[1]
      upperCI95[k] = t.test(vS)$conf.int[2]
      binaries1[k] = round(mean(binaries2))
    }
    cols = colourScale[(((means1-0)/(1-0))*100)+1]; # print(max(stdvs1))
    if (plottingStDevs) cols = colourScale[(((stdvs1-0)/(0.25-0))*100)+1]
    if (plottingLower95CI) cols = colourScale[(((lowerCI95-0)/(1-0))*100)+1]
    if (plottingUpper95CI) cols = colourScale[(((upperCI95-0)/(1-0))*100)+1]
    if (plottingBinaryMap) cols = c("gray90","#de4327")[binaries1+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA)
    plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
    mtext(scenario_names[i], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")	
    if (plottingStDevs) rast = raster(as.matrix(c(0,0.25)))	
    if ((i == length(scenarios))&(j == length(periods))&(!plottingBinaryMap))
    {
      plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
           axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
    }
  }
}
if (!plottingDifference) dev.off()

# difference maps
if (plottingDifference)
{
  vS = c()
  for (i in 1:length(scenarios))
  {
    for (j in 1:length(periods))
    {
      for (k in 1:length(models))
      {
        brt2 = as.matrix(brts[[i]][[j]][[k]])
        for (l in 1:dim(brt1)[1])
        {
          vS = c(vS, mean(brt1[l,]-brt2[l,]))
        }
      }
    }
  }
  vS = vS[which(!is.na(vS))]; minVS = min(vS); maxVS = max(vS)
  if (abs(minVS) < abs(maxVS)) minVS = -maxVS
  if (abs(maxVS) < abs(minVS)) maxVS = -minVS
  maxVS = 0.75; minVS = -0.75
  for (i in 2:3)
  {
    plot.new(); plot.new()	
    for (j in 1:length(periods))
    {
      differences = rep(NA, dim(brt1)[1])
      for (k in 1:dim(brt1)[1])
      {
        buffer = rep(NA,length(models))
        for (l in 1:length(models))
        {
          brt2 = as.matrix(brts[[i]][[j]][[l]]); buffer[l] = mean(brt2[k,]-brt1[k,])
        }
        differences[k] = mean(buffer)
      }
      cols = colourScale2[(((differences-minVS)/(maxVS-minVS))*length(colourScale2))+1]
      plot(coasts, lwd=0.4, border="gray30", col=NA)
      plot(coasts, add=T, border=NA, col="#E5E5E5")
      plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
      mtext(paste0(scenario_names[i]," - t0"), side=3, line=-1.7, at=1, cex=0.50, col="gray30")
      mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
      if ((i == length(scenarios))&(j == length(periods)))
      {
        rast = raster(as.matrix(c(minVS,maxVS))); par(lwd=0.1, fg="white")
        plot(rast, legend.only=T, add=T, col="#E5E5E5", legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96),
             legend.args=list(text="", cex=0.7, line=0.0, col=NA, col.axis=NA), col.tick=NA, adj=3, alpha=1, side=3,
             axis.args=list(cex.axis=0.65, lwd=0, col=NA, lwd.tick=0, col.tick=NA, tck=-1.2, col.axis=NA, line=0, mgp=c(0,0.45,0)))
        par(lwd=0.2, fg="gray30")
        plot(rast, legend.only=T, add=T, col=colourScale3, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
             axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, 
                            mgp=c(0,0.45,0), at=c(-0.6,-0.3,0,0.3,0.6)), alpha=1, side=3)
      }
    }
  }
  dev.off()
}

  # 5. Quantifying the population at risk of exposure
  # -------------------------------------------------
nruns = 100
exposure_outputs = 10  # 1 baseline + 9 future combinations

# Current (2000-2024) suitability values per NUTS3 and run
suitability_present = as.matrix(present)

# outputs
pop_exposure = data.frame(matrix(NA, nrow = nruns, ncol = exposure_outputs))

# Column names
colnames(pop_exposure)[1] = "Baseline_2024"
colnames(pop_exposure)[2:exposure_outputs] = paste0(
  rep(periods, each = length(scenarios)), "_",
  rep(scenarios, length(periods)))

# Baseline calculation (2024 population)
for (r in 1:nruns) {
  
  # Regions above suitability threshold
  threshold_value = as.numeric(thresholds_SI[r, "threshold"])
  valid_regions = which(suitability_present[, r] >= threshold_value)
  
  # Sum population in those regions (baseline 2024)
  pop_exposure[r, "Baseline_2024"] = sum(pops[[1]][valid_regions, "X2024"])
}
pop_exposure_fixed_pop = pop_exposure  # population fixed at 2024
pop_exposure_fixed_suit = pop_exposure  # suitability fixed

for (p in seq_along(periods)) {
  for (s in seq_along(scenarios)) {
    for (r in 1:nruns) {
      
      # compute mean suitability across models for each region
      n_regions = nrow(nutsM@data)
      n_models  = length(models)
      
      suitability_models = matrix(NA, nrow = n_regions, ncol = n_models)
      
      for (m in 1:n_models) {
        suitability_models[, m] = brts[[s]][[p]][[m]][, r]
      }
      
      mean_suitability = rowMeans(suitability_models)
      
      # find regions above threshold ---
      threshold_value = thresholds_SI[r, "threshold"]
      valid_regions_future = which(mean_suitability >= threshold_value)
      
      period_col = paste0("X", gsub("-", ".", periods[p])) # column name
      col_index = 1 + ((p - 1) * length(scenarios)) + s
      
      # (1) Main calculation (both suitability and population change) 
      pop_exposure[r, col_index] = sum(pops[[s]][valid_regions_future, period_col])
      
      
      # (2) Fixed population (only suitability changes)
      pop_exposure_fixed_pop[r, col_index] = sum(pops[[s]][valid_regions_future, "X2024"])
      
      
      # (3) Fixed suitability (only population changes)
      valid_regions_present = which(suitability_present[, r] >= threshold_value)
      pop_exposure_fixed_suit[r, col_index] = sum(pops[[s]][valid_regions_present, period_col])
    }
  }
}

boxplot(cbind(pop_exposure,pop_exposure_fixed_pop,pop_exposure_fixed_suit))

for (h in 1:3)
{
  if (h == 1) pop_exposures = pop_exposure
  if (h == 2) pop_exposures = pop_exposure_fixed_pop
  if (h == 3) pop_exposures = pop_exposure_fixed_suit
  message("Shapiro test:")
  for (i in 1:dim(pop_exposures)[2])
  {
    print(shapiro.test(pop_exposures[,i])$p.value)
  }
  pValues = matrix(0, nrow=dim(pop_exposures)[2], ncol=4) # Benjamini-Hochberg correction
  colnames(pValues) = c("p-value","rank","new_threshold_B-H","significant_after_B-H")
  message("wilcox test:")
  for (i in 2:dim(pop_exposures)[2])
  {
    print(round(wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value,5))
    pValues[i-1,1] = wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value
  }
  for (i in 2:dim(pop_exposures)[2])
  {
    pValues[i-1,2] = min(which(pValues[order(pValues[,1]),1]==pValues[i-1,1]))
    pValues[i-1,3] = (pValues[i-1,2]/dim(pValues)[1])*0.05
    if (pValues[i-1,1] < pValues[i-1,3]) pValues[i-1,4] = 1
    pValues[i-1,1] = round(pValues[i-1,1],3)
    pValues[i-1,3] = round(pValues[i-1,3],3)
  }
}	# for "pop_exposure": the three first paired Wilcoxon lead to adjusted p-values > 0.05

# pop_exposure: Future scenarios + Future Population
# - The first 3 Wilcoxon tests (Period 1 across SSPs) show p-values ~1.0.
# - In the immediate future (Period 1), the combined risk is not yet statistically higher than the 2024 baseline.
# - Risk increases significantly in mid-to-late century.

# pop_exposure_fixed_pop: Future scenarios + Fixed 2024 Population
# - he first 3 Wilcoxon tests (Period 1 across SSPs) show p-values ~1.0.
# - Climate change alone is a massive driver of risk, as even if the population remained frozen at current levels, 
#   the expansion of suitability would still put significantly more people in at risk areas over time.

# pop_exposure_fixed_suit: Fixed 2024 scenario + Future Population: no significant increase in risk over time.



  # 6. Estimating the ecological suitability variation at the continent scale
  # -------------------------------------------------------------------------

# Choose continent or country to consider
 region_indices = 1:nrow(nutsM@data) # all European countries
 region_indices = which(nutsM@data$CountryNam %in% c("Belgium", "Netherlands")) # Belgium and The Netherlands

# GCM models
climate_models = c("CanESM5-2","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2")


# current suitability subset to chosen region (rows = NUTS3 regions, columns = brt suitability values per run)
current_suitability = read.csv("github_data/Brt_projections_outputs/Current_projections/current_predictions_130526.csv")[region_indices, ]

# future suitability for 2 climate change scenarios (SSP3 & SSP5)
future_suitability_ssp3 = list()
future_suitability_ssp5 = list()

# subset to chosen region and load future suitability for each model and scenario
for (m in seq_along(climate_models)) {
  future_suitability_ssp3[[m]] = read.csv(paste0("github_data/Brt_projections_outputs/Future_projections/","Future_preds_ssp370_2075-2100_", climate_models[m], "_120526.csv"))[region_indices, -1]
  future_suitability_ssp5[[m]] = read.csv(paste0("github_data/Brt_projections_outputs/Future_projections/","Future_preds_ssp585_2075-2100_", climate_models[m], "_120526.csv"))[region_indices, -1]
}

# percentage change in suitability function
percent_change = function(future, current) {((future - current) / current) * 100}

# function to compute mean percentage change in suitability across all models for each run
compute_mean_change = function(current, future_list) {
  
  nruns = ncol(current)
  n_models = length(future_list)
  results = numeric(nruns)
  
  for (r in 1:nruns) {
    model_means = numeric(n_models)
    for (m in 1:n_models) {
      change_matrix = percent_change(future_list[[m]][, r],current[, r])
      model_means[m] = mean(change_matrix, na.rm = TRUE)
    }
    results[r] = mean(model_means)
  }
  return(results)
}

change_ssp3 = compute_mean_change(current_suitability, future_suitability_ssp3)
change_ssp5 = compute_mean_change(current_suitability, future_suitability_ssp5)

# summarizing function (mean and confidence intervals)
summarize = function(x) {
  m  = round(mean(x))
  ci = round(t.test(x)$conf.int) 
  c(mean = m, lower = ci[1], upper = ci[2])
}

summarize(change_ssp3) # for all countries: 62%, 95% CI = [55, 69]; for Belgium and The Netherlands: 96%, 95% CI = [85, 107]
                       # suitability nearly 1.5  across Europe (change of 62%),  twice as suitable for Belgium and The Netherlands (change of 96%)

summarize(change_ssp5) # for all countries: 115%, 95% CI = [102, 127]; for Belgium and The Netherlands: 100%, 95% CI = [88, 112]
                       # suitability over twice as suitable across Europe (change of 115%), doubles for Belgium and The Netherlands (change of 100%)

  # 7. Estimating the population at risk variation at the continent scale
  # ----------------------------------------------------------------------
thresholds = thresholds_SI$threshold
# subset to chosen region and load population counts for each scenario
pop_ssp3 = read.csv("github_data/Environmental_data/Population_density_data/NUTS3_population_ssp370_010526.csv")[region_indices, ]
pop_ssp5 = read.csv("github_data/Environmental_data/Population_density_data/NUTS3_population_ssp585_010526.csv")[region_indices, ]

calculate_population_exposure = function(current_suit, future_suit_list, pop_data) {
  nruns = ncol(current_suit)
  n_models = length(future_suit_list)
  
  exposure = matrix(NA, nrow = nruns, ncol = 2)
  colnames(exposure) = c("baseline_2024", "future_2075_2100")
  
  for (r in 1:nruns) {
    threshold = thresholds[r, "thresholds"]
    regions_baseline = which(current_suit[, r] >= threshold)
    exposure[r, 1] = sum(pop_data[regions_baseline, "X2024"])
    
    future_values = numeric(n_models)
    
    for (m in 1:n_models) {
      regions_future = which(future_suit_list[[m]][, r] >= threshold)
      future_values[m] = sum(pop_data[regions_future, "X2075.2100"])
    }
    exposure[r, 2] = mean(future_values)
  }
  return(exposure)
}

# calculate population exposure for SSP3 and SSP5 scenarios
exposure_ssp3 = calculate_population_exposure(current_suitability, future_suitability_ssp3, pop_ssp3) # row = sum of all NUTS3 regions above threshold for each run, column 1 = population at risk in 2024, column 2 = population at risk in 2075-2100 ssp3
exposure_ssp5 = calculate_population_exposure(current_suitability, future_suitability_ssp5, pop_ssp5)

# current population at risk (2024 baseline)
print(summarize(exposure_ssp3[,1]))
# 116 484 623 CI95%[110 942 922,122 026 325]  people or approx. 116 million [CI95% 111,122 million] 

# future population at risk (2075-2100)
print(summarize(exposure_ssp3[,2])) # ssp3: 158 548 155 CI95%[148 774 697,168 321 612] people or approx. 159 million [CI95% 149,168 million]
print(summarize(exposure_ssp5[,2])) # ssp5: 246 607 066 CI95%[231 160 080,262 054 051] people or approx. 247 million [CI95% 231,262 million]

# number of addtional people at risk in the future compared to 2024 baseline
delta_ssp3 = exposure_ssp3[,2] - exposure_ssp3[,1]
delta_ssp5 = exposure_ssp5[,2] - exposure_ssp5[,1]

summarize(delta_ssp3) # 42 063 531 CI95%[35 969 981,48 157 082] or approx. 42 million [CI95% 36,48 million] additional people at risk 
summarize(delta_ssp5) # 130 122 442 CI95%[118 866 220,141 378 664] or approx. 130 million [CI95% 119,141 million] additional people at risk



