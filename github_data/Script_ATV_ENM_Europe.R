# Geospatial libraries
library(ncdf4)
library(exactextractr)
library(raster)
library(terra)
library(sf)
library(sp)

# Data manipulation libraries
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

# Statistical modeling libraries
library(gmodels)
library(devtools)
library(dismo)
library(gbm)
library(blockCV)

# Plotting libraries
library(RColorBrewer)
library(scales)

# Parallel processing libraries
library(furrr)
library(future.apply)



##################################  SUMMARY  ###################################  

    ## I. DATA PREPARATION
        # 1. Optimized admin. area European shapefile (VectorNet)  
        # 2. Aedes-borne viral (DENV, CHIKV, ZIKV) occurrence data curation

    ## II. BRT MODEL TRAINING AND EVALUATION
        # 1. Pseudo-absence sampling setup
        # 2. Create pseudo-absence shapefile
        # 3. Environmental data extraction
        # 4. BRT model training
        # 5. Model predictive performance evaluation

    ## III. VARIABLE RELATIVE INFLUENCES AND RESPONSE CURVES
        # 1. Relative influences 
        # 2. Response Curves

    ## IV. CURRENT RISK MAP (2024)
        # 1. 2024 data extraction
        # 2. 2024 risk map

    ## V. PAST PROJECTIONS
        # 1. Data extraction and historical predictions
        # 2. Saving past predictions

    ## VI. FUTURE PROJECTIONS
        # 1. Data extraction and future projections
        # 2. Saving future projections

    ## VII. VISUILASATION SCRIPT (code by Dr Simon Dellicour)
        # 1. Preparing the optimized NUTS3 shapefile
        # 2. Generating the figure with the past projections
        # 3. Generating the figure with the future projections
        # 4. Estimating the variation at the continent scale
        # 5. Estimating the variation for Belgium (Le Soir)


################################################################################
#                         I. DATA PREPARATION
################################################################################


# 1. Optimized admin. area European shapefile (VectorNet)  
#########################################################

nutsM = crop(shapefile("Shapefiles/MOOD_EUROPE_NUTS3_WillyW/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # modified NUT3 shapefile
# Only keep EU areas
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
# remove the small islands where land-use data is not available
nutsM = subset(nutsM, !LocationCo%in%c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533","ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65")) 

nutsM_sf = st_as_sf(nutsM)
nutsM_sf$area = nutsM_sf %>% st_area %>% as.numeric/(1000*1000) ; #st_write(nutsM_sf, "Shapefiles/nutsM.gpkg")

# Correspondences from NUTS3 to optimized admin. areas (NUTSM)
correspondences = shapefile("Shapefiles/MOOD_EUROPE_NUTS3_WillyW/vnMOODdatamapcodejoinpoly.shp")@data[,c("DATLOCODE","MAPLOCODE")]

# Subset polygon IDs
nutsM_IDs = as.data.frame(matrix(nrow = dim(nutsM@data)[1], ncol = 3))
colnames(nutsM_IDs) = c("country","region_name", "nutM_code"); nutsM_IDs$nutM_code = nutsM@data[,"LocationCo"]
nutsM_IDs$country = nutsM@data[,"CountryNam"]; nutsM_IDs$region_name = nutsM@data[,"LocationNa"]



# 2. Aedes-borne viral (DENV, CHIKV, ZIKV) occurrence data curation
###################################################################
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
}else{Aedes_presences = read.csv("github_data/Occurrence_data/Aedes_viruses_European_presences.csv", header = T)}

# clean environment
keep = c("Aedes_presences", "nutsM", "contour","nutsM_sf","nutsM_IDs") 
rm(list = setdiff(ls(), keep))


################################################################################
#                   II. BRT MODEL TRAINING AND EVALUATION
################################################################################

# 1. Pseudo-absence sampling setup
##################################

# Group presence data by year ranges
# To avoid over weighting admin areas with repeated presences across years, 
# we keep only the year range per location. This allows averaging 
# environmental variables over that period.


# create the function
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

# 2. Create pseudo-absence shapefile
####################################
# For pseudo-absence (PA) sampling, remove polygons adjacent to presence sites 
# in order to avoid sampling pseudo-absences in areas too close to known presences.

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


# 3. Environmental data extraction
########################################
 

# Data extraction #
plan(multisession, workers = parallel::detectCores() - 1) # set parallele process
nruns = 100 # this will be used as the number of repetitions for the brt training
data_train = list()

for(n in 1:nruns) {
  indices = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices, 1*length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$Probabilities) 
  absences = Aedes_abs[selected_absences, ]
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign pseudo-years
  absences$year = sample(years, size = length(absences$year), replace = TRUE, prob = year_probs)
  
  # Add response
  absences$response = 0
  absences = absences[, !(colnames(absences) %in% c("surv", "Probabilities"))]
  Aedes_abs = Aedes_abs[, !(colnames(Aedes_abs) %in% c("surv", "Probabilities"))]
  data = rbind(ae.presences_agg, absences) # bind presences and pseudo-absences
  
  
  ## data extraction ##
  data_for_brt = as.data.frame(matrix(NA,nrow=nrow(data), ncol =28))
  data_for_brt[,c(1:5)] = data[,c("nutM_code","year","year_start","year_end", "response")]
  colnames(data_for_brt) = c("nutM_code","year","year_start","year_end", "response","temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                              "seasonal_variation","precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                              "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                              "croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonForest","secondaryForest","secondaryNonForest", "total_population","population_density")
  
  geometry_matches = nutsM_sf$geometry[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
  area_matches = nutsM_sf$area[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
  nutsM_train = st_sf(data_for_brt, geometry = geometry_matches); nutsM_train$area = area_matches
  names_env = names(data_for_brt[,c(6:27)])
    
  cat("Starting data extraction", n, "of", nruns)
  
  # parallel data extraction  
  results = future_lapply(1:nrow(data_for_brt), function(k) {
    if (!is.na(data_for_brt$year[k])) { 
      # For temporal training, split climate variables into seasonal subsets
      message("temporal training on single year")
      year = as.numeric(data_for_brt[k,"year"])
      year_index = year - 1900
      
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
        
        for (yr in year_start:year_end) {
          month_start=((yr - 1) - 1901) * 12 + 1
          winter_indices=c(winter_indices, month_start + 11, month_start + 12, month_start + 13) # December (previous year), January, February
          spring_indices=c(spring_indices, month_start + 14, month_start + 15, month_start + 16) # March, April, May
          summer_indices=c(summer_indices, month_start + 17, month_start + 18, month_start + 19) # June, July, August
          fall_indices  =c(fall_indices,   month_start + 20, month_start + 21, month_start + 22) # September, October, November
        }
      }
      
    # Extract environmental variables for the given year or year range
    covariates_temp = list()
    covariates_temp[["temp_winter"]] = mean(covariates[[1]][[winter_indices]]) - 273.15 # conversion to Celsius
    covariates_temp[["temp_spring"]] = mean(covariates[[1]][[spring_indices]]) - 273.15
    covariates_temp[["temp_summer"]] = mean(covariates[[1]][[summer_indices]]) - 273.15
    covariates_temp[["temp_fall"]]   = mean(covariates[[1]][[fall_indices]]) - 273.15
    
    covariates_temp[["seasonal_var"]] = covariates_temp[["temp_summer"]] - covariates_temp[["temp_winter"]]
    
    covariates_temp[["prec_winter"]] = mean(covariates[[2]][[winter_indices]]) * 60 * 60 * 24 # conversion to vers
    covariates_temp[["prec_spring"]] = mean(covariates[[2]][[spring_indices]]) * 60 * 60 * 24
    covariates_temp[["prec_summer"]] = mean(covariates[[2]][[summer_indices]]) * 60 * 60 * 24
    covariates_temp[["prec_fall"]]   = mean(covariates[[2]][[fall_indices]]) * 60 * 60 * 24
    
    covariates_temp[["relh_winter"]] = mean(covariates[[3]][[winter_indices]])
    covariates_temp[["relh_spring"]] = mean(covariates[[3]][[spring_indices]])
    covariates_temp[["relh_summer"]] = mean(covariates[[3]][[summer_indices]])
    covariates_temp[["relh_fall"]]   = mean(covariates[[3]][[fall_indices]])
    
    covariates_temp[["crops"]] = mean(covariates[[4]][[year_index]])
    covariates_temp[["pastures"]] = mean(covariates[[5]][[year_index]])
    covariates_temp[["urbanAreas"]] = mean(covariates[[6]][[year_index]])
    covariates_temp[["rangelands"]] = mean(covariates[[7]][[year_index]])
    covariates_temp[["primForest"]] = mean(covariates[[8]][[year_index]])
    covariates_temp[["primNonForest"]] = mean(covariates[[9]][[year_index]])
    covariates_temp[["secForest"]] = mean(covariates[[10]][[year_index]])
    covariates_temp[["secNonForest"]] = mean(covariates[[11]][[year_index]])
    covariates_temp[["totalPop"]] = mean(covariates[[12]][[year_index]])
    
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
  data_for_brt[,28] = data_for_brt$total_population/nutsM_train$area # population density by NUTS3
  data_for_brt$year = as.integer(unlist(data_for_brt$year))
  data_train[[n]] = data_for_brt
  write.csv(data_for_brt,paste0("github_data/Environmental_data/Brt_training_data/European_BRT_data/20crv3-era5_env_data_replicate",n,".csv"), row.names = F)  
  cat("Data extraction run", n, "of", nruns, "completed.\n")
} 
  
#saveRDS(data_train, file = "European_model/All_the_brt_models/training_data.rds")



# 4. BRT model training
#######################

# Calculate the block size for spatial cross-validation
data = data_train[[1]]

centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  


matching_centroids = as.data.frame(matrix(nrow = nrow(data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
matching_centroids$NUTM = data$nutM_code 
matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
correlogram = ncf::correlog(matching_centroids[,2], matching_centroids[,3], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)

plottingCorrelogram = FALSE
if (plottingCorrelogram == TRUE)
{
  dev.new(width=4.5, height=3); par(mar=c(2.5,2.2,1.0,1.5))
  plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,3500))
  abline(h=0, lwd=0.5, col="red", lty=2)
  lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
  points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray90", pch=16)
  points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.25, col="gray30", pch=1)
  axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
  axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
  title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
  title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
}# range size 680 km2

# brt model training
brt_model_scvs = list() # brt with spatial cross-validations (SCVs)
tabs_list1 = list()

AUCs = matrix(nrow = nruns, ncol = 1); colnames(AUCs) = "AUCs"
SIppcs = matrix(nrow = nruns, ncol = 1); colnames(SIppcs) = "SIpcc"
thresholds = matrix(nrow = nruns, ncol = 1); colnames(thresholds) = "thresholds"

for(i in 1:length(data_train)){
  current_data = data_train[[i]]
  
  # BRT model parameters
  gbm.x = 6:28 # covariate columns
  gbm.y = colnames(current_data)[5]
  offset = NULL
  tree.complexity = 5 # "tc" = number of nodes in the trees
  learning.rate = 0.001 # "lr" = contribution of each tree to the growing model
  bag.fraction = 0.80 # proportion of data used to train a given tree
  site.weights = rep(1, dim(current_data)[1])
  var.monotone = rep(0, length(gbm.x))
  prev.stratify = TRUE
  family = "bernoulli"
  n.trees = 100 # initial number of trees
  step.size = 10 # interval at which the predictive deviance is computed and logged
  # (at each interval, the folds are successively used as test data set
  # and the remaining folds as training data sets to compute the deviance)
  max.trees = 10000 # maximum number of trees that will be considered
  tolerance.method = "auto"
  tolerance = 0.001
  plot.main = TRUE
  plot.folds = FALSE
  verbose = TRUE
  silent = FALSE
  keep.fold.models = FALSE
  keep.fold.vector = FALSE
  keep.fold.fit = FALSE
  showingFoldsPlot = FALSE	
  n.folds = 5
  theRanges = c(680, 680)*1000 # transform block size distance to meters
      
  matching_centroids = as.data.frame(matrix(nrow = nrow(current_data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
  matching_centroids$NUTM = current_data$nutM_code 
  matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
  matching_centroids$response = current_data$response
      
  coordinates(matching_centroids) = ~longitude + latitude
  crs(matching_centroids) = crs(nutsM)
      
  spdf = SpatialPointsDataFrame(matching_centroids, current_data[,1:dim(current_data)[2]])
  myblocks = cv_spatial(spdf, column ="response", k=n.folds, size=theRanges, selection="random", plot = F)
  fold.vector = myblocks$folds_ids; n.trees = 100
      
  brt_model_scvs[[i]] =  gbm.step(current_data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
                                  var.monotone, n.folds, prev.stratify, family, n.trees, step.size,max.trees, tolerance.method, tolerance, plot.main = FALSE, 
                                  plot.folds = FALSE, verbose = TRUE, silent = FALSE, keep.fold.models = FALSE, keep.fold.vector = FALSE, keep.fold.fit = FALSE)
  
  
# 5. Model predictive performance evaluation
#############################################
  
  #### Calculation of the AUC ####
  AUCs[i,] = brt_model_scvs[[i]]$cv.statistics$discrimination.mean # Mean test AUC
  
  #### Calculation of the sorensen index to establish predictive performance of our models ####
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
    tmp[which(tmp[,1]==threshold),2] = sorensen_ppc_tmp
        
    if (sorensen_ppc < sorensen_ppc_tmp) {
      sorensen_ppc = sorensen_ppc_tmp
      optimised_threshold = threshold
      }
    }
  
  tabs_list1[[i]] = tmp
  SIppcs[i,] = sorensen_ppc
  thresholds[i,] = optimised_threshold
}

#saveRDS(brt_model_scvs, file = "European_model/All_the_brt_models/brt_trained_model.rds")
#saveRDS(tabs_list1, file = "European_model/All_the_brt_models/tabs_list1.rds")

# Predictive performance metrics
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

write.csv(metrics, "github_data/Brt_outputs/Brt_metrics/20crv3-era5_global_brt_performances.csv", row.names = FALSE)

################################################################################
#             III. VARIABLE RELATIVE INFLUENCES AND RESPONSE CURVES
################################################################################

# 1. Relative influences 
########################
brt_model_scvs = readRDS("European_model/All_the_brt_models/brt_trained_model.rds")
brt_model_scv = brt_model_scvs

envVariableNames = brt_model_scvs[[1]]$var.names

nruns = 100
relativeInfluences = matrix(0, nrow=length(envVariableNames), ncol=3) #nrow depends on the number of covariates
relativeInfluences_all = matrix(0, nrow=length(envVariableNames), ncol=nruns) #nrow depends on the number of covariates


# get relative influence of each covariate
for (j in 1:length(brt_model_scv)) {
  for (k in 1:length(envVariableNames)) {
      relativeInfluences[k,1] = relativeInfluences[k] + summary(brt_model_scv[[j]])[envVariableNames[k],"rel.inf"]
      relativeInfluences_all[k,j] = summary(brt_model_scv[[j]])[envVariableNames[k],"rel.inf"]
      
    }
  }
  
row.names(relativeInfluences) = envVariableNames
relativeInfluences[,1] = relativeInfluences[,1]/length(brt_model_scv)
relativeInfluences[,2] = as.data.frame(t(apply(as.matrix(relativeInfluences_all), 1, function(x) ci(x))))$`CI lower`
relativeInfluences[,3] = as.data.frame(t(apply(as.matrix(relativeInfluences_all), 1, function(x) ci(x))))$`CI upper`
write.csv(relativeInfluences, file = "European_model/All_the_brt_models/Relative_influences.csv", row.names = TRUE)
write.csv(relativeInfluences_all, file = "European_model/All_the_brt_models/All_relative_influences.csv", row.names = TRUE)

# 2. Response curves  
###############################
data_curv = brt_model_scv[[1]]$gbm.call$dataframe
data_curv = data_curv[data_curv[, "response"] == 1, ] 
envVariableValues = matrix(nrow=3, ncol=length(envVariableNames))
row.names(envVariableValues) = c("median","minV","maxV")
colnames(envVariableValues) = envVariableNames
  
  for (j in 1:length(envVariableNames))
  {
    minV = min(data_curv[,envVariableNames[j]], na.rm=T)
    maxV = max(data_curv[,envVariableNames[j]], na.rm=T)
    medianV = median(data_curv[,envVariableNames[j]], na.rm=T)
    envVariableValues[,j] = rbind(medianV, minV, maxV)
  }
  
  envVariableValues_list = envVariableValues
  
  
  
  # plot with Y axis scaled from 0 to 1
  pdf("Figures/Fig4_response_curves_eu.pdf", width = 4.75, height = 5.125) 
  #dev.new(width = 4.75, height = 5.125)
  par(mfrow=c(6,4), oma=c(0.5,0.5,0,0.5), mar=c(2.3,2,0.3,0), lwd=0.2, col="gray30")
 # oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
  ## prepare for plot (i used 6 rows, 4 columns, for 22 covariates)
  for (i in 1:length(envVariableNames))
  {
    predictions_list = list(); dfs = list()
    valuesInterval = 0.1; valuesInterval = (envVariableValues_list["maxV",i]-envVariableValues_list["minV",i])/100
    df = data.frame(matrix(nrow=length(seq(envVariableValues_list["minV",i],envVariableValues_list["maxV",i],valuesInterval)),ncol=length(envVariableNames)))
    colnames(df) = envVariableNames
    
    for (k in 1:length(envVariableNames))
    {
      valuesInterval = 0.1; valuesInterval = (envVariableValues_list["maxV",k]-envVariableValues_list["minV",k])/100
      if (i == k) df[,envVariableNames[k]] = seq(envVariableValues_list["minV",k],envVariableValues_list["maxV",k],valuesInterval)
      if (i != k) df[,envVariableNames[k]] = rep(envVariableValues_list["median",k],dim(df)[1])
    }
    dfs = df 
    predictions = list()

    for (j in 1:length(brt_model_scv))
    {
      n.trees = brt_model_scv[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      prediction = predict.gbm(brt_model_scv[[j]], newdata=df, n.trees, type, single.tree)
      if ((j == 1))
      {
        minX = min(df[,envVariableNames[i]]) 
        maxX = max(df[,envVariableNames[i]])
        minY = min(prediction); maxY = max(prediction)
      }	else	{
        if (minX > min(df[,envVariableNames[i]])) minX = min(df[,envVariableNames[i]])
        if (maxX < max(df[,envVariableNames[i]])) maxX = max(df[,envVariableNames[i]])
        if (minY > min(prediction)) minY = min(prediction)
        if (maxY < max(prediction)) maxY = max(prediction)
      }
      
      predictions[[j]] = prediction
    }
    predictions_list = predictions
    
    cols = "firebrick2"
    for (l in 1:length(brt_model_scv))
    {
      if ((l == 1))
      {
        plot(dfs[,envVariableNames[i]],predictions_list[[l]],col=cols,ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(0,1))
      }	else	{
        lines(dfs[,envVariableNames[i]],predictions_list[[l]],col=cols,lwd=0.2)
      }
      
    }
    
    
    envVariableNames1 = c("Winter temperature", "Spring temperature", "Summer temperature", "Fall temperature","Temp. seasonal variation",
                          "Winter precipitation", "Spring precipitation", "Summer precipitation", "Fall precipitation",
                          "Winter humidity", "Spring humidity", "Summer humidity", "Fall humidity","Croplands", "Pastures", 
                          "Rangelands", "Urban areas", "Primary forested areas", "Primary non-forested areas", "Secondary forested areas",
                          "Secondary non-for. areas", "Total Human population","Human population density")  # plot instructions
    env_label = envVariableNames1[i]
    ri_label = paste0("RI = ",round(relativeInfluences[i,1],1), "% [",round(relativeInfluences[i,2],1), "-", 
    round(relativeInfluences[i,3],1), "]")
                       
    
    axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
    axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.2,0))
    n_col = 4; panel_id = i
    if ((panel_id-1)%% n_col==0){title(ylab="Predicted value", cex.lab=0.65, mgp=c(1.3,0,0), col.lab="gray30")}
    title(xlab = env_label, cex.lab=0.65,mgp=c(0.7,0,0), col.lab="gray30")
    title(xlab = ri_label, cex.lab=0.65, mgp=c(1.3,0,0), col.lab="gray30")
    box(lwd=0.2, col="gray30")
  }
  
  dev.off()

  
################################################################################
#                             IV. CURRENT RISK MAP (2024)
################################################################################
  
  # 1. 2024 data extraction
  #########################
month_start = (2023 - 1901)*12+1; month_end = (2024 - 1901)*12 + 12
winter_indices = sort(c(month_start + 11, month_start + 12, month_start + 13)) # December (previous year), January, February
spring_indices = sort(c(month_start + 14, month_start + 15, month_start + 16)) # March, April, May
summer_indices = sort(c(month_start + 17, month_start + 18, month_start + 19)) # June, July, August
fall_indices   = sort(c(month_start + 20, month_start + 21, month_start + 22)) # September, October, November
year_index = 2024-1900

# Extract environmental variables for 2024
covariates_2024 = list()
covariates_2024[["temp_winter"]] = mean(covariates[[1]][[winter_indices]]) - 273.15 # conversion to Celsius
covariates_2024[["temp_spring"]] = mean(covariates[[1]][[spring_indices]]) - 273.15
covariates_2024[["temp_summer"]] = mean(covariates[[1]][[summer_indices]]) - 273.15
covariates_2024[["temp_fall"]]   = mean(covariates[[1]][[fall_indices]]) - 273.15

covariates_2024[["seasonal_var"]] = covariates_2024[["temp_summer"]] - covariates_2024[["temp_winter"]]

covariates_2024[["prec_winter"]] = mean(covariates[[2]][[winter_indices]]) * 60 * 60 * 24 # conversion to vers
covariates_2024[["prec_spring"]] = mean(covariates[[2]][[spring_indices]]) * 60 * 60 * 24
covariates_2024[["prec_summer"]] = mean(covariates[[2]][[summer_indices]]) * 60 * 60 * 24
covariates_2024[["prec_fall"]]   = mean(covariates[[2]][[fall_indices]]) * 60 * 60 * 24

covariates_2024[["relh_winter"]] = mean(covariates[[3]][[winter_indices]])
covariates_2024[["relh_spring"]] = mean(covariates[[3]][[spring_indices]])
covariates_2024[["relh_summer"]] = mean(covariates[[3]][[summer_indices]])
covariates_2024[["relh_fall"]]   = mean(covariates[[3]][[fall_indices]])

covariates_2024[["crops"]] = mean(covariates[[4]][[year_index]])
covariates_2024[["pastures"]] = mean(covariates[[5]][[year_index]])
covariates_2024[["urbanAreas"]] = mean(covariates[[6]][[year_index]])
covariates_2024[["rangelands"]] = mean(covariates[[7]][[year_index]])
covariates_2024[["primForest"]] = mean(covariates[[8]][[year_index]])
covariates_2024[["primNonForest"]] = mean(covariates[[9]][[year_index]])
covariates_2024[["secForest"]] = mean(covariates[[10]][[year_index]])
covariates_2024[["secNonForest"]] = mean(covariates[[11]][[year_index]])
covariates_2024[["totalPop"]] = mean(covariates[[12]][[year_index]])

envVariableNames = brt_model_scvs[[1]]$var.names
data2024 = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(covariates_2024) + 2))
colnames(data2024) = c("NUTM", envVariableNames); data2024$NUTM = nutsM_IDs$nutM_code

for (i in seq_along(covariates_2024)) {
  if (i == length(covariates_2024)) {
    extracted_values = exactextractr::exact_extract(covariates_2024[[i]], nutsM, fun = "sum")
  } else {
    extracted_values = exact_extract(covariates_2024[[i]], nutsM, fun = "mean")
    }
  data2024[, i + 1] = extracted_values  # + 1 because first column is NUTM
  }
data2024[,24] = data2024$total_population/nutsM_sf$area
write.csv(data2024,"github_data/Environmental_data/ISIMIP3a_variables_eu/20crv3-era5_europe_data2024.csv", row.names = F)

# Current predictions
eu_2024 = as.data.frame(matrix(nrow=nrow(data2024), ncol=100))
for(j in 1:length(brt_model_scvs)){
  object = brt_model_scvs[[j]]
  n.trees = brt_model_scvs[[j]]$gbm.call$best.trees 
  type = "response"; single.tree = FALSE
  prediction = predict.gbm(object, data2024[,2:24], n.trees, type, single.tree)
  eu_2024[,j] = prediction
}
colnames(eu_2024)[1:100] = paste0("rep_", 1:100); write.csv(eu_2024,"github_data/Brt_outputs/ISIMIP3a_2024_projections/2024_Europe_model.csv", row.names = F)

# 2. 2024 risk map
##################
plotting = FALSE
if (plotting == TRUE) {
eu_2024 = read.csv("github_data/Brt_outputs/ISIMIP3a_2024_projections/2024_Europe_model.csv", header = T)
eu_2024[, 102] = rowMeans(eu_2024[, 2:101])
pdf("Figures/current_riskmap.pdf", width=6, height=6)
legend1 = raster(as.matrix(c(0,1)))
par(oma=c(0,0.1,0.5,1.5), mar=c(0,2,0,0), lwd=0.1, col="gray30")
colourScale = rev(colorRampPalette(brewer.pal(11, "RdBu"))(121)[11:111])
cols = colourScale[(((eu_2024[, 102] - 0) / (1 - 0)) * 100) + 1]
plot(contour, lwd=0.4, border="gray30", col=NA)
plot(nutsM, col = cols, border = NA, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
centroids_coords <- st_coordinates(nutsM_FR_centroids)
points(centroids_coords, col = "black", pch = 16, cex= 0.5)
plot(legend1, col= colourScale, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3,
     smallplot=c(0.93,0.96,0.15,0.85), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2,
                                                             col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(2,0.5,0),
                                                             at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3)
dev.off()
}

  
################################################################################
#                             V. PAST PREDICTIONS
################################################################################

# 1. Data extraction and historical predictions
################################################

# Past time intervals
years_start = c(1901,1925,1950,1975,2000)
years_end = c(1924,1949,1974,1999,2021)

# start past data extraction and predictions
models = c("obsclim","counterclim")
brt_model_scvs = readRDS("European_model/All_the_brt_models/brt_trained_model.rds")
envVariableNames = brt_model_scvs[[1]]$var.names

options(future.globals.maxSize = 8 * 1024^3)
plan(multisession, workers = parallel::detectCores() - 4) # parallel across models

past_predictions_list = future_lapply(seq_along(models), function(m) {
  past_predictions = list()
  log_file = "European_model/models_progress_log.txt"
  
  temperature = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_tas_1901_2021","_monmean.nc")); temperature = mask(crop(temperature, nutsM), nutsM)
  precipitation = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_pr_1901_2021","_monmean.nc")); precipitation = mask(crop(precipitation, nutsM), nutsM)
  relativehumidity = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/", models[m], "/20CRv3-ERA5/20crv3-era5_", models[m],"_hurs_1901_2021","_monmean.nc")); relativehumidity = mask(crop(relativehumidity, nutsM), nutsM)
  
  for(y in 1:length((years_start))) {
    past_predictions[[y]] = list()
    start = years_start[[y]] ; end = years_end[[y]] 
    year_index = (start-1900):(end-1900)
    
    if (y==1) {month_start = (start - 1901)*12+1; month_end = (end - 1901)*12 + 12
    }else{month_start = ((start-1) - 1901)*12+1; month_end = (end - 1901)*12 + 12}
    
    temperature_temp = temperature[[month_start:month_end]]
    precipitation_temp = precipitation[[month_start:month_end]]
    relativehumidity_temp = relativehumidity[[month_start:month_end]]
    
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
    envVariables[[22]] =  mean(covariates[[12]][[year_index]]) # total population
    
    nutsM_data_past = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(envVariableNames) + 1))
    colnames(nutsM_data_past) = c("NUTM", envVariableNames)
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
    write.csv(nutsM_data_past, paste0("github_data/Environmental_data/ISIMIP3a_variables_eu/20crv3-era5_",models[[m]],"_",years_start[[y]],"-",years_end[[y]],"_europe.csv"), row.names = F)
    
    for(j in 1:length(brt_model_scvs)) {
      object = brt_model_scvs[[j]]
      n.trees = brt_model_scvs[[j]]$gbm.call$best.trees; 
      type = "response"; single.tree = FALSE
      prediction = predict.gbm(object, nutsM_data_past[, -1], n.trees, type, single.tree)
      past_predictions[[y]][[j]] = prediction 
    }
  }
  saveRDS(past_predictions, paste0("European_model/All_the_brt_models/past_predictions_", models[m], ".rds"))
  return(past_predictions)
})

# 2. Saving past predictions
#####################################
year_intervals = c("1901-1924","1925-1949","1950-1974","1975-1999","2000-2021")
path = "github_data/Brt_outputs/ISIMIP3a_past_projections/"

for (m in 1:length(models)) {
  for (y in 1:length(year_intervals)) {
    df = as.data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = 100))
    for (j in 1:100) {
      df[,j] = past_predictions_list[[m]][[y]][[j]]
    }
    df = cbind(NUTSM = nutsM_IDs$nutM_code, df); colnames(df)[2:101] = paste0("rep_",1:100)
    write.csv(df, paste0(path,"Past_preds_",models[[m]],"_",year_intervals[[y]],"_20CRv3-ERA5.csv"), row.names = F)
  }
}


  
################################################################################
#                             VI. FUTURE PREDICTIONS
################################################################################
  
# 1. Data extraction and future projections
################################################

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
  # crop to europe
  land_covers2 = lapply(land_covers2, function(x) crop(x, nutsM, snap="out"))
  land_covers2 = lapply(land_covers2, function(x) mask(x, nutsM))
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

# data extraction
areas = brick("Rasters/Environmental_rasters/clm45_area.nc4"); areas = areas/(1000*1000) # conversion to km2 # total_population = population_density*areas
years_start = c(2025,2050,2075)
years_end = c(2049,2074,2100)
year_intervals = c("2025-2049","2050-2074","2075-2100")

MODELS_isimip3b = c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2","CanESM5-2",
                    "CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3", "MIROC6")
models_isimip3b = c("gfdl-esm4","ipsl-cm6a-lr","mpi-esm1-2-hr","mri-esm2-0","ukesm1-0-ll", "canesm5",
                    "cnrm-cm6-1","cnrm-esm2-1","ec-earth3","miroc6")
nmodels = length(models_isimip3b)
brt_model_scvs = readRDS("European_model/All_the_brt_models/brt_trained_model.rds")


options(future.globals.maxSize = 8 * 1024^3)
plan(multisession, workers = parallel::detectCores() - 4) # parallel across scenarios

future_lapply(seq_along(scenarios), function(s) {
  
  # Load population raster per scenario
  Pop_D = brick(paste0("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_2026_2100_",scenarios2[[s]], ".nc"))
  Pop_D = crop(Pop_D, nutsM, snap="out"); Pop_D = mask(Pop_D, nutsM)
  areas_eu = crop(areas, nutsM, snap="out"); areas_eu = mask(areas_eu, nutsM)
  totalPopulation = Pop_D * areas_eu
  
  scenario_predictions = list()
  log_file = paste0("European_model/model_log_", scenarios[s], ".txt")
  
  for (y in seq_along(years_start)) {
    year_start = years_start[y]; year_end = years_end[y]
    year_predictions = list()
    
    # Month indicies 
    month_start = ((year_start - 1) - 2015) * 12 + 1
    month_end = (year_end - 2015) * 12 + 12
    
    for (m in seq_along(models_isimip3b)) {
      
      output_file <- paste0("github_data/Environmental_data/ISIMIP3b_variables_eu/data_",scenarios[[s]], "_", year_intervals[[y]], "_", MODELS_isimip3b[[m]], ".csv")
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
        
        # Crop/mask to Europe
        temperature = mask(crop(temperature, nutsM, snap="out"), nutsM)
        precipitation = mask(crop(precipitation, nutsM, snap="out"), nutsM)
        relative_humidity = mask(crop(relative_humidity, nutsM, snap="out"), nutsM)
        
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
          tempm_summer - tempm_winter,
          precp_winter, precp_spring, precp_summer, precp_fall,
          relh_winter, relh_spring, relh_summer, relh_fall,
          croplands_temp, pastures_temp, rangelands_temp, urbanAreas_temp,
          primaryForest_temp, primaryNonForest_temp, secondaryForest_temp, secondaryNonForest_temp,
          population
        )
        
        # Prepare data.frame
        nutsM_future = data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = length(envVariableNames)+1))
        colnames(nutsM_future) = c("NUTM", envVariableNames)
        nutsM_future$NUTM = nutsM_IDs$nutM_code
        
        # Extract values
        for (v in seq_along(envVariables)) {
          fun_type = if (v == length(envVariables)) "sum" else "mean"
          nutsM_future[, v+1] <- exactextractr::exact_extract(envVariables[[v]], nutsM, fun = fun_type)
        }
        nutsM_future$total_population_density = nutsM_future$total_population / nutsM_sf$area
        
        # Save CSV
        write.csv(nutsM_future, output_file, row.names = FALSE)
      } # end else (CSV does not exist)
      
      # Predict for all BRT models
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
  saveRDS(scenario_predictions, file = paste0("European_model/All_the_brt_models/future_predictions_", scenarios[s], ".rds"))
  return(paste("Scenario", scenarios[s], "done"))
})


# 2. Saving future projections
##############################

ssp1 = readRDS("European_model/All_the_brt_models/future_predictions_ssp126.rds")
ssp3 = readRDS("European_model/All_the_brt_models/future_predictions_ssp370.rds")
ssp5 = readRDS("European_model/All_the_brt_models/future_predictions_ssp585.rds")
predictions_future = list (ssp1,ssp3,ssp5)

year_intervals = c("2025-2049","2050-2074","2075-2100")
path = "github_data/Brt_outputs/ISIMIP3b_future_projections/"
for (s in 1:length(scenarios)) {
  for (y in 1:length(year_intervals)) {
    for (m in 1:length(models_isimip3b)){
      df = as.data.frame(matrix(nrow = nrow(nutsM_IDs), ncol = 100))
      for (j in 1:100) {
      df[,j] = predictions_future[[s]][[y]][[m]][[j]]
    }
    df = cbind(NUTSM = nutsM_IDs$nutM_code, df); colnames(df)[2:101] = paste0("rep_",1:100)
    write.csv(df, paste0(path,"Future_preds_",scenarios[[s]],"_",year_intervals[[y]],"_",MODELS_isimip3b[[m]],".csv"), row.names = F)
    }
  }
}


################################################################################
#               VII. VISUILASATION SCRIPT (code by Dr Simon Dellicour)
################################################################################

date = "10102025"

# 1. Preparing the optimised NUTS3 shapefile

nutsM = crop(shapefile("Europe_NUTS3_shapefile/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # to load the optimised NUTS3 shapefile
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ")) # to only keep EU areas
nutsM = subset(nutsM, !LocationCo%in%c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533","ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65"))
# to remove the small islands where land-use data is not available 
coasts = unionSpatialPolygons(nutsM, rep(1,length(nutsM))) # if "rgeos" is available
coasts = st_union(st_as_sf(nutsM)) # if the "rgeos" is not available for the R version

# 2. Generating the figure with the past projections

brts = list(); scenarios = c("obsclim","counterclim"); periods = c("1901-1924","1925-1949","1950-1974","1975-1999","2000-2021")
for (i in 1:length(scenarios))
{
  buffer = list()
  for (j in 1:length(periods))
  {
    buffer[[j]] = read.csv(paste0("BRT_projection_outputs/All_past_projections/Past_preds_",scenarios[i],"_",periods[j],"_20CRv3-ERA5.csv"), head=F)
    buffer[[j]] = buffer[[j]][,2:dim(buffer[[j]])[2]]
  }
  brts[[i]] = buffer
}

plottingStDevs = FALSE; suffix = "A" # splottingStDevs =TRUE; suffix = "B"
scenario_names = c("Historical","Counterfactual") # plottingStDevs = TRUE
pdf(paste0("Figure_1_",date,"_",suffix,"_NEW.pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # "D9" = 85%
colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") # "00" = 0%
colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5"
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18])
colourScale2[9] = paste0(colourScale2[9],"00")
colourScale3 = colourScale2; colourScale3[9] = "#E5E5E5"
for (i in 1:length(scenarios))
{
  for (j in 1:length(periods))
  {
    brt = as.matrix(brts[[i]][[j]])
    means = rep(NA, dim(brt)[1]); stdvs = rep(NA, dim(brt)[1])
    for (k in 1:dim(brt)[1])
    {
      means[k] = mean(brt[k,]); stdvs[k] = sd(brt[k,])
    }
    cols = colourScale1[(((means-0)/(1-0))*100)+1]
    if (plottingStDevs) cols = colourScale1[(((stdvs-0)/(1-0))*100)+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA)
    plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
    # if ((i == 1)&(j == 1)) mtext(expression(bold(A)), side=3, line=-1.45, at=-8.5, cex=0.70, col="gray30")
    # if ((i == 2)&(j == 1)) mtext(expression(bold(B)), side=3, line=-1.45, at=-8.5, cex=0.70, col="gray30")
    mtext(scenario_names[i], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")
    if ((i == length(scenarios))&(j == length(periods)))
    {
      plot(rast, legend.only=T, add=T, col=colourScale1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
           axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
    }								
  }
}
vS = c()
for (j in 1:length(periods))
{
  brt1 = as.matrix(brts[[1]][[j]]); brt2 = as.matrix(brts[[2]][[j]])
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
  brt1 = as.matrix(brts[[1]][[j]]); brt2 = as.matrix(brts[[2]][[j]])
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

# 3. Generating the figure with the future projections

present = read.csv(paste0("BRT_projection_outputs/Current_projections/2024_Europe_model.csv"), head=F)
brts = list(); pops = list(); scenarios = c("ssp126","ssp370","ssp585"); periods = c("2025-2049","2050-2074","2075-2100")
models = c("CanESM5-2","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2")
for (i in 1:length(scenarios))
{
  buffer1 = list()
  for (j in 1:length(periods))
  {
    buffer2 = list()
    for (k in 1:length(models))
    {
      buffer2[[k]] = read.csv(paste0("BRT_projection_outputs/Future_projections/Future_preds_",scenarios[i],"_",periods[j],"_",models[k],".csv"), head=F)
      buffer2[[k]] = buffer2[[k]][,2:dim(buffer2[[k]])[2]]
    }
    buffer1[[j]] = buffer2
  }
  brts[[i]] = buffer1
}
for (i in 1:length(scenarios))
{
  pops[[i]] = read.csv(paste0("Population_density_data/NUTS3_population_",scenarios[i],".csv"), head=T)
}

plottingDifference = TRUE; plottingStDevs = FALSE; suffix = "A"
# plottingDifference = FALSE; plottingStDevs = FALSE; suffix = "B"
scenario_names = c("SSP1-2.6","SSP3-7.0","SSP5-8.5"); # plottingStDevs = TRUE
pdf(paste0("Figure_2_",date,"_",suffix,"_NEW.pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")		
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # "D9" = 85%
colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") # "00" = 0%
colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5"
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18])
colourScale2[9] = paste0(colourScale2[9],"00")
colourScale3 = colourScale2; colourScale3[9] = "#E5E5E5"
I1 = 1; I2 = length(scenarios); c = 0
if (plottingDifference) { I1 = 2; I2 = 2 }
for (i in I1:I2)
{
  c = c+1
  plot.new()
  if (c == 1)
  {
    brt = as.matrix(present[,2:dim(present)[2]]); brt1 = brt
    means = rep(NA, dim(brt)[1]); stdvs = rep(NA, dim(brt)[1])
    for (j in 1:dim(brt)[1])
    {
      means[j] = mean(brt[j,]); stdvs[j] = sd(brt[j,])
    }
    cols = colourScale1[(((means-0)/(1-0))*100)+1]
    if (plottingStDevs) cols = colourScale1[(((stdvs-0)/(1-0))*100)+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
    mtext("Present time", side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext("(2024, t0)", side=3, line=-2.5, at=1, cex=0.50, col="gray30"); rast = raster(as.matrix(c(0,1)))
  }	else	{
    plot.new()
  }
  for (j in 1:length(periods))
  {
    means1 = rep(NA, length(means)); stdvs1 = rep(NA, length(means))
    for (k in 1:length(means))
    {
      means2 = rep(NA, length(models)); stdvs2 = rep(NA, length(models))
      for (l in 1:length(models))
      {
        brt = as.matrix(brts[[i]][[j]][[l]])
        means2[l] = mean(brt[k,]); stdvs2[l] = sd(brt[k,])
      }
      means1[k] = mean(means2); stdvs1[k] = mean(stdvs2)
    }
    cols = colourScale1[(((means1-0)/(1-0))*100)+1]
    if (plottingStDevs) cols = colourScale1[(((stdvs1-0)/(1-0))*100)+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA)
    plot(nutsM, col=cols, border=NA, lwd=0.1, add=T); rast = raster(as.matrix(c(0,1)))
    mtext(scenario_names[i], side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext(gsub("_","-",periods[j]), side=3, line=-2.5, at=1, cex=0.50, col="gray30")						
    if ((i == length(scenarios))&(j == length(periods)))
    {
      plot(rast, legend.only=T, add=T, col=colourScale1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.060,0.080,0.75,0.96), adj=3,
           axis.args=list(cex.axis=0.65, lwd=0, col="gray30", lwd.tick=0.2, col.tick="gray30", tck=-1.2, col.axis="gray30", line=0, mgp=c(0,0.45,0)), alpha=1, side=3)
    }
  }
}
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
        buffer = rep(NA,length(models) )
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
}
dev.off()

thresholds_SI = read.csv("BRT_projection_outputs/SIppc_threshold.csv", head=T)
pop_exposures_1 = matrix(nrow=100, ncol=10); temp = as.matrix(present[,2:dim(present)[2]])
for (i in 1:dim(pop_exposures_1)[1])
{
  indexes = which(temp[,i]>=as.numeric(thresholds_SI[i,"threshold"]))
  pop_exposures_1[i,1] = sum(pops[[1]][indexes,"X2024"])
}
pop_exposures_2 = pop_exposures_1 # for when we keep the population at its 2024 level
for (i in 1:length(periods))
{
  for (j in 1:length(scenarios))
  {	
    for (k in 1:dim(pop_exposures)[1])
    {
      temp1 = matrix(nrow=dim(nutsM@data)[1], ncol=10)
      temp2 = matrix(nrow=dim(nutsM@data)[1], ncol=1)
      for (l in 1:length(models)) temp1[,l] = brts[[j]][[i]][[l]][,k]
      for (l in 1:dim(temp2)[1]) temp2[l,1] = mean(temp1[l,])
      indexes = which(temp2>=thresholds_SI[k,"threshold"])
      colname = paste0("X",gsub("-",".",periods[i]))
      pop_exposures_1[k,1+((i-1)*length(scenarios))+j] = sum(pops[[j]][indexes,colname])
      pop_exposures_2[k,1+((i-1)*length(scenarios))+j] = sum(pops[[j]][indexes,"X2024"])
    }
  }
}
boxplot(pop_exposures_1)
boxplot(cbind(pop_exposures_1,pop_exposures_2))

pop_exposure = pop_exposures_1 # pop_exposure = pop_exposures_2
for (i in 1:dim(pop_exposures)[2])
{
  print(shapiro.test(pop_exposures[,i])$p.value)
}
pValues = matrix(0, nrow=dim(pop_exposures)[2]-1, ncol=4) # Benjamini-Hochberg correction
colnames(pValues) = c("p-value","rank","new_threshold_B-H","significant_after_B-H")
for (i in 2:dim(pop_exposures)[2])
{
  print(round(wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value,5))
  pValues[i-1,1] = wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value
}
for (i in 2:dim(pop_exposures)[2])
{
  pValues[i-1,2] = which(pValues[order(pValues[,1]),1]==pValues[i-1,1])
  pValues[i-1,3] = (pValues[i-1,2]/dim(pValues)[1])*0.05
  if (pValues[i-1,1] < pValues[i-1,3]) pValues[i-1,4] = 1
  pValues[i-1,1] = round(pValues[i-1,1],3)
  pValues[i-1,3] = round(pValues[i-1,3],3)
}


# 4. Estimating the variation at the continent scale

indices = 1:dim(nutsM@data)[1] # when considering all the countries
indices = which(nutsM@data[,"CountryNam"]=="Belgium") # estimating the variation for Belgium
indices = which(nutsM@data[,"CountryNam"]%in%c("Belgium","Netherlands")) # estimating the variation for Belgium and The Netherlands

buffer_ssp3 = list(); buffer_ssp5 = list(); brts_ssp3 = list(); brts_ssp5 = list(); pops_ssp3 = list(); pops_ssp5 = list()
models = c("CanESM5-2","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2")
brts_ssp3[[1]] = read.csv(paste0("BRT_projection_outputs/Current_projections/2024_Europe_model.csv"), head=F)
brts_ssp3[[1]] = brts_ssp3[[1]][indices,2:dim(brts_ssp3[[1]])[2]]; brts_ssp5[[1]] = brts_ssp3[[1]]
for (i in 1:length(models))
{
  buffer_ssp3[[i]] = read.csv(paste0("BRT_projection_outputs/Future_projections/Future_preds_ssp370_2075-2100_",models[i],".csv"), head=F)
  buffer_ssp5[[i]] = read.csv(paste0("BRT_projection_outputs/Future_projections/Future_preds_ssp585_2075-2100_",models[i],".csv"), head=F)
  buffer_ssp3[[i]] = buffer_ssp3[[i]][indices,2:dim(buffer_ssp3[[i]])[2]]; buffer_ssp5[[i]] = buffer_ssp5[[i]][indices,2:dim(buffer_ssp5[[i]])[2]]
}
brts_ssp3[[2]] = buffer_ssp3; brts_ssp5[[2]] = buffer_ssp5; difference_ES_ssp3 = list(); difference_ES_ssp5 = list()
for (i in 1:length(models))
{
  buffer1_ssp3 = list(); buffer1_ssp5 = list()
  for (j in 1:dim(brts_ssp3[[1]])[2])
  {
    buffer2_ssp3 = list()
    for (k in 1:dim(brts_ssp3[[1]])[1])
    {
      buffer2_ssp3[[k]] = ((brts_ssp3[[2]][[i]][k,j]-brts_ssp3[[1]][k,j])/brts_ssp3[[1]][k,j])*100
    }
    buffer1_ssp3[[j]] = buffer2_ssp3
  }
  difference_ES_ssp3[[i]] = buffer1_ssp3
  for (j in 1:dim(brts_ssp5[[1]])[2])
  {
    buffer2_ssp5 = list()
    for (k in 1:dim(brts_ssp5[[1]])[1])
    {
      buffer2_ssp5[[k]] = ((brts_ssp5[[2]][[i]][k,j]-brts_ssp5[[1]][k,j])/brts_ssp5[[1]][k,j])*100
    }
    buffer1_ssp5[[j]] = buffer2_ssp5
  }
  difference_ES_ssp5[[i]] = buffer1_ssp5
}
dS2_ssp3 = rep(NA, dim(brts_ssp3[[1]])[2]) # first averaging across models and the replicates
dS2_ssp5 = rep(NA, dim(brts_ssp5[[1]])[2]) # first averaging across models and the replicates
for (i in 1:dim(brts_ssp3[[1]])[2])
{
  dS_ssp3 = rep(NA, length(models))
  for (j in 1:length(models))
  {
    dS_ssp3[j] = mean(unlist(difference_ES_ssp3[[j]][[i]]))
  }
  dS2_ssp3[i] = mean(dS_ssp3)
}
for (i in 1:dim(brts_ssp5[[1]])[2])
{
  dS_ssp5 = rep(NA, length(models))
  for (j in 1:length(models))
  {
    dS_ssp5[j] = mean(unlist(difference_ES_ssp5[[j]][[i]]))
  }
  dS2_ssp5[i] = mean(dS_ssp5)
}
m = round(mean(dS2_ssp3)); ci = round(t.test(dS2_ssp3)$conf.int)
cat("dS2 (SSP3) = ",m,"%, 95% CI = [",ci[1],", ",ci[2],"]",sep="")
# for all countries: 68%, 95% CI = [58, 78]
# for Belgium only: 170%, 95% CI = [145, 194]
# --> the ecological suitability could more than double before the end of this century
# for Belgium and the Netherlands: 123%, 95% CI = [103, 144]
m = round(mean(dS2_ssp5)); ci = round(t.test(dS2_ssp5)$conf.int)
cat("dS2 (SSP5) = ",m,"%, 95% CI = [",ci[1],", ",ci[2],"]",sep="")
# for all countries: 107%, 95% CI = [89, 125]
# for Belgium only: 145%, 95% CI = [120, 169]
# for Belgium and the Netherlands: 127%, 95% CI = [104, 150]

thresholds_SI = read.csv("BRT_projection_outputs/SIppc_threshold.csv", head=T)
pops_ssp3 = read.csv("Population_density_data/NUTS3_population_ssp370.csv")[indices,]
pops_ssp5 = read.csv("Population_density_data/NUTS3_population_ssp585.csv")[indices,]
pop_exposures_ssp3 = matrix(nrow=100, ncol=2); pop_exposures_ssp5 = matrix(nrow=100, ncol=2)
for (i in 1:dim(pop_exposures_ssp3)[1])
{
  indexes = which(brts_ssp3[[1]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
  pop_exposures_ssp3[i,1] = sum(pops_ssp3[indexes,"X2024"]); buffer = rep(NA, length(brts_ssp3[[2]]))
  for (j in 1:length(brts_ssp3[[2]]))
  {
    indexes = which(brts_ssp3[[2]][[j]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
    buffer[j] = sum(pops_ssp3[indexes,"X2075.2100"])
  }
  pop_exposures_ssp3[i,2] = mean(buffer)
}
for (i in 1:dim(pop_exposures_ssp5)[1])
{
  indexes = which(brts_ssp5[[1]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
  pop_exposures_ssp5[i,1] = sum(pops_ssp5[indexes,"X2024"]); buffer = rep(NA, length(brts_ssp5[[2]]))
  for (j in 1:length(brts_ssp5[[2]]))
  {
    indexes = which(brts_ssp5[[2]][[j]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
    buffer[j] = sum(pops_ssp5[indexes,"X2075.2100"])
  }
  pop_exposures_ssp5[i,2] = mean(buffer)
}
dS_ssp3 = pop_exposures_ssp3[,2]-pop_exposures_ssp3[,1]; dS_ssp5 = pop_exposures_ssp5[,2]-pop_exposures_ssp5[,1]
m = round(mean(dS_ssp3)); qS = round(quantile(dS_ssp3,c(0.025,0.975)),1); ci = round(t.test(dS_ssp3)$conf.int)
cat("dS (SSP3) = ",m,", 95% CI = [",ci[1],", ",ci[2],"]",sep="") # 56275182, 95% CI = [48611073, 63939291]
m = round(mean(pop_exposures_ssp5[,2])); ci = round(t.test(pop_exposures_ssp5[,2])$conf.int)
cat("dS (SSP5) = ",m,", 95% CI = [",ci[1],", ",ci[2],"]",sep="") # 5906807, 95% CI = [4970787, 6842828]

m = round(mean(dS_ssp5)); qS = round(quantile(dS_ssp5,c(0.025,0.975)),1); ci = round(t.test(dS_ssp5)$conf.int)
cat("dS (SSP5) = ",m,", 95% CI = [",ci[1],", ",ci[2],"]",sep="") # 5906807, 95% CI = [4970787, 6842828]

# 5. Estimating the variation for Belgium (Le Soir)

indices = which(nutsM@data[,"CountryNam"]=="Belgium"); buffer = list(); brts = list(); pops = list()
models = c("CanESM5-2","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL-2")
brts[[1]] = read.csv(paste0("BRT_projection_outputs/Current_projections/2024_Europe_model.csv"), head=F); brts[[1]] = brts[[1]][indices,2:dim(brts[[1]])[2]]
for (i in 1:length(models))
{
  buffer[[i]] = read.csv(paste0("BRT_projection_outputs/Future_projections/Future_preds_ssp370_2075-2100_",models[i],".csv"), head=F)
  buffer[[i]] = buffer[[i]][indices,2:dim(buffer[[i]])[2]]
}
brts[[2]] = buffer; difference_ES = list()
for (i in 1:length(models))
{
  buffer1 = list()
  for (j in 1:dim(brts[[1]])[2])
  {
    buffer2 = list()
    for (k in 1:dim(brts[[1]])[1])
    {
      buffer2[[k]] = ((brts[[2]][[i]][k,j]-brts[[1]][k,j])/brts[[1]][k,j])*100
    }
    buffer1[[j]] = buffer2
  }
  difference_ES[[i]] = buffer1
}
dS1 = rep(NA, length(models)) # first averaging across replicates and the models
dS2 = rep(NA, dim(brts[[1]])[2]) # first averaging across models and the replicates
for (i in 1:length(models))
{
  dS = rep(NA, dim(brts[[1]])[2])
  for (j in 1:dim(brts[[1]])[2])
  {
    dS[j] = mean(unlist(difference_ES[[i]][[j]]))
  }
  dS1[i] = mean(dS)
}
for (i in 1:dim(brts[[1]])[2])
{
  dS = rep(NA, length(models))
  for (j in 1:length(models))
  {
    dS[j] = mean(unlist(difference_ES[[j]][[i]]))
  }
  dS2[i] = mean(dS)
}
m = round(mean(dS1),1); qS = round(quantile(dS1,c(0.025,0.975)),1)
cat("dS1 = ",m,"%, 95% CI = [",qS[1],", ",qS[2],"]",sep="") # dS1 = 169.6%, 95% CI = [109.1, 212.0]
m = round(mean(dS2),1); qS = round(quantile(dS2,c(0.025,0.975)),1)
cat("dS2 = ",m,"%, 95% CI = [",qS[1],", ",qS[2],"]",sep="") # dS2 = 169.6%, 95% CI = [29.5, 527.8]
# --> the ecological suitability could more than double before the end of this century
thresholds_SI = read.csv("BRT_projection_outputs/SIppc_threshold.csv", head=T)
pops = read.csv("Population_density_data/NUTS3_population_ssp370.csv")[indices,]
pop_exposures = matrix(nrow=100, ncol=2)
for (i in 1:dim(pop_exposures)[1])
{
  eco_suitability[i,1] = mean(brts[[1]][,i])
  indexes = which(brts[[1]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
  pop_exposures[i,1] = sum(pops[indexes,"X2024"]); buffer = rep(NA, length(brts[[2]]))
  for (j in 1:length(brts[[2]]))
  {
    indexes = which(brts[[2]][[j]][,i]>=as.numeric(thresholds_SI[i,"threshold"]))
    buffer[j] = sum(pops[indexes,"X2075.2100"])
  }
  pop_exposures[i,2] = mean(buffer)
}
dS = pop_exposures[,2]-pop_exposures[,1]; qS = quantile(dS, c(0.025,0.975))
mean(dS) # in 2071-2100, ~5.9 million people could leave in areas ecologically at risk of human exposure

