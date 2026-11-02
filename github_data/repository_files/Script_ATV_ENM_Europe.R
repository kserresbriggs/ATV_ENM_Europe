# ------------------------------------------------------------------------------
# Script name: Script_ATV_ENM_Europe
#
# Purpose: Increasing numbers of people living in European areas ecologically 
#          suitable for local circulation of Aedes mosquito-borne viruses
#
# Author(s): Kyla Serres
# ------------------------------------------------------------------------------

# Load necessary libraries -----------------------------------------------------
library(blockCV)
library(dismo)
library(exactextractr)
library(gbm)
library(raster)
library(RColorBrewer)
library(sf)
library(dplyr)
library(gmodels)
library(future.apply)
library(ade4)


# ------------------------------------------------------------------------------
#                                   SUMMARY
# ------------------------------------------------------------------------------

    ## I. DATA PREPARATION
        # 1. Load optimized European administrative area shapefile (VectorNet)  
        # 2. loading Aedes-borne (DENV, CHIKV, ZIKV) viral occurrence data

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
          # 4.3 Caracterising the environmental space associated with Ps and PAs
          # 4.4 Assessing model predictive performance from the cross-validation
          # 4.5 Assessing model predictive performances with a temporal validation
          # 4.6 Comparing model predictive performance associated with different PA sampling approaches
          # 4.7 Caracterising the composition (presences, PAs) of the spatial folds
          # 4.8 Retrieving the relative importance of each environmental variable
          # 4.9 Generating and plotting the response curves for each variable

    ## III. PAST PROJECTIONS
        # 1. Conducting past projections and generating the past projections figure


    ## IV. FUTURE PREDICTIONS
        # 1. Conducting future projections and generating the future projections figure

    ## V. Esrimating the total number of people living in area ecologically suitable
        # 1. Generating a figure gathering different set of boxplots
        # 2. Generating a global figure showing the evolution through time

    ## VI. Estimating the ecological suitability variations at the continent scale
setwd("~/Dropbox/PhD/Research/ATV_ENM_Europe/github_data/repository_files/")

# I. DATA PREPARATION ----------------------------------------------------------

  # 1. Load optimized European administrative area shapefile (VectorNet)  
  # --------------------------------------------------------------------

nutsM = crop(shapefile("europe_shapefile/VectornetMAPforMOODjan21.shp"), extent(-10,33,34.5,72)) # modified NUT3 shapefile
# Only keep EU areas
nutsM = subset(nutsM, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
# remove the small islands where land-use data is not available
nutsM = subset(nutsM, !LocationCo%in%c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533","ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65")) 

nutsM_sf = st_as_sf(nutsM)
coasts = st_union(nutsM_sf)
nutsM_sf$area = nutsM_sf %>% st_area %>% as.numeric/(1000*1000) 

# Subset polygon IDs
nutsM_IDs = as.data.frame(matrix(nrow = dim(nutsM@data)[1], ncol = 3))
colnames(nutsM_IDs) = c("country","region_name", "nutM_code"); nutsM_IDs$nutM_code = nutsM@data[,"LocationCo"]
nutsM_IDs$country = nutsM@data[,"CountryNam"]; nutsM_IDs$region_name = nutsM@data[,"LocationNa"]


  # 2. loading Aedes-borne (DENV, CHIKV, ZIKV) viral occurrence data 
  # ------------------------------------------------------------------

# Note: Occurrence data used in this study has has been obtained from:
#   - ECDC surveillance platform: https://atlas.ecdc.europa.eu/public/
#   - Santé Publique France: https://odisse.santepubliquefrance.fr/page/accueil/
#   - Literature: https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(25)00023-7/fulltext
#   - Padi-Web: https://padi-web.cirad.fr/en/

# all occurrence data used in this study:
Aedes_presences = read.csv("occurrence_data/Aedes_viruses_European_presences.csv", header = T)

aggregatePresences = FALSE
if (aggregatePresences==TRUE){
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
write.csv(ae.presences_agg,"occurrence_data/European_brt_training_presences.csv", row.names = FALSE)
} else{
  ae.presences_agg = read.csv("occurrence_data/European_brt_training_presences.csv", header = T)
}


# II. BRT MODEL TRAINING AND EVALUATION ----------------------------------------

  # 1. Pseudo-absence sampling 
  # ---------------------------

# Creating the pseudo-absence shapefile: for pseudo-absence (PA) sampling, 
# remove polygons adjacent to presence sites

# Identify adjacent polygons for each presence location
ae.pres_geoms = which(nutsM_sf$LocationCo %in% ae.presences_agg$nutM_code)
ae.pres_geoms = nutsM_sf[ae.pres_geoms, ] 
ae.adjacent = st_intersects(ae.pres_geoms$geometry, nutsM_sf$geometry)
ae.adj = unlist(ae.adjacent)
ae.adj_geoms = nutsM_sf[ae.adj, ] # presence and neighboring polygons


# Remove adjacent polygons from PA sampling shapefile
absences = setdiff(nutsM_IDs$nutM_code, ae.presences_agg$nutM_code) #ids of polygons without recorded presence
Aedes_abs = data.frame(nutM_code = absences, year = NA)
remove_rows = which(Aedes_abs$nutM_code %in% ae.adj_geoms$LocationCo) 
Aedes_abs = Aedes_abs[(-remove_rows), ] #absence polygons that are not adjacent to presence polygons


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
rast_surv = raster("brt_training_data/pseudo_abcence_strategies/Rasters/Lim2025_Surveillance_map_wmean.tif")

# extract surveillance capability values per Nuts area
surv_nutsM = exact_extract(rast_surv, nutsM, fun="mean")
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

# Note: The raster files are not included in the GitHub repo because of their size.
# The extracted raster values per NUTS3 region, for all climate scenarios and
# time periods, are available as CSV files in the `environmental_data/` folder.
# To download the original raster files please visit: https://data.isimip.org/
directory1 = "../../Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/"
directory2 = "../../Rasters/Environmental_rasters/Land-use/"
directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"

temperature = brick(paste0(directory1,"20crv3-era5_obsclim_tas_1901_2021_monmean.nc")); temperature2 = brick(paste0(directory1,"20crv3-era5_obsclim_tas_2022_2024_monmean.nc"))
temperature = addLayer(temperature,temperature2)

precipitation = brick(paste0(directory1,"20crv3-era5_obsclim_pr_1901_2021_monmean.nc")); precipitation2 = brick(paste0(directory1,"20crv3-era5_obsclim_pr_2022_2024_monmean.nc"))
precipitation = addLayer(precipitation, precipitation2)

relativehumidity = brick(paste0(directory1,"20crv3-era5_obsclim_hurs_1901_2021_monmean.nc")); relativehumidity2 = brick(paste0(directory1,"20crv3-era5_obsclim_hurs_2022_2024_monmean.nc"))
relativehumidity= addLayer(relativehumidity,relativehumidity2) 

croplands = brick(paste0(directory2,"croplands_LUH2-GCB2024_1901_2024.nc"))
pastures = brick(paste0(directory2,"pastures_LUH2-GCB2024_1901_2024.nc"))
rangelands = brick(paste0(directory2,"rangeland_LUH2-GCB2024_1901_2024.nc"))
urbanAreas = brick(paste0(directory2,"urbanAreas_LUH2-GCB2024_1901_2024.nc"))
primaryForest = brick(paste0(directory2,"primaryForest_LUH2-GCB2024_1901_2024.nc"))
primaryNonForest = brick(paste0(directory2,"primaryNonF_LUH2-GCB2024_1901_2024.nc"))
secondaryForest = brick(paste0(directory2,"secondaryForest_LUH2-GCB2024_1901_2024.nc"))
secondaryNonForest = brick(paste0(directory2,"secondaryNonF_LUH2-GCB2024_1901_2024.nc"))

population_density = brick(paste0(directory3,"Population_density_05deg_1900_2025_historical.nc"))
population_density = population_density[[2:125]] # 1901 to 2024

# cell grid areas 
areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")
areas = areas/(1000*1000) # conversion to km2

# convert population density to total population count 
total_population = population_density * areas

# add all variables
covariates = list(temperature, precipitation, relativehumidity, croplands, 
                  pastures, rangelands, urbanAreas, primaryForest, primaryNonForest, 
                  secondaryForest, secondaryNonForest, total_population)

names(covariates) = c("temperature", "precipitation", "relative_humidity", "croplands", 
                    "pastures", "rangelands", "urbanAreas", "primaryForest", "primaryNonForest", 
                    "secondaryForest", "secondaryNonForest", "total_population")

pres_codes = unique(ae.presences_agg$nutM_code)
pres_years = c(ae.presences_agg$year, ae.presences_agg$year_start, ae.presences_agg$year_end)
data_adj = data.frame(nutM_code = setdiff(unique(ae.adj_geoms$LocationCo), pres_codes), year = NA_integer_,
                      year_start = min(pres_years, na.rm = TRUE), year_end = max(pres_years, na.rm = TRUE),
                      response = 0.5)
adj_file = "brt_training_data/Dataframe_adj.csv"
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
    data = rbind(ae.presences_agg, absences) # bind presences and pseudo-absences
    if (n == 1) data = rbind(data, data_adj[, colnames(data)]) # adjacent areas 


    ## data extraction ##
    names_env = c("temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                 "seasonal_variation","precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                 "croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonForest","secondaryForest","secondaryNonForest", 
                 "total_population")
    
    data_for_brt = as.data.frame(matrix(NA,nrow=nrow(data), ncol = length(names_env) + 6))
    data_for_brt[,c(1:5)] = data[,c("nutM_code","year","year_start","year_end", "response")]
    colnames(data_for_brt) = c("nutM_code","year","year_start","year_end", "response", names_env, "population_density")
    
    geometry_matches = nutsM_sf$geometry[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
    area_matches = nutsM_sf$area[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
    nutsM_train = st_sf(data_for_brt, geometry = geometry_matches); nutsM_train$area = area_matches

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
    data_for_brt$population_density = data_for_brt$total_population/nutsM_train$area # population density by NUTS3
    if (n == 1) {
      write.csv(data_for_brt[data_for_brt$response == 0.5, ], adj_file, row.names = F)
      data_for_brt = data_for_brt[data_for_brt$response != 0.5, ]
    }
    data_for_brt$year = as.integer(unlist(data_for_brt$year))
    data_train[[n]] = data_for_brt
    write.csv(data_for_brt,paste0("brt_training_data/20crv3-era5_env_data_replicate",n,".csv"), row.names = F)  
    cat("Data extraction run", n, "of", nruns, "completed.\n")
  } 
  plan(sequential)  
} else {
  prefix = "20crv3-era5_env_data_replicate"   
  n_rep  = 100                                
  files = file.path("brt_training_data", paste0(prefix, seq_len(n_rep), ".csv"))
  stopifnot(all(file.exists(files)))          
  data_train = lapply(files,read.csv)
}

  # 3. Investigating the environmental variables to include in the model (shadow variables approach)
  # ------------------------------------------------------------------------------------------------
nutsM_centroids = coordinates(nutsM); dataframes = list()
EnvNames = names(data_train[[1]][,6:29])

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

rP_mean = matrix(nrow=length(EnvNames), ncol=length(EnvNames)); colnames(rP_mean) = EnvNames
rP_lowerCI = matrix(nrow=length(EnvNames), ncol=length(EnvNames)); colnames(rP_lowerCI) = EnvNames
rP_upperCI = matrix(nrow=length(EnvNames), ncol=length(EnvNames)); colnames(rP_upperCI) = EnvNames
for (i in 1:length(EnvNames))
{
  for (j in 1:length(EnvNames))
  {
    if (i != j)
    {
      vS = rep(NA, length(dataframes))
      for (k in 1:length(dataframes))
      {
        vS[k] = cor(dataframes[[k]][,EnvNames[i]], dataframes[[k]][,EnvNames[j]], method="pearson")
      }
      rP_mean[i,j] = mean(vS); rP_lowerCI[i,j] = t.test(vS)$conf.int[1]; rP_upperCI[i,j] = t.test(vS)$conf.int[2]
    }
  }
}
for (i in 2:length(EnvNames))
{
  for (j in 1:(i-1))
  {
    to_report = FALSE
    if ((rP_mean[i,j] < 0)&(rP_upperCI[i,j] <= -0.7)) to_report = TRUE
    if ((rP_mean[i,j] > 0)&(rP_lowerCI[i,j] >= 0.7)) to_report = TRUE
    if (to_report)
    {
      cat(EnvNames[i]," - ",EnvNames[j],": rP = ",round(rP_mean[i,j],2),", 95% CI = [",round(rP_lowerCI[i,j],2),", ",round(rP_upperCI[i,j],2),"]\n",sep="")
    }
  }
}
# temperature_spring - temperature_winter: rP = 0.81, 95% CI = [0.81, 0.82]
# temperature_summer - temperature_spring: rP = 0.85, 95% CI = [0.85, 0.85]
# temperature_fall - temperature_winter: rP = 0.86, 95% CI = [0.86, 0.87]
# temperature_fall - temperature_spring: rP = 0.89, 95% CI = [0.88, 0.89]
# temperature_fall - temperature_summer: rP = 0.81, 95% CI = [0.8, 0.81]
# relative_humidity_summer - temperature_summer: rP = -0.71, 95% CI = [-0.72, -0.7]
# relative_humidity_summer - relative_humidity_spring: rP = 0.73, 95% CI = [0.72, 0.74]
# relative_humidity_fall - temperature_summer: rP = -0.72, 95% CI = [-0.72, -0.71]
# relative_humidity_fall - relative_humidity_summer: rP = 0.72, 95% CI = [0.71, 0.72]
# population_density - urbanAreas: rP = 0.96, 95% CI = [0.96, 0.97]

    # 3.2 Investagating spatial autocorrelation and defining block size for spatial cross-validation
nberOfReplicates = 100
# nutsM polygon centroids
centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  

plottingCorrelogram = TRUE
if (plottingCorrelogram)
{
  dev.new(width=8, height=2.2); par(mar=c(2.2,1.5,1.0,1.2), oma=c(0.5,0.5,0,0), lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
  correlograms = list(); values = rep(NA, nberOfReplicates)
  for (i in 1:nberOfReplicates)
  {
    data = data_train[[i]]
    matching_centroids = as.data.frame(matrix(nrow = nrow(data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
    matching_centroids$NUTM = data$nutM_code 
    matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
    
    correlograms[[i]] = ncf::correlog(matching_centroids[,2], matching_centroids[,3], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)
  }
  plot(correlograms[[1]]$mean.of.class, correlograms[[1]]$correlation, ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-1,1))
  for (i in 1:nberOfReplicates)
  {
    values[i] = correlograms[[i]]$mean.of.class[which(correlograms[[i]]$correlation <= 0)[1]]
    lines(correlograms[[i]]$mean.of.class, correlograms[[i]]$correlation, lwd=0.1, col="gray70")
  }
  meanV = round(mean(values)); percentiles = round(quantile(values, c(0.025,0.975)))
  cat(meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],",",percentiles[2],"]\n",sep="") # 620, 95% CI = [415,755]
  abline(h=0, lwd=0.75, col="#de4327", lty=2); abline(v=meanV, lwd=0.75, col="#de4327", lty=1)
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=0.55, tck=-0.03, col.axis="gray30", mgp=c(0,-0.08,0), seq(0,4000,500))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=0.55, tck=-0.03, col.axis="gray30", mgp=c(0,0.20,0), pos=-50)
  title(xlab="Distance (km)", cex.lab=0.70, mgp=c(0.9,0,0), col.lab="gray30")
  title(ylab="Correlation", cex.lab=0.70, mgp=c(0.4,0,0), col.lab="gray30")
}


    # 3.3 Comparing the relative importance (RI) of each variable with the RI of its shadow version

brt_analysis = function(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
{
  theRanges = c(theRange,theRange)*1000
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
  plot.main = TRUE
  plot.folds = FALSE
  verbose = FALSE
  silent = TRUE
  keep.fold.models = FALSE
  keep.fold.vector = keep.fold.fit
  keep.fold.fit = keep.fold.fit
  showingFoldsPlot = FALSE
  spdf = SpatialPointsDataFrame(dataframe[c("x","y")], dataframe[,c("response",gbm.x)], proj4string=crs)
  myblocks = cv_spatial(spdf, column="response", k=n.folds, size=theRanges[1], selection="random", plot=F, progress=F, report=F)
  fold.vector = myblocks$folds_ids
  brt_model = gbm.step(dataframe, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
                       var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance,
                       plot.main, plot.folds, verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
  return(brt_model)
}

gbm.y = "response"
gbm.x = EnvNames
tree.complexity = 3
learning.rate = 0.001
bag.fraction = 0.75
n.trees = 10
theRange = 650 # km
n.folds = 5
crs = crs(nutsM)
keep.fold.fit = FALSE

tab = matrix(0, nrow=length(EnvNames), ncol=1); row.names(tab) = EnvNames
for (h in 1:10)
{
  for (i in 1:length(dataframes))
  {
    shadows = dataframes[[i]][,gbm.x]; colnames(shadows) = paste0(colnames(shadows),"_shadow")
    for (k in 1:dim(shadows)[2]) shadows[,k] = shadows[sample(1:dim(shadows)[1],dim(shadows)[1],replace=F),k]
    dataframe_shadow = cbind(dataframes[[i]], shadows); gbm.x_shadow = c(gbm.x, colnames(shadows))
    brt_model = brt_analysis(dataframe_shadow, gbm.y, gbm.x_shadow, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
    relativeImportances = matrix(nrow=length(gbm.x_shadow), ncol=1); row.names(relativeImportances) = gbm.x_shadow
    for (j in 1:length(gbm.x_shadow)) relativeImportances[j,1] = summary(brt_model, plotit=F)[gsub("-","\\.",gbm.x_shadow)[j],"rel.inf"]
    for (j in 1:length(gbm.x))
    {
      if (relativeImportances[gbm.x[j],1] < relativeImportances[paste0(gbm.x[j],"_shadow"),1]) tab[j,] = tab[j,]+1
    }
  }
}

tab = tab/(length(dataframes)*10); colnames(tab) = c("p-value") # to get p-values
tab = tab[!rownames(tab) %in% c("total_population", "population_density"), ,drop=F ] # remove population to avoid co linearity with urban areas
write.csv(round(tab,3), paste0("All_RI_shadow_pvalues.csv"), quote=F)


# Final variables to include, those with p-value < 0.05 
# (i.e., with RI higher than their shadow version in at least 95% of the runs): 
included_vars = which(tab < 0.05) # variables with p-value < 0.05
included_vars = tab[included_vars,]  
print(included_vars)
# temperature_winter: 0.012         
# temperature_fall:  0.022    
# precipitation_spring: 0.035      
# precipitation_fall:  0.004
# relative_humidity_winter: 0.005 
# relative_humidity_fall: 0.001               
# pastures: 0.000              
# urbanAreas: 0.031
# primaryNonForest: 0.022


  # 4. BRT model training
  # ----------------------

# final inlcuded variables
training_df = list()
for (i in 1:length(dataframes)) {
  df = dataframes[[i]]
  selected_vars = rownames(tab)[tab[, 1] < 0.05]
  df= df[, c("nutM_code","year","year_start","year_end","response", selected_vars, "x", "y")]
  training_df[[i]] = df
}

var_names = names(training_df[[1]][,selected_vars])

    # 4.1 hyperparameter tuning on selected set of variables

# fixed parameters
gbm.y = "response"
gbm.x = var_names
bag.fraction = 0.75
theRange = 650 
n.trees = 10
crs_obj = crs(nutsM) 
tuning_grid = expand.grid(learning_rate   = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005), tree_complexity = c(2, 3, 5, 6, 7, 10))

base_seed = 123
results_list = future_lapply(seq_along(training_df), function(j) {
  current_data = as.data.frame(training_df[[j]])
  rep_seed = base_seed + j   # one seed per replicate -> same folds for every lr/tc combination
  
  local_results = vector("list", nrow(tuning_grid))
  for (g in seq_len(nrow(tuning_grid))) {
    lr = tuning_grid$learning_rate[g]; tc = tuning_grid$tree_complexity[g]
    set.seed(rep_seed)
    model = tryCatch(brt_analysis(dataframe = current_data, gbm.y = gbm.y, gbm.x = gbm.x,
                                  tree.complexity = tc, learning.rate = lr,
                                  bag.fraction = bag.fraction, n.trees = n.trees,
                                  theRange = theRange, crs = crs_obj,
                                  keep.fold.fit = FALSE),
                     error = function(e) NULL)
    local_results[[g]] = data.frame(replicate       = j,
                                    learning_rate   = lr,
                                    tree_complexity = tc,
                                    bag_fraction    = bag.fraction,
                                    n_trees     = if (is.null(model)) NA else model$gbm.call$best.trees,
                                    cv_deviance = if (is.null(model)) NA else model$cv.statistics$deviance.mean,
                                    AUC         = if (is.null(model)) NA else model$cv.statistics$discrimination.mean)
  }
  do.call(rbind, local_results)
}, future.seed = TRUE, future.packages = c("dismo", "gbm", "sp", "blockCV"))


results_df = do.call(rbind, results_list)   
safe_ci = function(v) { v = v[!is.na(v)]
if (length(v) < 2 || length(unique(v)) < 2) return(c(NA_real_, NA_real_))
t.test(v)$conf.int[1:2]
}

summary_results = results_df %>%
  group_by(learning_rate, tree_complexity) %>%
  summarise(
    n_ok     = sum(!is.na(cv_deviance)),            
    mean_auc = mean(AUC, na.rm = TRUE),
    mean_dev = mean(cv_deviance, na.rm = TRUE),
    sd_dev   = sd(cv_deviance, na.rm = TRUE),
    ci_low   = safe_ci(cv_deviance)[1],
    ci_high  = safe_ci(cv_deviance)[2],
    mean_nt  = mean(n_trees, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(mean_dev)                                 

best = summary_results %>% filter(mean_nt > 1000)
# we chose lr = 0.001, tc = 3 and bag fraction = 0.75.


    # 4.2 BRT training, with selected variables and hyperparameters, 
    # and spatial cross-validation (SCV) ---------------------------------------

nberOfReplicates = 100
envVariableNames = var_names
envVariableTitles = c("Air temperature (winter, °C)","Air temperature (fall, °C)",expression("Precipitation (spring, kg/m"^2*"/day)"),expression("Precipitation (fall, kg/m"^2*"/day)"),
                      "Relative humidity (winter, %)","Relative humidity (fall, %)","Pastures","Urban areas","Primary non-forest areas")

if (!file.exists("brt_training_data/Complete_BRT_models.rds"))
{
  tree.complexity = 3 
  learning.rate = 0.001 
  bag.fraction = 0.75 # 
  keep.fold.fit = FALSE; brt_models = list()
  for (i in 1:nberOfReplicates)
  {
    dataframe = training_df[[i]]; n.trees = 10
    brt_models[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
  }
  saveRDS(brt_models, "brt_training_data/Complete_BRT_models.rds")
}
if (!file.exists("brt_training_data/FoldFit_BRT_models.rds"))
{
  tree.complexity = 3 
  learning.rate = 0.001 
  bag.fraction = 0.75 
  keep.fold.fit = TRUE; brt_models = list()
  for (i in 1:nberOfReplicates)
  {
    dataframe = training_df[[i]]; n.trees = 10
    brt_models[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
  }
  saveRDS(brt_models, "brt_training_data/FoldFit_BRT_models.rds")
}
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds")


  # 4.3. Caracterising the environmental space associated with Ps and PAs ------

dataframe_adj = read.csv("brt_training_data/Dataframe_adj.csv", head=T)
dataframe_adj = cbind(rep("0.5",dim(dataframe_adj)[1]),dataframe_adj); colnames(dataframe_adj)[1] = "response"

df = training_df[[1]][which(training_df[[1]][,"response"]==1),c("response","nutM_code",gbm.x)]
df = rbind(df, dataframe_adj); df_all = df; df1 = df; df2 = df
for (i in 1:length(training_df))
{
  tmp = training_df[[i]][which(training_df[[i]][,"response"]==0),c("response","nutM_code",gbm.x)]
  df_all = rbind(df_all, tmp[which(!tmp[,"nutM_code"]%in%df_all[,"nutM_code"]),])
  if (i == 1) df1 = rbind(df1, tmp[which(!tmp[,"nutM_code"]%in%df1[,"nutM_code"]),])
  if (i == 2) df2 = rbind(df2, tmp[which(!tmp[,"nutM_code"]%in%df2[,"nutM_code"]),])
}
pca_all = dudi.pca(df_all[,3:dim(df_all)[2]], scannf=F, nf=9); lis_all = pca_all$li[,1:2]; cos_all = pca_all$co
pca_df1 = dudi.pca(df1[,3:dim(df1)[2]], scannf=F, nf=9); lis_df1 = pca_df1$li[,1:2]; cos_df1 = pca_df1$co
pca_df2 = dudi.pca(df2[,3:dim(df2)[2]], scannf=F, nf=9); lis_df2 = pca_df2$li[,1:2]; cos_df2 = pca_df2$co
colours = c("#4676bb66","#faa521b3","#de4327b3"); cols_all = colours[(as.numeric(df_all[,"response"])*2)+1]
cols_df1 = colours[(as.numeric(df1[,"response"])*2)+1]; cols_df2 = colours[(as.numeric(df2[,"response"])*2)+1]
showingPlots = TRUE
if (showingPlots)
{
  dev.new(width=8*(2/3), height=2.7*2)
  par(mar=c(3,3,1.0,1.0), oma=c(0.75,0.75,0,0), lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
  plot(lis_all, col=NA, cex=0.7, pch=16, ann=F, axes=F)
  for (i in length(cols_all):1)
  {
    points(lis_all[i,1], lis_all[i,2], col=cols_all[i], cex=0.7, pch=16)
    points(lis_all[i,1], lis_all[i,2], col="gray30", cex=0.7, pch=1, lwd=0.3)
  }
  s.corcircle(4*cos_all, xax=1, yax=2, box=F, sub="", csub=1, clabel=0.7, possub="topleft", grid=T, cgrid=1, full=T, add.plot=T)
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=0.7, tck=-0.015, col.axis="gray30", mgp=c(0,0.2,0), at=seq(-8,8,2))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=0.7, tck=-0.015, col.axis="gray30", mgp=c(0,0.4,0), at=seq(-8,8,2))
  title(xlab=paste0("PCA axis 1 (",round((pca_all$eig[1]/sum(pca_all$eig))*100,1),"%)"), cex.lab=0.9, mgp=c(1.5,0,0), col.lab="gray30")
  title(ylab=paste0("PCA axis 2 (",round((pca_all$eig[2]/sum(pca_all$eig))*100,1),"%)"), cex.lab=0.9, mgp=c(1.7,0,0), col.lab="gray30")
  dev.new(width=8, height=2.6)
  par(mfrow=c(1,3), mar=c(3.1,3.6,1.0,1.0), oma=c(0.75,0.75,0,0), lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
  plot(lis_all, col=NA, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-4.5,4.5))
  for (i in length(cols_all):1)
  {
    points(lis_all[i,1], lis_all[i,2], col=cols_all[i], cex=0.9, pch=16)
    points(lis_all[i,1], lis_all[i,2], col="gray30", cex=0.9, pch=1, lwd=0.3)
  }
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  title(xlab=paste0("PCA axis 1 (",round((pca_all$eig[1]/sum(pca_all$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.0,0,0), col.lab="gray30")
  title(ylab=paste0("PCA axis 2 (",round((pca_all$eig[2]/sum(pca_all$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.1,0,0), col.lab="gray30")
  mtext(expression(bold(A)), side=3, line=-0.95, cex=0.9, at=-6.9)
  plot(lis_df1, col=NA, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-4.5,4.5))
  for (i in length(cols_df1):1)
  {
    points(lis_df1[i,1], lis_df1[i,2], col=cols_df1[i], cex=0.9, pch=16)
    points(lis_df1[i,1], lis_df1[i,2], col="gray30", cex=0.9, pch=1, lwd=0.3)
  }
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  title(xlab=paste0("PCA axis 1 (",round((pca_df1$eig[1]/sum(pca_df1$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.0,0,0), col.lab="gray30")
  title(ylab=paste0("PCA axis 2 (",round((pca_df1$eig[2]/sum(pca_df1$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.1,0,0), col.lab="gray30")
  mtext(expression(bold(B)), side=3, line=-0.95, cex=0.9, at=-6.9)
  plot(lis_df2, col=NA, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-4.5,4.5))
  for (i in length(cols_df2):1)
  {
    points(lis_df2[i,1], lis_df2[i,2], col=cols_df2[i], cex=0.9, pch=16)
    points(lis_df2[i,1], lis_df2[i,2], col="gray30", cex=0.9, pch=1, lwd=0.3)
  }
  axis(side=1, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  axis(side=2, lwd.tick=0.3, lwd=0.3, cex.axis=1.00, tck=-0.03, col.axis="gray30", mgp=c(0,0.5,0), at=seq(-8,8,2))
  title(xlab=paste0("PCA axis 1 (",round((pca_df2$eig[1]/sum(pca_df2$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.0,0,0), col.lab="gray30")
  title(ylab=paste0("PCA axis 2 (",round((pca_df2$eig[2]/sum(pca_df2$eig))*100,1),"%)"), cex.lab=1.1, mgp=c(2.1,0,0), col.lab="gray30")
  mtext(expression(bold(C)), side=3, line=-0.95, cex=0.9, at=-6.9)
}

  # 4.4 Assessing model predictive performance from the cross-validation

brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds"); AUCs = rep(NA, length(brt_models))
for (i in 1:length(brt_models))
{
  AUCs[i] = brt_models[[i]]$cv.statistics$discrimination.mean
}
meanV = round(mean(AUCs), 3); percentiles = round(quantile(AUCs,c(0.025,0.975)), 3)
cat("mean AUC = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# mean AUC = 0.878, 2.5th-97.5th percentiles range = [0.771, 0.956]
SIs = rep(NA, length(brt_models)); thresholds = rep(NA, length(brt_models))
for (i in 1:length(brt_models))
{
  responses = brt_models[[i]]$gbm.call$dataframe[,gbm.y]; foldFit = brt_models[[i]]$fold.fit
  if (length(foldFit) != length(responses))
  {
    stop(paste0("Replicate ",i,": fold.fit and responses do not have the same length"))
  }
  prediction = stats::plogis(foldFit) # back-transformation from the logit to the probability scale
  P = sum(responses==1); A = sum(responses==0); prev = P/(P+A); x = (P/A)*((1-prev)/prev)
  SI = 0; optimised_threshold = NA
  for (threshold in seq(0,1,0.01))
  {
    TP = length(which((responses==1)&(prediction>=threshold))) # true positives
    FN = length(which((responses==1)&(prediction<threshold))) # false negatives
    FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
    SI_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
    if (SI < SI_tmp)
    {
      SI = SI_tmp; optimised_threshold = threshold
    }
  }
  SIs[i] = SI; thresholds[i] = optimised_threshold
}
meanV = round(mean(SIs), 3); percentiles = round(quantile(SIs,c(0.025,0.975)),3)
cat("mean SI = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# mean SI = 0.808, 2.5th-97.5th percentiles range = [0.724, 0.893]
meanV = round(mean(thresholds), 3); percentiles = round(quantile(thresholds,c(0.025,0.975)),3)
cat("mean SI threshold = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# mean SI threshold = 0.37, 2.5th-97.5th percentiles range = [0.18, 0.575]
tab = cbind(round(AUCs,3), round(SIs,3), thresholds); colnames(tab) = c("AUC","SI","threshold")
write.csv(tab, "brt_training_data/SI_&_thresholds.csv", row.names=F, quote=F)

  # 4.5 Assessing model predictive performances with a temporal validation -----

training_proportion = 0.50; training_proportion = 0.75
if (!file.exists(paste0("brt_training_data/Dataframes_on_",training_proportion*100,"%.rds")))
{
  buffer = list(); dataframes_training = list(); dataframes_testing = list()
  for (i in 1:nberOfReplicates)
  {
    dataframe = training_df[[i]]; years = dataframe$year; n.trees = 10
    for (j in 1:length(years))
    {
      if (is.na(years[j])) years[j] = mean(c(dataframe[j,"year_start"],dataframe[j,"year_end"]))
    }
    dataframe = dataframe[order(years),]; years = years[order(years)]
    for (y in min(years):max(years))
    {
      indices = which(years==y)
      if (length(indices) > 1) dataframe[indices,] = dataframe[sample(indices),]
    }
    training_indices = 1:round((dim(dataframe)[1]*training_proportion)); all_indices = 1:dim(dataframe)[1]
    testing_indices = all_indices[which(!all_indices%in%training_indices)]
    dataframes_training[[i]] = dataframe[training_indices,]; dataframes_testing[[i]] = dataframe[testing_indices,]
  }
  buffer[[1]] = dataframes_training; buffer[[2]] = dataframes_testing
  saveRDS(buffer, paste0("brt_training_data/Dataframes_on_",training_proportion*100,"%.rds"))
}
buffer = readRDS(paste0("brt_training_data/Dataframes_on_",training_proportion*100,"%.rds"))
dataframes_training = buffer[[1]]; dataframes_testing = buffer[[2]]
if (!file.exists(paste0("brt_training_data/BRT_models_on_",training_proportion*100,"%.rds")))
{		
  tree.complexity = 3 # (cfr. hyperparameters tuning)
  learning.rate = 0.001 # (cfr. hyperparameters tuning)
  bag.fraction = 0.75 # (cfr. Elith et al. 2008)
  keep.fold.fit = FALSE; brt_models_validation = list()
  for (i in 1:nberOfReplicates)
  {
    dataframe = dataframes_training[[i]]; n.trees = 10
    brt_models_validation[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
  }
  saveRDS(brt_models_validation, paste0("brt_training_data/BRT_models_on_",training_proportion*100,"%.rds"))
}
brt_models_validation = readRDS(paste0("brt_training_data/BRT_models_on_",training_proportion*100,"%.rds"))

dataframe_2024 = read.csv(paste0("brt_training_data/Dataframe_2024.csv"), head=T)
BIs = rep(NA, nberOfReplicates); bins_min = seq(0.0,0.9,0.1); bins_max = seq(0.1,1.0,0.1) # BIs = "Boyce indexes"
AUCs = rep(NA, nberOfReplicates); SIs = rep(NA, nberOfReplicates); thresholds = rep(NA, nberOfReplicates)
for (i in 1:nberOfReplicates)
{
  n.trees = brt_models_validation[[i]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
  vS_all = predict.gbm(brt_models_validation[[i]], newdata=dataframe_2024, n.trees, type, single.tree)
  df_testing_presence = dataframes_testing[[i]][which(dataframes_testing[[i]][,"response"]==1),]
  vS_presence = predict.gbm(brt_models_validation[[i]], newdata=df_testing_presence, n.trees, type, single.tree)
  Fis = rep(NA, length(bins_min))
  for (j in 1:length(bins_min))
  {
    P = sum((vS_presence>=bins_min[j])&(vS_presence<bins_max[j]))/length(vS_presence)
    E = sum((vS_all>=bins_min[j])&(vS_all<bins_max[j]))/length(vS_all)
    Fis[j] = P/E
  }
  bins_mid = bins_min+0.05; bins_mid = bins_mid[which(!is.na(Fis))]; Fis = Fis[which(!is.na(Fis))]
  BIs[i] = cor(bins_mid, Fis, method="spearman")
  prediction = predict.gbm(brt_models_validation[[i]], newdata=dataframes_testing[[i]], n.trees, type, single.tree)
  observations = dataframes_testing[[i]]$response
  AUCs[i] = dismo::evaluate(p=prediction[which(observations==1)], a=prediction[which(observations==0)])@auc
  P = sum(observations==1); A = sum(observations==0); prev = P/(P+A); x = (P/A)*((1-prev)/prev) # here always = 1
  SI = 0; optimised_threshold = NA
  for (threshold in seq(0,1,0.01))
  {
    TP = length(which((observations==1)&(prediction>=threshold))) # true positives
    FN = length(which((observations==1)&(prediction<threshold))) # false negatives
    FP_pa = length(which((observations==0)&(prediction>=threshold))) # false positives
    SI_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
    if (SI < SI_tmp)
    {
      SI = SI_tmp; optimised_threshold = threshold
    }
  }
  SIs[i] = SI; thresholds[i] = optimised_threshold
  if ((showingPlots)&(i == 1))
  {
    dev.new(width=3.6, height=5); par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.4, col="gray30")
    cols = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])[(((vS_all-0)/(1-0))*100)+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
    presence_pts = dataframes_testing[[i]][which(dataframes_testing[[i]]$response==1),c("x","y")]
    pa_points = dataframes_testing[[i]][which(dataframes_testing[[i]]$response==0),c("x","y")]
    points(presence_pts, pch=16, cex=0.8, col="gray30"); points(presence_pts, pch=1, cex=0.8, col="gray30")
    points(pa_points, pch=16, cex=0.8, col="white"); points(pa_points, pch=1, cex=0.8, col="gray30")
  }
}
meanV = round(mean(BIs),3); percentiles = round(quantile(BIs,c(0.025,0.975)),3)
cat("\tmean BI = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# with a testing proportion of 50%: mean BI = 0.731, 2.5th-97.5th percentiles range = [0.394, 0.956]
# with a testing proportion of 25%: mean BI = 0.420, 2.5th-97.5th percentiles range = [0.058, 0.817]
# (the models showed a moderate positive association between predicted suitability and observed occurrence
# patterns, although the wide range indicates substantial variability in predictive performance)
meanV = round(mean(AUCs),3); percentiles = round(quantile(AUCs,c(0.025,0.975)),3)
cat("\tmean AUC = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# with a testing proportion of 50%: mean AUC = 0.875, 2.5th-97.5th percentiles range = [0.799, 0.946]
# with a testing proportion of 25%: mean AUC = 0.865, 2.5th-97.5th percentiles range = [0.750, 0.963]
meanV = round(mean(SIs),3); percentiles = round(quantile(SIs,c(0.025,0.975)),3)
cat("\tmean SI = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# with a testing proportion of 50%: mean SI = 0.838, 2.5th-97.5th percentiles range = [0.761, 0.899]
# with a testing proportion of 25%: mean SI = 0.851, 2.5th-97.5th percentiles range = [0.764, 0.941]

  # 4.6 Comparing model predictive performance associated with different PA sampling approaches

approaches = c("random","Ae_density","surveillance","espace"); ratios = c("1-1","1-3","1-6")
performance_table = matrix(nrow=2*length(approaches)*length(ratios), ncol=5); n = 30
colnames(performance_table) = c("PA_sampling_approach","P_PA_ratio","AUC","SI","SI_threshold")
for (g in 1:length(approaches))
{
  for (h in 1:length(ratios))
  {
    if (!file.exists(paste0("brt_training_data/pseudo_abcence_strategies/DFs_testing_approach_",g,"_ratio_",ratios[h],".rds")))
    {
      buffer = list()
      for (i in 1:n)
      {
        dataframe = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/Dataframe_",approaches[g],"_",ratios[h],"_",i,".csv"))
        years = dataframe$year; n.trees = 10
        for (j in 1:length(years))
        {
          if (is.na(years[j])) years[j] = mean(c(dataframe[j,"year_start"],dataframe[j,"year_end"]))
        }
        dataframe = dataframe[order(years),]; years = years[order(years)]
        for (y in min(years):max(years))
        {
          indices = which(years==y)
          if (length(indices) > 1) dataframe[indices,] = dataframe[sample(indices),]
        }
        training_indices = 1:round((dim(dataframe)[1]*training_proportion)); all_indices = 1:dim(dataframe)[1]
        testing_indices = all_indices[which(!all_indices%in%training_indices)]
        dataframes_training[[i]] = dataframe[training_indices,]; dataframes_testing[[i]] = dataframe[testing_indices,]
      }
      buffer[[1]] = dataframes_training; buffer[[2]] = dataframes_testing
      saveRDS(buffer, paste0("brt_training_data/pseudo_abcence_strategies/DFs_testing_approach_",g,"_ratio_",ratios[h],".rds"))		
    }
  }
}
for (g in 1:length(approaches))
{
  for (h in 1:length(ratios))
  {
    tree.complexity = 3 # (cfr. hyperparameters tuning)
    learning.rate = 0.001 # (cfr. hyperparameters tuning)
    bag.fraction = 0.75 # (cfr. Elith et al. 2008)
    if (!file.exists(paste0("brt_training_data/pseudo_abcence_strategies/BRT_models_approach_",g,"_ratio_",ratios[h],".rds")))
    {
      brt_models = list(); keep.fold.fit = FALSE
      for (i in 1:n)
      {
        n.trees = 10
        dataframe = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/Dataframe_",approaches[g],"_",ratios[h],"_",i,".csv"))
        dataframe$x = rep(NA, dim(dataframe)[1]); dataframe$y = rep(NA, dim(dataframe)[1])
        for (j in 1:dim(dataframe)[1])
        {
          centroid = nutsM_centroids[which(nutsM@data[,"LocationCo"]==dataframe[j,"nutM_code"]),]
          dataframe[j,"x"] = centroid[1]; dataframe[j,"y"] = centroid[2]
        }
        brt_models[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
      }
      saveRDS(brt_models, paste0("brt_training_data/pseudo_abcence_strategies/BRT_models_approach_",g,"_ratio_",ratios[h],".rds"))
    }
    if (!file.exists(paste0("brt_training_data/pseudo_abcence_strategies/BRT_foldFit_approach_",g,"_ratio_",ratios[h],".rds")))
    {
      brt_models = list(); keep.fold.fit = TRUE
      for (i in 1:n)
      {
        n.trees = 10
        dataframe = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/Dataframe_",approaches[g],"_",ratios[h],"_",i,".csv"))
        for (j in 1:dim(dataframe)[1])
        {
          centroid = nutsM_centroids[which(nutsM@data[,"LocationCo"]==dataframe[j,"nutM_code"]),]
          dataframe[j,"x"] = centroid[1]; dataframe[j,"y"] = centroid[2]
        }
        brt_models[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
      }
      saveRDS(brt_models, paste0("brt_training_data/pseudo_abcence_strategies/BRT_foldFit_approach_",g,"_ratio_",ratios[h],".rds"))
    }
    if (!file.exists(paste0("brt_training_data/pseudo_abcence_strategies/BRT_testing_approach_",g,"_ratio_",ratios[h],".rds")))
    {
      brt_models = list(); keep.fold.fit = FALSE; training_proportion = 0.75
      dataframes_training = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/DFs_testing_approach_",g,"_ratio_",ratios[h],".rds"))[[1]]
      for (i in 1:n)
      {
        dataframe = dataframes_training[[i]]; n.trees = 10
        for (j in 1:dim(dataframe)[1])
        {
          centroid = nutsM_centroids[which(nutsM@data[,"LocationCo"]==dataframe[j,"nutM_code"]),]
          dataframe[j,"x"] = centroid[1]; dataframe[j,"y"] = centroid[2]
        }
        brt_models[[i]] = brt_analysis(dataframe, gbm.y, gbm.x, tree.complexity, learning.rate, bag.fraction, n.trees, theRange, crs, keep.fold.fit)
      }
      saveRDS(brt_models, paste0("brt_training_data/pseudo_abcence_strategies/BRT_testing_approach_",g,"_ratio_",ratios[h],".rds"))
    }
  }
}
I = 0
for (g in 1:length(approaches))
{ 
  for (h in 1:length(ratios))
  {
    I = I+1; performance_table[I,"PA_sampling_approach"] = approaches[g]; performance_table[I,"P_PA_ratio"] = ratios[h]
    brt_models = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/BRT_foldFit_approach_",g,"_ratio_",ratios[h],".rds"))
    AUCs = rep(NA, length(brt_models)); SIs = rep(NA, length(brt_models)); thresholds = rep(NA, length(brt_models))
    for (i in 1:length(brt_models)) AUCs[i] = brt_models[[i]]$cv.statistics$discrimination.mean
    meanV = round(mean(AUCs), 3); percentiles = round(quantile(AUCs,c(0.025,0.975)), 3)
    performance_table[I,"AUC"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
    for (i in 1:length(brt_models))
    {
      responses = brt_models[[i]]$gbm.call$dataframe[,gbm.y]
      foldFit = brt_models[[i]]$fold.fit
      if (length(foldFit) != length(responses))
      {
        stop(paste0("Replicate ",i,": fold.fit and responses do not have the same length"))
      }
      prediction = stats::plogis(foldFit) # back-transformation from the logit to the probability scale
      P = sum(responses==1); A = sum(responses==0); prev = P/(P+A); x = (P/A)*((1-prev)/prev)
      SI = 0; optimised_threshold = NA
      for (threshold in seq(0,1,0.01))
      {
        TP = length(which((responses==1)&(prediction>=threshold))) # true positives
        FN = length(which((responses==1)&(prediction<threshold))) # false negatives
        FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
        SI_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
        if (SI < SI_tmp)
        {
          SI = SI_tmp; optimised_threshold = threshold
        }
      }
      SIs[i] = SI; thresholds[i] = optimised_threshold
    }
    meanV = round(mean(SIs), 3); percentiles = round(quantile(SIs,c(0.025,0.975)),3)
    performance_table[I,"SI"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
    meanV = round(mean(thresholds), 3); percentiles = round(quantile(thresholds,c(0.025,0.975)),3)
    performance_table[I,"SI_threshold"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
    tab = cbind(round(AUCs,3), round(SIs,3), thresholds); colnames(tab) = c("AUC","SI","threshold")
    dir.create(paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/"), showWarnings=F)
    write.csv(tab, paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/SI_&_thresholds.csv"), row.names=F, quote=F)
  }
}
for (g in 1:length(approaches))
{ 
  for (h in 1:length(ratios))
  {
    I = I+1; performance_table[I,"PA_sampling_approach"] = approaches[g]; performance_table[I,"P_PA_ratio"] = ratios[h]
    buffer = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/DFs_testing_approach_",g,"_ratio_",ratios[h],".rds"))
    dataframes_training = buffer[[1]]; dataframes_testing = buffer[[2]]
    brt_models_validation = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/BRT_testing_approach_",g,"_ratio_",ratios[h],".rds"))
    AUCs = rep(NA, length(brt_models_validation)); SIs = rep(NA, length(brt_models_validation)); thresholds = rep(NA, length(brt_models_validation))
    for (i in 1:n)
    {
      dataframe = dataframes_testing[[i]]; n.trees = brt_models_validation[[i]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      prediction = predict.gbm(brt_models_validation[[i]], newdata=dataframe, n.trees, type, single.tree)
      observations = dataframes_testing[[i]]$response
      AUCs[i] = dismo::evaluate(p=prediction[which(observations==1)], a=prediction[which(observations==0)])@auc
      P = sum(observations==1); A = sum(observations==0); prev = P/(P+A); x = (P/A)*((1-prev)/prev) # here always = 1
      SI = 0; optimised_threshold = NA
      for (threshold in seq(0,1,0.01))
      {
        TP = length(which((observations==1)&(prediction>=threshold))) # true positives
        FN = length(which((observations==1)&(prediction<threshold))) # false negatives
        FP_pa = length(which((observations==0)&(prediction>=threshold))) # false positives
        SI_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
        if (SI < SI_tmp)
        {
          SI = SI_tmp; optimised_threshold = threshold
        }
      }
      SIs[i] = SI; thresholds[i] = optimised_threshold
    }
    meanV = round(mean(AUCs), 3); percentiles = round(quantile(AUCs,c(0.025,0.975)), 3)
    performance_table[I,"AUC"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
    meanV = round(mean(SIs), 3); percentiles = round(quantile(SIs,c(0.025,0.975)),3)
    performance_table[I,"SI"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
    meanV = round(mean(thresholds), 3); percentiles = round(quantile(thresholds,c(0.025,0.975)),3)
    performance_table[I,"SI_threshold"] = paste0(meanV," [",percentiles[1],", ",percentiles[2],"]")
  }
}
write.table(performance_table, "brt_training_data/pseudo_abcence_strategies/PA_sampling_Table_S2.csv", row.names=F, quote=F, sep=";")


  # 4.7 Caracterising the composition (presences, PAs) of the spatial folds ----

nberOfPresences = c(); nberOfPseudoAbsences = c(); brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds")
for (i in 1:nberOfReplicates)
{
  dataframe = dataframes[[i]]; fold.vector = brt_models[[i]]$fold.vector
  for (j in 1:n.folds)
  {
    nberOfPresences = c(nberOfPresences, sum((dataframe[,"response"]==1)&(fold.vector==j)))
    nberOfPseudoAbsences = c(nberOfPseudoAbsences, sum((dataframe[,"response"]==0)&(fold.vector==j)))
  }
}
meanV = round(mean(nberOfPresences)); percentiles = round(quantile(nberOfPresences,c(0.025,0.975)))
cat("number of presence entries per fold = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# number of presence entries per fold = 12, (2.5th-97.5th percentiles range = [5, 19])
meanV = round(mean(nberOfPseudoAbsences)); percentiles = round(quantile(nberOfPseudoAbsences,c(0.025,0.975)))
cat("Number of pseudo-absence entries per fold = ",meanV,", 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# number of pseudo-absence entries per fold = 12, (2.5th-97.5th percentiles range = [7, 16])

  # 4.8 Retrieving the relative importance of each environmental variable
envVariableNames = var_names
if (!file.exists("brt_training_data/Relative_importances.csv"))
{
  brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds")
  relativeImportances = matrix(nrow=length(envVariableNames), ncol=1)
  row.names(relativeImportances) = envVariableNames
  for (i in 1:length(envVariableNames))
  {
    vS = rep(NA, length(brt_models))
    for (j in 1:length(brt_models))
    {
      tab = summary(brt_models[[j]])
      row.names(tab) = gsub("`","",row.names(tab))
      vS[j] = tab[envVariableNames[i],"rel.inf"]
    }
    meanV = round(mean(vS,na.rm=T),1)
    minV = round(min(vS,na.rm=T),1); maxV = round(max(vS,na.rm=T),1)
    percentiles = round(quantile(vS,c(0.025,0.975)),1)
    if (!grepl("\\.",meanV)) meanV = paste0(meanV,".0")
    if (!grepl("\\.",percentiles[1])) percentiles[1] = paste0(percentiles[1],".0")
    if (!grepl("\\.",percentiles[2])) percentiles[2] = paste0(percentiles[2],".0")
    relativeImportances[i,1] = paste0(meanV,"% [",percentiles[1],"-",percentiles[2],"]")
  }
  write.table(relativeImportances, "Relative_importances.csv", quote=F, col.names=F, sep=",")
}


  # 4.9. Generating and plotting the response curves for each variable ---------

pdf("All_response_curves.pdf", height=5, width=8) # dev.new(height=4.5, width=8)
par(mfrow=c(3,3), oma=c(0.3,0.3,0.5,0.5), mar=c(2.7,1.9,0.5,0.5), lwd=0.4, col="gray30", col.axis="gray30", fg="gray30")
plotting_all_curves = TRUE; plotting_all_curves = FALSE
relativeImportances = read.csv("Relative_importances.csv", head=F)
relativeImportanceValues = as.numeric(gsub("%","",relativeImportances[,2]))
envVariableValues = matrix(nrow=5, ncol=length(envVariableNames))
colnames(envVariableValues) = envVariableNames
row.names(envVariableValues) = c("quartile1","quartile2","quartile3","minV","maxV")
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds")
for (i in 1:length(envVariableNames))
{
  vS = c()
  for (j in 1:nberOfReplicates) vS = c(vS, training_df[[j]][,gsub("-","\\.",envVariableNames)[i]])
  qS = quantile(vS, c(0.25,0.50,0.75)); minV = min(vS); maxV = max(vS)
  envVariableValues[,i] = cbind(qS[1], qS[2], qS[3], minV, maxV)
}
for (i in order(relativeImportanceValues, decreasing=T))
{
  valuesInterval = (envVariableValues["maxV",i]-envVariableValues["minV",i])/100
  dataframe_q1 = data.frame(matrix(nrow=length(seq(envVariableValues["minV",i],envVariableValues["maxV",i],valuesInterval)), ncol=length(envVariableNames)))
  dataframe_q2 = data.frame(matrix(nrow=length(seq(envVariableValues["minV",i],envVariableValues["maxV",i],valuesInterval)), ncol=length(envVariableNames)))
  dataframe_q3 = data.frame(matrix(nrow=length(seq(envVariableValues["minV",i],envVariableValues["maxV",i],valuesInterval)), ncol=length(envVariableNames)))
  colnames(dataframe_q1) = envVariableNames; colnames(dataframe_q2) = envVariableNames; colnames(dataframe_q3) = envVariableNames
  for (j in 1:length(envVariableNames))
  {
    interval = (envVariableValues["maxV",j]-envVariableValues["minV",j])/100
    if (i == j)
    {
      dataframe_q1[,envVariableNames[j]] = seq(envVariableValues["minV",j], envVariableValues["maxV",j],interval)
      dataframe_q2[,envVariableNames[j]] = seq(envVariableValues["minV",j], envVariableValues["maxV",j],interval)
      dataframe_q3[,envVariableNames[j]] = seq(envVariableValues["minV",j], envVariableValues["maxV",j],interval)
    }	else	{
      dataframe_q1[,envVariableNames[j]] = rep(envVariableValues["quartile1",j],dim(dataframe_q1)[1])
      dataframe_q2[,envVariableNames[j]] = rep(envVariableValues["quartile2",j],dim(dataframe_q2)[1])
      dataframe_q3[,envVariableNames[j]] = rep(envVariableValues["quartile3",j],dim(dataframe_q3)[1])
    }
  }
  predictions_q1 = list(); predictions_q2 = list(); predictions_q3 = list()
  for (j in 1:length(brt_models))
  {
    n.trees = brt_models[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
    prediction_q1 = predict.gbm(brt_models[[j]], newdata=dataframe_q1, n.trees, type, single.tree)
    prediction_q2 = predict.gbm(brt_models[[j]], newdata=dataframe_q2, n.trees, type, single.tree)
    prediction_q3 = predict.gbm(brt_models[[j]], newdata=dataframe_q3, n.trees, type, single.tree)
    if (j == 1)
    {
      prediction_matrix_q1 = matrix(nrow=length(prediction_q1), ncol=length(brt_models))
      prediction_matrix_q2 = matrix(nrow=length(prediction_q2), ncol=length(brt_models))
      prediction_matrix_q3 = matrix(nrow=length(prediction_q3), ncol=length(brt_models))
      minX = min(dataframe_q1[,envVariableNames[i]]); maxX = max(dataframe_q1[,envVariableNames[i]])
      minY = min(prediction_q1); maxY = max(prediction_q1)
    }
    predictions_q1[[j]] = prediction_q1; prediction_matrix_q1[,j] = prediction_q1
    predictions_q2[[j]] = prediction_q2; prediction_matrix_q2[,j] = prediction_q2
    predictions_q3[[j]] = prediction_q3; prediction_matrix_q3[,j] = prediction_q3
    if (minX > min(dataframe_q1[,envVariableNames[i]])) minX = min(dataframe_q1[,envVariableNames[i]])
    if (maxX < max(dataframe_q1[,envVariableNames[i]])) maxX = max(dataframe_q1[,envVariableNames[i]])
    if (minX > min(dataframe_q2[,envVariableNames[i]])) minX = min(dataframe_q2[,envVariableNames[i]])
    if (maxX < max(dataframe_q2[,envVariableNames[i]])) maxX = max(dataframe_q2[,envVariableNames[i]])
    if (minX > min(dataframe_q3[,envVariableNames[i]])) minX = min(dataframe_q3[,envVariableNames[i]])
    if (maxX < max(dataframe_q3[,envVariableNames[i]])) maxX = max(dataframe_q3[,envVariableNames[i]])
    if (minY > min(prediction_q1)) minY = min(prediction_q1)
    if (maxY < max(prediction_q1)) maxY = max(prediction_q1)
    if (minY > min(prediction_q2)) minY = min(prediction_q2)
    if (maxY < max(prediction_q2)) maxY = max(prediction_q2)
    if (minY > min(prediction_q3)) minY = min(prediction_q3)
    if (maxY < max(prediction_q3)) maxY = max(prediction_q3)
  }
  minY = minY-((maxY-minY)*0.3); maxY = maxY+((maxY-minY)*0.04)
  cols = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)); cols[5] = "gray50"
  if (plotting_all_curves)
  {
    for (j in 1:length(brt_models))
    {
      if (j == 1)
      {
        plot(dataframe_q1[,envVariableNames[i]], predictions_q1[[j]], col=rgb(250,165,33,255,maxColorValue=255), ann=F, axes=F, lwd=0.1, type="l", xlim=c(minX,maxX), ylim=c(minY,maxY))
      }	else	{
        lines(dataframe_q1[,envVariableNames[i]], predictions_q1[[j]], col=rgb(250,165,33,255,maxColorValue=255), lwd=0.1)
      }
    }
    for (j in 1:length(brt_models)) lines(dataframe_q3[,envVariableNames[i]], predictions_q3[[j]], col=rgb(70,118,187,maxColorValue=255), lwd=0.1)
    for (j in 1:length(brt_models)) lines(dataframe_q2[,envVariableNames[i]], predictions_q2[[j]], col=rgb(222,67,39,255,maxColorValue=255), lwd=0.1)
  }	else	{
    xS_q1 = dataframe_q1[,envVariableNames[i]]; yS_median_q1 = apply(prediction_matrix_q1, 1, median)
    xS_q2 = dataframe_q2[,envVariableNames[i]]; yS_median_q2 = apply(prediction_matrix_q2, 1, median)
    xS_q3 = dataframe_q3[,envVariableNames[i]]; yS_median_q3 = apply(prediction_matrix_q3, 1, median)
    yS_min_q1 = apply(prediction_matrix_q1, 1, min); yS_max_q1 = apply(prediction_matrix_q1, 1, max)
    yS_min_q2 = apply(prediction_matrix_q2, 1, min); yS_max_q2 = apply(prediction_matrix_q2, 1, max)
    yS_min_q3 = apply(prediction_matrix_q3, 1, min); yS_max_q3 = apply(prediction_matrix_q3, 1, max)
    plot(xS_q2, yS_median_q2, type="n", ann=F, axes=F, lwd=0.1, xlim=c(minX,maxX), ylim=c(minY,maxY))
    # polygon(c(xS_q1, rev(xS_q1)), c(yS_min_q1, rev(yS_max_q1)), col=rgb(250,165,33,20,maxColorValue=255), border=NA)
    # polygon(c(xS_q3, rev(xS_q3)), c(yS_min_q3, rev(yS_max_q3)), col=rgb(70,118,187,20,maxColorValue=255), border=NA)
    polygon(c(xS_q2, rev(xS_q2)), c(yS_min_q2, rev(yS_max_q2)), col=rgb(222,67,39,20,maxColorValue=255), border=NA)
    lines(xS_q1, yS_min_q1, col=rgb(250,165,33,255,maxColorValue=255), lwd=0.7, lty=2)
    lines(xS_q1, yS_max_q1, col=rgb(250,165,33,255,maxColorValue=255), lwd=0.7, lty=2)
    lines(xS_q1, yS_median_q1, col=rgb(250,165,33,255,maxColorValue=255), lwd=0.9)
    lines(xS_q3, yS_min_q3, col=rgb(70,118,187,255,maxColorValue=255), lwd=0.7, lty=2)
    lines(xS_q3, yS_max_q3, col=rgb(70,118,187,255,maxColorValue=255), lwd=0.7, lty=2)
    lines(xS_q3, yS_median_q3, col=rgb(70,118,187,255,maxColorValue=255), lwd=0.9)
    lines(xS_q2, yS_median_q2, col=rgb(222,67,39,255,maxColorValue=255), lwd=0.9)
  }		
  axis(side=1, lwd.tick=0.4, cex.axis=0.85, lwd=0, tck=-0.035, col.axis="gray30", mgp=c(0,0.25,0))
  axis(side=2, lwd.tick=0.4, cex.axis=0.85, lwd=0, tck=-0.035, col.axis="gray30", mgp=c(0,0.40,0))
  title(xlab=envVariableTitles[i], cex.lab=0.95, mgp=c(1.5,0,0), col.lab="gray30")
  txt = relativeImportances[which(relativeImportances[,1]==envVariableNames[i]),2]
  mtext(paste0("RI = ",txt), side=3, line=-8.5, cex=0.6, col="gray30")
  box(lwd=0.4, col="gray30")
}
dev.off()


# ------------------------------------------------------------------------------
#                           III. PAST PROJECTIONS
# ------------------------------------------------------------------------------

# 1. Conducting past projections and generating the past projections figure ----
# ------------------------------------------------------------------------------

# note: original climate and land use variables are available at:
seasons = c("winter","spring","summer","fall"); season_months = list(); season_years = list()
season_months[[1]] = c("12","01","02"); season_months[[2]] = c("03","04","05"); season_months[[3]] = c("06","07","08"); season_months[[4]] = c("09","10","11")
season_years[[1]] = c(-1,0,0); season_years[[2]] = c(0,0,0); season_years[[3]] = c(0,0,0); season_years[[4]] = c(0,0,0)

scenarios = c("Obsclim","Counterclim"); periods = c("1901-1924","1925-1949","1950-1974","1975-1999","2000-2021"); envValues_list = list()
allEnvVariableNames = c(names_env, "total_population_log","population_density","population_density_log")
for (i in 1:length(scenarios))
{
  directory1 = paste0("../../Rasters/Environmental_rasters/ISIMIP3a/",scenarios[i],"/20CRv3-ERA5")
  directory2 = "../../Rasters/Environmental_rasters/Land-use/"
  directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"
  temperature = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_tas_1901_2021_monmean.nc"))-273.15 # conversion to Celcius degrees
  precipitation = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_pr_1901_2021_monmean.nc"))*60*60*24 # conversion to kg/m2/day
  relative_humidity = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_hurs_1901_2021_monmean.nc"))
  croplands = brick(paste0(directory2,"croplands_LUH2-GCB2024_1901_2024.nc"))
  pastures = brick(paste0(directory2,"pastures_LUH2-GCB2024_1901_2024.nc"))
  rangelands = brick(paste0(directory2,"rangeland_LUH2-GCB2024_1901_2024.nc"))
  urbanAreas = brick(paste0(directory2,"urbanAreas_LUH2-GCB2024_1901_2024.nc"))
  primaryForest = brick(paste0(directory2,"primaryForest_LUH2-GCB2024_1901_2024.nc"))
  primaryNonForest = brick(paste0(directory2,"primaryNonF_LUH2-GCB2024_1901_2024.nc"))
  secondaryForest = brick(paste0(directory2,"secondaryForest_LUH2-GCB2024_1901_2024.nc"))
  secondaryNonForest = brick(paste0(directory2,"secondaryNonF_LUH2-GCB2024_1901_2024.nc"))
  
  population_density = brick(paste0(directory3,"Population_density_05deg_1900_2025_historical.nc"), varname="Population_density")
  #population_density = population_density[[2:125]] # 1901 to 2024
  population_density_log = population_density; population_density_log[] = log(population_density[]+1); pop_years = names(population_density)
  for (j in 1:length(pop_years)) pop_years[j] = gsub("X","",unlist(strsplit(pop_years[j],"\\."))[1])
  cell_areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")/(1000*1000)
  population = population_density*cell_areas # sum(population[[126]][], na.rm=T) # 8236298503 --> ok
  buffer = list()
  for (j in 1:length(periods))
  {		
    year1 = as.numeric(unlist(strsplit(periods[[j]],"-"))[1]); year2 = as.numeric(unlist(strsplit(periods[[j]],"-"))[2]); years = year1:year2
    envValues = matrix(nrow=length(nutsM), ncol=length(allEnvVariableNames)); colnames(envValues) = allEnvVariableNames
    for (k in 1:length(seasons))
    {
      indexes = c()
      for (l in 1:length(years))
      {
        for (m in 1:length(season_months[[k]]))
        {
          indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(temperature))))
        }
      }
      rast = mean(temperature[[indexes]])
      envValues[,paste0("temperature_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      indexes = c()
      for (l in 1:length(years))
      {
        for (m in 1:length(season_months[[k]]))
        {
          indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(precipitation))))
        }
      }
      rast = mean(precipitation[[indexes]])
      envValues[,paste0("precipitation_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      indexes = c()
      for (l in 1:length(years))
      {
        for (m in 1:length(season_months[[k]]))
        {
          indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(relative_humidity))))
        }
      }
      rast = mean(relative_humidity[[indexes]])
      envValues[,paste0("relative_humidity_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
    }	
    rasts = list()
    for (k in 1:length(years))
    {
      indexes = c()
      for (l in 1:length(season_months[[3]])) # summer
      {
        indexes = c(indexes, which(grepl(paste0(years[k]+season_years[[3]][l],".",season_months[[3]][l]),names(temperature))))
      }
      rast_summer = mean(temperature[[indexes]])
      indexes = c()
      for (l in 1:length(season_months[[1]])) # winter
      {
        indexes = c(indexes, which(grepl(paste0(years[k]+season_years[[1]][l],".",season_months[[1]][l]),names(temperature))))
      }
      rast_winter = mean(temperature[[indexes]])
      rasts[[k]] = rast_summer-rast_winter
    }
    envValues[,"seasonal_variation"] = exactextractr::exact_extract(mean(stack(rasts)),nutsM,fun="mean")
    land_cover_rasters = list(
      pastures = pastures,
      rangelands = rangelands,
      croplands = croplands,
      urbanAreas = urbanAreas,
      primaryForest = primaryForest,
      primaryNonForest = primaryNonForest,
      secondaryForest = secondaryForest,
      secondaryNonForest = secondaryNonForest
    )
    
    for (k in names(land_cover_rasters))
    {
      indexes = which((1900 + as.numeric(gsub("X", "", names(land_cover_rasters[[k]])))) %in% years)
      rast = mean(land_cover_rasters[[k]][[indexes]])
      envValues[, k] = exactextractr::exact_extract(rast, nutsM, fun = "mean")
    }    
    rast = mean(population[[which((1899+as.numeric(gsub("layer\\.","",names(population))))%in%years)]])
    envValues[,"total_population"] = exactextractr::exact_extract(rast,nutsM,fun="sum")
    envValues[,"total_population_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="sum")+1)
    rast[] = mean(population_density[[which(as.numeric(pop_years)%in%years)]])[]
    envValues[,"population_density"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
    envValues[,"population_density_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="mean")+1)
    buffer[[j]] = envValues	
  }
  envValues_list[[i]] = buffer
}
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds"); thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
envVariableNames = c("temperature_winter","temperature_fall","precipitation_spring","precipitation_fall",
                     "relative_humidity_winter","relative_humidity_fall","pastures","urbanAreas","primaryNonForest")
for (i in 1:length(scenarios))
{
  for (j in 1:length(periods))
  {
    mat = matrix(nrow=length(nutsM), ncol=nberOfReplicates)
    row.names(mat) = nutsM@data[,"LocationCo"]
    dataframe = as.data.frame(envValues_list[[i]][[j]][,envVariableNames])
    for (k in 1:length(brt_models))
    {
      n.trees = brt_models[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      mat[,k] = predict.gbm(brt_models[[k]], newdata=dataframe, n.trees, type, single.tree)
    }
    write.table(mat, paste0("brt_projections_outputs/Past_projections/Past_",tolower(scenarios[i]),"_",periods[j],"_20CRv3-ERA5.csv"), col.names=F, quote=F, sep=",")
  }
}
brts = list()
for (i in 1:length(scenarios))
{
  buffer = list()
  for (j in 1:length(periods))
  {
    buffer[[j]] = read.csv(paste0("brt_projections_outputs/Past_projections/Past_",tolower(scenarios[i]),"_",periods[j],"_20CRv3-ERA5.csv"), head=F)
    buffer[[j]] = buffer[[j]][,2:dim(buffer[[j]])[2]]
  }
  brts[[i]] = buffer
}

suffix = "A"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "B"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "C"; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE

pdf(paste0("Figure_1",suffix,".pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")
scenario_names = c("Historical","Counterfactual")
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # "D9" = 85%
colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") # "00" = 0%
colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5"
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18])
colourScale2[9] = paste0(colourScale2[9],"00")
colourScale3 = colourScale2; colourScale3[9] = "#E5E5E5"
colourScale4 = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[11:111]
colourScale = colourScale1
if (plottingStDevs) { colourScale = colourScale4 }
for (i in 1:length(scenarios))
{
  for (j in 1:length(periods))
  {
    brt = as.matrix(brts[[i]][[j]])
    means = rep(NA, dim(brt)[1]); stdvs = rep(NA, dim(brt)[1])
    binaries = rep(0, dim(brt)[1])
    for (k in 1:dim(brt)[1])
    {
      means[k] = mean(brt[k,]); stdvs[k] = sd(brt[k,])
      temp = rep(0, dim(brt)[2])
      temp[which(brt[k,]>=thresholds_SI[,"threshold"])] = 1
      binaries[k] = round(mean(temp))
    }
    cols = colourScale[(((means-0)/(1-0))*100)+1]; # print(max(stdvs))
    if (plottingStDevs) cols = colourScale[(((stdvs-0)/(0.25-0))*100)+1]
    if (plottingBinaryMap) cols = c("gray90","#de4327")[binaries+1]
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
  if (j == length(periods))
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

# ------------------------------------------------------------------------------
#                           IV. FUTURE PREDICTIONS
# ------------------------------------------------------------------------------

# 1. Conducting future projections and generating the future projections figure
# ------------------------------------------------------------------------------

scenarios = c("SSP1-2.6","SSP3-7.0","SSP5-8.5"); scenarios2 = c("ssp126","ssp370","ssp585");
scenarios3= c("SSP1","SSP3","SSP5"); periods = c("2020-2024","2025-2049","2050-2074","2075-2100")
models = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
approaches = c("random","Ae_density","surveillance","espace"); ratios = c("1-1","1-3","1-6")

if (!file.exists("EnvValues_2020-2024.rds"))
{
  directory1 = "../../Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5"
  directory2 = "../../Rasters/Environmental_rasters/Land-use/"
  directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"		
  temperature_1 = brick(paste0(directory1,"/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))-273.15 # conversion to Celcius degrees
  precipitation_1 = brick(paste0(directory1,"/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"))*60*60*24 # conversion to kg/m2/day
  relative_humidity_1 = brick(paste0(directory1,"/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"))
  temperature_2 = brick(paste0(directory1,"/20crv3-era5_obsclim_tas_2022_2024_monmean.nc"))-273.15 # conversion to Celcius degrees
  precipitation_2 = brick(paste0(directory1,"/20crv3-era5_obsclim_pr_2022_2024_monmean.nc"))*60*60*24 # conversion to kg/m2/day
  relative_humidity_2 = brick(paste0(directory1,"/20crv3-era5_obsclim_hurs_2022_2024_monmean.nc"))
  temperature = brick(stack(temperature_1, temperature_2))
  precipitation = brick(stack(precipitation_1, precipitation_2))
  relative_humidity = brick(stack(relative_humidity_1, relative_humidity_2))
  
  croplands = brick(paste0(directory2,"croplands_LUH2-GCB2024_1901_2024.nc"))
  pastures = brick(paste0(directory2,"pastures_LUH2-GCB2024_1901_2024.nc"))
  rangelands = brick(paste0(directory2,"rangeland_LUH2-GCB2024_1901_2024.nc"))
  urbanAreas = brick(paste0(directory2,"urbanAreas_LUH2-GCB2024_1901_2024.nc"))
  primaryForest = brick(paste0(directory2,"primaryForest_LUH2-GCB2024_1901_2024.nc"))
  primaryNonForest = brick(paste0(directory2,"primaryNonF_LUH2-GCB2024_1901_2024.nc"))
  secondaryForest = brick(paste0(directory2,"secondaryForest_LUH2-GCB2024_1901_2024.nc"))
  secondaryNonForest = brick(paste0(directory2,"secondaryNonF_LUH2-GCB2024_1901_2024.nc"))
  
  population_density = brick(paste0(directory3,"/Population_density_05deg_1900_2025_historical.nc"), varname="Population_density")
  population_density_log = population_density; population_density_log[] = log(population_density[]+1); pop_years = names(population_density)
  
  for (i in 1:length(pop_years)) pop_years[i] = gsub("X","",unlist(strsplit(pop_years[i],"\\."))[1])
  cell_areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")/(1000*1000)
  population = population_density*cell_areas # sum(population[[126]][], na.rm=T) # 8236298503 --> ok
  year1 = 2020; year2 = 2024; years = year1:year2
  envValues = matrix(nrow=length(nutsM), ncol=length(allEnvVariableNames)); colnames(envValues) = allEnvVariableNames
  for (l in 1:length(seasons))
  {
    indexes = c()
    for (m in 1:length(years))
    {
      for (n in 1:length(season_months[[l]]))
      {
        indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(temperature))))
      }
    }
    rast = mean(temperature[[indexes]])
    envValues[,paste0("temperature_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
    indexes = c()
    for (m in 1:length(years))
    {
      for (n in 1:length(season_months[[l]]))
      {
        indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(precipitation))))
      }
    }
    rast = mean(precipitation[[indexes]])
    envValues[,paste0("precipitation_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
    indexes = c()
    for (m in 1:length(years))
    {
      for (n in 1:length(season_months[[l]]))
      {
        indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(relative_humidity))))
      }
    }
    rast = mean(relative_humidity[[indexes]])
    envValues[,paste0("relative_humidity_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
  }
  rasts = list()
  for (l in 1:length(years))
  {
    indexes = c()
    for (m in 1:length(season_months[[3]])) # summer
    {
      indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[3]][m],".",season_months[[3]][m]),names(temperature))))
    }
    rast_summer = mean(temperature[[indexes]])
    indexes = c()
    for (m in 1:length(season_months[[1]])) # winter
    {
      indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[1]][m],".",season_months[[1]][m]),names(temperature))))
    }
    rast_winter = mean(temperature[[indexes]])
    rasts[[l]] = rast_summer-rast_winter
  }
  envValues[,"seasonal_variation"] = exactextractr::exact_extract(mean(stack(rasts)),nutsM,fun="mean")
  land_cover_rasters = list(
    pastures = pastures,
    rangelands = rangelands,
    croplands = croplands,
    urbanAreas = urbanAreas,
    primaryForest = primaryForest,
    primaryNonForest = primaryNonForest,
    secondaryForest = secondaryForest,
    secondaryNonForest = secondaryNonForest
  )
  
  for (k in names(land_cover_rasters))
  {
    indexes = which((1900 + as.numeric(gsub("X", "", names(land_cover_rasters[[k]])))) %in% years)
    rast = mean(land_cover_rasters[[k]][[indexes]])
    envValues[, k] = exactextractr::exact_extract(rast, nutsM, fun = "mean")
    }    
  
  rast = mean(population[[which((1899+as.numeric(gsub("layer\\.","",names(population))))%in%years)]])
  envValues[,"total_population"] = exactextractr::exact_extract(rast,nutsM,fun="sum")
  envValues[,"total_population_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="sum")+1)
  rast[] = mean(population_density[[which(as.numeric(pop_years)%in%years)]])[]
  envValues[,"population_density"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
  envValues[,"population_density_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="mean")+1)
  write.csv(envValues, "EnvValues_2020-2024.csv", row.names=F, quote=F)
  saveRDS(envValues, "EnvValues_2020-2024.rds")
}	else	{
  envValues = readRDS("EnvValues_2020-2024.rds")
}
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds"); dataframe = as.data.frame(envValues[,envVariableNames])
mat = matrix(nrow=length(nutsM), ncol=nberOfReplicates); row.names(mat) = nutsM@data[,"LocationCo"]
for (l in 1:length(brt_models))
{
  n.trees = brt_models[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
  mat[,l] = predict.gbm(brt_models[[l]], newdata=dataframe, n.trees, type, single.tree)
}
write.table(mat, paste0("brt_projections_outputs/Current_projections/Proj_2020-2024.csv"), col.names=F, quote=F, sep=",")
for (g in 1:length(approaches))
{
  for (h in 1:length(ratios))
  {
    brt_models = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/BRT_foldFit_approach_",g,"_ratio_",ratios[h],".rds"))
    mat = matrix(nrow=length(nutsM), ncol=length(brt_models)); row.names(mat) = nutsM@data[,"LocationCo"]
    for (l in 1:length(brt_models))
    {
      n.trees = brt_models[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      mat[,l] = predict.gbm(brt_models[[l]], newdata=dataframe, n.trees, type, single.tree)
    }
    write.table(mat, paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/Proj_2020-2024.csv"), col.names=F, quote=F, sep=",")
  }
}

if (!file.exists("EnvValues_future_proj.rds"))
{
  directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"; envValues_list = list()
  population_density_1 = brick(paste0(directory3,"Population_density_05deg_1900_2025_historical.nc"), varname="Population_density")
  for (i in 1:length(scenarios))
  {
    buffer1 = list()
    directory2 = paste0("../../Rasters/Environmental_rasters/Land-use/Future_landuse_annual/",scenarios2[i])
    landCoverVariableNames = as.character(read.csv("../../Rasters/Environmental_rasters/Land-use/LC_vars.csv")[1:12,2])
    landCoverVariableIDs = as.character(read.csv("../../Rasters/Environmental_rasters/Land-use/LC_vars.csv")[1:12,1])
    land_covers1 = list(); land_covers2 = list()
    for (j in 1:length(landCoverVariableIDs))
    {
      file = list.files(directory2); land_covers1[[j]] = brick(paste0(directory2,"/",file), varname=landCoverVariableIDs[j])
    }
    population_density_2 = brick(paste0(directory3,"/Population_density_05deg_2026_2100_",unlist(strsplit(scenarios3[i],"-"))[1],".nc"), varname="Population_density")
    population_density = brick(stack(population_density_1, population_density_2))
    population_density_log = population_density; population_density_log[] = log(population_density[]+1); pop_years = names(population_density)
    for (j in 1:length(pop_years)) pop_years[j] = gsub("X","",unlist(strsplit(pop_years[j],"\\."))[1])
    cell_areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")/(1000*1000)
    population = population_density*cell_areas # sum(population[[126]][], na.rm=T) # 8236298503 --> ok
    for (j in 1:length(models))
    {
      print(c(i,j))
      buffer2 = list()
      directory1 = paste0("../../Rasters/Environmental_rasters/ISIMIP3b/",scenarios2[i],"/",models[j])
      files = list.files(directory1); files = files[which(grepl("tas",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) temperature = brick(paste0(directory1,"/",files[k]))-273.15 # conversion to Celcius degrees
        if (k != 1) temperature = brick(stack(temperature, brick(paste0(directory1,"/",files[k]))-273.15))
      }
      files = list.files(directory1); files = files[which(grepl("pr",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) precipitation = brick(paste0(directory1,"/",files[k]))*60*60*24 # conversion to kg/m2/day
        if (k != 1) precipitation = brick(stack(precipitation, brick(paste0(directory1,"/",files[k]))*60*60*24))
      }
      files = list.files(directory1); files = files[which(grepl("hurs",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) relative_humidity = brick(paste0(directory1,"/",files[k]))
        if (k != 1) relative_humidity = brick(stack(relative_humidity, brick(paste0(directory1,"/",files[k]))))
      }
      for (k in 1:length(periods))
      {
        year1 = as.numeric(unlist(strsplit(periods[[k]],"-"))[1]); year2 = as.numeric(unlist(strsplit(periods[[k]],"-"))[2]); years = year1:year2
        envValues = matrix(nrow=length(nutsM), ncol=length(allEnvVariableNames)); colnames(envValues) = allEnvVariableNames
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(temperature))))
        rast = mean(temperature[[indexes]])
        envValues[,"temperature"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(precipitation))))
        rast = mean(precipitation[[indexes]])
        envValues[,"precipitation"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(relative_humidity))))
        rast = mean(relative_humidity[[indexes]])
        envValues[,"relative_humidity"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        for (l in 1:length(seasons))
        {
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(temperature))))
            }
          }
          rast = mean(temperature[[indexes]])
          envValues[,paste0("temperature_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(precipitation))))
            }
          }
          rast = mean(precipitation[[indexes]])
          envValues[,paste0("precipitation_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(relative_humidity))))
            }
          }
          rast = mean(relative_humidity[[indexes]])
          envValues[,paste0("relative_humidity_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        }
        rasts = list()
        for (l in 1:length(years))
        {
          indexes = c()
          for (m in 1:length(season_months[[3]])) # summer
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[3]][m],".",season_months[[3]][m]),names(temperature))))
          }
          rast_summer = mean(temperature[[indexes]])
          indexes = c()
          for (m in 1:length(season_months[[1]])) # winter
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[1]][m],".",season_months[[1]][m]),names(temperature))))
          }
          rast_winter = mean(temperature[[indexes]])
          rasts[[l]] = rast_summer-rast_winter
        }
        envValues[,"seasonal_temperature_variation"] = exactextractr::exact_extract(mean(stack(rasts)),nutsM,fun="mean")
        variable_codes = c("pastures","rangelands","croplands","urbanAreas","primaryForest","primaryNonF","secondary","secondaryForest","secondaryNonF")
        variable_names = c("pasture","rangeland","crops","urban land","forested primary land","non-forested primary land","secondary",
                           "potentially forested secondary land","potentially non-forested secondary land"); land_covers2 = list()
        for (l in 1:length(variable_names))
        {
          names = gsub("\\."," ",landCoverVariableNames); indexes1 = which(landCoverVariableNames==variable_names[l])
          if (length(indexes1) == 0) indexes1 = which(grepl(variable_names[l],names))
          # if (variable_names[l] == "pasture") indexes1 = c(indexes1, which(grepl("rangeland",names)))
          indexes2 = which((2015+as.numeric(gsub("X","",names(land_covers1[[indexes1[1]]]))))%in%years); rasts = list()
          for (m in 1:length(indexes2))
          {
            rast = land_covers1[[indexes1[1]]][[indexes2[m]]]
            if (length(indexes1) > 1)
            {
              for (n in 2:length(indexes1)) rast[] = rast[]+land_covers1[[indexes1[n]]][[indexes2[m]]][]
            }
            rasts[[m]] = rast
          }
          land_cover = mean(stack(rasts)); names(land_cover) = variable_codes[l]; land_covers2[[l]] = raster::aggregate(land_cover,2)
        }
        envValues[,"primaryForest"] = exactextractr::exact_extract(land_covers2[[5]],nutsM,fun="mean")
        envValues[,"primaryNonForest"] = exactextractr::exact_extract(land_covers2[[6]],nutsM,fun="mean")								
        envValues[,"secondaryForest"] = exactextractr::exact_extract(land_covers2[[7]],nutsM,fun="mean")
        envValues[,"secondary_forest_areas"] = exactextractr::exact_extract(land_covers2[[8]],nutsM,fun="mean")
        envValues[,"secondaryNonForest"] = exactextractr::exact_extract(land_covers2[[9]],nutsM,fun="mean")								
        envValues[,"pastures"] = exactextractr::exact_extract(land_covers2[[1]],nutsM,fun="mean")
        envValues[,"rangelands"] = exactextractr::exact_extract(land_covers2[[2]],nutsM,fun="mean")
        envValues[,"croplands"] = exactextractr::exact_extract(land_covers2[[3]],nutsM,fun="mean")
        envValues[,"urbanAreas"] = exactextractr::exact_extract(land_covers2[[4]],nutsM,fun="mean")
        rast = mean(population[[which((1899+as.numeric(gsub("layer\\.","",names(population))))%in%years)]])
        envValues[,"total_population"] = exactextractr::exact_extract(rast,nutsM,fun="sum")
        envValues[,"total_population_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="sum")+1)
        rast[] = mean(population_density[[which(as.numeric(pop_years)%in%years)]])[]
        envValues[,"population_density"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        envValues[,"population_density_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="mean")+1)
        buffer2[[k]] = envValues
      }
      buffer1[[j]] = buffer2
    }
    envValues_list[[i]] = buffer1; saveRDS(envValues_list, "brt_projections_outputs/EnvValues_future_proj.rds")	
  }
  saveRDS(envValues_list, "brt_projections_outputs/EnvValues_future_proj.rds")
}	else	{
  envValues_list = readRDS("brt_projections_outputs/EnvValues_future_proj.rds")
}
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds"); thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
for (i in 1:length(scenarios))
{
  for (j in 1:length(models))
  {
    for (k in 1:length(periods))
    {
      mat = matrix(nrow=length(nutsM), ncol=nberOfReplicates)
      row.names(mat) = nutsM@data[,"LocationCo"]
      dataframe = as.data.frame(envValues_list[[i]][[j]][[k]][,envVariableNames])
      for (l in 1:length(brt_models))
      {
        n.trees = brt_models[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
        mat[,l] = predict.gbm(brt_models[[l]], newdata=dataframe, n.trees, type, single.tree)
      }
      write.table(mat, paste0("brt_projections_outputs/Future_projections/Future_",tolower(scenarios[i]),"_",periods[k],"_",models[j],".csv"), col.names=F, quote=F, sep=",")
    }
  }
}
for (g in 1:length(approaches))
{
  for (h in 1:length(ratios))
  {
    brt_models = readRDS(paste0("brt_training_data/pseudo_abcence_strategies/BRT_foldFit_approach_",g,"_ratio_",ratios[h],".rds"))
    thresholds_SI = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/SI_&_thresholds.csv"), head=T)
    directory =  paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/Future_projections/")
    dir.create(directory, showWarnings=F)
    for (i in 1:length(scenarios))
    {
      for (j in 1:length(models))
      {
        for (k in 1:length(periods))
        {
          mat = matrix(nrow=length(nutsM), ncol=30)
          row.names(mat) = nutsM@data[,"LocationCo"]
          dataframe = as.data.frame(envValues_list[[i]][[j]][[k]][,envVariableNames])
          for (l in 1:length(brt_models))
          {
            n.trees = brt_models[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
            mat[,l] = predict.gbm(brt_models[[l]], newdata=dataframe, n.trees, type, single.tree)
          }
          fileName = paste0("Future_",tolower(scenarios[i]),"_",periods[k],"_",models[j],".csv")
          write.table(mat, paste0(directory,fileName), col.names=F, quote=F, sep=",")
        }
      }
    }
  }
}
thresholds_SI = read.csv(paste0("brt_training_data/SI_&_thresholds.csv"), head=T)
present = read.csv(paste0("brt_projections_outputs/Current_projections/Proj_2020-2024.csv"), head=F)
brts = list(); pops_scenarios = list()
for (i in 1:length(scenarios))
{
  buffer1 = list()
  for (j in 2:length(periods))
  {
    buffer2 = list()
    for (k in 1:length(models))
    {
      buffer2[[k]] = read.csv(paste0("brt_projections_outputs/Future_projections/Future_",tolower(scenarios[i]),"_",periods[j],"_",models[k],".csv"), head=F)
      buffer2[[k]] = buffer2[[k]][,2:dim(buffer2[[k]])[2]]
    }
    buffer1[[j]] = buffer2
  }
  brts[[i]] = buffer1
}

suffix = "A"; plottingDifference = TRUE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "B"; plottingDifference = FALSE; plottingStDevs = TRUE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "C"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE
suffix = "D"; plottingDifference = FALSE; plottingStDevs = FALSE; plottingLower95CI = FALSE; plottingUpper95CI = FALSE; plottingBinaryMap = FALSE

pdf(paste0("Figure_2",suffix,".pdf"), width=(8/6)*5, height=5.8)
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30", col.axis="gray30", fg="gray30")		
scenario_names = c("SSP1-2.6","SSP3-7.0","SSP5-8.5")
colourScale1 = rev(colorRampPalette(brewer.pal(11,"RdBu"))(121)[11:111])
colourScale2 = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(171)[11:161])
colourScale2[c(1:70,82:151)] = paste0(colourScale2[c(1:70,82:151)],"D9") # "D9" = 85%
colourScale2[(76-5):(76+5)] = paste0(colourScale2[(76-5):(76+5)],"00") # "00" = 0%
colourScale3 = colourScale2; colourScale3[(76-5):(76+5)] = "#E5E5E5"
colourScale2 = rev(colorRampPalette(brewer.pal(11,"PuOr"))(19)[2:18])
colourScale2[9] = paste0(colourScale2[9],"00")
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
    binaries = rep(0, dim(brt)[1])
    for (j in 1:dim(brt)[1])
    {
      means[j] = mean(brt[j,]); stdvs[j] = sd(brt[j,])
      temp = rep(0, dim(brt)[2])
      temp[which(brt[j,]>=thresholds_SI[,"threshold"])] = 1
      binaries[j] = round(mean(temp))
    }
    cols = colourScale[(((means-0)/(1-0))*100)+1]; # print(max(stdvs))
    if (plottingStDevs) cols = colourScale[(((stdvs-0)/(0.25-0))*100)+1]
    if (plottingBinaryMap) cols = c("gray90","#de4327")[binaries+1]
    plot(coasts, lwd=0.4, border="gray30", col=NA); plot(nutsM, col=cols, border=NA, lwd=0.1, add=T)
    mtext("Present time", side=3, line=-1.7, at=1, cex=0.50, col="gray30")
    mtext("(2020-2024)", side=3, line=-2.5, at=1, cex=0.50, col="gray30"); rast = raster(as.matrix(c(0,1)))
  }	else	{
    plot.new()
  }
  for (j in 2:length(periods))
  {
    means1 = rep(NA, length(means)); binaries1 = rep(NA, dim(brt)[1])
    stdvs1 = rep(NA, length(means)) # average of standard deviations by GCM
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
      binaries1[k] = round(mean(binaries2))
    }
    cols = colourScale[(((means1-0)/(1-0))*100)+1]; # print(max(stdvs1))
    if (plottingStDevs) cols = colourScale[(((stdvs1-0)/(0.25-0))*100)+1]
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
if (plottingDifference)
{
  vS = c()
  for (i in 1:length(scenarios))
  {
    for (j in 2:length(periods))
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
    for (j in 2:length(periods))
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

# ------------------------------------------------------------------------------
# V. Esrimating the total number of people living in area ecologically suitable
# ------------------------------------------------------------------------------

scenarios = c("SSP1-2.6","SSP3-7.0","SSP5-8.5"); periods = c("2020-2024","2025-2049","2050-2074","2075-2100")
models = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
approaches = c("random","Ae_density","surveillance","espace"); ratios = c("1-1","1-3","1-6")

  # 1. Generating a figure gathering different set of boxplots -----------------


present = read.csv(paste0("brt_projections_outputs/Current_projections/Proj_2020-2024.csv"), head=F)
brts = list(); pops_scenarios = list()
for (i in 1:length(scenarios))
{
  buffer1 = list()
  for (j in 2:length(periods))
  {
    buffer2 = list()
    for (k in 1:length(models))
    {
      buffer2[[k]] = read.csv(paste0("brt_projections_outputs/Future_projections/Future_",tolower(scenarios[i]),"_",periods[j],"_",models[k],".csv"), head=F)
      buffer2[[k]] = buffer2[[k]][,2:dim(buffer2[[k]])[2]]
    }
    buffer1[[j]] = buffer2
  }
  brts[[i]] = buffer1
}
for (i in 1:length(scenarios))
{
  pops_scenarios[[i]] = read.csv(paste0("environmental_data/Population_density_data/NUTS3_population_",tolower(scenarios2[i]),".csv"), head=T)
}
thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
pop_exposures_1 = matrix(nrow=100, ncol=10); temp = as.matrix(present[,2:dim(present)[2]])
for (i in 1:dim(pop_exposures_1)[1])
{
  indexes = which(temp[,i]>=as.numeric(thresholds_SI[i,"threshold"]))
  pop_exposures_1[i,1] = sum(pops_scenarios[[1]][indexes,"X2024"])
}
pop_exposures_2 = pop_exposures_1 # for when we keep the population at its 2024 level
pop_exposures_3 = pop_exposures_1 # for when we only consider the population changes
for (i in 2:length(periods))
{
  for (j in 1:length(scenarios))
  {
    for (k in 1:dim(pop_exposures_1)[1])
    {
      temp1 = matrix(nrow=dim(nutsM@data)[1], ncol=10)
      temp2 = matrix(nrow=dim(nutsM@data)[1], ncol=1)
      for (l in 1:length(models)) temp1[,l] = brts[[j]][[i]][[l]][,k]
      for (l in 1:dim(temp2)[1]) temp2[l,1] = mean(temp1[l,])
      indexes = which(temp2>=thresholds_SI[k,"threshold"])
      colname = paste0("X",gsub("-",".",periods[i]))
      pop_exposures_1[k,1+((i-2)*length(scenarios))+j] = sum(pops_scenarios[[j]][indexes,colname])
      pop_exposures_2[k,1+((i-2)*length(scenarios))+j] = sum(pops_scenarios[[j]][indexes,"X2024"])
      indexes = which(temp[,k]>=as.numeric(thresholds_SI[k,"threshold"]))
      pop_exposures_3[k,1+((i-2)*length(scenarios))+j] = sum(pops_scenarios[[j]][indexes,colname])
    }
  }
}
# boxplot(cbind(pop_exposures_1,pop_exposures_2,pop_exposures_3))
pop_exposures_list1 = list()
for (j in 1:length(scenarios))
{
  pop_exposures_4 = matrix(nrow=100, ncol=30)
  for (i in 2:length(periods))
  {
    for (l in 1:length(models))
    {
      for (k in 1:dim(pop_exposures_4)[1])
      {
        temp = brts[[j]][[i]][[l]][,k]
        indexes = which(temp>=thresholds_SI[k,"threshold"])
        colname = paste0("X",gsub("-",".",periods[i]))
        pop_exposures_4[k,((i-2)*length(models))+l] = sum(pops_scenarios[[j]][indexes,colname])
      }
    }
  }
  pop_exposures_list1[[j]] = pop_exposures_4
}
# boxplot(pop_exposures_list1[[2]]) # only for SSP3
pop_exposures_list2 = list()
for (g in c(1,2,4))
{
  for (h in 1:1) # only for the first ratio (1:1)
  {
    directory =  paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/Future_projections/")
    thresholds_SI = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/SI_&_thresholds.csv"), head=T)
    present = read.csv(paste0("brt_training_data/pseudo_abcence_strategies/BRT_projections_approach_",g,"_ratio_",ratios[h],"/Proj_2020-2024.csv"), head=F)
    pop_exposures = matrix(nrow=30, ncol=10); temp = as.matrix(present[,2:dim(present)[2]])
    for (i in 1:dim(pop_exposures)[1])
    {
      indexes = which(temp[,i]>=as.numeric(thresholds_SI[i,"threshold"]))
      pop_exposures[i,1] = sum(pops_scenarios[[1]][indexes,"X2024"])
    }
    brts = list(); pops_scenarios = list()
    for (i in 1:length(scenarios))
    {
      buffer1 = list()
      for (j in 2:length(periods))
      {
        buffer2 = list()
        for (k in 1:length(models))
        {
          fileName = paste0("Future_",tolower(scenarios[i]),"_",periods[j],"_",models[k],".csv")
          buffer2[[k]] = read.csv(paste0(directory,fileName), head=F)
          buffer2[[k]] = buffer2[[k]][,2:dim(buffer2[[k]])[2]]
        }
        buffer1[[j]] = buffer2
      }
      brts[[i]] = buffer1
    }
    for (i in 1:length(scenarios))
    {
      pops_scenarios[[i]] = read.csv(paste0("environmental_data/Population_density_data/NUTS3_population_",tolower(scenarios2[i]),".csv"), head=T)
    }
    for (i in 2:length(periods))
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
          pop_exposures[k,1+((i-2)*length(scenarios))+j] = sum(pops_scenarios[[j]][indexes,colname])
        }
      }
    }				
  }
  pop_exposures_list2[[g]] = pop_exposures
}
# boxplot(cbind(pop_exposures_list2[[1]],pop_exposures_list2[[2]],pop_exposures_list2[[3]]))

dev.new(height=8, width=8); ats = c(0:6)*(10^8); labels = paste0(seq(0,600,100),"M")
par(mfrow=c(3,1), oma=c(0.3,0.3,0.5,0.5), mar=c(2.7,1.9,0.5,0.5), lwd=0.4, col="gray30", col.axis="gray30", fg="gray30")
colours1 = c("#4d4d4d","#565f8a","#aa7a9b","#455962"); colours2 = c("#4d4d4d1A","#565f8a4D","#aa7a9b4D","#4559624D") # 30% transparency
cols1 = rep(c(colours1[1],rep(colours1[2:4],3)),3); cols2 = rep(c(colours2[1],rep(colours2[2:4],3)),3)
boxplot(cbind(pop_exposures_1,pop_exposures_2,pop_exposures_3), axes=F, xlim=c(1.3,29.3), ylim=c(0.8*10^8,5.2*10^8), boxwex=0.7, lwd=0.5, medlwd=2, outcex=1, border=cols1, col=cols2)
axis(side=2, lwd=0.4, lwd.tick=0.4, cex.axis=0.90, tck=-0.02, col.axis="gray30", mgp=c(0,0.45,0), at=ats, label=labels)
cols1 = rep(c(colours1[1],rep(colours1[2:4],3)),3); cols2 = rep(c(colours2[1],rep(colours2[2:4],3)),3)
boxplot(cbind(pop_exposures_list2[[1]],pop_exposures_list2[[2]],pop_exposures_list2[[4]]), axes=F, xlim=c(1.3,29.3), ylim=c(0.8*10^8,5.2*10^8), boxwex=0.7, lwd=0.5, medlwd=2, outcex=1, border=cols1, col=cols2)
axis(side=2, lwd=0.4, lwd.tick=0.4, cex.axis=0.90, tck=-0.02, col.axis="gray30", mgp=c(0,0.45,0), at=ats, label=labels)
cols1 = c(rep(colours1[2],10),rep(colours1[3],10),rep(colours1[4],10)); cols2 = c(rep(colours2[2],10),rep(colours2[3],10),rep(colours2[4],10))
boxplot(pop_exposures_list1[[2]], axes=F, xlim=c(1.3,29.3), ylim=c(0.8*10^8,5.2*10^8), boxwex=0.7, lwd=0.5, medlwd=2, outcex=1, border=cols1, col=cols2) # only for SSP3
axis(side=2, lwd=0.4, lwd.tick=0.4, cex.axis=0.90, tck=-0.02, col.axis="gray30", mgp=c(0,0.45,0), at=ats, label=labels)

for (h in 1:1)
{
  if (h == 1) pop_exposures = pop_exposures_1
  if (h == 2) pop_exposures = pop_exposures_2
  if (h == 3) pop_exposures = pop_exposures_3
  for (i in 1:dim(pop_exposures)[2])
  {
    # print(shapiro.test(pop_exposures[,i])$p.value)
  }
  pValues = matrix(0, nrow=(dim(pop_exposures)[2]-1), ncol=4) # Benjamini-Hochberg correction
  colnames(pValues) = c("p-value","rank","new_threshold_B-H","significant_after_B-H")
  for (i in 2:dim(pop_exposures)[2])
  {
    print(round(wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value,5))
    pValues[i-1,1] = wilcox.test(pop_exposures[,1], pop_exposures[,i], alternative="less", paired=T)$p.value
  }
  pValues_adjusted = round(p.adjust(pValues[,"p-value"],method="BH"),4)
  for (i in 2:dim(pop_exposures)[2])
  {
    pValues[i-1,2] = which(pValues[order(pValues[,1]),1]==pValues[i-1,1])
    pValues[i-1,3] = (pValues[i-1,2]/dim(pValues)[1])*0.05
    if (pValues[i-1,1] < pValues[i-1,3]) pValues[i-1,4] = 1
    pValues[i-1,1] = round(pValues[i-1,1],3)
    pValues[i-1,3] = round(pValues[i-1,3],3)
  }
}	# for "pop_exposures_1": only the second paired Wilcoxon lead to an adjusted p-value > 0.05

  
  # 2. Generating a global figure showing the evolution through time -----------

years_1 = 1901:2024; years_2 = 2025:2100; envValues_list1 = list(); envValues_list2 = list()
scenarios_1 = c("obsclim","counterclim"); scenarios_2 = c("ssp126","ssp370","ssp585")
if (!file.exists("EnvValues_1901-2024.rds"))
{
  scenarios = scenarios_1
  for (i in 1:length(scenarios))
  {
    directory1 = paste0("../../Rasters/Environmental_rasters/ISIMIP3a/",scenarios[i],"/20CRv3-ERA5")
    directory2 = "../../Rasters/Environmental_rasters/Land-use/"
    directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"		
    if (scenarios[i] == "obsclim")
    {
      years_1 = 1901:2024
      temperature_1 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_tas_1901_2021_monmean.nc"))-273.15 # conversion to Celcius degrees
      precipitation_1 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_pr_1901_2021_monmean.nc"))*60*60*24 # conversion to kg/m2/day
      relative_humidity_1 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_hurs_1901_2021_monmean.nc"))
      temperature_2 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_tas_2022_2024_monmean.nc"))-273.15 # conversion to Celcius degrees
      precipitation_2 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_pr_2022_2024_monmean.nc"))*60*60*24 # conversion to kg/m2/day
      relative_humidity_2 = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_hurs_2022_2024_monmean.nc"))
      temperature = brick(stack(temperature_1, temperature_2))
      precipitation = brick(stack(precipitation_1, precipitation_2))
      relative_humidity = brick(stack(relative_humidity_1, relative_humidity_2))				
    }
    if (scenarios[i] == "counterclim")
    {
      years_1 = 1901:2021
      temperature = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_tas_1901_2021_monmean.nc"))-273.15 # conversion to Celcius degrees
      precipitation = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_pr_1901_2021_monmean.nc"))*60*60*24 # conversion to kg/m2/day
      relative_humidity = brick(paste0(directory1,"/20crv3-era5_",tolower(scenarios[i]),"_hurs_1901_2021_monmean.nc"))
    }
    landCoverVariableNames = as.character(read.csv(paste0(directory2,"/Luse.csv"))[1:12,2])
    landCoverVariableIDs = as.character(read.csv(paste0(directory2,"/Luse.csv"))[1:12,1])
    land_covers1 = list(); land_covers2 = list()
    for (j in 1:length(landCoverVariableIDs))
    {
      land_covers1[[j]] = brick(paste0(directory2,"LUH2_GCB2024_annual_850_2024.nc"), varname=landCoverVariableIDs[j])
    }
    population_density = brick(paste0(directory3,"/Population_density_05deg_1900_2025_historical.nc"), varname="Population_density")
    population_density_log = population_density; population_density_log[] = log(population_density[]+1); pop_years = names(population_density)
    for (j in 1:length(pop_years)) pop_years[j] = gsub("X","",unlist(strsplit(pop_years[j],"\\."))[1])
    cell_areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")/(1000*1000)
    population = population_density*cell_areas # sum(population[[126]][], na.rm=T) # 8236298503 --> ok
    buffer = list()
    for (j in 1:length(years_1))
    {		
      year1 = years_1[j]; year2 = years_1[j]; years = year1:year2
      envValues = matrix(nrow=length(nutsM), ncol=length(allEnvVariableNames)); colnames(envValues) = allEnvVariableNames
      indexes = c()
      for (k in 1:length(years)) indexes = c(indexes, which(grepl(years[k],names(temperature))))
      rast = mean(temperature[[indexes]])
      envValues[,"temperature"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      indexes = c()
      for (k in 1:length(years)) indexes = c(indexes, which(grepl(years[k],names(precipitation))))
      rast = mean(precipitation[[indexes]])
      envValues[,"precipitation"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      indexes = c()
      for (k in 1:length(years)) indexes = c(indexes, which(grepl(years[k],names(relative_humidity))))
      rast = mean(relative_humidity[[indexes]])
      envValues[,"relative_humidity"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      for (k in 1:length(seasons))
      {
        indexes = c()
        for (l in 1:length(years))
        {
          for (m in 1:length(season_months[[k]]))
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(temperature))))
          }
        }
        rast = mean(temperature[[indexes]])
        envValues[,paste0("temperature_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years))
        {
          for (m in 1:length(season_months[[k]]))
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(precipitation))))
          }
        }
        rast = mean(precipitation[[indexes]])
        envValues[,paste0("precipitation_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years))
        {
          for (m in 1:length(season_months[[k]]))
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[k]][m],".",season_months[[k]][m]),names(relative_humidity))))
          }
        }
        rast = mean(relative_humidity[[indexes]])
        envValues[,paste0("relative_humidity_",seasons[k])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      }	
      rasts = list()
      for (k in 1:length(years))
      {
        indexes = c()
        for (l in 1:length(season_months[[3]])) # summer
        {
          indexes = c(indexes, which(grepl(paste0(years[k]+season_years[[3]][l],".",season_months[[3]][l]),names(temperature))))
        }
        rast_summer = mean(temperature[[indexes]])
        indexes = c()
        for (l in 1:length(season_months[[1]])) # winter
        {
          indexes = c(indexes, which(grepl(paste0(years[k]+season_years[[1]][l],".",season_months[[1]][l]),names(temperature))))
        }
        rast_winter = mean(temperature[[indexes]])
        rasts[[k]] = rast_summer-rast_winter
      }
      envValues[,"seasonal_temperature_variation"] = exactextractr::exact_extract(mean(stack(rasts)),nutsM,fun="mean")
      variable_codes = c("pastures","rangelands","croplands","urbanAreas","primaryForest","primaryNonF","secondary","secondaryForest","secondaryNonF")
      variable_names = c("pasture","rangeland","crops","urban land","forested primary land","non-forested primary land","secondary",
                         "potentially forested secondary land","potentially non-forested secondary land"); land_covers2 = list()
      for (k in 1:length(variable_names))
      {
        names = gsub("\\."," ",landCoverVariableNames); indexes1 = which(landCoverVariableNames==variable_names[k])
        if (length(indexes1) == 0) indexes1 = which(grepl(variable_names[k],names))
        # if (variable_names[k] == "pasture") indexes1 = c(indexes1, which(grepl("rangeland",names)))
        indexes2 = which((850+as.numeric(gsub("X","",names(land_covers1[[indexes1[1]]]))))%in%years); rasts = list()
        for (l in 1:length(indexes2))
        {
          rast = land_covers1[[indexes1[1]]][[indexes2[l]]]
          if (length(indexes1) > 1)
          {
            for (m in 2:length(indexes1)) rast[] = rast[]+land_covers1[[indexes1[m]]][[indexes2[l]]][]
          }
          rasts[[l]] = rast
        }
        land_cover = mean(stack(rasts)); names(land_cover) = variable_codes[k]; land_covers2[[k]] = raster::aggregate(land_cover,2)
      }
      envValues[,"primaryForest"] = exactextractr::exact_extract(land_covers2[[5]],nutsM,fun="mean")
      envValues[,"primaryNonForest"] = exactextractr::exact_extract(land_covers2[[6]],nutsM,fun="mean")								
      envValues[,"secondary_areas"] = exactextractr::exact_extract(land_covers2[[7]],nutsM,fun="mean")
      envValues[,"secondaryForest"] = exactextractr::exact_extract(land_covers2[[8]],nutsM,fun="mean")
      envValues[,"secondaryNonForest"] = exactextractr::exact_extract(land_covers2[[9]],nutsM,fun="mean")								
      envValues[,"pastures"] = exactextractr::exact_extract(land_covers2[[1]],nutsM,fun="mean")
      envValues[,"rangelands"] = exactextractr::exact_extract(land_covers2[[2]],nutsM,fun="mean")
      envValues[,"croplands"] = exactextractr::exact_extract(land_covers2[[3]],nutsM,fun="mean")
      envValues[,"urbanAreas"] = exactextractr::exact_extract(land_covers2[[4]],nutsM,fun="mean")
      rast = mean(population[[which((1899+as.numeric(gsub("layer\\.","",names(population))))%in%years)]])
      envValues[,"total_population"] = exactextractr::exact_extract(rast,nutsM,fun="sum")
      envValues[,"total_population_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="sum")+1)
      rast[] = mean(population_density[[which(as.numeric(pop_years)%in%years)]])[]
      envValues[,"population_density"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
      envValues[,"population_density_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="mean")+1)
      buffer[[j]] = envValues	
    }
    envValues_list1[[i]] = buffer
  }
  saveRDS(envValues_list1, "EnvValues_1901-2024.rds")
}	else	{
  envValues_list1 = readRDS("EnvValues_1901-2024.rds")
}
models = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
if (!file.exists("EnvValues_2025-2100.rds"))
{
  scenarios = scenarios_2
  directory3 = "../../Rasters/Environmental_rasters/Population/Population_density_1900_2100/"; envValues_list = list()
  population_density_1 = brick(paste0(directory3,"/Population_density_05deg_1900_2025_historical.nc"), varname="Population_density")
  for (i in 1:length(scenarios))
  {
    buffer1 = list()
    directory2 = paste0("../../Rasters/Environmental_rasters/Land-use/Future_landuse_annual/",scenarios[i])
    landCoverVariableNames = as.character(read.csv("../../Rasters/Environmental_rasters/Land-use/LC_vars.csv")[1:12,2])
    landCoverVariableIDs = as.character(read.csv("../../Rasters/Environmental_rasters/Land-use/LC_vars.csv")[1:12,1])
    land_covers1 = list(); land_covers2 = list()
    for (j in 1:length(landCoverVariableIDs))
    {
      file = list.files(directory2); land_covers1[[j]] = brick(paste0(directory2,"/",file), varname=landCoverVariableIDs[j])
    }
    population_density_2 = brick(paste0(directory3,"/Population_density_05deg_2026_2100_",unlist(strsplit(scenarios[i],"-"))[1],".nc"), varname="Population_density")
    population_density = brick(stack(population_density_1, population_density_2))
    population_density_log = population_density; population_density_log[] = log(population_density[]+1); pop_years = names(population_density)
    for (j in 1:length(pop_years)) pop_years[j] = gsub("X","",unlist(strsplit(pop_years[j],"\\."))[1])
    cell_areas = brick("../../Rasters/Environmental_rasters/clm45_area.nc4")/(1000*1000)
    population = population_density*cell_areas # sum(population[[126]][], na.rm=T) # 8236298503 --> ok
    for (j in 1:length(models))
    {
      print(c(i,j))
      buffer2 = list()
      directory1 = paste0("../../Rasters/Environmental_rasters/ISIMIP3b/",scenarios[i],"/",models[j])
      files = list.files(directory1); files = files[which(grepl("tas",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) temperature = brick(paste0(directory1,"/",files[k]))-273.15 # conversion to Celcius degrees
        if (k != 1) temperature = brick(stack(temperature, brick(paste0(directory1,"/",files[k]))-273.15))
      }
      files = list.files(directory1); files = files[which(grepl("pr",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) precipitation = brick(paste0(directory1,"/",files[k]))*60*60*24 # conversion to kg/m2/day
        if (k != 1) precipitation = brick(stack(precipitation, brick(paste0(directory1,"/",files[k]))*60*60*24))
      }
      files = list.files(directory1); files = files[which(grepl("hurs",files))]; r = NULL
      for (k in 1:length(files))
      {
        if (k == 1) relative_humidity = brick(paste0(directory1,"/",files[k]))
        if (k != 1) relative_humidity = brick(stack(relative_humidity, brick(paste0(directory1,"/",files[k]))))
      }
      for (k in 1:length(years_2))
      {
        year1 = years_2[k]; year2 = years_2[k]; years = year1:year2
        envValues = matrix(nrow=length(nutsM), ncol=length(allEnvVariableNames)); colnames(envValues) = allEnvVariableNames
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(temperature))))
        rast = mean(temperature[[indexes]])
        envValues[,"temperature"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(precipitation))))
        rast = mean(precipitation[[indexes]])
        envValues[,"precipitation"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        indexes = c()
        for (l in 1:length(years)) indexes = c(indexes, which(grepl(years[l],names(relative_humidity))))
        rast = mean(relative_humidity[[indexes]])
        envValues[,"relative_humidity"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        for (l in 1:length(seasons))
        {
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(temperature))))
            }
          }
          rast = mean(temperature[[indexes]])
          envValues[,paste0("temperature_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(precipitation))))
            }
          }
          rast = mean(precipitation[[indexes]])
          envValues[,paste0("precipitation_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
          indexes = c()
          for (m in 1:length(years))
          {
            for (n in 1:length(season_months[[l]]))
            {
              indexes = c(indexes, which(grepl(paste0(years[m]+season_years[[l]][n],".",season_months[[l]][n]),names(relative_humidity))))
            }
          }
          rast = mean(relative_humidity[[indexes]])
          envValues[,paste0("relative_humidity_",seasons[l])] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        }
        rasts = list()
        for (l in 1:length(years))
        {
          indexes = c()
          for (m in 1:length(season_months[[3]])) # summer
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[3]][m],".",season_months[[3]][m]),names(temperature))))
          }
          rast_summer = mean(temperature[[indexes]])
          indexes = c()
          for (m in 1:length(season_months[[1]])) # winter
          {
            indexes = c(indexes, which(grepl(paste0(years[l]+season_years[[1]][m],".",season_months[[1]][m]),names(temperature))))
          }
          rast_winter = mean(temperature[[indexes]])
          rasts[[l]] = rast_summer-rast_winter
        }
        envValues[,"seasonal_variation"] = exactextractr::exact_extract(mean(stack(rasts)),nutsM,fun="mean")
        variable_codes = c("pastures","rangelands","croplands","urbanAreas","primaryForest","primaryNonF","secondary","secondaryForest","secondaryNonF")
        variable_names = c("pasture","rangeland","crops","urban land","forested primary land","non-forested primary land","secondary",
                           "potentially forested secondary land","potentially non-forested secondary land"); land_covers2 = list()
        for (l in 1:length(variable_names))
        {
          names = gsub("\\."," ",landCoverVariableNames); indexes1 = which(landCoverVariableNames==variable_names[l])
          if (length(indexes1) == 0) indexes1 = which(grepl(variable_names[l],names))
          # if (variable_names[l] == "pasture") indexes1 = c(indexes1, which(grepl("rangeland",names)))
          indexes2 = which((2015+as.numeric(gsub("X","",names(land_covers1[[indexes1[1]]]))))%in%years); rasts = list()
          for (m in 1:length(indexes2))
          {
            rast = land_covers1[[indexes1[1]]][[indexes2[m]]]
            if (length(indexes1) > 1)
            {
              for (n in 2:length(indexes1)) rast[] = rast[]+land_covers1[[indexes1[n]]][[indexes2[m]]][]
            }
            rasts[[m]] = rast
          }
          land_cover = mean(stack(rasts)); names(land_cover) = variable_codes[l]; land_covers2[[l]] = raster::aggregate(land_cover,2)
        }
        envValues[,"primaryForest"] = exactextractr::exact_extract(land_covers2[[5]],nutsM,fun="mean")
        envValues[,"primaryNonForest"] = exactextractr::exact_extract(land_covers2[[6]],nutsM,fun="mean")								
        envValues[,"secondaryForest"] = exactextractr::exact_extract(land_covers2[[7]],nutsM,fun="mean")
        envValues[,"secondary_forest_areas"] = exactextractr::exact_extract(land_covers2[[8]],nutsM,fun="mean")
        envValues[,"secondaryNonForest"] = exactextractr::exact_extract(land_covers2[[9]],nutsM,fun="mean")								
        envValues[,"pastures"] = exactextractr::exact_extract(land_covers2[[1]],nutsM,fun="mean")
        envValues[,"rangelands"] = exactextractr::exact_extract(land_covers2[[2]],nutsM,fun="mean")
        envValues[,"croplands"] = exactextractr::exact_extract(land_covers2[[3]],nutsM,fun="mean")
        envValues[,"urbanAreas"] = exactextractr::exact_extract(land_covers2[[4]],nutsM,fun="mean")
        rast = mean(population[[which((1899+as.numeric(gsub("layer\\.","",names(population))))%in%years)]])
        envValues[,"total_population"] = exactextractr::exact_extract(rast,nutsM,fun="sum")
        envValues[,"total_population_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="sum")+1)
        rast[] = mean(population_density[[which(as.numeric(pop_years)%in%years)]])[]
        envValues[,"population_density"] = exactextractr::exact_extract(rast,nutsM,fun="mean")
        envValues[,"population_density_log"] = log(exactextractr::exact_extract(rast,nutsM,fun="mean")+1)
        buffer2[[k]] = envValues
      }
      buffer1[[j]] = buffer2
    }
    envValues_list2[[i]] = buffer1; saveRDS(envValues_list2, "EnvValues_2025-2100.rds")
  }
  saveRDS(envValues_list2, "EnvValues_2025-2100.rds")
}	else	{
  envValues_list2 = readRDS("EnvValues_2025-2100.rds")
}
brt_models = readRDS("brt_training_data/FoldFit_BRT_models.rds"); thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
envVariableNames = c("temperature_winter","temperature_fall","precipitation_spring","precipitation_fall",
                     "relative_humidity_winter","relative_humidity_fall","pastures","urbanAreas","primaryNonForest")
population_at_risk_1 = list(); population_at_risk_2 = list()
if (!file.exists("brt_projections_outputs/PopAtRisk_1901-2024.rds"))
{
  for (i in 1:length(scenarios_1))
  {
    if (scenarios_1[i] == "obsclim") years_1 = 1901:2024
    if (scenarios_1[i] == "counterclim") years_1 = 1901:2021
    mat = matrix(nrow=length(years_1), ncol=length(brt_models))
    for (k in 1:length(brt_models))
    {
      for (j in 1:length(years_1))
      {
        dataframe = as.data.frame(envValues_list1[[i]][[j]][,envVariableNames])
        n.trees = brt_models[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
        temp = predict.gbm(brt_models[[k]], newdata=dataframe, n.trees, type, single.tree)
        indexes = which(temp>=as.numeric(thresholds_SI[k,"threshold"]))
        mat[j,k] = sum(envValues_list1[[i]][[j]][indexes,"total_population"])
      }				
    }
    row.names(mat) = years_1; population_at_risk_1[[i]] = mat
  }
  saveRDS(population_at_risk_1, "brt_projections_outputs/PopAtRisk_1901-2024.rds")
}	else	{
  population_at_risk_1 = readRDS("brt_projections_outputs/PopAtRisk_1901-2024.rds")
}
if (!file.exists("brt_projections_outputs/PopAtRisk_2025-2100.rds"))
{
  for (i in 1:length(scenarios_2))
  {
    mat = matrix(nrow=length(years_2), ncol=length(models)*length(brt_models))
    for (k in 1:length(models))
    {
      for (l in 1:length(brt_models))
      {
        for (j in 1:length(years_2))
        {
          dataframe = as.data.frame(envValues_list2[[i]][[k]][[j]][,envVariableNames])
          n.trees = brt_models[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
          temp = predict.gbm(brt_models[[l]], newdata=dataframe, n.trees, type, single.tree)
          indexes = which(temp>=as.numeric(thresholds_SI[l,"threshold"]))
          mat[j,(((k-1)*length(brt_models))+l)] = sum(envValues_list2[[i]][[k]][[j]][indexes,"total_population"])
        }				
      }
    }
    row.names(mat) = years_2; population_at_risk_2[[i]] = mat
  }
  saveRDS(population_at_risk_2, "brt_projections_outputs/PopAtRisk_2025-2100.rds")
}	else	{
  population_at_risk_2 = readRDS("brt_projections_outputs/PopAtRisk_2025-2100.rds")
}
population_at_risk_1_sd3y = population_at_risk_1; population_at_risk_2_sd3y = population_at_risk_2
for (i in 1:length(population_at_risk_1))
{
  for (j in 2:(dim(population_at_risk_1[[i]])[1]-1))
  {
    for (k in 1:dim(population_at_risk_1[[i]])[2])
    {
      population_at_risk_1_sd3y[[i]][j,k] = mean(population_at_risk_1[[i]][(j-1):(j+1),k])
    }
  }
}
for (i in 1:length(population_at_risk_2))
{
  for (j in 2:(dim(population_at_risk_2[[i]])[1]-1))
  {
    for (k in 1:dim(population_at_risk_2[[i]])[2])
    {
      population_at_risk_2_sd3y[[i]][j,k] = mean(population_at_risk_2[[i]][(j-1):(j+1),k])
    }
  }
}
population_at_risk_1_sd5y = population_at_risk_1; population_at_risk_2_sd5y = population_at_risk_2
for (i in 1:length(population_at_risk_1))
{
  for (j in 3:(dim(population_at_risk_1[[i]])[1]-2))
  {
    for (k in 1:dim(population_at_risk_1[[i]])[2])
    {
      population_at_risk_1_sd5y[[i]][j,k] = mean(population_at_risk_1[[i]][(j-2):(j+2),k])
    }
  }
}
for (i in 1:length(population_at_risk_2))
{
  for (j in 3:(dim(population_at_risk_2[[i]])[1]-2))
  {
    for (k in 1:dim(population_at_risk_2[[i]])[2])
    {
      population_at_risk_2_sd5y[[i]][j,k] = mean(population_at_risk_2[[i]][(j-2):(j+2),k])
    }
  }
}

dev.new(height=2.5, width=8)
par(oma=c(0.0,0.3,0.5,0.5), mar=c(1.6,4.1,0.5,0.5), lwd=0.4, col="gray30", col.axis="gray30", fg="gray30")
xS_list = list(); yS_list = list(); cols1 = list(); cols2 = list()
cols1[[1]] = c("#447ba2","#747c80"); cols1[[2]] = c("#565f8a","#aa7a9b","#455962")
cols2[[1]] = c("#447ba21A","#747c801A"); cols2[[2]] = c("#565f8a1A","#aa7a9b1A","#4559621A") # 10% transparency
for (i in 1:length(scenarios_1))
{
  xS = row.names(population_at_risk_1_sd3y[[i]]); yS = matrix(nrow=length(xS), ncol=3)
  for (j in 1:dim(yS)[1])
  {
    yS[j,1] = median(population_at_risk_1_sd3y[[i]][j,]); yS[j,2:3] = quantile(population_at_risk_1_sd3y[[i]][j,], c(0.05,0.95))
  }
  xS_list[[i]] = xS; yS_list[[i]] = yS
}
for (i in 1:length(scenarios_1))
{
  indexes = c(2:(length(xS_list[[i]])-1))
  if (i == 1)
  {
    plot(xS_list[[i]][indexes], yS_list[[i]][indexes,1], type="n", ann=F, axes=F, lwd=0.1, xlim=c(1908,2094), ylim=c(11530000,500000000))
    axis(side=1, lwd.tick=0.4, cex.axis=0.6, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,-0.02,0), at=seq(1925,2075,25))
    axis(side=2, lwd.tick=0.4, cex.axis=0.6, lwd=0, tck=-0.025, col.axis="gray30", mgp=c(0,0.45,0), las=2, padj=0.4, label=c("0","100M","200M","300M","400M","500M"), at=seq(0,5*10^8,10^8))
    title(ylab="Estimated European population living in areas ecologically", cex.lab=0.7, mgp=c(2.8,0,0), col.lab="gray30")
    title(ylab=expression("for the local circulation of "*italic(paste("Aedes"))*"-transmitted viruses"), cex.lab=0.7, mgp=c(2.0,0,0), col.lab="gray30")
    grid(lwd=0.75, lty=3)
  }
  polygon(c(xS_list[[i]][indexes], rev(xS_list[[i]][indexes])), c(yS_list[[i]][indexes,2], rev(yS_list[[i]][indexes,3])), col=cols2[[1]][i], border=NA)
}
for (i in 1:length(scenarios_1))
{
  indexes = c(2:(length(xS_list[[i]])-1))
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,2], col=cols1[[1]][i], lwd=0.3)
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,3], col=cols1[[1]][i], lwd=0.3)
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,1], col=cols1[[1]][i], lwd=2.0)
}
for (i in 1:length(scenarios_2))
{
  xS = row.names(population_at_risk_2_sd3y[[i]]); yS = matrix(nrow=length(xS), ncol=3)
  for (j in 1:dim(yS)[1])
  {
    yS[j,1] = median(population_at_risk_2_sd3y[[i]][j,]); yS[j,2:3] = quantile(population_at_risk_2_sd3y[[i]][j,], c(0.05,0.95))
  }
  xS_list[[i]] = xS; yS_list[[i]] = yS
}
indexes = c(2:(length(xS_list[[i]])-1))
for (i in 1:length(scenarios_2))
{
  polygon(c(xS_list[[i]][indexes], rev(xS_list[[i]][indexes])), c(yS_list[[i]][indexes,2], rev(yS_list[[i]][indexes,3])), col=cols2[[2]][i], border=NA)
}
for (i in 1:length(scenarios_2))
{
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,2], col=cols1[[2]][i], lwd=0.3)
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,3], col=cols1[[2]][i], lwd=0.3)
  lines(xS_list[[i]][indexes], yS_list[[i]][indexes,1], col=cols1[[2]][i], lwd=2.0)
}
box(lwd=0.4, col="gray30")

# ------------------------------------------------------------------------------
# VI. Estimating the ecological suitability variations at the continent scale
# ------------------------------------------------------------------------------

indices = which(nutsM@data[,"CountryNam"]=="Belgium") # estimating the variation for Belgium
indices = which(nutsM@data[,"CountryNam"]%in%c("Belgium","Netherlands")) # estimating the variation for Belgium and The Netherlands
indices = 1:dim(nutsM@data)[1] # when considering all the countries

models = c("CanESM5","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
buffer_ssp3 = list(); buffer_ssp5 = list(); brts_ssp3 = list(); brts_ssp5 = list(); pops_ssp3 = list(); pops_ssp5 = list()
brts_ssp3[[1]] = read.csv(paste0("brt_projections_outputs/Proj_2020-2024.csv"), head=F)
brts_ssp3[[1]] = brts_ssp3[[1]][indices,2:dim(brts_ssp3[[1]])[2]]; brts_ssp5[[1]] = brts_ssp3[[1]]
for (i in 1:length(models))
{
  buffer_ssp3[[i]] = read.csv(paste0("brt_projections_outputs/Future_projections/Future_ssp3-7.0_2075-2100_",models[i],".csv"), head=F)
  buffer_ssp5[[i]] = read.csv(paste0("brt_projections_outputs/Future_projections/Future_ssp5-8.5_2075-2100_",models[i],".csv"), head=F)
  buffer_ssp3[[i]] = buffer_ssp3[[i]][indices,2:dim(buffer_ssp3[[i]])[2]]; buffer_ssp5[[i]] = buffer_ssp5[[i]][indices,2:dim(buffer_ssp5[[i]])[2]]
}
brts_ssp3[[2]] = buffer_ssp3; brts_ssp5[[2]] = buffer_ssp5
dS2_ssp3 = rep(NA, dim(brts_ssp3[[1]])[2]) # first averaging across models and then replicates
dS2_ssp5 = rep(NA, dim(brts_ssp5[[1]])[2]) # first averaging across models and then replicates
for (i in 1:dim(brts_ssp3[[1]])[2])
{
  dS_ssp3 = rep(NA, length(models))
  for (j in 1:length(models))
  {
    dS_ssp3[j] = ((sum(brts_ssp3[[2]][[j]][,i])-sum(brts_ssp3[[1]][,i]))/sum(brts_ssp3[[1]][,i]))*100
  }
  dS2_ssp3[i] = mean(dS_ssp3)
}
for (i in 1:dim(brts_ssp5[[1]])[2])
{
  dS_ssp5 = rep(NA, length(models))
  for (j in 1:length(models))
  {
    dS_ssp5[j] = ((sum(brts_ssp5[[2]][[j]][,i])-sum(brts_ssp5[[1]][,i]))/sum(brts_ssp5[[1]][,i]))*100
  }
  dS2_ssp5[i] = mean(dS_ssp5)
}
meanV = round(mean(dS2_ssp3)); percentiles = round(quantile(dS2_ssp3, c(0.025,0.975)))
cat("\tdS2 (SSP3) = ",meanV,"%, 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# for all countries: dS2 (SSP3) = 14%, 2.5th-97.5th percentiles range = [0, 28]
# for Belgium only: dS2 (SSP3) = 125%, 2.5th-97.5th percentiles range = [45, 221]
# --> the ecological suitability could more than double before the end of this century
# for Belgium and the Netherlands: dS2 (SSP3) = 88%, 2.5th-97.5th percentiles range = [32, 169]
meanV = round(mean(dS2_ssp5)); percentiles = round(quantile(dS2_ssp5, c(0.025,0.975)))
cat("\tdS2 (SSP5) = ",meanV,"%, 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# for all countries: dS2 (SSP5) = 38%, 2.5th-97.5th percentiles range = [19, 59]
# for Belgium only: dS2 (SSP5) = 86%, 2.5th-97.5th percentiles range = [23, 175]
# for Belgium and the Netherlands: dS2 (SSP5) = 81%, 2.5th-97.5th percentiles range = [30, 157]

thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
nuts_areas = raster::area(nutsM)[indices]/(1000*1000) # in km2
suitable_areas_ssp3 = matrix(nrow=100, ncol=2); suitable_areas_ssp5 = matrix(nrow=100, ncol=2)
for (i in 1:dim(suitable_areas_ssp3)[1])
{
  threshold = as.numeric(thresholds_SI[i,"threshold"])
  suitable_areas_ssp3[i,1] = sum(nuts_areas[which(brts_ssp3[[1]][,i]>=threshold)]); buffer_ssp3 = rep(NA, length(models))
  suitable_areas_ssp5[i,1] = sum(nuts_areas[which(brts_ssp5[[1]][,i]>=threshold)]); buffer_ssp5 = rep(NA, length(models))
  for (j in 1:length(models))
  {
    buffer_ssp3[j] = sum(nuts_areas[which(brts_ssp3[[2]][[j]][,i]>=threshold)])
    buffer_ssp5[j] = sum(nuts_areas[which(brts_ssp5[[2]][[j]][,i]>=threshold)])
  }
  suitable_areas_ssp3[i,2] = mean(buffer_ssp3); suitable_areas_ssp5[i,2] = mean(buffer_ssp5)
}
dA_ssp3 = ((suitable_areas_ssp3[,2]-suitable_areas_ssp3[,1])/suitable_areas_ssp3[,1])*100
dA_ssp5 = ((suitable_areas_ssp5[,2]-suitable_areas_ssp5[,1])/suitable_areas_ssp5[,1])*100
meanV = round(mean(dA_ssp3)); percentiles = round(quantile(dA_ssp3, c(0.025,0.975)))
cat("\tdA (SSP3) = ",meanV,"%, 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# for all countries: dA (SSP3) = 27%, 2.5th-97.5th percentiles range = [-11, 57]
# for Belgium and the Netherlands (based on mean across replicates): 20211 km2 --> 45206 km2
# for Belgium (based on mean across replicates): 7863 km2 --> 24067 km2
meanV = round(mean(dA_ssp5)); percentiles = round(quantile(dA_ssp5, c(0.025,0.975)))
cat("\tdA (SSP5) = ",meanV,"%, 2.5th-97.5th percentiles range = [",percentiles[1],", ",percentiles[2],"]\n",sep="")
# for all countries: dA (SSP5) = 82%, 2.5th-97.5th percentiles range = [28, 155]
# for Belgium and the Netherlands (based on mean across replicates): 20211 km2 --> 47550 km2
# for Belgium (based on mean across replicates): 7863 km2 --> 21766 km2

thresholds_SI = read.csv("brt_training_data/SI_&_thresholds.csv", head=T)
pops_ssp3 = read.csv("environmental_data/Population_density_data/NUTS3_population_ssp370.csv")[indices,]
pops_ssp5 = read.csv("environmental_data/Population_densities_data/NUTS3_population_ssp585.csv")[indices,]
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
meanV = round(mean(pop_exposures_ssp3[,1])); percentiles = round(quantile(pop_exposures_ssp3[,1], c(0.025,0.975)))
cat("n (2024) = ",meanV,"; 2.5th-97.5th percentiles range = [",percentiles[1],"; ",percentiles[2],"]",sep="")
# n (2024) = 175,546,357; 2.5th-97.5th percentiles range = [93,212,610; 322,608,194]
meanV = round(mean(pop_exposures_ssp3[,2])); percentiles = round(quantile(pop_exposures_ssp3[,2], c(0.025,0.975)))
cat("n (SSP3) = ",meanV,"; 2.5th-97.5th percentiles range = [",percentiles[1],"; ",percentiles[2],"]",sep="")
# n (SSP3) = 244,661,541; 2.5th-97.5th percentiles range = [105,485,111; 420,120,632]
dS_ssp3 = pop_exposures_ssp3[,2]-pop_exposures_ssp3[,1]; dS_ssp5 = pop_exposures_ssp5[,2]-pop_exposures_ssp5[,1]
meanV = round(mean(dS_ssp3)); qS = round(quantile(dS_ssp3,c(0.025,0.975)),1); percentiles = round(quantile(dS_ssp3, c(0.025,0.975)))
cat("dS (SSP3) = ",meanV,"; 2.5th-97.5th percentiles range = [",percentiles[1],"; ",percentiles[2],"]",sep="")
# dS (SSP3) = 69,115,185; 2.5th-97.5th percentiles range = [2,113,062; 131,278,667]
meanV = round(mean(pop_exposures_ssp5[,2])); percentiles = round(quantile(pop_exposures_ssp5[,2], c(0.025,0.975)))
cat("n (SSP5) = ",meanV,"; 2.5th-97.5th percentiles range = [",percentiles[1],"; ",percentiles[2],"]",sep="")
# n (SSP5) = 349,885,959; 2.5th-97.5th percentiles range = [151,620,632; 486,896,100]
meanV = round(mean(dS_ssp5)); qS = round(quantile(dS_ssp5,c(0.025,0.975)),1); percentiles = round(quantile(dS_ssp5, c(0.025,0.975)))
cat("dS (SSP5) = ",meanV,"; 2.5th-97.5th percentiles range = [",percentiles[1],"; ",percentiles[2],"]",sep="")
# dS (SSP5) = 174,339,602; 2.5th-97.5th percentiles range = [48,960,089; 269,808,390]


