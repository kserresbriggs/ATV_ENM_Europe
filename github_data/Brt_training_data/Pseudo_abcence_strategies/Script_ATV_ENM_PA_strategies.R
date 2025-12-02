# Geospatial libraries
library(ncdf4)
library(exactextractr)
library(raster)
library(terra)
library(sf)
library(sp)

# Data manipulation libraries
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
library(openxlsx)
library(overlapping)

# Statistical modeling libraries
library(gmodels)
library(devtools)
library(USE)
library(dismo)
library(gbm)
library(blockCV)
library(ecospat)

# Plotting libraries
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(hrbrthemes)
library(scales)

# Parallel processing libraries
library(furrr)
library(future.apply)


################################################################################
#                    I. LOADING ALL DATA
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

# 2. Load all presence records in Europe
########################################
Aedes_presences = read.csv("github_data/Occurrence_data/Aedes_viruses_European_presences.csv", header = T)


################################################################################
#                    II. PSEUDO-ABSENCE SAMPLING METHODS
################################################################################

# 1. Group presence data by years
##########################################

# In order to not give more importance to same areas where presences have been 
# recorded accross multiple years, we will only keep the year range so that the 
# environmental variables can be averaged across that range

# create the funtion
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

# 2. Create pseudo-absence shapefile
####################################

# For each sampling strategy, the adjacent polygons to presence sites are removed

# find adjacent polygons to each disease presence location
ae.pres_geoms = which(nutsM_sf$LocationCo %in% ae.presences_agg$nutM_code)
ae.pres_geoms = nutsM_sf[ae.pres_geoms, ]
ae.adjacent = st_intersects(ae.pres_geoms$geometry, nutsM_sf$geometry)
ae.adj = unlist(ae.adjacent)
ae.adj_geoms = nutsM_sf[ae.adj, ]


# removing adjacent polygons from the possible pseudo-absence areas
absences = setdiff(nutsM_IDs$nutM_code, ae.presences_agg) # select all NUTM areas that are not in presences
Aedes_abs = data.frame(
  nutM_code = absences, 
  year = NA
)

remove_rows = which(Aedes_abs$nutM_code %in% ae.adj_geoms$LocationCo) 
Aedes_abs = Aedes_abs[(-remove_rows), ] # remove adjacent polygons from the possible pseudo-absences


# Pseudo-absences must be assigned to a specific year,
# as the model is trained temporally. Sampled years are
# weighted according to the frequency of presence years.

# get presence outbreak dates (year of infection)
Aedes_presences.unique = Aedes_presences[!duplicated(Aedes_presences[c("year", "nutM_code")]), ] # remove duplicated presences
presence_years = Aedes_presences.unique$year

year_counts = table(presence_years)
year_probs = year_counts / sum(year_counts)
years = as.integer(names(year_probs))



# 3. All pseudo-absence sampling methods
########################################

## 3A. Random PA sampling strategy ##
# Pseudo-absence sampling strategies to be evaluated within the model training loop
random_sampling = quote({
  indices1 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices1, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F) # set ratio (1:1, 1:3, 1:6)
  absences = Aedes_abs[selected_absences, ]
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign pseudo-years
  absences$year = sample(years, size = length(absences$year), replace = TRUE, prob = year_probs)
  
  # Add response
  absences$response = 0
  
  # Combine with presences
  data = rbind(ae.presences_agg, absences)
})


## 3B. Density-weighted geographic PA sampling strategy ##
# the Aedes Albopictus ecological suitability raster is used to 
# sample PAs in geographic regions where there the vector is abundant

# Load GeoTIFF aedes albopictus ecological suitability (M.Kraemer)
geo_density_sampling = quote({
  aedes = raster("Rasters/albopictus_suitability_Kraemer2015.tif")
  aedes_cropped = crop(aedes, nutsM) # crop to Europe
  
  # extract suitability values per nutM area
  aedes_nutsM = exact_extract(aedes_cropped, nutsM, fun ="sum")
  nutsM_data = nutsM_IDs
  nutsM_data$aedes_suitability = aedes_nutsM
  
  # Compute weighted probabilities
  Aedes_abs$aedes_suitability = NA
  Aedes_abs$probabilities = NA
  for (a in seq_len(nrow(Aedes_abs))) {
    index = which(nutsM_data$nutM_code == Aedes_abs[a, "nutM_code"])
    Aedes_abs[a, "aedes_suitability"] = nutsM_data[index, "aedes_suitability"]
  }
  Aedes_abs$probabilities = Aedes_abs$aedes_suitability / sum(Aedes_abs$aedes_suitability)
  
  # the following code is to be in introduced in large loop to sample different PA
  # for each replicate.
  indices2 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices2, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$probabilities) # set ratio (1:1, 1:3, 1:6)
  absences = Aedes_abs[selected_absences, ]
  
  absences$year = NA_integer_
  absences$year_start = NA_integer_
  absences$year_end = NA_integer_
  
  # Assign pseudo-years
  absences$year = sample(years, size = length(absences$year), replace = TRUE, prob = year_probs)
  
  # Add response
  absences$response = 0
  absences = absences[, !(colnames(absences) %in% c("aedes_suitability", "probabilities"))]
  Aedes_abs = Aedes_abs[, !(colnames(Aedes_abs) %in% c("aedes_suitability", "probabilities"))]
  
  # Combine with presences
  data = rbind(ae.presences_agg, absences)  
})

## 3C. Surveillance raster PA sampling strategy ##
# the surveillance suitability raster is used to 
# sample PAs in geographic regions where there arbovirus surveillance is robust and where presences have not been reported

# Load GeoTIFF surveillance raster Lim et al 2025
surveillance_sampling = quote({
  rast_surv = raster("Rasters/Lim_et_al_2025_data/outputs/Rasters/Surveillance_map_wmean.tif")
  rast_surv_cropped = crop(rast_surv, nutsM) # crop to Europe
  
  # extract suitability values per nutM area
  surv_nutsM = exact_extract(rast_surv_cropped, nutsM, fun="mean")
  surv.nutsM_data = nutsM_IDs
  surv.nutsM_data$surv = surv_nutsM 
  surv.nutsM_data$surv[is.na(surv.nutsM_data$surv)] = 0
  
  # Compute weighted probabilities
  Aedes_abs$surv = NA
  Aedes_abs$Probabilities = NA
  for (m in seq_len(nrow(Aedes_abs))) {
    index = which(surv.nutsM_data$nutM_code == Aedes_abs[m, "nutM_code"])
    Aedes_abs[m, "surv"] = surv.nutsM_data[index, "surv"]
  }
  Aedes_abs$Probabilities = Aedes_abs$surv / sum(Aedes_abs$surv)
  
  # the following code is to be in introduced in large loop to sample different PA
  # for each replicate.
  indices3 = seq_len(nrow(Aedes_abs))
  selected_absences = sample(indices3, ratios[[r]] * length(ae.presences_agg$nutM_code), replace = F, prob = Aedes_abs$Probabilities) # set ratio (1:1, 1:3, 1:6)
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
  
  data = rbind(ae.presences_agg, absences)
})


## 3D. Environmental space PA sampling strategy ##
ES_setup = quote({
  # Load the reference raster (e.g., one of the ISIMIP climatic layers)
  reference = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  reference = reference[[1]] # select the first layer
  reference = crop(reference, nutsM, snap="out") # crop to europe
  reference = mask(reference, nutsM)
  values(reference) = NA # Reference raster becomes NULL raster
  
  # duplicate reference raster
  rast_remove = reference
  
  # Use exact_extract to get cell indices for all presence geometries
  # All cells that are within a presence polygon are assigned a value of 1
  presence_indices = exact_extract(reference, ae.pres_geoms, include_cell = TRUE)
  indices = list()
  for (i in 1:length(presence_indices)){Idx = presence_indices[[i]][[2]];indices[[i]] = Idx}
  indices = unlist(indices, use.names = FALSE)
  # set presence cells of null raster to 1
  reference[indices] = 1
  
  # Use exact_extract to get cell indices for all adjacent polygons to presence geometries
  # All cells that are within an adjacent polygon are assigned a value of 1
  # this raster will be used to remove the adjacent polygons
  
  # NA_ids = c("EL304","EL307","EL413","EL421","EL422","EL624","ES531","ES533",
  #            "ES630","ES640","GG","GI","IM","JE","MEG125354","MEG125368","UKM65") # these IDs correspond to small islands in Europe, 
  # # unfortunately the landcover data is unavailable for these areas (remove from possible absences)
  # na_geoms = which(nutsM_sf$LocationCo %in% NA_ids)
  # na_geoms = nutsM_sf[na_geoms, ]
  Ae.adj_geoms = ae.adj_geoms
  
  adj_indices = exact_extract(reference, Ae.adj_geoms, include_cell = TRUE)
  Adj_indices = list()
  for (i in 1:length(adj_indices)){
    idx = adj_indices[[i]][[2]]; Adj_indices[[i]] = idx}
  Adj_indices = unlist(Adj_indices, use.names = FALSE)
  # set adjacent cells of second null raster to 2
  rast_remove[Adj_indices] = 2
  
  # Find cells where both presence and adjacent cells overlap (1 + 2)
  # remove the overlap, the final raster will only have cells corresponding to 
  # adjacent polygons
  overlap_mask = reference == 1 & rast_remove == 2
  rast_remove[overlap_mask] = NA
  
  # Extract data from the reference presence raster
  ref_data = as.data.frame(matrix(nrow=ncell(reference), ncol=4)); colnames(ref_data) = c("indices", "longitude", "latitude", "response")
  ref_data$indices = 1:ncell(reference) # get presence cell indices
  coords = xyFromCell(reference, 1:ncell(reference))
  ref_data$longitude = coords[, 1]; ref_data$latitude = coords[, 2] # get presence cell coordinates
  ref_data$response = getValues(reference) # presence cells (value = 1)
  
  myPres = ref_data[!is.na(ref_data$response), ]  # presences (value=1)
  myPres = st_as_sf(myPres, coords = c("longitude", "latitude")) # convert to sf object
  
  
  
  ## Load climatic rasters for the environmental space pseudo-absence sampling ##
  # load raster layers, subset over the past 20 years and by season 
  year_start = (1999 - 1901) * 12 + 1; year_end = (2024 - 1901 + 1) * 12
  temp1 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"); temp2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_2022_2024_monmean.nc")
  temp = addLayer(temp1, temp2);
  temp = subset(temp, year_start:year_end)
  
  winter_indices = sort(c(12, 1, 2)) + rep((1999:2024 - 1999) * 12, each = 3); winter_indices = winter_indices[-c(1, 2, length(winter_indices))] # December (previous year), January, February
  spring_indices = sort(c(3,4,5) + rep((2000:2024-1999)*12, each= 3))
  summer_indices = sort(c(6,7,8) + rep((2000:2024-1999)*12, each= 3))
  fall_indices = sort(c(9,10,11) + rep((2000:2024-1999)*12, each= 3))
  
  tepmm_winter = mean(temp[[winter_indices]]) - 273.15
  tepmm_spring = mean(temp[[spring_indices]]) - 273.15
  tepmm_summer = mean(temp[[summer_indices]]) - 273.15
  tepmm_fall = mean(temp[[fall_indices]]) - 273.15
  
  precp1 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"); precp2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_2022_2024_monmean.nc")
  precp = addLayer(precp1, precp2);
  precp = subset(precp, year_start:year_end)
  
  precp_winter = mean(precp[[winter_indices]]) * 60 * 60 * 24
  precp_spring = mean(precp[[spring_indices]]) * 60 * 60 * 24
  precp_summer = mean(precp[[summer_indices]]) * 60 * 60 * 24
  precp_fall = mean(precp[[fall_indices]]) * 60 * 60 * 24
  
  relh1 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"); relh2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_2022_2024_monmean.nc")
  relh = addLayer(relh1,relh2)
  relh = subset(relh, year_start:year_end)
  
  relh_winter = mean(relh[[winter_indices]])
  relh_spring = mean(relh[[spring_indices]]) 
  relh_summer = mean(relh[[summer_indices]])
  relh_fall = mean(relh[[fall_indices]])
  
  year_start_temp = 2000 - 1900;  year_end_temp = 2024 - 1900
  crops_temp = brick("Rasters/Environmental_rasters/Land-use/croplands_LUH2-GCB2024_1901_2024.nc"); crops_temp = mean(crops_temp[[year_start_temp:year_end_temp]])
  past_temp = brick("Rasters/Environmental_rasters/Land-use/pastures_LUH2-GCB2024_1901_2024.nc"); past_temp = mean(past_temp[[year_start_temp:year_end_temp]])
  ranges_temp = brick("Rasters/Environmental_rasters/Land-use/rangeland_LUH2-GCB2024_1901_2024.nc"); ranges_temp = mean(ranges_temp[[year_start_temp:year_end_temp]])
  urban_temp = brick("Rasters/Environmental_rasters/Land-use/urbanAreas_LUH2-GCB2024_1901_2024.nc"); urban_temp = mean(urban_temp[[year_start_temp:year_end_temp]])
  primForest_temp = brick("Rasters/Environmental_rasters/Land-use/primaryForest_LUH2-GCB2024_1901_2024.nc"); primForest_temp = mean(primForest_temp[[year_start_temp:year_end_temp]])
  primNonForest_temp = brick("Rasters/Environmental_rasters/Land-use/primaryNonF_LUH2-GCB2024_1901_2024.nc"); primNonForest_temp = mean(primNonForest_temp[[year_start_temp:year_end_temp]])
  secForest_temp = brick("Rasters/Environmental_rasters/Land-use/secondaryForest_LUH2-GCB2024_1901_2024.nc"); secForest_temp = mean(secForest_temp[[year_start_temp:year_end_temp]])
  secNonForest_temp = brick("Rasters/Environmental_rasters/Land-use/secondaryNonF_LUH2-GCB2024_1901_2024.nc"); secNonForest_temp = mean(secNonForest_temp[[year_start_temp:year_end_temp]])
  
  pop = brick("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_1900_2025_historical.nc")
  year_start_temp = 2000 - 1900;  year_end_temp = 2024 - 1900
  pop_temp = mean(pop[[year_start_temp:year_end_temp]]); pop_temp = log10(pop_temp + 1)
  
  
  # group rasters
  envVariables = list()
  envVariables[[1]] = tepmm_winter
  envVariables[[2]] = tepmm_spring
  envVariables[[3]] =  tepmm_summer
  envVariables[[4]] =  tepmm_fall
  envVariables[[5]] =  precp_winter
  envVariables[[6]] =  precp_spring
  envVariables[[7]] =  precp_summer
  envVariables[[8]] =  precp_fall
  envVariables[[9]] =  relh_winter
  envVariables[[10]] =  relh_spring
  envVariables[[11]] =  relh_summer
  envVariables[[12]] =  relh_fall
  envVariables[[13]] = crops_temp
  envVariables[[14]] = past_temp
  envVariables[[15]] = ranges_temp
  envVariables[[16]] = urban_temp
  envVariables[[17]] = primForest_temp
  envVariables[[18]] = primNonForest_temp
  envVariables[[19]] = secForest_temp
  envVariables[[20]] = secNonForest_temp
  envVariables[[21]] = pop_temp
  
  # match extent to europe ans stack
  for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], nutsM, snap="out")
  for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], nutsM)
  variable_stack = stack(envVariables)
  variable_stack[rast_remove] = NA # remove cells corresponding to adjacent polygons from the rasters and thus the environmental space
  
  ## 3. Run the PCA
  PCA = USE::rastPCA(env.rast = variable_stack, nPC = 2, stand = TRUE)
  pca_result = PCA$pca ; pca_summary = summary(pca_result); print(pca_summary) 
  PCA_stack = c(PCA$PCs$PC1, PCA$PCs$PC2)
  
  
  ## 4. Find the optimal grid resolution for pseudo-absence sampling
  PCAstack_df = as.data.frame(PCA_stack, xy = T, na.rm = T)
  PCAstack_sp = st_as_sf(PCAstack_df, coords = c("PC1", "PC2"))
  opt_res = optimRes(sdf = PCAstack_sp, grid.res = seq(1, 15), showOpt = T)
  
  # 5. rasterize the Europe polygon shapefile to associate raster cells to a polygon ID
  rast = brick(paste0("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"))
  rast = rast[[1]] # select the first layer to make the polygon file the same resolution and extent
  rast = crop(rast, nutsM, snap="out")
  rast = mask(rast, nutsM)
  
  nutsM$LocationCo_numeric = as.numeric(as.factor(nutsM$LocationCo)) # Convert polygon IDs to numeric indices
  rasterized_nutsM = rasterize(nutsM, rast, field="LocationCo_numeric")
  
  cell_ids = which(!is.na(values(rasterized_nutsM)))  # cell numbers with data
  polygon_indices = rasterized_nutsM[cell_ids]        # the numeric values (polygon IDs)
  lookup_table = unique(as.data.frame(nutsM)[, c("LocationCo_numeric", "LocationCo")])
  
})

# sampling strategy for Environmental space sampling
ES_sampling = quote({
  ps.abs = USE::paSampling(variable_stack,
                           pres = myPres,
                           thres = 0.75,
                           H = NULL,
                           grid.res = opt_res$Opt_res, 
                           prev = ratios[[r]], 
                           sub.ts = FALSE,
                           plot_proc = F)
  
  
  ABS = ps.abs[,c("x","y")]
  raster_extent = extent(variable_stack[[1]])
  
  # Create an empty raster with the same extent and resolution
  raster_cells = raster(raster_extent, nrows = nrow(variable_stack[[1]]), ncols = ncol(variable_stack[[1]]))
  
  # Initialize the raster with NAs (no data)
  raster_cells[] = NA
  
  # Set pseudo-absence points (ABS) as 1 in the raster grid
  # 'ABS' contains the coordinates of pseudo-absence points
  for(j in 1:nrow(ABS)) {
    x = ABS$x[j]
    y = ABS$y[j]
    
    # Use the raster cell coordinates to set the value to 1
    cell = cellFromXY(raster_cells, c(x, y))
    raster_cells[cell] = 1
  }
  
  df_cell.nutsM = as.data.frame(matrix(nrow = length(values(raster_cells)), ncol = 3))
  colnames(df_cell.nutsM) = c("longitude", "latitude", "PA_value")
  df_cell.nutsM$longitude = xFromCol(raster_cells, 1:ncol(raster_cells))
  df_cell.nutsM$latitude = yFromRow(raster_cells, 1:nrow(raster_cells))
  df_cell.nutsM$PA_value = getValues(raster_cells) # pseudo absence cells (value=1)
  df_cell.nutsM = na.omit(df_cell.nutsM)
  
  # Keep only the overlapping cells in both rasters (nutsM ids and pseudoabsences)
  overlap = mask(raster_cells, rasterized_nutsM)
  overlap_idx = which(!is.na(values(overlap)))
  
  # Extract the polygon IDs from rasterized_nutsM at those indices
  polygon_ids_numeric = values(rasterized_nutsM)[overlap_idx]
  
  # Match numeric IDs to NUTS3 ID codes
  matched_ids = lookup_table$LocationCo[match(polygon_ids_numeric, lookup_table$LocationCo_numeric)]
  df_cell.nutsM$polygon_ID = matched_ids 
  
  final_ps.abs_df = as.data.frame(matrix(nrow = dim(df_cell.nutsM[1]), ncol = 8))
  colnames(final_ps.abs_df) = c("disease", "nutM_code", "year", "year_start","year_end","region_name","case_numbers", "response")
  final_ps.abs_df$year = NA;final_ps.abs_df$year_start = NA; final_ps.abs_df$year_end = NA;
  final_ps.abs_df$nutM_code = df_cell.nutsM$polygon_ID;
  final_ps.abs_df$disease = "aedes transmitted virus"
  final_ps.abs_df$response = 0 ; final_ps.abs_df$case_numbers = NA
  final_ps.abs_df$region_name = nutsM_IDs$region_name[match(final_ps.abs_df$nutM_code, nutsM_IDs$nutM_code)]
  
  n_abs = nrow(final_ps.abs_df)
  n_range_abs = round(n_abs * p_ranges)
  range_indices = sample(seq_len(n_abs), size = n_range_abs, replace = FALSE)
  single_indices = setdiff(seq_len(n_abs), range_indices)
  
  
  # Assign pseudo-years
  final_ps.abs_df$year[single_indices] = sample(
    years, size = length(single_indices), replace = TRUE, prob = year_probs
  )
  
  final_ps.abs_df = final_ps.abs_df[, !(colnames(final_ps.abs_df) %in% c("disease", "case_numbers"))]
  
  # Combine presences and absences
  data = rbind(ae.presences_agg, final_ps.abs_df)
})



################################################################################
#                   III. BRT MODEL TRAINING AND EVALUATION
################################################################################


# 1. Environmental data extraction
########################################

nruns = 30 # this will be used as the number of repetitions for the brt training
pa.strategy = c("random", "geo_ae.density", "surveillance","espace") # pseudo-absence sampling strategy
namesRatios = c("1.1", "1.3", "1.6") # ratio names
ratios = c(1,3,6) # corresponds to the ratio of desired sampled PAs (1:1, 1:3, 1:6)
ratios2 =c(3.7, 1.3, 0.64); # interpolated prevalence values for Environmental space sampling (to get correct ratio values)

croplands = brick("Rasters/Environmental_rasters/Land-use/croplands_LUH2-GCB2024_1901_2024.nc")
pastures = brick("Rasters/Environmental_rasters/Land-use/pastures_LUH2-GCB2024_1901_2024.nc")
rangelands = brick("Rasters/Environmental_rasters/Land-use/rangeland_LUH2-GCB2024_1901_2024.nc")
urbanAreas = brick("Rasters/Environmental_rasters/Land-use/urbanAreas_LUH2-GCB2024_1901_2024.nc")
primaryForest = brick("Rasters/Environmental_rasters/Land-use/primaryForest_LUH2-GCB2024_1901_2024.nc")
primaryNonForest = brick("Rasters/Environmental_rasters/Land-use/primaryNonF_LUH2-GCB2024_1901_2024.nc")
secondaryForest = brick("Rasters/Environmental_rasters/Land-use/secondaryForest_LUH2-GCB2024_1901_2024.nc")
secondaryNonForest = brick("Rasters/Environmental_rasters/Land-use/secondaryNonF_LUH2-GCB2024_1901_2024.nc")
population = brick("Rasters/Environmental_rasters/Population/Population_density_1900_2100/Population_density_05deg_1900_2025_historical.nc")

# run PCA for the environmental space sampling
eval(ES_setup)
# set parallele process
plan(multisession, workers = parallel::detectCores() - 1)

# Loop through all pseudo-absence strategies and ratios
for(s in 1:length(pa.strategy)){
  data_train = list()
  # load climate rasters
  temperature = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_1901_2021_monmean.nc"); temperature2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_tas_2022_2024_monmean.nc")
  precipitation = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_1901_2021_monmean.nc"); precipitation2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_pr_2022_2024_monmean.nc")
  relativehumidity = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_1901_2021_monmean.nc"); relativehumidity2 = brick("Rasters/Environmental_rasters/ISIMIP3a/obsclim/20CRv3-ERA5/20crv3-era5_obsclim_hurs_2022_2024_monmean.nc")
  temperature = addLayer(temperature, temperature2); precipitation = addLayer(precipitation, precipitation2); relativehumidity = addLayer(relativehumidity,relativehumidity2)
  
  # ratio loop
  for (r in 1:length(ratios)){
    data_train[[r]] = list()

    # repition loop 
    for(n in 1:nruns){
      if(s==1){eval(random_sampling)}
      if(s==2){eval(geo_density_sampling)}
      if(s==3){eval(surveillance_sampling)}
      if(s==4){ratios=ratios2; eval(ES_sampling)}
      
      ## data extraction ##
      
      data_for_brt = as.data.frame(matrix(NA,nrow=nrow(data), ncol =27))
      data_for_brt[,c(1:5)] = data[,c("nutM_code","year","year_start","year_end", "response")]
      colnames(data_for_brt) = c("nutM_code","year","year_start","year_end", "response",
                                 "temperature_winter","temperature_spring","temperature_summer","temperature_fall",
                                 "seasonal_variation","precipitation_winter","precipitation_spring","precipitation_summer","precipitation_fall",
                                 "relative_humidity_winter","relative_humidity_spring","relative_humidity_summer","relative_humidity_fall",
                                 "croplands","pastures","rangelands","urbanAreas","primaryForest","primaryNonForest","secondaryForest",
                                 "secondaryNonForest", "population")
      
      geometry_matches = nutsM_sf$geometry[match(data_for_brt$nutM_code, nutsM_sf$LocationCo)]
      nutsM_train = st_sf(data_for_brt, geometry = geometry_matches)
      names_env = names(data_for_brt[,c(6:27)])
      
      cat("Starting data extraction", n, "of", nruns, "for PA sampling method", pa.strategy[s], "ratio", namesRatios[r])
      # parallel data extraction block 
      results = future_lapply(1:nrow(data_for_brt), function(k) {
        if (!is.na(data_for_brt$year[k])) { 
          # Temporal training - split climate variables into seasons
          message("temporal training on single year")
          year = as.numeric(data_for_brt[k,"year"]);  year2 = year
          start = (year - 1901)*12+1; end = start + 11
          month_start = ((year-1) - 1901)*12+1
          winter_indices = sort(c(month_start + 11, month_start + 12, month_start + 13)) # December (previous year), January, February
          spring_indices = sort(c(month_start + 14, month_start + 15, month_start + 16)) # March, April, May
          summer_indices = sort(c(month_start + 17, month_start + 18, month_start + 19)) # June, July, August
          fall_indices   = sort(c(month_start + 20, month_start + 21, month_start + 22)) # September, October, November
          
          year_index = year2 - 1900 # landuse is updated and goes until 2024
          
          # seasonal variation
          yearly_temp = temperature[[start:end]]-273.15
          temp_raster = crop(yearly_temp, nutsM, snap = "out")
          temp_raster = mask(temp_raster, nutsM)
          
          max = cellStats(temp_raster, stat = "max", na.rm = TRUE)
          min = cellStats(temp_raster, stat = "min", na.rm = TRUE)
          
          hotest_indice = which.max(max)
          hotest_indice = start + hotest_indice - 1
          
          coldest_indice = which.min(min)
          coldest_indice = start + coldest_indice - 1
          
          
        } else {
          message("temporal training on average year range")
          year_start = as.numeric(data_for_brt[k,"year_start"]); year_end = as.numeric(data_for_brt[k,"year_end"])
          year_start_temp = year_start - 1900;  year_end_temp = year_end - 1900
          year_index = year_start_temp:year_end_temp
          
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
          
          # seasonal variation
          allyears = year_start:year_end
          hotest_indice = list()
          coldest_indice = list()
          for (i in 1:length(allyears)){
            yr = allyears[[i]]
            start = (yr - 1901)*12+1; end = start + 11
            yearly_temp = temperature[[start:end]]-273.15
            temp_raster = crop(yearly_temp, nutsM, snap = "out")
            temp_raster = mask(temp_raster, nutsM)
            
            max = cellStats(temp_raster, stat = "max", na.rm = TRUE)
            min = cellStats(temp_raster, stat = "min", na.rm = TRUE)
            
            Hotest_indice = which.max(max)
            hotest_indice[[i]] = start + Hotest_indice - 1
            Coldest_indice = which.min(min)
            coldest_indice[[i]] = start + Coldest_indice - 1
            
          }
        }
        
        # Extract environmental variables for the given year or year range
        temperature_temp_winter = mean(temperature[[winter_indices]]) - 273.15 # conversion to Celsius
        temperature_temp_spring = mean(temperature[[spring_indices]]) - 273.15
        temperature_temp_summer = mean(temperature[[summer_indices]]) - 273.15
        temperature_temp_fall   = mean(temperature[[fall_indices]]) - 273.15
        
        hotmonth_temp = mean(temperature[[hotest_indice]]) - 273.15
        coldmonth_temp = mean(temperature[[coldest_indice]]) - 273.15
        seasonal_var_temp = hotmonth_temp - coldmonth_temp
        
        precipitation_temp_winter = mean(precipitation[[winter_indices]]) * 60 * 60 * 24 # conversion to kg/m2/day
        precipitation_temp_spring = mean(precipitation[[spring_indices]]) * 60 * 60 * 24
        precipitation_temp_summer = mean(precipitation[[summer_indices]]) * 60 * 60 * 24
        precipitation_temp_fall   = mean(precipitation[[fall_indices]]) * 60 * 60 * 24
        
        relhumidity_temp_winter = mean(relativehumidity[[winter_indices]])
        relhumidity_temp_spring = mean(relativehumidity[[spring_indices]])
        relhumidity_temp_summer = mean(relativehumidity[[summer_indices]])
        relhumidity_temp_fall   = mean(relativehumidity[[fall_indices]])
        
        croplands_temp = mean(croplands[[year_index]])
        pastures_temp = mean(pastures[[year_index]])
        urbanAreas_temp = mean(urbanAreas[[year_index]])
        rangelands_temp = mean(rangelands[[year_index]])
        primaryForest_temp = mean(primaryForest[[year_index]])
        primaryNonForest_temp = mean(primaryNonForest[[year_index]])
        secondaryForest_temp = mean(secondaryForest[[year_index]])
        secondaryNonForest_temp = mean(secondaryNonForest[[year_index]])
        population_temp = mean(population[[year_index]])
        
        envVariables = list(
          temperature_temp_winter, temperature_temp_spring, temperature_temp_summer, temperature_temp_fall,
          seasonal_var_temp, precipitation_temp_winter, precipitation_temp_spring, precipitation_temp_summer, precipitation_temp_fall,
          relhumidity_temp_winter, relhumidity_temp_spring, relhumidity_temp_summer, relhumidity_temp_fall,
          croplands_temp, pastures_temp, urbanAreas_temp, rangelands_temp,
          primaryForest_temp, primaryNonForest_temp, secondaryForest_temp, secondaryNonForest_temp,
          population_temp
        )
        
        envVariables = lapply(envVariables, crop, y = nutsM, snap = "out")
        envVariables = lapply(envVariables, mask, mask = nutsM)
        
        vals = numeric(length(envVariables))
        for (i in seq_along(envVariables)) {
          if (i == length(envVariables)) {
            vals[i] = log(exact_extract(envVariables[[i]], nutsM_train[k,], fun = "sum") + 1)
          } else {
            vals[i] = exact_extract(envVariables[[i]], nutsM_train[k,], fun = "mean")
          }
        }
        
        return(vals)
      })
      data_for_brt[, names_env] = do.call(rbind, results)
      data_train[[r]][[n]] = data_for_brt
      
      cat("Data extraction run", n, "of", nruns, "for PA sampling method", 
          pa.strategy[s], "ratio", namesRatios[r], "completed.\n")
    } # end of repitiions loop for environmental data extraction
  }
  saveRDS(data_train, paste0(file = "European_model/All_the_brt_models/PA_strategies_brts/training_data", pa.strategy[s],".rds"))
}


# 2. BRT model training
#######################

ratios = c(1,3,6) # corresponds to the ratio of desired sampled PAs (1:1, 1:3, 1:6)

# Calculate the mean block size for spatial cross-validation per PA strategy and ratio
file_paths = c(
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_random.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_geo_ae.density.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_surveillance.rds",
  "./European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds"
)
#data_train4 = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds")

Data = lapply(file_paths, readRDS)
# keep one training dataset per PA sampling method
Data = lapply(Data, function(three_list) {
  lapply(three_list, function(sublist) {
    sublist[[1]]
  })
})

centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  

plottingCorrelogram = FALSE

for (i in 1:length(Data)){
  for (j in 1:length(Data[[i]])){
    data = Data[[i]][[j]]
    matching_centroids = as.data.frame(matrix(nrow = nrow(data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
    matching_centroids$NUTM = data$nutM_code 
    matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
    correlogram = ncf::correlog(matching_centroids[,2], matching_centroids[,3], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)

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
    }
  }
}

# chose block size per strategy and ratio based on plotted correlogram (when curve intercepts 0)
blocks = list(Random = list(600, 800, 800),
     Ae = list(500,700,700),
     surv = list(600, 500,800),
     es = list(1000, 1100, 1000))

for(s in 1:length(pa.strategy)){
  brt_model_scvs = list() # brt with spatial cross-validations (SCVs)
  AUCs = matrix(nrow=nruns, ncol= length(ratios)); colnames(AUCs) = paste0(pa.strategy[[s]], "_", namesRatios)
  SIppcs = matrix(nrow=nruns, ncol=length(ratios)); colnames(SIppcs) = paste0(pa.strategy[[s]], "_", namesRatios)
  thresholds = matrix(nrow=nruns, ncol=length(ratios)); colnames(thresholds) = paste0(pa.strategy[[s]], "_", namesRatios)
  tabs_list1 = list()
  
  if(s==1){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_random.rds")}
  if(s==2){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_geo_ae.density.rds")}
  if(s==3){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_surveillance.rds")}
  if(s==4){data_train = readRDS("European_model/All_the_brt_models/PA_strategies_brts/training_data_espace.rds")}
  
  for (r in 1:length(ratios)){
    brt_model_scvs[[r]] = list()
    tabs_list1[[r]] = list()
    current_block = blocks[[s]][[r]]
    for(b in 1:nruns){
      current_data = data_train[[r]][[b]]
      # BRT model parameters
      gbm.x = 6:27 # environmental predictors
      gbm.y = colnames(current_data)[5] # response variable
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
      theRanges = c(current_block, current_block)*1000 # transform block size distance to meters
      
      centroids = as.data.frame(matrix(nrow=length(nutsM), ncol=3)); colnames(centroids) = c("NUTM","longitude","latitude")  
      centroids$NUTM = nutsM@data$LocationCo; coords = coordinates(nutsM)  
      centroids$longitude = coords[,1]; centroids$latitude = coords[,2]  
      
      matching_centroids = as.data.frame(matrix(nrow = nrow(current_data), ncol = 3)); colnames(matching_centroids) = c("NUTM","longitude", "latitude")
      matching_centroids$NUTM = current_data$nutM_code 
      matching_centroids$longitude = centroids$longitude[match(matching_centroids$NUTM, centroids$NUTM)]; matching_centroids$latitude = centroids$latitude[match(matching_centroids$NUTM, centroids$NUTM)]
      matching_centroids$response = current_data$response
      
      coordinates(matching_centroids) = ~longitude + latitude
      crs(matching_centroids) = crs(nutsM)
      
      spdf = SpatialPointsDataFrame(matching_centroids, current_data[,1:dim(current_data)[2]])
      myblocks = cv_spatial(spdf, column ="response", k=n.folds, size=theRanges, selection="random", plot = F)
      fold.vector = myblocks$folds_ids; n.trees = 100
      
      brt_model_scvs[[r]][[b]] =  gbm.step(current_data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
                                           var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
                                           verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit, plot.main = F)
      
      
      # 3. Model predictive performance evaluation (3 different metrics)
      ##################################################################
      
      #### Calculation of the AUC ####
      AUCs[b,r] = brt_model_scvs[[r]][[b]]$cv.statistics$discrimination.mean # Mean test AUC
      
      #### Calculation of the sorensen index to establish predictive performance of our models ####
      # Sources:
      # - computation performed according to the formulas of Leroi et al. (2018, J. Biogeography)
      # - optimisation of the threshold with a 0.01 step increment according to Li & Guo (2013, Ecography)
      
      tmp = matrix(nrow=101, ncol=2); tmp[,1] = seq(0,1,0.01)
      df = brt_model_scvs[[r]][[b]]$gbm.call$dataframe
      responses = df$response;
      data = df[,6:dim(df)[2]]
      n.trees = brt_model_scvs[[r]][[b]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
      prediction = predict.gbm(brt_model_scvs[[r]][[b]], data, n.trees, type, single.tree)
      N = dim(data)[1]; P = sum(responses==1); A = sum(responses==0)
      prev = P/(P+A) # proportion of recorded sites where the species is present
      x = (P/A)*((1-prev)/prev);
      sorensen_ppc = 0
      for (threshold in seq(0,1,0.01))
      {
        TP = length(which((responses==1)&(prediction>=threshold))) # true positives
        FN = length(which((responses==1)&(prediction<threshold))) # false negatives
        FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
        sorensen_ppc_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
        tmp[which(tmp[,1]==threshold),2] = sorensen_ppc_tmp
        if (sorensen_ppc < sorensen_ppc_tmp)
        {
          sorensen_ppc = sorensen_ppc_tmp
          optimised_threshold = threshold
        }
      }
      tabs_list1[[r]][[b]] = tmp
      SIppcs[b,r] = sorensen_ppc
      thresholds[b,r] = optimised_threshold
      
     }
  }
  
  print(paste0("Saving results for ", pa.strategy[s], " sampling"))
  saveRDS(brt_model_scvs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/brt_trained_model_", pa.strategy[s],".rds"))
  saveRDS(tabs_list1, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/tabs_list1_", pa.strategy[s],".rds"))
  write.csv(SIppcs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/sippc_index_", pa.strategy[s],".csv"))
  write.csv(thresholds, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/sippc_thresholds_", pa.strategy[s],".csv"))
  write.csv(AUCs, file = paste0("European_model/All_the_brt_models/PA_strategies_brts/AUCs_", pa.strategy[s], ".csv"))
}


# 4. Model performance figures 
##############################

# average model performance metrics table 
# Define statistical function
compute_stats = function(df, metric_name, ratios = c(1, 3, 6)) {
  summary_df = data.frame()
  
  for (i in seq_along(ratios)) {
    x = df[[i]]
    mean_val = mean(x, na.rm = TRUE)
    se = sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
    ci = 1.96 * se  # 95% CI
    
    CI_lower = mean_val - ci
    CI_upper = mean_val + ci
    
    formatted = sprintf("%.3f (%.3f – %.3f)", mean_val, CI_lower, CI_upper)
    
    summary_df = rbind(summary_df, data.frame(
      ratio = ratios[i],
      mean = mean_val,
      CI_lower = CI_lower,
      CI_upper = CI_upper,
      metric = metric_name,
      formatted = formatted
    ))
  }
  return(summary_df)
}

# Set file path
path = "European_model/All_the_brt_models/PA_strategies_brts/"

# Load all data files 
aucs = list(
  random = read.csv(paste0(path, "AUCs_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "AUCs_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "AUCs_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "AUCs_espace.csv"), header = TRUE)
)

sippc = list(
  random = read.csv(paste0(path, "SIppcs_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "SIppcs_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "sippc_index_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "sippc_index_espace.csv"), header = TRUE)
)

thresh = list(
  random = read.csv(paste0(path, "thresholds_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "thresholds_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "sippc_thresholds_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "sippc_thresholds_espace.csv"), header = TRUE)
)

bi = list(
  random = read.csv(paste0(path, "boyce_index_random.csv"), header = TRUE),
  geo    = read.csv(paste0(path, "boyce_index_geo_ae.density.csv"), header = TRUE),
  surv   = read.csv(paste0(path, "boyce_index_surveillance.csv"), header = TRUE),
  es     = read.csv(paste0(path, "boyce_index_espace.csv"), header = TRUE)
)

# labels
pa_ratios = c("ratio 1:1", "ratio 1:3", "ratio 1:6")
strategies = c("random", "geo", "surv", "es")
strategy_names = c("Random sampling",
                   "Weigthed vector distribution sampling",
                   "Weighted surveillance capability sampling",
                   "Environmental space sampling")

# Initialize table
metrics.table = data.frame(
  `Pseudo-Absence (PA) sampling strategy` = rep(strategy_names, each = 3),
  `PA ratio` = rep(pa_ratios, times = 4),
  stringsAsFactors = FALSE
)

# Fill the table 
row_counter = 1
for (strategy in strategies) {
  auc_data = aucs[[strategy]][, -1]
  sippc_data = sippc[[strategy]][, -1]
  thresh_data = thresh[[strategy]][, -1]
  bi_data = bi[[strategy]][, -1]
  
  auc_stats = compute_stats(auc_data, "AUC")
  sippc_stats = compute_stats(sippc_data, "SIppc")
  thresh_stats = compute_stats(thresh_data, "Threshold")
  bi_stats = compute_stats(bi_data, "BI")
  
  for (i in 1:3) {
    metrics.table[row_counter, "Area Under the Curve (AUC)"] = auc_stats$formatted[i]
    metrics.table[row_counter, "Prevalence-pseudo-absence calibrated Sorensen index (SIppc)"] = sippc_stats$formatted[i]
    metrics.table[row_counter, "Threshold value maximising the SIppc"] = thresh_stats$formatted[i]
    metrics.table[row_counter, "Boyce Index (BI)"] = bi_stats$formatted[i]
    row_counter = row_counter + 1
  }
}

# write to Excel
write.xlsx(metrics.table, paste0(path,"summary_metrics_table.xlsx"), rowNames = FALSE)

################################################################################
#                       IV. SUPPLEMENTARY FIGURES
################################################################################

# 1. pseudo-absence maps
##########################################
path = "all_data/all_the_brt_models/parallele/all_pa_method_results/"

# Load files 
datas = list(
  random = readRDS(paste0(path, "env_training_data_random.rds")),
  geo    = readRDS(paste0(path, "env_training_data_geo_ae.density.rds")),
  surv   = readRDS(paste0(path, "env_training_data_surveillance.rds")),
  es     = readRDS(paste0(path, "env_training_data_espace.rds"))
)
ae = expression(italic("Aedes albopictus"))
names = c("Random PA sampling"," " ,"Disease surveillance capability \n weighted PA sampling", "Environmental space \n PA sampling")
ratios = c("Ratio 1:1", "Ratio 1:3", "Ratio 1:6")

pres = datas$random[[1]][[1]]$nutM_code[datas$random[[1]][[1]]$response == 1]
pres = subset(nutsM, LocationCo %in% pres)

pdf("github/Figures/FigS1_pa_maps.pdf", width=6, height=10)
par(mfrow=c(4,3), oma=c(1,2.5,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
for (i in 1:length(datas)) {
  for (j in 1:length(ratios)) {
    plot(contour, lwd=0.4, border="gray30", col=NA)
    if (j==1){
      if(i==2){
        mtext(expression(italic("Aedes albopictus")~"suitability"), side = 2, line = 1.0, cex = 0.6)
        mtext("weighted PA sampling", side = 2, line = 0.2, cex = 0.6) }else{mtext(names[[i]], side = 2, line= 0.5,cex = 0.6)}
      }
      
    if (i==1) {mtext(ratios[[j]], side = 3, line= -1.5, cex = 0.6)}
    abs = datas[[i]][[j]][[1]]$nutM_code[datas[[i]][[j]][[1]]$response == 0]
    abs = subset(nutsM, LocationCo %in% abs)
    
    plot(nutsM, col = "gray95", border = T, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
    plot(ae.adj_geoms, col = "gray75", border = T, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
    plot(pres, col = "gray50", border = T, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
    plot(abs, col = "lightskyblue3", border = T, lwd = 0.1, ann = FALSE, legend = FALSE, axes = FALSE, box = FALSE, add=T)
  }
}

dev.off()







