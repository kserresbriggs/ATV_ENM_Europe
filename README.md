<h1> Escalating human exposure to <i>Aedes</i> mosquito-borne viruses in Europe </h1>


This repo gathers the input files and scripts related to our study entitled "**Escalating human exposure to <i>Aedes</i> mosquito-borne viruses in Europe**" (*submitted*). R scripts used to conduct the ecological niche modelling analyses described in this study are all available in `Script_ATV_ENM_Europe.r`.

**Abstract:** Climate change has multiple public health impacts, including the expansion of the geographical range of vector-borne diseases. Pathogens such as dengue, chikungunya and Zika viruses, transmitted by <i>Aedes</i> mosquitoes, can cause severe health outcomes ranging from acute febrile illness, chronic joint pain, to birth defects and even death. Evaluating the future risk of human population exposure is therefore crucial as large outbreaks could overwhelm healthcare systems. Europe, one of the fastest warming regions globally, harbours the competent mosquito vector <i>Aedes</i> albopictus in over 20 countries, making <i>Aedes</i>-borne viruses an increasing threat to the continent, which has already experienced local outbreaks over the past two decades. Here we use an ecological niche modelling approach to assess past, present, and future risk of human population exposure to <i>Aedes</i>-borne viruses in Europe. Our results show that recent climate change has already increased the potential exposure to these viruses, particularly across the Mediterranean basin, which is a current hotspot for local outbreaks. Major metropolitan areas in Spain, France, Italy, and Croatia are by now located in at-risk areas, and this risk is projected to intensify and expand northward by mid-century. Under a high greenhouse gas emissions scenario, European areas ecologically suitable for <i>Aedes</i>-borne virus circulation could nearly double, leading to an additional ~50 million people living in areas at risk by the end of the century. These findings underscore the urgent need for strengthened vector and epidemiological surveillance, as well as preparedness strategies across newly suitable regions to anticipate future public health threats associated with these arboviral diseases.

** NEW FIGURE TO BE INSERTED **

**Future projections of the European areas ecologically suitable for the local circulation of <i>Aedes</i>-borne viruses and for the associated human population at risk of exposure under three different shared socio-economic pathways (SSPs).** The first row of maps displays the spatial distribution of the presence data (administrative areas with confirmed local infection cases) from 2007 to 2024 used to train the ecological niche models, the predicted ecological suitability for the local circulation of <i>Aedes</i>-borne viruses at the present time (2020-2024, referred to as “t0” in the map headings), and the future projections for three successive periods (2025-2049, 2050-2074, 2075-2100) under the high emission scenario SSP3-7.0 (see Figure S4 for the projections under the low and very high emission scenarios, SSP1-2.6 and SSP5-8.5, as well as Figure S5 for predictive presence-absence maps based on a binary classification of ecological suitability estimations). The second and third rows of maps display the difference in ecological suitability between estimates obtained for specific future projections for scenarios SSP3-7.0 or SSP5-8.5 and estimates obtained for the present-day (“t0”). Finally, the embedded boxplot graphic reports the number of people estimated to live in European areas at risk of exposure to local <i>Aedes</i>-borne virus circulation for the different periods and SSPs considered (see the Methods section for further detail on the approach used to define the areas at risk). (*) refers to a significant increase in the estimated number of people living in at-risk areas as compared to 2024 (as evaluated with a paired Wilcoxon test).

### Installation and system requirements ###
All available code is written in R (v4.4.1 or higher) using Rstudio (2024.04.2+764 or higher).
To run the scripts locally, adjust the root directory and file paths as needed for your system. The full workflow in Script_ATV_ENM_Europe.r was run on a MacBook Pro equipped with an Apple M3 Pro chip, 12 CPU cores (6 performance and 6 efficiency cores), and 36 GB RAM. On this system, the complete script runtime is approximately 3 hours. Runtime may vary depending on hardware configuration, available memory and storage speed.

- R installation: https://cran.rstudio.com/
- RStudio installation: https://posit.co/download/rstudio-desktop

### Analyses summary ###
`Script_ATV_ENM_Europe.R` performs the complete ecological niche modelling workflow:

1. **Data preparation** — administrative boundaries and Aedes-borne viral occurrence data curation.  
2. **Model training and evaluation** — pseudo-absence generation, environmental predictor selection, spatial cross-validation, hyperparameter tuning, and BRT model fitting.  
3. **Model interpretation** — variable importance assessment and response curves.  
4. **Current risk mapping** — ecological suitability predictions for 2020–2024.  
5. **Historical projections** — past suitability projections and visualisation.  
6. **Future projections** — future suitability projections and population exposure estimation.

`Script_ATV_ENM_PA_strategies.R` performs the different pseudo-absence sampling strategies workflows evaluated in the study.

### Files and structure ###

| Folder | Contents | Purpose |
|----------|----------|----------|
| `brt_projection_outputs/` | CSV files containing estimated ecological suitability values for each European NUTS3 region across current, historical (observed and counterfactual), and future (all SSP scenarios) periods for each BRT replicate | Stores model projection outputs|
| `brt_training_data/` | Training datasets containing presence records and sampled pseudo-absences linked to selected environmental predictors; includes 100 CSV files used to train independent BRT replicates | Stores model training inputs |
| `brt_training_data/pseudo_absence_strategies/` | Pseudo-absence sampling script together with datasets generated for each pseudo-absence strategy and sampling ratio | Stores pseudo-absence workflow and sensitivity analyses |
| `environmental_data/` | CSV files containing environmental predictor values for each European NUTS3 region across current, historical, and future periods | Stores environmental covariates used in modelling |
| `Europe_NUTS3_shapefile/` | Optimised European administrative boundary shapefiles | Provides spatial boundaries for analyses and mapping |
| `occurrence_data/` | CSV files containing viral occurrence records and aggregated presence datasets used for BRT training | Stores occurrence inputs for the model |


