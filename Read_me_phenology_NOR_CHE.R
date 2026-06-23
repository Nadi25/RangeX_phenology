

# READ ME -----------------------------------------------------------------


# Data preparation --------------------------------------------------------

# phenology
source("Data_preparation_phenology_NOR.R")
source("Data_preparation_phenology_NOR_CHE_combined.R")

# climate stations
source("Data_preparation_climate_station_NOR.R")
source("Data_preparation_climate_station_CHE.R")

# tomst and calculate GDD
# predict tms data with weather station data in NOR
source("TMS4_Weatherstation_correlation_predictions.R")
#source("Data_preparation_TMS4_CHE.R")
source("Data_preparation_TMS4_CHE.R")

# biomass
source("Data_praparation_biomass_NOR.R")
source("Data_preparation_traits_biomass_NOR_24.R")

# traits
source("RangeX_data_paper_cleaning_demographic_traits_23.R")



# Onset -------------------------------------------------------------------

# DOY - julian days
source("Transplantation_warming_onset_predictions_DOY_NOR_CHE.R")

# overall model
source("Onset_DOY_overal_model_across_stages.R")

# GDD
# only ambient low vs high
source("GDD_Transplantation_onset_predictions_NOR_CHE.R")

# overall model
source("Onset_GDD_overal_model_across_stages.R")

# GDD vs DOY
source("DOY_GDD_TMS4_NOR_CHE.R")

# 

# Duration ----------------------------------------------------------------

source("Duration_bud_flower_fruit_predictions_NOR_CHE.R")



# Temperature sensitivity -------------------------------------------------

source("Temperature_sensitivity_functions.R")

source("Temperature_sensitivity_analysis.R")

source("Temperature_sensitivity_plots.R")




# Flower number -----------------------------------------------------------------

# flower number as it is
source("Flower_number_predictions_transplantation_warming.R")

# predicting biomass for 23
source("Biomass_traits_correlation_per_species.R")
source("Biomass_prediction_24_compare_with_real_24.R")
source("Biomass_prediction_23_per_species.R")

# combine biomass with phenology
source("Biomass_phenology_combine_species_models.R")

# flower number adjusted for biomass
source("Biomass_flower_number_predictions.R")












