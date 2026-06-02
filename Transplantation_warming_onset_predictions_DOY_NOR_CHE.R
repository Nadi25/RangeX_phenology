

# 01 Onset julian days -------------------------------------------------------


# Onset bud, flower, fruit model with site_treat_warming combination --------
# julian days

# load library ---------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(performance)
library(see)




# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology2)



# function to get first onset per individual ----------------------------------------------
get_onset <- function(data, stage_name, onset_type) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, treatment_site_temp, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID
    ) |>
    summarise(
      onset = min(.data[[onset_type]]),
      .groups = "drop"
    ) |> 
    # remove groups where budding never occurred
    filter(is.finite(onset))
}


# calculate first onset per species and plot for all stages ----------------
onset_bud    <- get_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_onset(phenology2, "No_Seeds", "jday")



# function of model per stage and region ---------------------------------------------














# function for predictions ------------------------------------------------




















