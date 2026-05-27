

# Temperature sensitivity analysis -------------------------------------------------


# RangeX phenology effect of transplantation on temperature sensitivity ------------

## Data used: RangeX_clean_phenology_2023_NOR.csv
##            
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
##            RangeX_clean_climate_station_NOR_2021-2025.csv
## Date:      27.05.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on temperature sensitivity


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lme4)
library(performance)
library(see)
library(emmeans)


# set theme for plots for presentation ------------------------------------
# theme_set(theme_bw(base_size = 20))


# source script with functions --------------------------------------------
source("Temperature_sensitivity_functions.R")


# calculate delta T -------------------------------------------------------
# use growing season mean because this is not specific to species
# we have species with different onset timings

# source climate scripts --------------------------------------------------
source("Data_preparation_climate_station_NOR.R")
climate_23 
climate_23_daily

source("Data_preparation_climate_station_CHE.R")
climate_che_combined
# which is daily per site



# mean temperature growing season mean ------------------------------------------
# NOR
climate_nor_23_gs <- climate_23_daily |>
  filter(
    date >= as.Date("2023-04-01"),
    date <= as.Date("2023-09-30")
  ) |> 
  mutate(region = "Norway")


# CHE
climate_che_22_gs <- climate_che_combined |>
  filter(
    date >= as.Date("2022-04-01"),
    date <= as.Date("2022-09-30")
  ) |> 
  mutate(region = "Switzerland")


# combine
climate_all_gs <- bind_rows(climate_che_22_gs, climate_nor_23_gs)


# calculate mean per site in this time period
pre_climate_gs <-  climate_all_gs |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate_gs




# source clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)



# rename infructescence stage of NOR to FlowWithrd ------------------------
# this will be the fruiting stage
# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology <- phenology |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "No_Infructescences" = "No_FloWithrd"))

# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology3 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional 



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology3 <- phenology3 |> 
  filter(treat_warming == "ambi")



phenology3 <- phenology3 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))

































