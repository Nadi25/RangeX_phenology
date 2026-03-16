


# Data preparation - phenology NOR 23 ----------------------------------------





# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(purrr)

# load clean phenology data -----------------------------------------------
pheno_23_NOR <- read.csv("Data/Clean/RangeX_clean_Phenology_2023_NOR.csv")



# import metadata NOR -----------------------------------------------------
meta_NOR <- read.csv("Data/RangeX_clean_MetadataFocal_NOR.csv")


phenology <- left_join(meta_NOR, pheno_23_NOR, by = c("unique_plant_ID", "species"))


# Filter buds, flowers, infructescences -----------------------------------
## don't use seeds_collected
phenology <- phenology |> 
  filter(phenology$phenology_stage != "No_Seeds")

# combined treatment column -----------------------------------------------
phenology$treatment <- paste(phenology$site, phenology$treat_warming, phenology$treat_competition, sep = "_")

# change region and treatment names  --------------------------------------
phenology <- phenology |>
  mutate(region = case_when(
    region == "NOR" ~ "Norway",
    TRUE ~ region
  ))

phenology <- phenology |>
  mutate(treat_competition = case_when(
    treat_competition == "bare" ~ "without",
    treat_competition == "vege" ~ "with",
    TRUE ~ treat_competition
  ))

phenology <- phenology |>
  mutate(treat_warming = case_when(
    treat_warming == "ambi" ~ "ambient",
    treat_warming == "warm" ~ "warmed",
    TRUE ~ treat_warming
  ))

phenology <- phenology |>
  mutate(site = case_when(
    site == "lo" ~ "low",
    site == "hi" ~ "high",
    TRUE ~ site
  ))




# and add julian days --------------------------------------
phenology <- phenology |> 
  mutate(jday = yday(date_measurement),   # Julian day (1–365)
         jday_scaled = scale(jday))   



# rename pheno stages to match regions ------------------------------------
# "number_infructescences" in NOR correpsonfs to "No_FloWithrd" in CHE?
phenology <- phenology |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "number_buds" = "No_Buds",
                                  "number_flowers" = "No_FloOpen",
                                  #"number_infructescences" = "No_FloWithrd",
                                  "seeds_collected" = "No_Seeds"))


# only flowers ----------------------------------------------------------
phenology <- phenology |> 
  filter(phenology_stage == "No_FloOpen")


# change reference to be low site ----------------------------------------------
# by factor
phenology <- phenology |>
  mutate(
    site = factor(site),
    site = relevel(site, ref = "low")
  )


# do the filtering in the model











