


# RangeX phenology data preparation NOR and CHE ------------

## Data used: RangeX_clean_Phenology_2022_CHE.csv
##            RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
## Date:      14.08.25
## Author:    Nadine Arzt
## Purpose:   Combine NOR and CHE data and prepare


# comment -----------------------------------------------------------------
# CHE Flowers withered corresponds to infructescences in NOR

# CHE.lo.ambi.bare.wf.10.22.1 silvul 2022-06-02 EI No_Buds 11 --> should be 1

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(purrr)


# import phenology data CHE ---------------------------------------
# Read as a single column of text
raw <- read_lines("Data/RangeX_clean_Phenology_2022_CHE.csv")

# Remove the outer wrapping quotes and fix doubled quotes
raw <- gsub('^"|"$', '', raw)   # remove first and last quote in each line
raw <- gsub('""', '"', raw)     # replace doubled quotes with single quotes

# Write a cleaned temporary CSV
write_lines(raw, "Data/RangeX_clean_Phenology_2022_CHE_clean.csv")

# Now read normally
pheno_che <- read_csv("Data/RangeX_clean_Phenology_2022_CHE_clean.csv")

# import phenology data NOR ---------------------------------------
# pheno_nor <- read.csv("Data/RangeX_clean_phenology_2023_NOR.csv")

pheno_nor <- read.csv("Data/Clean/RangeX_clean_Phenology_2023_NOR.csv")

# import metadata CHE -----------------------------------------------------
meta_CHE <- read_csv("Data/RangeX_clean_MetadataFocal_CHE.csv")

# import metadata NOR -----------------------------------------------------
meta_NOR <- read.csv("Data/RangeX_clean_MetadataFocal_NOR.csv")


# merge metadata with phenology -------------------------------------------
# CHE
pheno_22_CHE <- left_join(meta_CHE, pheno_che, by = c("unique_plant_ID", "species"))

# NOR
pheno_23_NOR <- left_join(meta_NOR, pheno_nor, by = c("unique_plant_ID", "species"))
# nor has column comment
# delete for now?

# pheno_23_NOR <- pheno_23_NOR |> 
#   select(-comment)


# one phenology data set -------------------------------------------------
# combine CHE and NOR
phenology <- rbind(pheno_22_CHE, pheno_23_NOR)


# rename pheno stages to match regions ------------------------------------
# "number_infructescences" in NOR correpsonds to "No_FloWithrd" in CHE?
phenology <- phenology |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "number_buds" = "No_Buds",
                                  "number_flowers" = "No_FloOpen",
                                  #"number_infructescences" = "No_FloWithrd",
                                  "seeds_collected" = "No_Seeds"))


# fix typo in che ---------------------------------------------------------
# CHE.lo.ambi.bare.wf.10.22.1 silvul 2022-06-02 EI No_Buds 11 --> should be 1
phenology <- phenology |>
  mutate(value = replace(value,
                         unique_plant_ID == "CHE.lo.ambi.bare.wf.10.22.1" &
                           value == 11, 1))



# save joint data set -----------------------------------------------------
#write.csv(phenology, file = "Data/RangeX_clean_phenology_NOR_CHE.csv")




# Filter buds, flowers, infructescences -----------------------------------
## don't use seeds_collected
# phenology <- phenology |> 
#   filter(phenology$phenology_stage != "No_Seeds")

# combined treatment column -----------------------------------------------
phenology$treatment <- paste(phenology$site, phenology$treat_warming, phenology$treat_competition, sep = "_")

# change region and treatment names  --------------------------------------
phenology <- phenology |>
  mutate(region = case_when(
    region == "NOR" ~ "Norway",
    region == "CHE" ~ "Switzerland",
    TRUE ~ region
  ))


# will change to biotic interactions - with and without
phenology <- phenology |>
  mutate(treat_competition = case_when(
    treat_competition == "bare" ~ "without",
    treat_competition == "vege" ~ "with",
    TRUE ~ treat_competition
  ))




# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology2 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional scaling if needed




# add site_warming treatment ----------------------------------------------
phenology2$treatment_site_temp <- paste(phenology2$site, phenology2$treat_warming, sep = "_")


phenology2 <- phenology2 |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_ambi",
                                                "hi_warm")))



phenology2 <- phenology2 |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))
































