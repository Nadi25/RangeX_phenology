

# RangeX phenology 2023 data cleaning -------------------------------------

## Data used:      RangeX_raw_phenology_low_2023.xlsx, 
##                 RangeX_raw_phenology_high_2023.xlsx,
##                 RangeX_clean_MetadataFocal_NOR.csv,
## Script sourced: RangeX_data_paper_cleaning_demographic_traits_23.R
## Date:      09.01.2025
## Author:    Nadine Arzt
## Purpose:   Cleaning of the complete raw data files of phenology 2023


# load library ------------------------------------------------------------

library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(readxl)
library(janitor)
library(here)

# Import data low + high --------------------------------------------------
## use read_delim from reader package
phenology_low_raw <- read_excel(here("Data/Data_phenology/RangeX_raw_phenology_low_2023.xlsx"),sheet = 2)
head(phenology_low_raw)
dput(colnames(phenology_low_raw))
str(phenology_low_raw)

phenology_high_raw <- read_excel(here("Data/Data_phenology/RangeX_raw_phenology_high_2023.xlsx"), sheet = 2)
head(phenology_high_raw)
dput(colnames(phenology_high_raw))
glimpse(phenology_high_raw)


# combine low and high site -----------------------------------------------
rangex_phenology_raw <- bind_rows(phenology_low_raw, phenology_high_raw)

head(rangex_phenology_raw)
str(rangex_phenology_raw)

length(rangex_phenology_raw)
dput(colnames(rangex_phenology_raw))

# check for rows that only have NAs
na_rows <- rangex_phenology_raw[rowSums(is.na(rangex_phenology_raw)) == ncol(rangex_phenology_raw), ]
na_rows # 0


# clean column names ------------------------------------------------------
rangex_phenology_raw <- clean_names(rangex_phenology_raw)
colnames(rangex_phenology_raw)

## rename column names to match with RangeX metadata file
rangex_phenology_raw <- rangex_phenology_raw |> 
  rename("block_ID_original" = "block_id",
         "plot_ID_original" = "plot_id",
         "position_ID_original" = "position_id")


# change data format ------------------------------------------------------
rangex_phenology_raw <- rangex_phenology_raw |> 
  mutate(date = as.Date(date))


# replace NAs with 0 ------------------------------------------------------
rangex_phenology_raw <- rangex_phenology_raw |> 
  mutate_at(vars(number_buds, number_flowers, number_infructescences, seeds_collected), ~replace(., is.na(.), 0))


# import meta data file ---------------------------------------------------
metadata <- read.csv("Data/Metadata/RangeX_clean_MetadataFocal_NOR.csv")
head(metadata)
colnames(metadata)

## clean column names
# makes ID to id
#metadata <- clean_names(metadata)


# merge metadata and phenology --------------------------------------------
rangex_phenology_raw <- left_join(rangex_phenology_raw, metadata,
                                  by = c("region", "site", "block_ID_original", "plot_ID_original", "position_ID_original", "species"))


# combined treatment column -----------------------------------------------
rangex_phenology_raw$treatment <- paste(rangex_phenology_raw$site, rangex_phenology_raw$treat_warming, rangex_phenology_raw$treat_competition, sep = "_")

colnames(rangex_phenology_raw)


# Problems to fix ---------------------------------------------------------
# 1. For cennig, we only counted stems and not flowers. So this is different to the rest.
# --> We know how many flowers we have per stem at peak season, when we did the demographic trait measurements.
# 
# 2. Cennig: can have more seeds collected then infructescences, because one stem has multiple flowers
# 
# 3. Seeds collected:
# --> When seeds have been collected, they obviously can't be counted as infructescences anymore next time.
# --> So we have to add the number of collected seeds to the number of infructescences for the next date.
# --> If 3 seeds have been collected on 27.09.23, we have to add 3 to infructescences on 10.10.23
# --> Can we specify that when the last seed has been collected, the next date does not get more seeds added to the infructescences? Because then the plant is completely done
# 
# 4. Date: Sometimes we have several data points for e.g. buds very close, when not all plots have been observed at the same date
# --> create new date column with same date for one block of observations



# Fix seeds collected -----------------------------------------------------
# Arrange the data by unique plant ID and date
rangex_phenology_clean <- rangex_phenology_raw |> 
  arrange(unique_plant_ID, date) |> 
  group_by(unique_plant_ID)  |> 
  mutate(
    # Use cumsum to sum up seeds collected in previous dates
    cumsum_seeds_collected = cumsum(seeds_collected),
    # Adjust the number of infructescences by adding the cumsum of seeds collected
    number_infructescences_adj = number_infructescences + lag(cumsum_seeds_collected, default = 0)
  ) |> 
  ungroup()

## rename number infructescences and delete old column
rangex_phenology_clean <- rangex_phenology_clean |> 
  rename("number_infructescences_old" = "number_infructescences",
         "number_infructescences" = "number_infructescences_adj")

## delete number_infructescences_old, cumsum_seeds_colleted
rangex_phenology_clean <- rangex_phenology_clean |> 
  select(-c(number_infructescences_old, cumsum_seeds_collected))



# Fix date ----------------------------------------------------------------
## change date name to date_old
rangex_phenology_clean <- rangex_phenology_clean %>%
  rename("date_old" = "date")

## count how many measurements per day to see which date needs fixing
count <- count(rangex_phenology_clean, date_old, site)
count

## these dates need fixing
## low site: 11.09 and 13.09 need to be together --> 13.09.2023
## high site: 15./17./18./19.07 --> 18.07.2023
## hi: 17./18.08 --> 18.08.2023
## hi: 29./30.08 --> 29.08.2023

rangex_phenology_clean <- rangex_phenology_clean |> 
  mutate(date = case_when(
    # For the low site where date is 2023-09-11, change it to 2023-09-13
    site == "lo" & date_old == "2023-09-11" ~ as.Date("2023-09-13"),
    # For the high site where date is 15./17./18./19.07, change it to 18.07.2023
    site == "hi" & date_old %in% c("2023-07-15", "2023-07-17", "2023-07-19") ~ as.Date("2023-07-18"),
    # For the high site where date is 17./18.08, change it to 18.08.2023
    site == "hi" & date_old == "2023-08-17" ~ as.Date("2023-08-18"),
    # For the high site where date is 29./30.08, change it to 29.08.2023
    site == "hi" & date_old == "2023-08-30" ~ as.Date("2023-08-29"),
    # For all other cases, keep the original date
    .default = date_old
  ))

count_new <- count(rangex_phenology_clean, date, site)
count_new

## delete date_old
rangex_phenology_clean <- rangex_phenology_clean |> 
  select(-date_old)



# Fix Centaurea problem ---------------------------------------------------

# Source demographic traits low and high 23 -------------------------------
source("RangeX_data_paper_cleaning_demographic_traits_23.R")

demographic_traits_23 <- demo_traits_2023


#demographic_traits_23 <- read.csv("Data/RangeX_clean_traits_2023.csv")

# Merge phenology with demographic traits ---------------------------------
# rename unique_plant_ID to match
# rangex_phenology_clean <- rangex_phenology_clean |> 
#   rename("unique_plant_ID" = "unique_plant_id")

rangex_phenology_demo_traits <- left_join(rangex_phenology_clean, demographic_traits_23,
                                                by = c("unique_plant_ID", "species"))



# Calculate number of flowers per stem ------------------------------------

rangex_phenology_demo_traits <- rangex_phenology_demo_traits |> 
  mutate(number_flowers_per_stem = number_flowers.y / no_rep_stems)



# Replace number_flowers with number_flowers_per_stem for cennig ----------

## rename column number of flowers to old
rangex_phenology_demo_traits <- rangex_phenology_demo_traits |> 
  rename("number_flowers_old" = "number_flowers.x")

## for cennig I want to have number_flowers_old (= no_stems) * number_flower_per_stem to get the actual number of flowers
rangex_phenology_demo_traits <- rangex_phenology_demo_traits |> 
  mutate(number_flowers.y = case_when( species == "cennig"  & number_flowers_old != 0 ~ number_flowers_per_stem*number_flowers_old,
                                     # For all other cases, keep the original no flowers
                                     .default = number_flowers_old
  ))


# keep only necessary columns for phenology ----------------------------
dput(colnames(rangex_phenology_demo_traits))
rangex_phenology_clean <- rangex_phenology_demo_traits |>
  select(unique_plant_ID, species, date.x, recorder,
         number_buds, number_flowers.y,
         number_infructescences, seeds_collected, comment)



# rename columns ----------------------------------------------------------
rangex_phenology_clean <- rangex_phenology_clean |> 
  rename("date_measurement" = "date.x",
         "collector" = "recorder",
         "number_flowers" =  "number_flowers.y")



# Pivot longer the data ---------------------------------------------------
rangex_phenology_clean_long <- rangex_phenology_clean |> 
  pivot_longer(cols = c(number_buds, number_flowers, number_infructescences, seeds_collected),
               names_to = "phenology_stage",
               values_to = "value")



# make names to initials --------------------------------------------------
rangex_phenology_clean_long <- rangex_phenology_clean_long |>
  mutate(collector = case_when(
    collector %in% c("Dagmar") ~ "DE",
    collector %in% c("Nadine") ~ "NMA",
    collector %in% c("Lucas") ~ "LP", # Lucas Parra
    collector %in% c("Ingrid") ~ "IE", # Ingrid Espeland
    collector %in% c("Malo") ~ "MF", # Malo Le Fur
    collector %in% c("Aleksandra", "Ola") ~ "AP", # Aleksandra Posłuszny
    collector %in% c("Julia Filetti") ~ "JF",
    collector %in% c("Julia") ~ "JS", # Julia Schlick-Steiner
    collector %in% c("Lizzy") ~ "ED", # Elizabeth Duke-Moe
    collector %in% c("Nadine / Ingrid") ~ "NMA/IE",
    collector %in% c("Dagmar / Julia Filetti") ~ "DE/JF",
    collector %in% c("Malo / Lucas") ~ "MF/LP",
    collector %in% c("Dagmar / Lizzy") ~ "DE/ED",
    collector %in% c("Ingrid / Aleksandra") ~ "IE/AP",
    collector %in% c("Dagmar / Ingrid") ~ "DE/IE",
    collector %in% c("Dagmar /Alexandra") ~ "DE/AP",
      TRUE ~ collector
    ))


# no comment column -------------------------------------------------------
rangex_phenology_clean_long <- rangex_phenology_clean_long |>
  select(-comment)


# rename number_infructescences -------------------------------------------
rangex_phenology_clean_long <- rangex_phenology_clean_long |> 
  mutate(phenology_stage = if_else(
    phenology_stage == "number_infructescences",
    "No_Infructescences",
    phenology_stage
  ))


# save the clean data -----------------------------------------------------
write.csv(rangex_phenology_clean_long, file = "Data/Data_phenology/Clean_Phenology/RangeX_clean_Phenology_2023_NOR.csv", row.names = FALSE)

phenology <- read.csv("Data/Data_phenology/Clean_Phenology/RangeX_clean_Phenology_2023_NOR.csv")












