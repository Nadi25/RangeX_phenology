

# BIOMASS 1 ---------------------------------------------------------------


# Data preparation - biomass 24 -------------------------------------------


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)


# comments ----------------------------------------------------------------
# low site: fix #Value and ?

# import metadata ---------------------------------------------------------
meta_NOR <- read.csv("Data/RangeX_clean_MetadataFocal_NOR.csv")


# import biomass data -------------------------------------------------------------
biomass_high <- read.csv("Data/Raw/RangeX_raw_NOR_biomass_high_2024_new.csv")

biomass_low <- read.csv("Data/Raw/RangeX_raw_NOR_biomass_low_2024_new.csv")

# remove empty row at end -------------------------------------------------
biomass_low <- biomass_low |>
  mutate(across(where(is.character), ~na_if(.x, "")))

biomass_low <- biomass_low |> 
  filter(!if_all(everything(), is.na))

# add column site ---------------------------------------------------------
biomass_high <- biomass_high |> 
  mutate(site = "hi")

biomass_low <- biomass_low |> 
  mutate(site = "lo")

# change to same as in high to merge
biomass_low <- biomass_low |> 
  rename("block" = "f")

biomass_high <- biomass_high |> 
  rename("block" = "X")

# delete superfluous columns -----------------------------------------------
biomass_high <- biomass_high |>
  select(-starts_with("X"))

biomass_low <- biomass_low |>
  select(-starts_with("X"))


# make weight cols numeric ----------------------------------------------------------
weight_cols <- c(
  "dry_weight_stem_g",
  "dry_weight_leaves_g",
  "dry_weight_flowers_g",
  "dry_weight_total_g")

biomass_high[weight_cols] <- lapply(biomass_high[weight_cols], \(x) as.numeric(x))
biomass_low[weight_cols]  <- lapply(biomass_low[weight_cols],  \(x) as.numeric(x))

# combine low and high ----------------------------------------------------
biomass <- bind_rows(biomass_high, biomass_low)



# rename columns to match meta data --------------------------------------------
biomass <- biomass |> 
  rename("block_ID_original" = "block",
         "plot_ID_original" = "treat",
         "position_ID_original" = "coord")


# merge biomass with metadata ---------------------------------------------
biomass_NOR <- left_join(meta_NOR, biomass, by = c("site", "block_ID_original",
                                                   "plot_ID_original",
                                                   "position_ID_original",
                                                   "species"))


# combined treatment column -----------------------------------------------
biomass_NOR$treatment <- paste(biomass_NOR$site, biomass_NOR$treat_warming, 
                               biomass_NOR$treat_competition, sep = "_")



# check structure ---------------------------------------------------------
glimpse(biomass_NOR)































