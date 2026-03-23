

# BIOMASS 5 ---------------------------------------------------------------


# Combine traits 23 with predicted biomass 23 --------------------------------

# source("RangeX_data_paper_cleaning_demographic_traits_23.R")
# # to use 
# rangex_traits_23_nor

# or import csv file


# import trait data 23 ----------------------------------------------------

rx_traits_23_nor <- read.csv("Data/RangeX_clean_YearlyDemographics_2023_NOR.csv")


# rename no_stems ---------------------------------------------------------
demo_traits_2023_bio2 <- rx_traits_23_nor |> 
  rename(no_stems = no_rep_stems)



# predict biomass 23 per species ------------------------------------------

# sucpra ------------------------------------------------------------------
df_2023_sucpra <- demo_traits_2023_bio2 |> filter(species == "sucpra")

df_2023_sucpra <- df_2023_sucpra |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_leaves)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_sucpra$pred_log_biomass <- predict(species_models$sucpra, newdata = df_2023_sucpra, re.form = NA)


# cennig ------------------------------------------------------------------
df_2023_cennig <- demo_traits_2023_bio2 |> filter(species == "cennig")

df_2023_cennig <- df_2023_cennig |> 
  mutate(
    log_height_reproductive_str = log1p(height_reproductive_str),
    log_no_stems                = log1p(no_stems)
  ) |> drop_na(log_height_reproductive_str, log_no_stems)

df_2023_cennig$pred_log_biomass <- predict(species_models$cennig, newdata = df_2023_cennig, re.form = NA)


# pimsax ------------------------------------------------------------------
df_2023_pimsax <- demo_traits_2023_bio2 |> filter(species == "pimsax")

df_2023_pimsax <- df_2023_pimsax |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_leaves)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_pimsax$pred_log_biomass <- predict(species_models$pimsax, newdata = df_2023_pimsax, re.form = NA)


# luzmul ------------------------------------------------------------------
# no leaves is actually tillers (bunch of leaves)

df_2023_luzmul <- demo_traits_2023_bio2 |> filter(species == "luzmul")

df_2023_luzmul <- df_2023_luzmul |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_tillers)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_luzmul$pred_log_biomass <- predict(species_models$luzmul, newdata = df_2023_luzmul, re.form = NA)


# leuvul ------------------------------------------------------------------
df_2023_leuvul <- demo_traits_2023_bio2 |> filter(species == "leuvul")

df_2023_leuvul <- df_2023_leuvul |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_leaves)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_leuvul$pred_log_biomass <- predict(species_models$leuvul, newdata = df_2023_leuvul, re.form = NA)


# tripra ------------------------------------------------------------------
df_2023_tripra <- demo_traits_2023_bio2 |> filter(species == "tripra")

df_2023_tripra <- df_2023_tripra |> 
  mutate(
    log_height_reproductive_str = log1p(height_reproductive_str),
    log_no_stems                = log1p(no_stems)
  ) |> drop_na(log_height_reproductive_str, log_no_stems)

df_2023_tripra$pred_log_biomass <- predict(species_models$tripra, newdata = df_2023_tripra, re.form = NA)


# hypmac ------------------------------------------------------------------
df_2023_hypmac <- demo_traits_2023_bio2 |> filter(species == "hypmac")

df_2023_hypmac <- df_2023_hypmac |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_leaves)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_hypmac$pred_log_biomass <- predict(species_models$hypmac, newdata = df_2023_hypmac, re.form = NA)


# plalan ------------------------------------------------------------------
df_2023_plalan <- demo_traits_2023_bio2 |> filter(species == "plalan")

df_2023_plalan <- df_2023_plalan |> 
  mutate(
    log_height_reproductive_str = log1p(height_reproductive_str),
    log_number_leaves           = log1p(number_leaves)
  ) |> drop_na(log_height_reproductive_str, log_number_leaves)

df_2023_plalan$pred_log_biomass <- predict(species_models$plalan, newdata = df_2023_plalan, re.form = NA)


# cyncri ------------------------------------------------------------------
# no leaves is actually tillers (bunch of leaves)

df_2023_cyncri <- demo_traits_2023_bio2 |> filter(species == "cyncri")

df_2023_cyncri <- df_2023_cyncri |> 
  mutate(
    log_no_stems      = log1p(no_stems),
    log_number_leaves = log1p(number_tillers)
  ) |> drop_na(log_no_stems, log_number_leaves)

df_2023_cyncri$pred_log_biomass <- predict(species_models$cyncri, newdata = df_2023_cyncri, re.form = NA)




# sildio ------------------------------------------------------------------
df_2023_sildio <- demo_traits_2023_bio2 |> filter(species == "sildio")

df_2023_sildio <- df_2023_sildio |> 
  mutate(
    log_number_leaves = log1p(number_leaves),
    log_no_stems      = log1p(no_stems)
  ) |> drop_na(log_number_leaves, log_no_stems)

df_2023_sildio$pred_log_biomass <- predict(species_models$sildio, newdata = df_2023_sildio, re.form = NA)


# combine all datasets ----------------------------------------------------

bio_pred_23_species <- bind_rows(
  df_2023_sucpra,
  df_2023_cennig,
  df_2023_pimsax,
  df_2023_luzmul,
  df_2023_leuvul,
  df_2023_tripra,
  df_2023_hypmac,
  df_2023_plalan,
  df_2023_cyncri,
  df_2023_sildio
)

bio_pred_23_species
# 1212 observations

















