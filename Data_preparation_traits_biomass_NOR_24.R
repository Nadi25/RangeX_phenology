
# BIOMASS 2 ---------------------------------------------------------------

# Data preparation - traits 24 --------------------------------------------

# Combine traits 24 with biomass 24 data ----------------------------------

# comments ----------------------------------------------------------------
# some leaves have been collected for functional traits (column leaf.collected)
# they should be included in the biomass
# are weighed but not digitized yet

# lo 3a tripra i9: probably seedlings from i8 --> deleted

# import metadata NOR -----------------------------------------------------
meta_NOR <- read.csv("Data/RangeX_clean_MetadataFocal_NOR.csv")


# import trait data 2024 --------------------------------------------------
trait_hi <- read.csv("Data/Raw/RangeX_raw_demographic_traits_high_2024.csv")

trait_lo <- read.csv("Data/Raw/RangeX_raw_demographic_traits_low_2024.csv")


# make heights numeric ----------------------------------------------------
trait_hi <- trait_hi |> 
  mutate(
    across(
      starts_with("height"),
      ~ .x |>
        as.numeric()
    )
  )

trait_lo <- trait_lo |> 
  mutate(
    across(
      starts_with("height"),
      ~ .x |>
        as.numeric()
    )
  )

trait_lo <- trait_lo |> 
  rename("block" = "X")

# remove empty row at end -------------------------------------------------
trait_lo <- trait_lo |>
  mutate(across(where(is.character), ~na_if(.x, "")))

trait_lo <- trait_lo |> 
  filter(!if_all(everything(), is.na))



# add column site ---------------------------------------------------------
trait_hi <- trait_hi |> 
  mutate(site = "hi")

trait_lo <- trait_lo |> 
  mutate(site = "lo")



# fix typos ---------------------------------------------------------------
trait_hi <- trait_hi |> 
  mutate(coord = case_when(
    coord == "f7-" ~ "f7",
    TRUE ~ coord))


# delete lo 3a tripra i9 --------------------------------------------------
# seedlings from i8?
trait_lo <- trait_lo |>
  filter(!(site == "lo" &
             block == 3 &
             treat == "a" &
             species == "tripra" &
             coord == "i9"))



# combine hi and lo -------------------------------------------------------
traits <- bind_rows(trait_hi, trait_lo)


# rename column names -----------------------------------------------------
traits <- traits |> 
  rename("block_ID_original" = "block",
         "plot_ID_original" = "treat",
         "position_ID_original" = "coord")


# merge traits with metadata ---------------------------------------------
traits_NOR <- left_join(meta_NOR, traits, by = c("site", "block_ID_original",
                                                 "plot_ID_original",
                                                 "position_ID_original",
                                                 "species"))


# combined treatment column -----------------------------------------------
traits_NOR$treatment <- paste(traits_NOR$site, traits_NOR$treat_warming, 
                              traits_NOR$treat_competition, sep = "_")




# replace # leaves and stems with the harvested counts -----------------------------------------
traits_NOR <- traits_NOR |>
  mutate(
    no..leaves.tillers.if.different...harvest = 
      as.numeric(no..leaves.tillers.if.different...harvest),
    no..leaves.tillers = if_else(
      !is.na(no..leaves.tillers.if.different...harvest),
      no..leaves.tillers.if.different...harvest,
      no..leaves.tillers
    )
  )


traits_NOR <- traits_NOR |>
  mutate(
    no..stems.if.different...harvest = 
      as.numeric(no..stems.if.different...harvest),
    no..stems = if_else(
      !is.na(no..stems.if.different...harvest),
      no..stems.if.different...harvest,
      no..stems
    )
  )



# rename columns ----------------------------------------------------------
traits_NOR <- traits_NOR |>
  rename(
    # Heights
    height_vegetative       = height.veg..cm.,
    height_vegetative_str   = height.veg.stretch..cm.,
    height_reproductive     = height.rep..cm.,
    height_reproductive_str = height.rep.stretch..cm.,
    height_Nathan_stretch_cm = height.Nathan.stretch..cm.,
    
    # Leaves / tillers
    number_leaves_tillers = no..leaves.tillers,
    
    # Stems
    no_stems = no..stems,
    
    # Herbivory
    herbivory_leaf = herbivory.leaf,
    herbivory_flower  = herbivory.flower,
    herbivory  = herbivory.if.different...harvest,
    
    sampled_quarter = sampled.quarter,
    sampled_half = sampled.half)



# sampled quarter/half ---------------------------------------------------------
# where sampled quarter, multiply by 4
# half by 2
traits_NOR$number_leaves <- ifelse(
  traits_NOR$sampled_quarter == "yes",
  traits_NOR$number_leaves_tillers * 4,
  ifelse(
    traits_NOR$sampled_half == "yes",
    traits_NOR$number_leaves_tillers * 2,
    traits_NOR$number_leaves_tillers))

# so number_leaves is the actual # leaves to use



# "clean" traits data --------------------------------------------------------------

traits_NOR_24 <- traits_NOR



# combine biomass with traits ---------------------------------------------
key <- "unique_plant_ID"

traits_vars <- traits_NOR_24 |>
  select(
    unique_plant_ID,
    height_vegetative,
    height_vegetative_str,
    height_reproductive,
    height_reproductive_str,
    height_Nathan_stretch_cm,
    number_leaves_tillers,
    number_leaves,
    no..leaves.tillers.if.different...harvest,
    no_stems,
    no..stems.if.different...harvest,
    herbivory_leaf,
    herbivory_flower
  )

# join biomass with traits
biomass_traits_NOR <- biomass_NOR |>
  left_join(traits_vars, by = key)


# filter out plants where weight is 0 -------------------------------------
biomass_traits_NOR_24 <- biomass_traits_NOR |>
  filter(dry_weight_total_g > 0)

# test distribution of biomass total
hist(biomass_traits_NOR_24$dry_weight_total_g)
hist(log(biomass_traits_NOR_24$dry_weight_total_g))


qqnorm(biomass_traits_NOR_24$dry_weight_total_g); qqline(biomass_traits_NOR_24$dry_weight_total_g)
qqnorm(log(biomass_traits_NOR_24$dry_weight_total_g)); qqline(log(biomass_traits_NOR_24$dry_weight_total_g))


# log biomass -------------------------------------------------------------
biomass_traits_NOR_24 <- biomass_traits_NOR_24 |>
  mutate(log_biomass = log(dry_weight_total_g))


## that leaves 1643 plants with biomass






