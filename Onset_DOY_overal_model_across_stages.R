

# 01.2 Onset julian days overall model -------------------------------------------------------

# overall model across stages
# Onset bud, flower, fruit model with site_treat_warming combination --------
# julian days


library(emmeans)


# source script with final tms data ---------------------------------------
# this one sources the NOR and CHE preparation scripts
source("DOY_GDD_TMS4_NOR_CHE.R")


# add region and combine the datasets -------------------------------------

# > season_start_nor_tms
# # A tibble: 6 × 4
# site  treat_warming treat_competition season_start
# <chr> <chr>         <chr>             <date>      
# 1 hi    ambi          bare              2023-04-12  
# 2 hi    ambi          vege              2023-04-12  
# 3 hi    warm          bare              2023-04-11  
# 4 hi    warm          vege              2023-04-11  
# 5 lo    ambi          bare              2023-04-11  
# 6 lo    ambi          vege              2023-04-11 
season_start_nor_tms$region <- "Norway"

# season_start_che_tms
# # A tibble: 6 × 4
# site  treat_warming treat_competition season_start
# <chr> <chr>         <chr>             <date>      
# 1 hi    ambi          bare              2022-04-16  
# 2 hi    ambi          vege              2022-04-14  
# 3 hi    warm          bare              2022-04-16  
# 4 hi    warm          vege              2022-04-13  
# 5 lo    ambi          bare              2022-03-27  
# 6 lo    ambi          vege              2022-03-27  
season_start_che_tms$region <- "Switzerland"



# combine nor and che tms data --------------------------------------------
# and change names of treat_competition to match
# this is the season start dataset for
season_start_all <- bind_rows(
  season_start_nor_tms,
  season_start_che_tms) |> 
  mutate(treat_competition = recode(treat_competition,
                               "bare" = "without",
                               "vege" = "with"))


tms_final_nor_che <- tms_final_nor_che |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))


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

phenology4 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional 



phenology4 <- phenology4 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))

phenology4$treatment_site_temp <- paste(phenology4$site, phenology4$treat_warming, sep = "_")




# Get first onset per individual ------------------------------------------
onset_ind <- phenology4 |>
  filter(value > 0) |>
  group_by(
    region,
    site,
    treat_competition,
    treat_warming,
    phenology_stage,
    unique_plant_ID
  ) |>
  summarise(
    onset = min(date_measurement),
    .groups = "drop"
  )


# get mean onset per treatment --------------------------------------------
onset_treat <- onset_ind |>
  group_by(
    region, site,
    treat_competition, treat_warming,
    phenology_stage
  ) |>
  summarise(
    onset = as.Date(mean(onset)),
    .groups = "drop"
  )

# make wide table
onset_wide <- onset_treat |>
  select(region, site, treat_competition,
         treat_warming, phenology_stage, onset) |>
  pivot_wider(
    names_from = phenology_stage,
    values_from = onset)



# Get all onsets together -------------------------------------------------
# combine season start tables and also with onset table
onset_wide2 <- onset_wide |>
  left_join(
    season_start_all,
    by = c("region", "site", "treat_warming", "treat_competition")
  )

stage_windows <- onset_wide2 |>
  mutate(
    bud_start = season_start,
    bud_end = No_Buds,
    
    flower_start = No_Buds,
    flower_end = No_FloOpen,
    
    fruit_start = No_FloOpen,
    fruit_end = No_FloWithrd,
    
    seed_start = No_FloWithrd,
    seed_end = No_Seeds
  )



# Function to calculate mean temp per stage -------------------------------
calc_mean_temp <- function(data, windows, start_col, end_col, label) {
  
  data |>
    inner_join(
      windows,
      by = c("region", "site", "treat_warming", "treat_competition")
    ) |>
    filter(
      date >= .data[[start_col]],
      date <= .data[[end_col]]
    ) |>
    group_by(region, site, treat_warming, treat_competition) |>
    summarise(
      mean_temp = mean(temp_mean, na.rm = TRUE),
      stage = label,
      .groups = "drop"
    )
}


# apply the function per stage -------------------------------------------
# bud
temp_bud <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "bud_start",
  "bud_end",
  "bud"
)

# flower
temp_flower <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "flower_start",
  "flower_end",
  "flower"
)

# fruit
temp_fruit <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "fruit_start",
  "fruit_end",
  "fruit"
)

# seed
temp_seed <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "seed_start",
  "seed_end",
  "seed"
)


# mean temp per stage and treatment - combine all stages ------------------------------------------------------
temp_all <- bind_rows(temp_bud, temp_flower, temp_fruit, temp_seed)
temp_all

# add combined treatment
temp_all$treatment_site_temp <- paste(temp_all$site, temp_all$treat_warming, sep = "_")

# plot
ggplot(temp_all,
       aes(x = stage,
           y = mean_temp,
           color = treat_warming,
           group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(region + site ~ treat_competition) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature during stage",
    color = "Warming"
  )






# first onset for all individuals -----------------------------------------------
onset_all <- phenology4 |>
  filter(value > 0) |>
  group_by(
    region,
    site,
    treat_warming,
    treat_competition,
    treatment_site_temp,
    species,
    block_ID,
    unique_plot_ID,
    unique_plant_ID,
    phenology_stage
  ) |>
  summarise(
    onset = min(jday),
    .groups = "drop"
  )


# make stage names the same as in temp ------------------------------------
onset_all <- onset_all |>
  mutate(
    stage = recode(
      phenology_stage,
      "No_Buds" = "bud",
      "No_FloOpen" = "flower",
      "No_FloWithrd" = "fruit",
      "No_Seeds" = "seed"))


# join temperature data with onset data -----------------------------------
onset_temp <- onset_all |>
  left_join(
    temp_all,
    by = c("region", "site", "treatment_site_temp",
           "treat_warming", "treat_competition", "stage"))



# make a model per region -----------------------------------------------

# Norway ------------------------------------------------------------------
onset_temp_nor <- onset_temp |> 
  filter(region == "Norway")

# are mean temperature per stage affect the onset across stages?
m_onset_temp_nor <- lmerTest::lmer(
  onset ~ stage * treat_competition * mean_temp +
    (1 | species) + (1 | block_ID),
  data = onset_temp_nor)

summary(m_onset_temp_nor)



# does site, warming and biotic interaction shift onset across stages?
m_onset_temp_nor2 <- lmerTest::lmer(
  onset ~ stage * treatment_site_temp * treat_competition +
    (1 | species) + (1 | block_ID),
  data = onset_temp_nor)
summary(m_onset_temp_nor2)



# Switzerland -------------------------------------------------------------
onset_temp_che <- onset_temp |> 
  filter(region == "Switzerland")

# are mean temperature per stage affect the onset across stages?
m_onset_temp_che <- lmerTest::lmer(
  onset ~ stage * treat_competition * mean_temp +
    (1 | species) + (1 | block_ID),
  data = onset_temp_che)

summary(m_onset_temp_che)



# does site, warming and biotic interaction shift onset across stages?
m_onset_temp_che2 <- lmerTest::lmer(
  onset ~ stage * treatment_site_temp * treat_competition +
    (1 | species) + (1 | block_ID),
  data = onset_temp_che)
summary(m_onset_temp_che2)



