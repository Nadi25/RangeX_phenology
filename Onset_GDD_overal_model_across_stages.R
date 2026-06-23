





source("Onset_DOY_overal_model_across_stages.R")



# GDD ---------------------------------------------------------------------


# first onset for all individuals -----------------------------------------------
onset_all_gdd <- phenology_gdd_nor_che |>
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
    onset_gdd = min(GDD_cum),
    .groups = "drop"
  )


# make stage names the same as in temp ------------------------------------
onset_all_gdd <- onset_all_gdd |>
  mutate(
    stage = recode(
      phenology_stage,
      "No_Buds" = "bud",
      "No_FloOpen" = "flower",
      "No_FloWithrd" = "fruit",
      "No_Seeds" = "seed"))


# join temperature data with onset data -----------------------------------
onset_temp_gdd <- onset_all_gdd |>
  left_join(
    temp_all,
    by = c("region", "site", "treatment_site_temp",
           "treat_warming", "treat_competition", "stage"))



# make a model per region -----------------------------------------------

# Norway ------------------------------------------------------------------
onset_temp_gdd_nor <- onset_temp_gdd |> 
  filter(region == "Norway")

# are mean temperature per stage affect the onset across stages?
m_onset_temp_gdd_nor <- lmerTest::lmer(
  onset_gdd ~ stage * treat_competition * mean_temp +
    (1 | species) + (1 | block_ID),
  data = onset_temp_gdd_nor)

summary(m_onset_temp_gdd_nor)
# Plot random effects
dotplot(ranef(m_onset_temp_gdd_nor, condVar=TRUE))



# does site, warming and biotic interaction shift onset across stages?
m_onset_temp_gdd_nor2 <- lmerTest::lmer(
  onset_gdd ~ stage * treatment_site_temp * treat_competition +
    (1 | species) + (1 | block_ID),
  data = onset_temp_gdd_nor)
summary(m_onset_temp_gdd_nor2)

dotplot(ranef(m_onset_temp_gdd_nor2, condVar=TRUE))



# Switzerland -------------------------------------------------------------
onset_temp_gdd_che <- onset_temp_gdd |> 
  filter(region == "Switzerland")

# are mean temperature per stage affect the onset across stages?
m_onset_temp_gdd_che <- lmerTest::lmer(
  onset_gdd ~ stage * treat_competition * mean_temp +
    (1 | species) + (1 | block_ID),
  data = onset_temp_gdd_che)

summary(m_onset_temp_gdd_che)



# does site, warming and biotic interaction shift onset across stages?
m_onset_temp_gdd_che2 <- lmerTest::lmer(
  onset_gdd ~ stage * treatment_site_temp * treat_competition +
    (1 | species) + (1 | block_ID),
  data = onset_temp_gdd_che)
summary(m_onset_temp_gdd_che2)

# Plot random effects
dotplot(ranef(m_onset_temp_gdd_che, condVar=TRUE))





