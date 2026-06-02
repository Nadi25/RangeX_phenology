

# 02 Onset GDD ------------------------------------------------------------

# Effect of upslope transplantation on onsets GDD  -------------------------------------------------------------

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(performance)
library(see)
library(emmeans)


source("Data_preparation_phenology_NOR_CHE_combined.R")

theme_set(theme_bw())

# filter only nor ------------------------------------------------------
phenology_nor <- phenology |> 
  filter(region == "Norway")


unique(phenology_nor$date_measurement)
# start: 2023-05-12
# end: 2023-09-13


source("Data_preparation_climate_station_NOR.R")

# use climate_gdd instead of climate_gdd_pt because it is without filtering
climate_gdd

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd <- phenology_nor |> 
  left_join(climate_gdd |> 
              select(site, date_measurement, GDD_cum, Tavg),
            by = c("site", "date_measurement"))



# filter only che ------------------------------------------------------
phenology_che <- phenology |> 
  filter(region == "Switzerland")


unique(phenology_che$date_measurement)
# start: "2022-05-04"
# end: "2022-09-27"


source("Data_preparation_climate_station_CHE.R")

# use 
# climate_gdd_che
# uses tms data for beginning until May and then climate station
climate_gdd_che_comb

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd_che <- phenology_che |> 
  left_join(climate_gdd_che_comb |> 
              select(site, date_measurement, GDD_cum, Tavg),
            by = c("site", "date_measurement"))




# combine nor and che -----------------------------------------------------
phenology_gdd_nor_che <- rbind(phenology_with_gdd, phenology_with_gdd_che)



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology_gdd_nor_che <- phenology_gdd_nor_che |> 
  filter(treat_warming == "ambi")

# and get julian days --------------------------------------
phenology_gdd_nor_che <- phenology_gdd_nor_che |> 
  mutate(jday = yday(date_measurement))   # Julian day (1–365)       


# one model per stage -----------------------------------------------------

# budding -----------------------------------------------------------------

# calculate budding onset ------------------------------------------------
budding_onset_gdd_ambi <- phenology_gdd_nor_che |> 
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))


m_gdd_bud <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                          (1|species) + (1|block_ID),
                        data = budding_onset_gdd_ambi)
summary(m_gdd_bud)
car::Anova(m_gdd_bud)

model_performance(m_gdd_bud)
check_collinearity(m_gdd_bud)
#check_model(m_gdd_bud)



# flowering ---------------------------------------------------------------

# calculate flowering onset ------------------------------------------------
flowering_onset_gdd_ambi <- phenology_gdd_nor_che |> 
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))



# model flowering onset nor che ------------------------------------------
m_gdd_flower <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                            (1|species) + (1|block_ID),
                                          data = flowering_onset_gdd_ambi)
summary(m_gdd_flower)


model_performance(m_gdd_flower)
check_collinearity(m_gdd_flower)
#check_model(m_gdd_flower)



# fruiting ----------------------------------------------------------------

# combine the fruiting stages nor and che to comapre the onset ------------
phenology_gdd_nor_che <- phenology_gdd_nor_che |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))


# calculate fruiting onset ------------------------------------------------
fruiting_onset_gdd_ambi <- phenology_gdd_nor_che |> 
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))



# model flowering onset nor che ------------------------------------------
m_gdd_fruit <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                           (1|species) + (1|block_ID),
                                         data = fruiting_onset_gdd_ambi)
summary(m_gdd_fruit)


model_performance(m_gdd_fruit)
check_collinearity(m_gdd_fruit)
#check_model(m_gdd_fruit)


# seeds ----------------------------------------------------------------

# calculate seed onset ------------------------------------------------
seed_onset_gdd_ambi <- phenology_gdd_nor_che |> 
  filter(phenology_stage == "No_Seeds", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))



# model flowering onset nor che ------------------------------------------
m_gdd_seed <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                (1|species) + (1|block_ID),
                              data = seed_onset_gdd_ambi)
summary(m_gdd_seed)


model_performance(m_gdd_seed)
check_collinearity(m_gdd_seed)
#check_model(m_gdd_seed)




# make predictions for each stage -----------------------------------------

# bud ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_gdd_bud <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("lo", "hi"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_gdd_bud



# predict  ----------------------------------------------------------------
newdat_gdd_bud$onset <- predict(
  m_gdd_bud,
  newdata = newdat_gdd_bud,
  re.form = NA 
)
newdat_gdd_bud


# make model matrix -------------------------------------------------------
mm_gdd_bud <- model.matrix(terms(m_gdd_bud), newdat_gdd_bud)
mm_gdd_bud


pvar1_gdd_bud <- diag(mm_gdd_bud %*% tcrossprod(vcov(m_gdd_bud), mm_gdd_bud))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_gdd_bud <- as.numeric(VarCorr(m_gdd_bud)$species)
var_block_gdd_bud <- as.numeric(VarCorr(m_gdd_bud)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_gdd_bud <- pvar1_gdd_bud + var_species_gdd_bud + var_block_gdd_bud + sigma(m_gdd_bud)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_gdd_bud <- data.frame(
  newdat_gdd_bud
  , plo = newdat_gdd_bud$onset-cmult*sqrt(pvar1_gdd_bud) # fixed effects only, confidence interval
  , phi = newdat_gdd_bud$onset+cmult*sqrt(pvar1_gdd_bud)
  , tlo = newdat_gdd_bud$onset-cmult*sqrt(tvar1_gdd_bud) # takes fixed and random effects, prediction intervall
  , thi = newdat_gdd_bud$onset+cmult*sqrt(tvar1_gdd_bud)
)
newdat_gdd_bud


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

b <- ggplot(newdat_gdd_bud, aes(x=treat_competition, y=onset, 
                           color=treat_competition, shape = site )) +
  geom_point(position = pd)+
  facet_wrap(~ region)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset GDD (budding)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(b)



# flower ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_gdd_flower <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("lo", "hi"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_gdd_flower



# predict  ----------------------------------------------------------------
newdat_gdd_flower$onset <- predict(
  m_gdd_flower,
  newdata = newdat_gdd_flower,
  re.form = NA 
)
newdat_gdd_flower


# make model matrix -------------------------------------------------------
mm_gdd_flower <- model.matrix(terms(m_gdd_flower), newdat_gdd_flower)
mm_gdd_flower


pvar1_gdd_flower <- diag(mm_gdd_flower %*% tcrossprod(vcov(m_gdd_flower), mm_gdd_flower))
# tvar1 <- pvar1+VarCorr(m_onset_flower_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_gdd_flower <- as.numeric(VarCorr(m_gdd_flower)$species)
var_block_gdd_flower <- as.numeric(VarCorr(m_gdd_flower)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_gdd_flower <- pvar1_gdd_flower + var_species_gdd_flower + var_block_gdd_flower + sigma(m_gdd_flower)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_gdd_flower <- data.frame(
  newdat_gdd_flower
  , plo = newdat_gdd_flower$onset-cmult*sqrt(pvar1_gdd_flower) # fixed effects only, confidence interval
  , phi = newdat_gdd_flower$onset+cmult*sqrt(pvar1_gdd_flower)
  , tlo = newdat_gdd_flower$onset-cmult*sqrt(tvar1_gdd_flower) # takes fixed and random effects, prediction intervall
  , thi = newdat_gdd_flower$onset+cmult*sqrt(tvar1_gdd_flower)
)
newdat_gdd_flower


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

f <- ggplot(newdat_gdd_flower, aes(x=treat_competition, y=onset, 
                                color=treat_competition, shape = site )) +
  geom_point(position = pd)+
  facet_wrap(~ region)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset GDD (flowering)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(f)




# fruit ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_gdd_fruit <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("lo", "hi"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_gdd_fruit



# predict  ----------------------------------------------------------------
newdat_gdd_fruit$onset <- predict(
  m_gdd_fruit,
  newdata = newdat_gdd_fruit,
  re.form = NA 
)
newdat_gdd_fruit


# make model matrix -------------------------------------------------------
mm_gdd_fruit <- model.matrix(terms(m_gdd_fruit), newdat_gdd_fruit)
mm_gdd_fruit


pvar1_gdd_fruit <- diag(mm_gdd_fruit %*% tcrossprod(vcov(m_gdd_fruit), mm_gdd_fruit))
# tvar1 <- pvar1+VarCorr(m_onset_fruit_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_gdd_fruit <- as.numeric(VarCorr(m_gdd_fruit)$species)
var_block_gdd_fruit <- as.numeric(VarCorr(m_gdd_fruit)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_gdd_fruit <- pvar1_gdd_fruit + var_species_gdd_fruit + var_block_gdd_fruit + sigma(m_gdd_fruit)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_gdd_fruit <- data.frame(
  newdat_gdd_fruit
  , plo = newdat_gdd_fruit$onset-cmult*sqrt(pvar1_gdd_fruit) # fixed effects only, confidence interval
  , phi = newdat_gdd_fruit$onset+cmult*sqrt(pvar1_gdd_fruit)
  , tlo = newdat_gdd_fruit$onset-cmult*sqrt(tvar1_gdd_fruit) # takes fixed and random effects, prediction intervall
  , thi = newdat_gdd_fruit$onset+cmult*sqrt(tvar1_gdd_fruit)
)
newdat_gdd_fruit


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

fr <- ggplot(newdat_gdd_fruit, aes(x=treat_competition, y=onset, 
                                   color=treat_competition, shape = site )) +
  geom_point(position = pd)+
  facet_wrap(~ region)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset GDD (fruiting)",
    title = "Effect of transplantation on fruiting onset across regions"
  )
print(fr)



# seed ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_gdd_seed <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("lo", "hi"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_gdd_seed



# predict  ----------------------------------------------------------------
newdat_gdd_seed$onset <- predict(
  m_gdd_seed,
  newdata = newdat_gdd_seed,
  re.form = NA 
)
newdat_gdd_seed


# make model matrix -------------------------------------------------------
mm_gdd_seed <- model.matrix(terms(m_gdd_seed), newdat_gdd_seed)
mm_gdd_seed


pvar1_gdd_seed <- diag(mm_gdd_seed %*% tcrossprod(vcov(m_gdd_seed), mm_gdd_seed))
# tvar1 <- pvar1+VarCorr(m_onset_fruit_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_gdd_seed <- as.numeric(VarCorr(m_gdd_seed)$species)
var_block_gdd_seed <- as.numeric(VarCorr(m_gdd_seed)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_gdd_seed <- pvar1_gdd_seed + var_species_gdd_seed + var_block_gdd_seed + sigma(m_gdd_seed)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_gdd_seed <- data.frame(
  newdat_gdd_seed
  , plo = newdat_gdd_seed$onset-cmult*sqrt(pvar1_gdd_fruit) # fixed effects only, confidence interval
  , phi = newdat_gdd_seed$onset+cmult*sqrt(pvar1_gdd_fruit)
  , tlo = newdat_gdd_seed$onset-cmult*sqrt(tvar1_gdd_fruit) # takes fixed and random effects, prediction intervall
  , thi = newdat_gdd_seed$onset+cmult*sqrt(tvar1_gdd_fruit)
)
newdat_gdd_seed


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

s <- ggplot(newdat_gdd_seed, aes(x=treat_competition, y=onset, 
                                   color=treat_competition, shape = site )) +
  geom_point(position = pd)+
  facet_wrap(~ region)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset GDD (seed)",
    title = "Effect of transplantation on seed onset across regions"
  )
print(s)



# combine all stages into one plot ----------------------------------------

stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636",
  "Seeds" = "grey50"  
)

plot_df_gdd_bud  <- newdat_gdd_bud   |> mutate(stage = "Budding")

plot_df_gdd_flower <- newdat_gdd_flower |> mutate(stage = "Flowering")

plot_df_gdd_fruit <- newdat_gdd_fruit |> mutate(stage = "Fruiting")

plot_df_gdd_seed <- newdat_gdd_seed |> mutate(stage = "Seeds")


plot_df_gdd_all <- bind_rows(
  plot_df_gdd_bud,
  plot_df_gdd_flower,
  plot_df_gdd_fruit,
  plot_df_gdd_seed
)
plot_df_gdd_all





# plot facet by stage --------------------------------------------------
region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

b_f_fr_gdd <- ggplot(plot_df_gdd_all, aes(
  x = treat_competition,
  y = onset,
  color = region,              
  shape = site  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_grid(~ stage) +     
  
  scale_color_manual(values = region_colors) +   
  
  scale_shape_manual(
    values = c(
      "lo" = 16,   # circle
      "hi" = 17    # triangle
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (GDD)",
    title = "Effect of transplantation on onset across regions",
    shape = "Site",
    color = "Region"
  )

print(b_f_fr_gdd)

# ggsave(filename = "Output/Onset/GDD_Transplantation_Warming_Onset_bud_flower_fruit_seed_predictions_tbase2.png", 
#        plot = b_f_fr_gdd, width = 18, height = 10, units = "in")




b_f_fr_gdd2 <- ggplot(plot_df_gdd_all, aes(
  x = treat_competition,
  y = onset,
  color = treat_competition,              
  shape = site  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_grid(region ~ stage) +
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  #scale_color_manual(values = region_colors) +   
  
  scale_shape_manual(
    values = c(
      "lo" = 16,   # circle
      "hi" = 17    # triangle
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (GDD2)",
    title = "Effect of transplantation on onset GDD",
    shape = "Site",
    color = "Region"
  )

print(b_f_fr_gdd2)


# ggsave(filename = "Output/Onset/GDD_Transplantation_Warming_Onset_bud_flower_fruit_seed_predictions_tbase2_separate_region.png", 
#        plot = b_f_fr_gdd2, width = 18, height = 10, units = "in")



