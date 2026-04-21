


# Temperature sensitivity -------------------------------------------------

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(lme4)
library(performance)
library(see)
library(emmeans)

# calcualte delta T -------------------------------------------------------

# mean temperature per site
# timeframe before the first phenological event

# can only do from mid april because che starts from 13.04.22

# che hi:	13/05/2022	12/10/2022
# che lo:	12/04/2022	18/10/2022


# when does onset start at each site
# first_onset_bud
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         158 = 07.06
# 2 Norway      lo     2023         132 = 12.05
# 3 Switzerland hi     2022         134 = 14.05
# 4 Switzerland lo     2022         124 = 04.05

# first_onset_flower
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         158 = 07.06
# 2 Norway      lo     2023         144 = 24.05
# 3 Switzerland hi     2022         134 = 14.05
# 4 Switzerland lo     2022         124 = 04.05

# first_onset_fruit
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         171
# 2 Norway      lo     2023         170
# 3 Switzerland hi     2022         150
# 4 Switzerland lo     2022         135 = 15.05


# source climate scripts --------------------------------------------------
source("Data_preparation_climate_station_CHE.R")

source("Data_preparation_climate_station_NOR.R")

# combine nor and che -----------------------------------------------------
climate_23
climate_che_22

# filter time period
climate_nor_23_pre <- climate_23 |>
  filter(
    date_time >= as.Date("2023-05-01"),
    date_time <= as.Date("2023-05-31")
  )

climate_che_22_pre <- climate_che_22 |>
  filter(
    date_time >= as.Date("2022-05-14"),
    date_time <= as.Date("2022-06-13")
  )

# climate_che_22_pre <- climate_che_22 |>
#   filter(
#     date_time >= as.Date("2022-05-14"),
#     date_time <= as.Date("2022-05-31")
#   )


# 23 day time frame before the first budding onset
# except for che hi
climate_che_22_pre <- climate_che_22 |>
  filter(
    (site == "lo" & date_time >= as.Date("2022-04-12") & date_time <= as.Date("2022-05-04")) |
      (site == "hi" & date_time >= as.Date("2022-05-13") & date_time <= as.Date("2022-06-05"))
  )

climate_nor_23_pre <- climate_23 |>
  filter(
    (site == "lo" & date_time >= as.Date("2023-04-19") & date_time <= as.Date("2023-05-12")) |
      (site == "hi" & date_time >= as.Date("2023-05-15") & date_time <= as.Date("2023-06-07"))
  )


# 14 day time frame before the first budding onset
# except for che hi
climate_che_22_pre <- climate_che_22 |>
  filter(
    (site == "lo" & date_time >= as.Date("2022-04-20") & date_time <= as.Date("2022-05-04")) |
      (site == "hi" & date_time >= as.Date("2022-04-30") & date_time <= as.Date("2022-06-05"))
  )

climate_nor_23_pre <- climate_23 |>
  filter(
    (site == "lo" & date_time >= as.Date("2023-04-28") & date_time <= as.Date("2023-05-12")) |
      (site == "hi" & date_time >= as.Date("2023-05-24") & date_time <= as.Date("2023-06-07"))
  )


# combine
climate_all <- bind_rows(climate_che_22_pre, climate_nor_23_pre)

climate_all <- climate_all |> 
  select(region, site, year, date_time, AirTemp_Avg, Humidity_Avg, WindDir_Avg, WindSpd_Avg, Radiation_Avg, Rainfall)


# calculate mean per site in this time period
pre_climate <-  climate_all |>
  mutate(date = as.Date(date_time)) |>
  group_by(region, site, year, date) |>
  summarise(T_day = mean(AirTemp_Avg, na.rm = TRUE)) |>
  group_by(region, site, year) |>
  summarise(Tmean = mean(T_day))
pre_climate

pre_climate <- pre_climate |>
  mutate(region = case_when(
    region == "NOR" ~ "Norway",
    region == "CHE" ~ "Switzerland",
    TRUE ~ region
  ))
pre_climate



# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))


# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology3 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional scaling if needed



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology3 <- phenology3 |> 
  filter(treat_warming == "ambi")



phenology3 <- phenology3 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))


# calculate onsets ---------------------------------------------------------


# bud ---------------------------------------------------------------------
onset_bud <- phenology3 |>
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))

first_onset_bud <- onset_bud |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_bud


# flower ---------------------------------------------------------------------
onset_flower <- phenology3 |>
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))

first_onset_flower <- onset_flower |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_flower

# fruit ---------------------------------------------------------------------
# combine the fruiting stages nor and che to comapre the onset ------------
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology3 <- phenology3 |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))

onset_fruit <- phenology3 |>
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where fruiting never occurred
  filter(is.finite(onset))

first_onset_fruit <- onset_fruit |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_fruit

# seed ---------------------------------------------------------------------
onset_seed <- phenology3 |>
  filter(phenology_stage == "No_Seeds", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where seeds never occurred
  filter(is.finite(onset))




# Average onset per site, species and treatment --------------------------------------------------
# could do per plot as well
onset_bud_mean <- onset_bud |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_bud_mean

onset_flower_mean <- onset_flower |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_flower_mean

onset_fruit_mean <- onset_fruit |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_fruit_mean

onset_seed_mean <- onset_seed |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_seed_mean



# Calculate temperature sensitivity ---------------------------------------


# bud ---------------------------------------------------------------------
sens_bud <- onset_bud_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_bud


sens_flower <- onset_flower_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_flower

sens_fruit <- onset_fruit_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_fruit

sens_seed <- onset_seed_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_seed



sens_all <- bind_rows(
  sens_bud   |> mutate(stage = "bud"),
  sens_flower|> mutate(stage = "flower"),
  sens_fruit |> mutate(stage = "fruit"),
  sens_seed  |> mutate(stage = "seed")
)
sens_all





# plot  -------------------------------------------------------------------

ggplot(sens_all,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 

sens_summary <- sens_all |>
  group_by(region, stage, treat_competition) |>
  summarise(
    mean = mean(temp_sens, na.rm = TRUE),
    se   = sd(temp_sens, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
sens_summary


region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

ggplot(sens_summary,
       aes(x = treat_competition,
           y = mean,
           color = region)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se,
        ymax = mean + 1.96 * se),
    width = 0.2,
    position = position_dodge(0.4)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid( ~ stage)+
  labs(y = "Temperature sensitivity")+
  scale_color_manual(values = region_colors)





# sensitivity models ------------------------------------------------------------------

# bud
m_sens_bud <- lmerTest::lmer(temp_sens ~ region * treat_competition + 
                            (1|species),
                          data = sens_bud)
summary(m_sens_bud)
anova(m_sens_bud)

model_performance(m_sens_bud)

# flower
m_sens_flower <- lmerTest::lmer(temp_sens ~ region * treat_competition + 
                               (1|species),
                             data = sens_flower)
summary(m_sens_flower)
anova(m_sens_flower)

model_performance(m_sens_flower)

# fruit
m_sens_fruit <- lmerTest::lmer(temp_sens ~ region * treat_competition + 
                                  (1|species),
                                data = sens_fruit)
summary(m_sens_fruit)
anova(m_sens_fruit)

model_performance(m_sens_fruit)

# seeds
m_sens_seed <- lmerTest::lmer(temp_sens ~ region * treat_competition + 
                                 (1|species),
                               data = sens_seed)
summary(m_sens_seed)
anova(m_sens_seed)

model_performance(m_sens_seed)


# make predictions for each stage -----------------------------------------

# bud ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_sens_bud <- expand.grid(
  region = c("Norway", "Switzerland"),
  treat_competition = c("with", "without"),
  temp_sens = 0
)
newdat_sens_bud



# predict  ----------------------------------------------------------------
newdat_sens_bud$temp_sens <- predict(
  m_sens_bud,
  newdata = newdat_sens_bud,
  re.form = NA 
)
newdat_sens_bud


# make model matrix -------------------------------------------------------
mm_sens_bud <- model.matrix(terms(m_sens_bud), newdat_sens_bud)
mm_sens_bud


pvar1_sens_bud <- diag(mm_sens_bud %*% tcrossprod(vcov(m_sens_bud), mm_sens_bud))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
tvar1_sens_bud <- as.numeric(VarCorr(m_sens_bud)$species)
#var_block_sens_bud <- as.numeric(VarCorr(m_sens_bud)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
#tvar1_sens_bud <- pvar1_sens_bud + var_species_sens_bud + var_block_sens_bud + sigma(m_sens_bud)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_sens_bud <- data.frame(
  newdat_sens_bud
  , plo = newdat_sens_bud$temp_sens-cmult*sqrt(pvar1_sens_bud) # fixed effects only, confidence interval
  , phi = newdat_sens_bud$temp_sens+cmult*sqrt(pvar1_sens_bud)
  , tlo = newdat_sens_bud$temp_sens-cmult*sqrt(tvar1_sens_bud) # takes fixed and random effects, prediction intervall
  , thi = newdat_sens_bud$temp_sens+cmult*sqrt(tvar1_sens_bud)
)
newdat_sens_bud


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

b<- ggplot(newdat_sens_bud, aes(x=treat_competition, y= temp_sens, 
                           color=treat_competition)) +
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
    y = "Predicted temperature sensitivity (budding)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(b)


# flower ---------------------------------------------------------------------
newdat_sens_flower <- expand.grid(
  region = c("Norway", "Switzerland"),
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_flower$temp_sens <- predict(
  m_sens_flower,
  newdata = newdat_sens_flower,
  re.form = NA 
)
newdat_sens_flower

mm_sens_flower <- model.matrix(terms(m_sens_flower), newdat_sens_flower)

pvar1_sens_flower <- diag(mm_sens_flower %*% tcrossprod(vcov(m_sens_flower), mm_sens_flower))

tvar1_sens_flower <- as.numeric(VarCorr(m_sens_flower)$species)
#var_block_sens_flower <- as.numeric(VarCorr(m_sens_flower)$block_ID)

#tvar1_sens_flower <- pvar1_sens_flower + var_species_sens_flower + var_block_sens_flower + sigma(m_sens_flower)^2

cmult <- 2

newdat_sens_flower <- data.frame(
  newdat_sens_flower,
  plo = newdat_sens_flower$temp_sens - cmult * sqrt(pvar1_sens_flower),
  phi = newdat_sens_flower$temp_sens + cmult * sqrt(pvar1_sens_flower),
  tlo = newdat_sens_flower$temp_sens - cmult * sqrt(tvar1_sens_flower),
  thi = newdat_sens_flower$temp_sens + cmult * sqrt(tvar1_sens_flower)
)
newdat_sens_flower

pd <- position_dodge(width = 0.4)

f <- ggplot(newdat_sens_flower, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  facet_wrap(~ region) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (flowering)",
    title = "Effect of transplantation on flowering onset across regions"
  )
print(f)

# fruit ---------------------------------------------------------------------
newdat_sens_fruit <- expand.grid(
  region = c("Norway", "Switzerland"),
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_fruit$temp_sens <- predict(
  m_sens_fruit,
  newdata = newdat_sens_fruit,
  re.form = NA 
)

mm_sens_fruit <- model.matrix(terms(m_sens_fruit), newdat_sens_fruit)

pvar1_sens_fruit <- diag(mm_sens_fruit %*% tcrossprod(vcov(m_sens_fruit), mm_sens_fruit))

tvar1_sens_fruit <- as.numeric(VarCorr(m_sens_fruit)$species)
#var_block_sens_fruit <- as.numeric(VarCorr(m_sens_fruit)$block_ID)

#tvar1_sens_fruit <- pvar1_sens_fruit + var_species_sens_fruit + var_block_sens_fruit + sigma(m_sens_fruit)^2

cmult <- 2

newdat_sens_fruit <- data.frame(
  newdat_sens_fruit,
  plo = newdat_sens_fruit$temp_sens - cmult * sqrt(pvar1_sens_fruit),
  phi = newdat_sens_fruit$temp_sens + cmult * sqrt(pvar1_sens_fruit),
  tlo = newdat_sens_fruit$temp_sens - cmult * sqrt(tvar1_sens_fruit),
  thi = newdat_sens_fruit$temp_sens + cmult * sqrt(tvar1_sens_fruit)
)

pd <- position_dodge(width = 0.4)

fr <- ggplot(newdat_sens_fruit, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  facet_wrap(~ region) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (fruiting)",
    title = "Effect of transplantation on fruiting onset across regions"
  )
print(fr)


# seed ---------------------------------------------------------------------
newdat_sens_seed <- expand.grid(
  region = c("Norway", "Switzerland"),
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_seed$temp_sens <- predict(
  m_sens_seed,
  newdata = newdat_sens_seed,
  re.form = NA 
)
newdat_sens_seed

mm_sens_seed <- model.matrix(terms(m_sens_seed), newdat_sens_seed)

pvar1_sens_seed <- diag(mm_sens_seed %*% tcrossprod(vcov(m_sens_seed), mm_sens_seed))

tvar1_sens_seed <- as.numeric(VarCorr(m_sens_seed)$species)
#var_block_sens_seed <- as.numeric(VarCorr(m_sens_seed)$block_ID)

#tvar1_sens_seed <- pvar1_sens_seed + var_species_sens_seed + var_block_sens_seed + sigma(m_sens_seed)^2

cmult <- 2

newdat_sens_seed <- data.frame(
  newdat_sens_seed,
  plo = newdat_sens_seed$temp_sens - cmult * sqrt(pvar1_sens_seed),
  phi = newdat_sens_seed$temp_sens + cmult * sqrt(pvar1_sens_seed),
  tlo = newdat_sens_seed$temp_sens - cmult * sqrt(tvar1_sens_seed),
  thi = newdat_sens_seed$temp_sens + cmult * sqrt(tvar1_sens_seed)
)
newdat_sens_seed


pd <- position_dodge(width = 0.4)

s <- ggplot(newdat_sens_seed, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  facet_wrap(~ region) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (seed set)",
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

plot_df_sens_bud  <- newdat_sens_bud   |> mutate(stage = "Budding")

plot_df_sens_flower <- newdat_sens_flower |> mutate(stage = "Flowering")

plot_df_sens_fruit <- newdat_sens_fruit |> mutate(stage = "Fruiting")

plot_df_sens_seed <- newdat_sens_seed |> mutate(stage = "Seeds")

plot_df_sens_all <- bind_rows(
  plot_df_sens_bud,
  plot_df_sens_flower,
  plot_df_sens_fruit,
  plot_df_sens_seed
)
plot_df_sens_all




# plot --------------------------------------------------------------------
pd <- position_dodge(width = 0.6) 

b_f_fr_s <- ggplot(plot_df_sens_all, aes(
  x = treat_competition,
  y = temp_sens,
  color = stage
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (days per °C)",
    title = "Effect of transplantation on onset across regions",
    color = "Phenological stage"
  )
print(b_f_fr_s)



# or plot facet by stage --------------------------------------------------
region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

b_f_fr_s2 <- ggplot(plot_df_sens_all, aes(
  x = treat_competition,
  y = temp_sens,
  color = region
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  #facet_wrap(~ stage) + 
  facet_grid(~ stage) +
  
  scale_color_manual(values = region_colors) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (days per °C)",
    title = "Effect of transplantation on onset across regions",
    color = "Region"
  )
print(b_f_fr_s2)


# ggsave(filename = "Output/Onset/Temperature_sensitivity_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#        plot = b_f_fr_s2,
#        width = 15, height = 10, units = "in")




