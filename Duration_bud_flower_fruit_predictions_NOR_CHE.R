

# 04 duration -------------------------------------------------------------

# Duration bud, flower, fruit model with site_treat_warming combination --------


# load library ---------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(performance)
library(see)




# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))


# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology2 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional scaling if needed




# add site_warming treatment ----------------------------------------------
phenology2$treatment_site_temp <- paste(phenology$site, phenology$treat_warming, sep = "_")


phenology2 <- phenology2 |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_ambi",
                                                "hi_warm")))





# onsets bud, flower, fruit, seed -----------------------------------------

onsets <- phenology2 |>
  filter(value > 0) |>
  group_by(region, treatment_site_temp, treat_competition,
           species, block_ID, unique_plot_ID, unique_plant_ID,
           phenology_stage) |>
  summarise(
    onset = min(jday),
    .groups = "drop"
  )


onsets_wide <- onsets |>
  pivot_wider(
    names_from = phenology_stage,
    values_from = onset
  )
onsets_wide



# duration of budding, flowering, fruiting --------------------------------
durations <- onsets_wide |>
  mutate(
    duration_bud = No_FloOpen - No_Buds,
    duration_flower = No_FloWithrd - No_FloOpen,
    duration_fruit = No_Seeds - No_FloWithrd,
  )
durations



# model for each stage ----------------------------------------------------

# model budding duration ------------------------------------------
m_duration_bud <- lmerTest::lmer(duration_bud ~ region * treatment_site_temp * treat_competition 
                                 + (1|species) + (1|block_ID),
                                 data = durations)
summary(m_duration_bud)
anova(m_duration_bud)

model_performance(m_duration_bud)


emm_bud <- emmeans(m_duration_bud, ~ treatment_site_temp * treat_competition)
pairs(emm_bud)
contrast(emm_bud, interaction = "pairwise")



# model flowering duration ------------------------------------------
m_duration_flowering <- lmerTest::lmer(duration_flower ~ region * treatment_site_temp * treat_competition 
                                 + (1|species) + (1|block_ID),
                                 data = durations)
summary(m_duration_flowering)
anova(m_duration_flowering)

model_performance(m_duration_flowering)

emm_flower <- emmeans(m_duration_flowering, ~ treatment_site_temp * treat_competition)
pairs(emm_flower)
contrast(emm_flower, interaction = "pairwise")


# model fruiting duration ------------------------------------------
m_duration_fruiting<- lmerTest::lmer(duration_fruit ~ region * treatment_site_temp * treat_competition 
                                       + (1|species) + (1|block_ID),
                                       data = durations)
summary(m_duration_fruiting)

model_performance(m_duration_fruiting)

emm_fruit <- emmeans(m_duration_fruiting, ~ treatment_site_temp * treat_competition)
pairs(emm_fruit)
contrast(emm_fruit, interaction = "pairwise")



# predictions for each stage ----------------------------------------------

# bud ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_duration_bud <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  duration_bud = 0
)
newdat_duration_bud



# predict  ----------------------------------------------------------------
newdat_duration_bud$duration_bud <- predict(
  m_duration_bud,
  newdata = newdat_duration_bud,
  re.form = NA 
)
newdat_duration_bud


# make model matrix -------------------------------------------------------
mm_duration_bud <- model.matrix(terms(m_duration_bud), newdat_duration_bud)
mm_duration_bud


pvar1_duration_bud <- diag(mm_duration_bud %*% tcrossprod(vcov(m_duration_bud), mm_duration_bud))
# tvar1 <- pvar1+VarCorr(m_onset_duration_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_duration_bud <- as.numeric(VarCorr(m_duration_bud)$species)
var_block_duration_bud <- as.numeric(VarCorr(m_duration_bud)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_duration_bud <- pvar1_duration_bud + var_species_duration_bud + var_block_duration_bud + sigma(m_duration_bud)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_duration_bud <- data.frame(
  newdat_duration_bud
  , plo = newdat_duration_bud$duration_bud-cmult*sqrt(pvar1_duration_bud) # fixed effects only, confidence interval
  , phi = newdat_duration_bud$duration_bud+cmult*sqrt(pvar1_duration_bud)
  , tlo = newdat_duration_bud$duration_bud-cmult*sqrt(tvar1_duration_bud) # takes fixed and random effects, prediction intervall
  , thi = newdat_duration_bud$duration_bud+cmult*sqrt(tvar1_duration_bud)
)
newdat_duration_bud


# plot predicted duration ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

db <- ggplot(newdat_duration_bud, aes(x=treat_competition, y=duration_bud, 
                           color=treat_competition, shape = treatment_site_temp )) +
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
    y = "Predicted duration (jdays)",
    title = "Effect of transplantation on budding duration across regions"
  )
print(db)


# flowering duration ------------------------------------------------------
# create new matrix for predicted data ------------------------------------
newdat_duration_flower <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  duration_flower = 0
)

# predict -----------------------------------------------------------------
newdat_duration_flower$duration_flower <- predict(
  m_duration_flowering,
  newdata = newdat_duration_flower,
  re.form = NA 
)

# model matrix ------------------------------------------------------------
mm_duration_flower <- model.matrix(terms(m_duration_flowering), newdat_duration_flower)

# variance (fixed effects) ------------------------------------------------
pvar1_duration_flower <- diag(mm_duration_flower %*% tcrossprod(vcov(m_duration_flowering), mm_duration_flower))

# random effects variance -------------------------------------------------
var_species_duration_flower <- as.numeric(VarCorr(m_duration_flowering)$species)
var_block_duration_flower   <- as.numeric(VarCorr(m_duration_flowering)$block_ID)

# total variance (prediction interval) ------------------------------------
tvar1_duration_flower <- pvar1_duration_flower + 
  var_species_duration_flower + 
  var_block_duration_flower + 
  sigma(m_duration_flowering)^2

# intervals ---------------------------------------------------------------
cmult <- 2

newdat_duration_flower <- data.frame(
  newdat_duration_flower,
  plo = newdat_duration_flower$duration_flower - cmult * sqrt(pvar1_duration_flower),
  phi = newdat_duration_flower$duration_flower + cmult * sqrt(pvar1_duration_flower),
  tlo = newdat_duration_flower$duration_flower - cmult * sqrt(tvar1_duration_flower),
  thi = newdat_duration_flower$duration_flower + cmult * sqrt(tvar1_duration_flower)
)

ggplot(newdat_duration_flower, aes(x=treat_competition, y=duration_flower, 
                                color=treat_competition, shape = treatment_site_temp )) +
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
    y = "Predicted duration flowering (jdays)",
    title = "Effect of transplantation on flowering duration across regions"
  )

# fruiting ----------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_duration_fruit <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  duration_fruit = 0
)

# predict -----------------------------------------------------------------
newdat_duration_fruit$duration_fruit <- predict(
  m_duration_fruiting,
  newdata = newdat_duration_fruit,
  re.form = NA 
)

# model matrix ------------------------------------------------------------
mm_duration_fruit <- model.matrix(terms(m_duration_fruiting), newdat_duration_fruit)

# variance (fixed effects) ------------------------------------------------
pvar1_duration_fruit <- diag(mm_duration_fruit %*% tcrossprod(vcov(m_duration_fruiting), mm_duration_fruit))

# random effects variance -------------------------------------------------
var_species_duration_fruit <- as.numeric(VarCorr(m_duration_fruiting)$species)
var_block_duration_fruit   <- as.numeric(VarCorr(m_duration_fruiting)$block_ID)

# total variance (prediction interval) ------------------------------------
tvar1_duration_fruit <- pvar1_duration_fruit + 
  var_species_duration_fruit + 
  var_block_duration_fruit + 
  sigma(m_duration_fruiting)^2

# intervals ---------------------------------------------------------------
cmult <- 2

newdat_duration_fruit <- data.frame(
  newdat_duration_fruit,
  plo = newdat_duration_fruit$duration_fruit - cmult * sqrt(pvar1_duration_fruit),
  phi = newdat_duration_fruit$duration_fruit + cmult * sqrt(pvar1_duration_fruit),
  tlo = newdat_duration_fruit$duration_fruit - cmult * sqrt(tvar1_duration_fruit),
  thi = newdat_duration_fruit$duration_fruit + cmult * sqrt(tvar1_duration_fruit)
)

ggplot(newdat_duration_fruit, aes(x=treat_competition, y=duration_fruit, 
                                color=treat_competition, shape = treatment_site_temp )) +
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
    y = "Predicted duration fruiting (jdays)",
    title = "Effect of transplantation on fruiting duration across regions"
  )




# combined plot -----------------------------------------------------------

newdat_duration_bud$stage <- "Budding"
newdat_duration_flower$stage <- "Flowering"
newdat_duration_fruit$stage <- "Fruiting"


pred_bud <- newdat_duration_bud |>
  rename(duration = duration_bud)

pred_flower <- newdat_duration_flower |>
  rename(duration = duration_flower)

pred_fruit <- newdat_duration_fruit |>
  rename(duration = duration_fruit)




pred_all <- bind_rows(pred_bud, pred_flower, pred_fruit)
pred_all

pred_all$stage <- factor(pred_all$stage,
                         levels = c("Budding", "Flowering", "Fruiting"))

pd <- position_dodge(width = 0.4)


region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

duration <- ggplot(pred_all,
       aes(x = treat_competition,
           y = duration,
           color = region,
           shape = treatment_site_temp)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = 0.2,
                position = pd) +
  facet_grid( ~ stage) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
    )
  ) +
  # scale_color_manual(values = c(
  #   "with" = "#528B8B",
  #   "without" = "#CD950C"
  # )) +
  labs(
    x = "Biotic interactions",
    y = "Predicted duration (days)",
    title = "Effect of transplantation on phenological durations across stages"
  )
duration

# ggsave(filename = "Output/Onset/Transplantation_Warming_duration_bud_flower_fruit_seeds_predictions2.png", 
#        plot = duration, width = 15, height = 12, units = "in")






