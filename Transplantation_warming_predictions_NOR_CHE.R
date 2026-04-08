

# Onset bud, flower, fruit model with site_treat_warming combination --------


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

# one model per stage -----------------------------------------------------

# bud ---------------------------------------------------------------------
onset_bud <- phenology2 |>
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, treatment_site_temp, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))

m_bud <- lmerTest::lmer(onset ~ region * treatment_site_temp * treat_competition + 
                (1|species) + (1|block_ID),
              data = onset_bud)
summary(m_bud)

model_performance(m_bud)
check_collinearity(m_bud)
#check_model(m_bud)



# flower ---------------------------------------------------------------------
onset_flower <- phenology2 |>
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, treatment_site_temp, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))

m_flower <- lmerTest::lmer(onset ~ region * treatment_site_temp * treat_competition + 
                          (1|species) + (1|block_ID),
                        data = onset_flower)
summary(m_flower)

model_performance(m_flower)
check_collinearity(m_flower)
#check_model(m_flower)


# fruit ---------------------------------------------------------------------
# combine the fruiting stages nor and che to comapre the onset ------------
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology2 <- phenology2 |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))

onset_fruit <- phenology2 |>
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, treatment_site_temp, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))

m_fruit <- lmerTest::lmer(onset ~ region * treatment_site_temp * treat_competition + 
                             (1|species) + (1|block_ID),
                           data = onset_fruit)
summary(m_fruit)

model_performance(m_fruit)
check_collinearity(m_fruit)
#check_model(m_fruit)




# make predictions for each stage -----------------------------------------

# bud ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_bud <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_bud



# predict  ----------------------------------------------------------------
newdat_bud$onset <- predict(
  m_bud,
  newdata = newdat_bud,
  re.form = NA 
)
newdat_bud


# make model matrix -------------------------------------------------------
mm_bud <- model.matrix(terms(m_bud), newdat_bud)
mm_bud


pvar1_bud <- diag(mm_bud %*% tcrossprod(vcov(m_bud), mm_bud))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_bud <- as.numeric(VarCorr(m_bud)$species)
var_block_bud <- as.numeric(VarCorr(m_bud)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_bud <- pvar1_bud + var_species_bud + var_block_bud + sigma(m_bud)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_bud <- data.frame(
  newdat_bud
  , plo = newdat_bud$onset-cmult*sqrt(pvar1_bud) # fixed effects only, confidence interval
  , phi = newdat_bud$onset+cmult*sqrt(pvar1_bud)
  , tlo = newdat_bud$onset-cmult*sqrt(tvar1_bud) # takes fixed and random effects, prediction intervall
  , thi = newdat_bud$onset+cmult*sqrt(tvar1_bud)
)
newdat_bud


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

p<- ggplot(newdat_bud, aes(x=treat_competition, y=onset, 
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
    y = "Predicted onset (budding)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(p)



# flower ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_flower <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_flower



# predict  ----------------------------------------------------------------
newdat_flower$onset <- predict(
  m_flower,
  newdata = newdat_flower,
  re.form = NA 
)
newdat_flower


# make model matrix -------------------------------------------------------
mm_flower <- model.matrix(terms(m_flower), newdat_flower)
mm_flower


pvar1_flower <- diag(mm_flower %*% tcrossprod(vcov(m_flower), mm_flower))
# tvar1 <- pvar1+VarCorr(m_onset_flower_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_flower <- as.numeric(VarCorr(m_flower)$species)
var_block_flower <- as.numeric(VarCorr(m_flower)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_flower <- pvar1_flower + var_species_flower + var_block_flower + sigma(m_flower)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_flower <- data.frame(
  newdat_flower
  , plo = newdat_flower$onset-cmult*sqrt(pvar1_flower) # fixed effects only, confidence interval
  , phi = newdat_flower$onset+cmult*sqrt(pvar1_flower)
  , tlo = newdat_flower$onset-cmult*sqrt(tvar1_flower) # takes fixed and random effects, prediction intervall
  , thi = newdat_flower$onset+cmult*sqrt(tvar1_flower)
)
newdat_flower


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

f <- ggplot(newdat_flower, aes(x=treat_competition, y=onset, 
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
    y = "Predicted onset (budding)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(f)



# fruit ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_fruit <- expand.grid(
  region = c("Norway", "Switzerland"),
  treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_fruit



# predict  ----------------------------------------------------------------
newdat_fruit$onset <- predict(
  m_fruit,
  newdata = newdat_fruit,
  re.form = NA 
)
newdat_fruit


# make model matrix -------------------------------------------------------
mm_fruit <- model.matrix(terms(m_fruit), newdat_fruit)
mm_fruit


pvar1_fruit <- diag(mm_fruit %*% tcrossprod(vcov(m_fruit), mm_fruit))
# tvar1 <- pvar1+VarCorr(m_onset_fruit_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species_fruit <- as.numeric(VarCorr(m_fruit)$species)
var_block_fruit <- as.numeric(VarCorr(m_fruit)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_fruit <- pvar1_fruit + var_species_fruit + var_block_fruit + sigma(m_fruit)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_fruit <- data.frame(
  newdat_fruit
  , plo = newdat_fruit$onset-cmult*sqrt(pvar1_fruit) # fixed effects only, confidence interval
  , phi = newdat_fruit$onset+cmult*sqrt(pvar1_fruit)
  , tlo = newdat_fruit$onset-cmult*sqrt(tvar1_fruit) # takes fixed and random effects, prediction intervall
  , thi = newdat_fruit$onset+cmult*sqrt(tvar1_fruit)
)
newdat_fruit


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

fr <- ggplot(newdat_fruit, aes(x=treat_competition, y=onset, 
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
    y = "Predicted onset (budding)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(fr)





# combine all stages into one plot ----------------------------------------

stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636"   
)

plot_df_bud  <- newdat_bud   |> mutate(stage = "Budding")

plot_df_flower <- newdat_flower |> mutate(stage = "Flowering")

plot_df_fruit <- newdat_fruit |> mutate(stage = "Fruiting")


plot_df_all <- bind_rows(
  plot_df_bud,
  plot_df_flower,
  plot_df_fruit
)
plot_df_all





# plot with filled and unfilled shapes ------------------------------------
plot_df_all <- plot_df_all |>
  mutate(
    shape_code = case_when(
      treatment_site_temp == "lo_ambi" & treat_competition == "with"    ~ 16,
      treatment_site_temp == "lo_ambi" & treat_competition == "without" ~ 1,
      
      treatment_site_temp == "hi_ambi" & treat_competition == "with"    ~ 17,
      treatment_site_temp == "hi_ambi" & treat_competition == "without" ~ 2,
      
      treatment_site_temp == "hi_warm" & treat_competition == "with"    ~ 15,
      treatment_site_temp == "hi_warm" & treat_competition == "without" ~ 0
    )
  )
plot_df_all


# define order of dots ----------------------------------------------------
plot_df_all <- plot_df_all |>
  mutate(shape_code = factor(shape_code,
                             levels = c(1,16, 2,17, 0,15)))



# plot --------------------------------------------------------------------
pd <- position_dodge(width = 0.6) 

b_f_fr <- ggplot(plot_df_all, aes(
  x = treat_competition,
  y = onset,
  color = stage,
  shape = factor(shape_code)
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  
  scale_shape_manual(
    values = c("0" = 0, "1" = 1, "2" = 2, "16" = 16, "17" = 17, "15" = 15),
    labels = c(
      "1"  = "lo_ambi, without",
      "16" = "lo_ambi, with",
      "2"  = "hi_ambi, without",
      "17" = "hi_ambi, with",
      "0"  = "hi_warm, without",
      "15" = "hi_warm, with"
    )
  ) +
  #guides(shape = "none")+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Site x treat warming x biotic interactions",
    color = "Phenological stage"
  )

print(b_f_fr)


# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_bud_flower_fruit_predictions.png", 
#        plot = b_f_fr, width = 18, height = 10, units = "in")


#



#######
plot_df_all2 <- plot_df_all |>
  mutate(
    shape_code = case_when(
      treatment_site_temp == "lo_ambi" & treat_competition == "with"    ~ 16,
      treatment_site_temp == "lo_ambi" & treat_competition == "without" ~ 16,
      
      treatment_site_temp == "hi_ambi" & treat_competition == "with"    ~ 2,
      treatment_site_temp == "hi_ambi" & treat_competition == "without" ~ 2,
      
      treatment_site_temp == "hi_warm" & treat_competition == "with"    ~ 17,
      treatment_site_temp == "hi_warm" & treat_competition == "without" ~ 17
    )
  )
plot_df_all2


# define order of dots ----------------------------------------------------
plot_df_all2 <- plot_df_all2 |>
  mutate(shape_code = factor(shape_code,
                             levels = c(16, 2, 17)))



# plot --------------------------------------------------------------------
pd <- position_dodge(width = 0.6) 

b_f_fr2 <- ggplot(plot_df_all2, aes(
  x = treat_competition,
  y = onset,
  color = stage,
  shape = factor(shape_code)
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  
  scale_shape_manual(
    values = c("2" = 2, "16" = 16, "17" = 17),
    labels = c(
      "16"  = "lo_ambi, without",
      "16" = "lo_ambi, with",
      "2"  = "hi_ambi, without",
      "2" = "hi_ambi, with",
      "17"  = "hi_warm, without",
      "17" = "hi_warm, with"
    )
  ) +
  #guides(shape = "none")+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Site x treat warming",
    color = "Phenological stage"
  )

print(b_f_fr2)




b_f_fr2 <- ggplot(plot_df_all2, aes(
  x = treat_competition,
  y = onset,
  color = stage,
  shape = treatment_site_temp 
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 2,    # triangle
      "hi_warm" = 17    # square or triangle
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Treatment site × warming",   # ← clean legend title
    color = "Phenological stage"
  )

print(b_f_fr2)


# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_bud_flower_fruit_predictions.png", 
#        plot = b_f_fr2, width = 18, height = 10, units = "in")




# or plot facet by stage --------------------------------------------------
region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

b_f_fr3 <- ggplot(plot_df_all2, aes(
  x = treat_competition,
  y = onset,
  color = region,              # ← color by region
  shape = treatment_site_temp  # ← shape by treatment
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ stage) +        # ← facet by stage
  
  scale_color_manual(values = region_colors) +   # ← correct palette
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 2,    # triangle
      "hi_warm" = 17    # square or whatever you prefer
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Treatment site × warming",
    color = "Region"
  )

print(b_f_fr3)

# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_bud_flower_fruit_predictions2.png", 
#        plot = b_f_fr3, width = 18, height = 10, units = "in")




