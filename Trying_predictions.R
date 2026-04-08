

# Predictions for onset ---------------------------------------------------
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions


theme_set(theme_bw(base_size = 20))

# Transplantation - budding -----------------------------------------------

source("Effect_cooling_on_budding_onset_NOR_CHE.R")


onset_bud_cool <- phenology_cool |> 
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, site, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))


m_onset_bud_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                        (1|species) + (1|block_ID), 
                                      data = onset_bud_cool)

summary(m_onset_bud_cooling)



# create new matrix for predicted data ------------------------------------
newdat <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("hi", "lo"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat



# predict  ----------------------------------------------------------------
newdat$onset <- predict(
  m_onset_bud_cooling,
  newdata = newdat,
  re.form = NA 
)
newdat


# make model matrix -------------------------------------------------------
mm <- model.matrix(terms(m_onset_bud_cooling),newdat)
mm


pvar1 <- diag(mm %*% tcrossprod(vcov(m_onset_bud_cooling),mm))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_species <- as.numeric(VarCorr(m_onset_bud_cooling)$species)
var_block <- as.numeric(VarCorr(m_onset_bud_cooling)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1 <- pvar1 + var_species + var_block + sigma(m_onset_bud_cooling)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat <- data.frame(
  newdat
  , plo = newdat$onset-cmult*sqrt(pvar1) # fixed effects only, confidence interval
  , phi = newdat$onset+cmult*sqrt(pvar1)
  , tlo = newdat$onset-cmult*sqrt(tvar1) # takes fixed and random effects, prediction intervall
  , thi = newdat$onset+cmult*sqrt(tvar1)
)
newdat



# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

p<- ggplot(newdat, aes(x=treat_competition, y=onset, 
                       color=treat_competition, shape = site)) +
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

# so we use the CI with fixed effects only




# Transplantation - flowering -----------------------------------------------

source("Effect_cooling_on_flowering_onset_NOR_CHE.R")

# calculate flowering onset ------------------------------------------------
flowering_onset_n_c_cool <- phenology_cool |> 
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))


m_onset_n_c_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                        (1|species) + (1|block_ID), 
                                      data = flowering_onset_n_c_cool)
summary(m_onset_n_c_cooling)



# create new matrix for predicted data ------------------------------------
newdat_cool_flower <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("hi", "lo"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_cool_flower



# predict  ----------------------------------------------------------------
newdat_cool_flower$onset <- predict(
  m_onset_n_c_cooling,
  newdata = newdat_cool_flower,
  re.form = NA 
)
newdat_cool_flower


# make model matrix -------------------------------------------------------
mm_c_f <- model.matrix(terms(m_onset_n_c_cooling), newdat_cool_flower)
mm_c_f


pvar1_c_f <- diag(mm_c_f %*% tcrossprod(vcov(m_onset_n_c_cooling),mm_c_f))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_c_f_species <- as.numeric(VarCorr(m_onset_n_c_cooling)$species)
var_c_f_block <- as.numeric(VarCorr(m_onset_n_c_cooling)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_c_f <- pvar1_c_f + var_c_f_species + var_c_f_block + sigma(m_onset_n_c_cooling)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_cool_flower <- data.frame(
  newdat_cool_flower
  , plo = newdat_cool_flower$onset-cmult*sqrt(pvar1) # fixed effects only, confidence interval
  , phi = newdat_cool_flower$onset+cmult*sqrt(pvar1)
  , tlo = newdat_cool_flower$onset-cmult*sqrt(tvar1) # takes fixed and random effects, prediction intervall
  , thi = newdat_cool_flower$onset+cmult*sqrt(tvar1)
)
newdat_cool_flower



# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

f<- ggplot(newdat_cool_flower, aes(x=treat_competition, y=onset, 
                       color=treat_competition, shape = site)) +
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
    y = "Predicted onset (flowering)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(f)



# Transplantation - fruiting -----------------------------------------------

source("Effect_cooling_on_fruiting_onset_NOR_CHE.R")

# calculate fruiting onset ------------------------------------------------
onset_fruit_cool <- phenology_cool |> 
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, site, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where fruits never occurred
  filter(is.finite(onset))

# model with region for fruiting onset lmer ----------------------------------
m_onset_fruit_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                          (1|species) + (1|block_ID), 
                                        data = onset_fruit_cool)

summary(m_onset_fruit_cooling)



# create new matrix for predicted data ------------------------------------
newdat_cool_fruit <- expand.grid(
  region = c("Norway", "Switzerland"),
  site = c("hi", "lo"),
  treat_competition = c("with", "without"),
  onset = 0
)
newdat_cool_fruit



# predict  ----------------------------------------------------------------
newdat_cool_fruit$onset <- predict(
  m_onset_fruit_cooling,
  newdata = newdat_cool_fruit,
  re.form = NA 
)
newdat_cool_fruit


# make model matrix -------------------------------------------------------
mm_c_fr <- model.matrix(terms(m_onset_fruit_cooling), newdat_cool_fruit)
mm_c_fr


pvar1_c_fr <- diag(mm_c_fr %*% tcrossprod(vcov(m_onset_fruit_cooling),mm_c_fr))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
var_c_fr_species <- as.numeric(VarCorr(m_onset_fruit_cooling)$species)
var_c_fr_block <- as.numeric(VarCorr(m_onset_fruit_cooling)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
tvar1_c_fr <- pvar1_c_fr + var_c_fr_species + var_c_fr_block + sigma(m_onset_fruit_cooling)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_cool_fruit <- data.frame(
  newdat_cool_fruit
  , plo = newdat_cool_fruit$onset-cmult*sqrt(pvar1) # fixed effects only, confidence interval
  , phi = newdat_cool_fruit$onset+cmult*sqrt(pvar1)
  , tlo = newdat_cool_fruit$onset-cmult*sqrt(tvar1) # takes fixed and random effects, prediction intervall
  , thi = newdat_cool_fruit$onset+cmult*sqrt(tvar1)
)
newdat_cool_fruit



# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.5) 

fr <- ggplot(newdat_cool_fruit, aes(x=treat_competition, y=onset, 
                                   color=treat_competition, shape = site)) +
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
    y = "Predicted onset (fruiting)",
    title = "Effect of transplantation on budding onset across regions"
  )
print(fr)






# combine stages into one figure ------------------------------------------

stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636"   
)

plot_df_bud  <- newdat   |> mutate(stage = "Budding")

plot_df_flower <- newdat_cool_flower |> mutate(stage = "Flowering")

plot_df_fruit <- newdat_cool_fruit |> mutate(stage = "Fruiting")


plot_df_all <- bind_rows(
  plot_df_bud,
  plot_df_flower,
  plot_df_fruit
)
plot_df_all


b_f_fr <- ggplot(plot_df_all, aes(x=treat_competition, y=onset, 
                                    color = stage, shape = site)) +
  geom_point(position = pd, size = 4, stroke = 1.2)+
  facet_wrap(~ region)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  #geom_line(aes(group = interaction(stage, site)), position = pd)+
  scale_color_manual(values = stage_colors)+
  scale_shape_manual(values = c("hi" = 24, "lo" = 21))+
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation on onset across regions"
  )#+
  #coord_flip()
print(b_f_fr)




plot_df_all <- plot_df_all |>
  mutate(
    shape_code = case_when(
      site == "hi" & treat_competition == "with" ~ 17,
      site == "hi" & treat_competition == "without" ~ 2,
      site == "lo" & treat_competition == "with" ~ 16,
      site == "lo" & treat_competition == "without" ~ 1
    )
  )
plot_df_all

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
    values = c("1" = 1, "2" = 2, "16" = 16, "17" = 17),
    labels = c(
      "1" = "lo site, without",
      "16" = "lo site, with",
      "2" = "hi site, without",
      "17" = "hi site, with"
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation on onset across regions",
    shape = "Site × biotic interactions",
    color = "Phenological stage")
print(b_f_fr)


# ggsave(filename = "Output/Onset/Transplantation_onset_bud_flower_fruit_predictions.png", 
#        plot = b_f_fr, width = 13, height = 10, units = "in")





# trying to get significance stars -------------------------------------------------

library(emmeans)

emm_bud <- emmeans(m_onset_bud_cooling, ~ treat_competition | region * site)

contrast_comp_bud <- contrast(emm_bud, method = "pairwise")
summary(contrast_comp_bud)


emm_flower <- emmeans(m_onset_n_c_cooling, ~ treat_competition | region * site)

contrast_comp_flower <- contrast(emm_flower, method = "pairwise")
summary(contrast_comp_flower)


emm_fruit <- emmeans(m_onset_fruit_cooling, ~ treat_competition | region * site)

contrast_comp_fruit <- contrast(emm_fruit, method = "pairwise")
summary(contrast_comp_fruit)




clean_contrasts <- function(contrast_df, stage_name) {
  contrast_df |>
    as.data.frame() |>
    mutate(
      stage = stage_name,
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) |>
    # emmeans names the comparison "with - without"
    separate(contrast, into = c("with", "without"), sep = " - ") |>
    select(region, site, stage, stars)
}


stars_bud    <- clean_contrasts(contrast_comp_bud,    "Budding")
stars_flower <- clean_contrasts(contrast_comp_flower, "Flowering")
stars_fruit  <- clean_contrasts(contrast_comp_fruit,  "Fruiting")

stars_all <- bind_rows(stars_bud, stars_flower, stars_fruit)
stars_all


plot_df_all <- plot_df_all |>
  left_join(stars_all, by = c("region", "site", "stage"))
plot_df_all

stars_for_plot <- plot_df_all |>
  filter(treat_competition == "with")
stars_for_plot



ggplot(plot_df_all, aes(
  x = treat_competition,
  y = onset,
  color = stage,
  shape = factor(shape_code)
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  # geom_line(aes(group = interaction(stage, site)),
  #           position = pd,
  #           alpha = 0.6)+
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  
  scale_shape_manual(
    values = c("1" = 1, "2" = 2, "16" = 16, "17" = 17),
    labels = c(
      "1" = "lo site, without",
      "16" = "lo site, with",
      "2" = "hi site, without",
      "17" = "hi site, with"
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation on onset across regions",
    shape = "Site × biotic interactions",
    color = "Phenological stage")+
  guides(shape = "none")+
  geom_text(
  data = stars_for_plot,
  aes(
    x = treat_competition,
    y = onset + 1,   # adjust upward offset
    label = stars
  ),
  color = "black",
  size = 6,
  vjust = 0,
  position = pd
)


##
plot_df_all2 <- plot_df_all |>
  mutate(
    x_num = as.numeric(factor(treat_competition)),
    dodge_id = as.numeric(factor(stage))  # dodge by stage
  )
plot_df_all2

dodge_width <- 0.2

ggplot(plot_df_all2, aes(
  y = onset,
  color = stage,
  shape = factor(shape_code)
)) +
  
  # lines (manually dodged)
  geom_line(aes(
    x = x_num + (dodge_id - 2) * dodge_width,
    group = interaction(stage, site)
  ),
  linewidth = 0.8,
  alpha = 0.7) +
  
  # points
  geom_point(aes(
    x = x_num + (dodge_id - 2) * dodge_width
  ),
  size = 4, stroke = 1.2) +
  
  # errorbars
  geom_errorbar(aes(
    x = x_num + (dodge_id - 2) * dodge_width,
    ymin = plo, ymax = phi
  ),
  width = .1) +
  
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("without", "with")
  ) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = stage_colors) +
  
  scale_shape_manual(
    values = c("1" = 1, "2" = 2, "16" = 16, "17" = 17)
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)"
  )









#####




