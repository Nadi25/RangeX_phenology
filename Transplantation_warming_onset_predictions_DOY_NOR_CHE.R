

# 01 Onset julian days -------------------------------------------------------


# Onset bud, flower, fruit model with site_treat_warming combination --------
# julian days

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
names(phenology2)



# function to get first onset per individual ----------------------------------------------
get_onset <- function(data, stage_name, onset_type) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, treatment_site_temp, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID
    ) |>
    summarise(
      onset = min(.data[[onset_type]]),
      .groups = "drop"
    ) |> 
    # remove groups where budding never occurred
    filter(is.finite(onset))
}


# calculate first onset per species and plot for all stages ----------------
onset_bud    <- get_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_onset(phenology2, "No_Seeds", "jday")



# function of model per stage and region ---------------------------------------------
fit_onset_model <- function(onset_data, region_name) {
  region <- onset_data |> 
    filter(region == region_name)
  
  model <- lmerTest::lmer(
    onset ~ treatment_site_temp + treat_competition + (1 | species) + (1 | block_ID),
    data = region
  )
  
  return(model)
}



# fit onset model ---------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_onset_bud_nor    <- fit_onset_model(onset_bud, "Norway")
m_onset_flower_nor <- fit_onset_model(onset_flower, "Norway")
m_onset_fruit_nor  <- fit_onset_model(onset_fruit, "Norway")
m_onset_seed_nor   <- fit_onset_model(onset_seed, "Norway")

# check model output
# bud
summary(m_onset_bud_nor)
anova(m_onset_bud_nor)
model_performance(m_onset_bud_nor)
#check_model(m_onset_bud_nor)

emmeans(m_onset_bud_nor,
        pairwise ~ treatment_site_temp * treat_competition)


# flower
summary(m_onset_flower_nor)
anova(m_onset_flower_nor)
model_performance(m_onset_flower_nor)
#check_model(m_onset_flower_nor)

emmeans(m_onset_flower_nor,
        pairwise ~ treatment_site_temp * treat_competition)

# fruit
summary(m_onset_fruit_nor)
anova(m_onset_fruit_nor)
model_performance(m_onset_fruit_nor)
#check_model(m_onset_fruit_nor)

emmeans(m_onset_fruit_nor,
        pairwise ~ treatment_site_temp * treat_competition)

# seed
summary(m_onset_seed_nor)
anova(m_onset_seed_nor)
model_performance(m_onset_seed_nor)
#check_model(m_onset_seed_nor)

emmeans(m_onset_seed_nor,
        pairwise ~ treatment_site_temp * treat_competition)




# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_onset_bud_che    <- fit_onset_model(onset_bud, "Switzerland")
m_onset_flower_che <- fit_onset_model(onset_flower, "Switzerland")
m_onset_fruit_che  <- fit_onset_model(onset_fruit, "Switzerland")
m_onset_seed_che   <- fit_onset_model(onset_seed, "Switzerland")

# check model output
# bud
summary(m_onset_bud_che)
anova(m_onset_bud_che)
model_performance(m_onset_bud_che)
#check_model(m_onset_bud_che)

emmeans(m_onset_bud_che,
        pairwise ~ treatment_site_temp * treat_competition)


# flower
summary(m_onset_flower_che)
anova(m_onset_flower_che)
model_performance(m_onset_flower_che)
#check_model(m_onset_flower_che)

emmeans(m_onset_flower_che,
        pairwise ~ treatment_site_temp * treat_competition)

# fruit
summary(m_onset_fruit_che)
anova(m_onset_fruit_che)
model_performance(m_onset_fruit_che)
#check_model(m_onset_fruit_che)

emmeans(m_onset_fruit_che,
        pairwise ~ treatment_site_temp * treat_competition)

# seed
summary(m_onset_seed_che)
anova(m_onset_seed_che)
model_performance(m_onset_seed_che)
#check_model(m_onset_seed_che)

emmeans(m_onset_seed_che,
        pairwise ~ treatment_site_temp * treat_competition)


# function for predictions ------------------------------------------------


make_onset_predictions <- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without")
  ) |>
    as_tibble()
  
  pred <- predict(
    model,
    newdata = newdata,
    re.form = NA,
    se.fit = TRUE) |> 
    
    as_tibble() |>
    mutate(upper = fit + 1.96 * se.fit,
           lower = fit - 1.96 * se.fit) |>
    bind_cols(newdata)
  
  return(pred)
}


# NOR ---------------------------------------------------------------------
# 
pred_onset_bud_nor    <- make_onset_predictions(m_onset_bud_nor)
pred_onset_flower_nor <- make_onset_predictions(m_onset_flower_nor)
pred_onset_fruit_nor  <- make_onset_predictions(m_onset_fruit_nor)
pred_onset_seed_nor   <- make_onset_predictions(m_onset_seed_nor)


# CHE ---------------------------------------------------------------------
# 
pred_onset_bud_che    <- make_onset_predictions(m_onset_bud_che)
pred_onset_flower_che <- make_onset_predictions(m_onset_flower_che)
pred_onset_fruit_che  <- make_onset_predictions(m_onset_fruit_che)
pred_onset_seed_che   <- make_onset_predictions(m_onset_seed_che)




# combine predictions into one dataframe ----------------------------------

# nor
plot_df_bud_nor  <- pred_onset_bud_nor   |> 
  mutate(stage = "Budding",
         region = "Norway")

plot_df_flower_nor  <- pred_onset_flower_nor   |> 
  mutate(stage = "Flowering",
         region = "Norway")

plot_df_fruit_nor  <- pred_onset_fruit_nor   |> 
  mutate(stage = "Fruiting",
         region = "Norway")

plot_df_seed_nor  <- pred_onset_seed_nor   |> 
  mutate(stage = "Seeds",
         region = "Norway")


# che
plot_df_bud_che  <- pred_onset_bud_che   |> 
  mutate(stage = "Budding",
         region = "Switzerland")

plot_df_flower_che  <- pred_onset_flower_che   |> 
  mutate(stage = "Flowering",
         region = "Switzerland")

plot_df_fruit_che  <- pred_onset_fruit_che   |> 
  mutate(stage = "Fruiting",
         region = "Switzerland")

plot_df_seed_che  <- pred_onset_seed_che   |> 
  mutate(stage = "Seeds",
         region = "Switzerland")


plot_df_all <- bind_rows(
  plot_df_bud_nor,
  plot_df_flower_nor,
  plot_df_fruit_nor,
  plot_df_seed_nor,
  plot_df_bud_che,
  plot_df_flower_che,
  plot_df_fruit_che,
  plot_df_seed_che
)
plot_df_all



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




# combine predictions into one dataframe ----------------------------------


raw_bud  <- onset_bud   |> 
  mutate(stage = "Budding")

raw_flower  <- onset_flower  |> 
  mutate(stage = "Flowering")

raw_fruit  <- onset_fruit   |> 
  mutate(stage = "Fruiting")

raw_seed  <- onset_seed   |> 
  mutate(stage = "Seeds")


plot_df_raw_all <- bind_rows(
  raw_bud,
  raw_flower,
  raw_fruit,
  raw_seed
)
plot_df_raw_all





# or plot facet by stage --------------------------------------------------
region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

b_f_fr4 <- ggplot(plot_df_all2, aes(
  x = treat_competition,
  y = fit,
  color = region,              
  shape = treatment_site_temp  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = .2,
                position = pd) +
  
  #facet_wrap(~ stage) + 
  facet_grid(~ stage) +
  
  scale_color_manual(values = region_colors) +   
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Treatment site × warming",
    color = "Region"
  )

print(b_f_fr4)



b_f_fr5 <- ggplot(plot_df_all2, aes(
  x = treat_competition,
  y = fit,
  color = treat_competition,              
  shape = treatment_site_temp  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = .2,
                position = pd) +
  
  #facet_wrap(~ stage) + 
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +   
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
    )
  ) +
  
  labs(
    x = "Biotic interactions",
    y = "Predicted onset (julian days)",
    title = "Effect of transplantation and warming on onset across regions",
    shape = "Treatment site × warming",
    color = "Biotic interactions"
  )+
  geom_violin(
    data = plot_df_raw_all,
    aes(
      x = treat_competition,
      y = onset,
      fill = treat_competition
    ),
    alpha = 0.25,
    position = pd,
    color = NA,
    trim = TRUE 
  )+
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))
print(b_f_fr5)



b_f_fr6 <- ggplot(plot_df_all2, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,              
  shape = treatment_site_temp  
)) +
  geom_violin(
    data = plot_df_raw_all,
    aes(
      x = treatment_site_temp,
      y = onset,
      fill = treat_competition
    ),
    alpha = 0.25,
    position = pd,
    color = NA,
    trim = TRUE)+
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"))+
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = .2,
                position = pd) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +   
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions"
  )+
  
  guides(fill = "none")
print(b_f_fr6)


b_f_fr6 <- ggplot(plot_df_all2, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # raw:
  geom_violin(
    data = plot_df_raw_all,
    aes(
      x = treatment_site_temp,
      y = onset,
      fill = treat_competition
      #group = interaction(treatment_site_temp, treat_competition)
    ),
    position = pd,
    alpha = 0.25,
    color = NA,
    trim = FALSE,
    adjust = 1.5 # makes it smoother so it looks less bubbly
  ) +
  
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  guides(shape = "none")
b_f_fr6


# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_bud_flower_fruit_seeds_predictions_violin3.png", 
#        plot = b_f_fr6, width = 18, height = 10, units = "in")



