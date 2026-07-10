
# Onset julian days region -------------------------------------------------------


# Onset bud, flower, fruit model with site_treat_warming combination --------
# julian days

# load library ---------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(performance)
library(see)
library(sjPlot)


# source script with functions --------------------------------------------
source("Functions_onset.R")


# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology2)

theme_set(theme_bw())


# calculate first onset per species and plot for all stages ----------------
onset_bud    <- get_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_onset(phenology2, "No_Seeds", "jday")





# fit onset model ---------------------------------------------------------

m_onset_bud    <- fit_onset_model_region(onset_bud)
m_onset_flower <- fit_onset_model_region(onset_flower)
m_onset_fruit  <- fit_onset_model_region(onset_fruit)
m_onset_seed   <- fit_onset_model_region(onset_seed)

# check model output
# bud
summary(m_onset_bud)

anova(m_onset_bud)
anova_bud <- as.data.frame(anova(m_onset_bud))
sjPlot::tab_df(
  anova_bud,
  title = "Type III ANOVA Bud Onset")

library(gt)

anova_gt <- anova(m_onset_bud) |>
  as.data.frame() |>
  tibble::rownames_to_column("Effect") |>
  gt()
anova_gt
gtsave(anova_gt, "Output/anova_bud.png")
       

model_performance(m_onset_bud)
onset_region_bud_anova <- tab_model(m_onset_bud)
onset_region_bud_anova
sjPlot::save_tab(onset_region_bud_anova, file = "Output/m_onset_bud.png")


emmeans(m_onset_bud, pairwise ~ region)

emmeans(m_onset_bud,
  pairwise ~ region | treatment_site_temp * treat_competition)

plot_model(m_onset_bud, type = "est", show.values = TRUE, value.size = 4,
           vline.color = "red")

plot_model(m_onset_bud,
           type = "pred",
           terms = c("treatment_site_temp",
            "treat_competition",
            "region"),
           show.values = TRUE)



# flower
summary(m_onset_flower)
anova(m_onset_flower)
model_performance(m_onset_flower)

emmeans(m_onset_flower, pairwise ~ region)

emmeans(m_onset_flower,
        pairwise ~ region | treatment_site_temp * treat_competition)

plot_model(m_onset_flower, type = "est", show.values = TRUE, value.size = 4,
           vline.color = "red")

plot_model(m_onset_flower,
           type = "pred",
           terms = c("treatment_site_temp",
                     "treat_competition",
                     "region"))


















