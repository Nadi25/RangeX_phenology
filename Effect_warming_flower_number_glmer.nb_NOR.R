

# Effect of warming through transplantation on flower number NOR glmer.nb ---------------------------------

## Data used: 
## Date:      11.03.2026
## Author:    Nadine Arzt
## Purpose:   Effect of warming on flower number NOR with glmer.nb



# load library ------------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))


# source the phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR.R")


# fit the model glmer.nb  ------------------------------------
m_flower_number_warm_nb <- glmer.nb(
  value ~ treat_warming * treat_competition +
    (1|species) + (1|block_ID),
  data = phenology,
  subset = site == "high",
  control = glmerControl(optimizer = "bobyqa")
)
summary(m_flower_number_warm_nb)



# predict number of flowers ----------------------------------------------------------------------
emm_fl_num_warm_nb <- emmeans(
  m_flower_number_warm_nb,
  ~ treat_warming | treat_competition,
  type = "response"
)

emm_fl_num_df_warm_nb <- as.data.frame(
  summary(emm_fl_num_warm_nb, infer = TRUE)
) |>
  rename(
    emmean = response,
    lower.CL = asymp.LCL,
    upper.CL = asymp.UCL
  )


# contrasts ---------------------------------------------------------------
contr_warm_nb <- contrast(emm_fl_num_warm_nb, method = "pairwise") |>
  as.data.frame() |>
  mutate(type = "warming",
         sig = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ "ns"
         ))

contr_comp_nb <- contrast(
  emmeans(m_flower_number_warm_nb, ~ treat_competition | treat_warming),
  method = "pairwise"
) |>
  as.data.frame() |>
  mutate(type = "competition",
         sig = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ "ns"
         ))



# annotations -------------------------------------------------------------
warm_levels <- levels(emm_fl_num_df_warm_nb$treat_warming)

ann_warm_nb <- contr_warm_nb |> 
  mutate(
    w1 = warm_levels[1],
    w2 = warm_levels[2],
    x_shift = ifelse(treat_competition == "with", -offset, +offset),
    xmin = as.numeric(factor(w1, levels = warm_levels)) + x_shift,
    xmax = as.numeric(factor(w2, levels = warm_levels)) + x_shift,
    y = y_max + row_number() * spacing,
    label = sig
  )

ann_comp_nb <- contr_comp_nb |> 
  mutate(
    warm_x = as.numeric(factor(treat_warming, levels = warm_levels)),
    xmin = warm_x - offset,
    xmax = warm_x + offset,
    y = y_max + (nrow(ann_warm_nb) + row_number()) * spacing,
    label = sig
  )



# raw data ----------------------------------------------------------------
raw_means <- phenology |>
  filter(site == "high") |> 
  group_by(treat_warming, treat_competition, species, block_ID) |> 
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")



# final plot ---------------------------------------------------------
p_nb_warm <- ggplot() +
  
  geom_jitter(
    data = raw_means,
    aes(x = treat_warming, y = mean_value, color = treat_competition),
    width = 0.1, alpha = 0.3, size = 1.5
  ) +
  
  geom_point(
    data = emm_fl_num_df_warm_nb,
    aes(x = treat_warming, y = emmean, color = treat_competition),
    size = 5,
    position = dodge
  ) +
  
  geom_errorbar(
    data = emm_fl_num_df_warm_nb,
    aes(x = treat_warming, ymin = lower.CL, ymax = upper.CL,
        color = treat_competition),
    width = 0.2,
    linewidth = 0.8,
    position = dodge
  ) +
  
  geom_line(
    data = emm_fl_num_df_warm_nb,
    aes(x = treat_warming, y = emmean, color = treat_competition,
        group = treat_competition),
    position = dodge
  ) +
  
  # warming effect
  bracket_geoms(ann_warm_nb) +
  
  # competition effect
  bracket_geoms(ann_comp_nb) +
  
  labs(
    x = "Warming treatment",
    y = "Predicted mean number of flowers",
    color = "Biotic interactions",
    title = "Effect of warming and competition on flower number\n(glmer.nb)"
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))

p_nb_warm

# ggsave(filename = "Output/Biomass/Warming_flower_number_NOR_glmer.nb.png", 
#        plot = p_nb_warm, width = 10, height = 8, units = "in")











