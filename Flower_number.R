


# Flower number predictions -----------------------------------------------

# Effect of transplantation and warming on flower number NOR glmer.nb ---------------------------------

## Data used: 
## Date:      11.03.2026
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation and warming on flower number NOR with glmer.nb


# comments ----------------------------------------------------------------

# NOR.hi.ambi.vege.wf.09.13.1 has NA in value but why

# CHE: assuming that NA means plant not found
# because there are many 0


# load library ------------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(ggplot2)
library(gt)
library(multcomp)
library(multcompView)
library(gtsummary)


# set theme for plots ------------------------------------
theme_set(theme_bw())


# source the phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use phenology2 which has combined site and temperature treatment

# compare low site ambient with high site ambient = site effect



# factor treatment --------------------------------------------------------
phenology2$treat_competition <- factor(phenology2$treat_competition)


# filter only Norway ------------------------------------------------------
phenology_flower_count_nor <- phenology2 |>
  filter(region == "Norway")


# exclude cennig and sildio -----------------------------------------------
# cennig is a calculation
# sildio is not correct counts because we also counted stems
# 
phenology_flower_count_nor <- phenology_flower_count_nor |> 
  filter(!species %in% c("cennig", "sildio"))



# filter flowers ----------------------------------------------------------
phenology_flower_count_nor <- phenology_flower_count_nor |>
  filter(phenology_stage == "No_FloOpen")








# presence/absence of flowers ---------------------------------------------

flower_pre_abs <- phenology_flower_count_nor |> 
  group_by(site, treatment_site_temp, species, treat_competition, block_ID, 
           unique_plant_ID) |>
  summarise(
    flower_presence = as.integer(any(value > 0)),
    .groups = "drop")

m_flowering <- glmer(flower_presence ~ treatment_site_temp * treat_competition
                     + (1|species) + (1|block_ID), 
                     family = binomial, data = flower_pre_abs)


summary(m_flowering)

car::Anova(m_flowering)

emmeans(m_flowering,
        pairwise ~ treatment_site_temp * treat_competition)

emm_flower <- emmeans(
  m_flowering,
  pairwise ~ treatment_site_temp * treat_competition,
  type = "response")
emm_flower


make_predictions_flowering <- function(model, data) {
  
  biomass_mean <- mean(data$pred_log_biomass_species, na.rm = TRUE)
  
  newdat <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_warm", "hi_ambi"),
    treat_competition = c("with", "without")) |>
    as_tibble()
  
  pred <- predict(
    model,
    newdata = newdat,
    re.form = NA,
    type = "link",
    se.fit = TRUE
  )
  
  newdat$fit <- pred$fit
  newdat$se.fit <- pred$se.fit
  
  newdat |>
    mutate(
      probability = plogis(fit),
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit)
    )
}

pred_flowering <-
  make_predictions_flowering(
    m_flowering,
    flower_pre_abs)
pred_flowering




flower_obs <- flower_pre_abs |>
  group_by(treatment_site_temp, treat_competition) |>
  summarise(
    obs_prop = mean(flower_presence),
    .groups = "drop"
  )
flower_obs

flower_obs_species <- flower_pre_abs |>
  group_by(treatment_site_temp, treat_competition, species) |>
  summarise(
    obs_prop = mean(flower_presence),
    .groups = "drop"
  )
flower_obs_species

pd <- position_dodge(width = 0.7)

pred_flowering <- pred_flowering |>
  mutate(
    point_fill = case_when(
      treatment_site_temp == "hi_warm" ~ "white",
      treat_competition == "with" ~ "#528B8B",
      treat_competition == "without" ~ "#CD950C"))

flower_presence_plot <- ggplot(
  pred_flowering,
  aes(
    x = treatment_site_temp,
    y = probability,
    colour = treat_competition,
    shape = treatment_site_temp)) +
  #raw 
  geom_col(
    data = flower_obs,
    aes(y = obs_prop,
        fill = treat_competition),
    alpha = 0.1,
    position = pd,
    width = 0.3, linewidth = 0.01,
    show.legend = F)+
  # confidence intervals
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.1,
    position = pd) +
  
  geom_point(position = pd, size = 2.5, stroke = 0.4,
             aes(fill = point_fill)) +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C")) +
  
  scale_fill_manual(values = c("#528B8B", "#CD950C", "white", "#528B8B", "white")) +
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 24)) +
  
  scale_x_discrete(
    labels = c(
      "lo_ambi" = "low ambient",
      "hi_ambi" = "high ambient",
      "hi_warm" = "high warmed")) +
  
  labs(
    x = "Site temperature treatment",
    y = "Probability of flowering",
    colour = "Biotic interactions") +
  
  scale_y_continuous(
    breaks = seq(0,1,0.2)) +
  coord_cartesian(
    ylim = c(0,1.25),
    clip = "off")+
  geom_jitter(
    data = flower_obs_species,
    aes(
      x = treatment_site_temp,
      y = obs_prop,
      color = treat_competition),
    position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
    alpha = 0.3)+
  
  theme(legend.position = "bottom")+
  guides(
    shape = "none",
    fill = "none")
flower_presence_plot




# plot with stars as significances ----------------------------------------

emm_flower_contr <- emm_flower$contrasts |> 
  as.data.frame()

emm_flower_contr <- emm_flower_contr |> mutate(stars = case_when(p.value < 0.001 ~ "***",
                                                                         p.value < 0.01 ~ "**",
                                                                         p.value < 0.05 ~ "*",
                                                                         TRUE ~ ""))
emm_flower_contr

emm_flower_contr2 <- emm_flower_contr |> 
  filter(contrast %in% c("lo_ambi with / lo_ambi without",
                         "hi_warm with / hi_warm without",
                         "hi_ambi with / hi_ambi without",
                         "lo_ambi with / hi_ambi with",
                         "lo_ambi without / hi_ambi without",
                         "hi_warm with / hi_ambi with",
                         "hi_warm without / hi_ambi without"))
emm_flower_contr2

pos_df <- tibble(
  group = c(
    "lo_ambi with",
    "lo_ambi without",
    "hi_warm with",
    "hi_warm without",
    "hi_ambi with",
    "hi_ambi without"),
  x = c(
    0.775, 1.225,
    1.775, 2.225,
    2.775, 3.225))



brackets <- emm_flower_contr2 |>
  separate(
    contrast,
    into = c("group1", "group2"),
    sep = " / "
  )
brackets <- brackets |>
  left_join(pos_df, by = c("group1" = "group")) |>
  rename(xmin = x) |>
  left_join(pos_df, by = c("group2" = "group")) |>
  rename(xmax = x)
brackets


brackets |>
  select(group1, group2, xmin, xmax, stars)

top_y <- max(pred_flowering$upper)


brackets <- brackets |>
  mutate(
    y = top_y + c(0.09, 0.03, 0.2, 0.03, 0.03, 0.15, 0.25))
brackets


pd <- position_dodge(width = 0.7) 

theme_set(theme_bw(base_size = 8))

flower_presence_plot2 <- flower_presence_plot +
  geom_segment(
    data = brackets, linewidth = 0.25,
    aes(x = xmin, xend = xmax,
        y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,linewidth = 0.25,
    aes(x = xmin, xend = xmin,
        y = y, yend = y - 0.02),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,linewidth = 0.25,
    aes(x = xmax, xend = xmax,
        y = y, yend = y - 0.02),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = brackets,
    aes(
      x = (xmin + xmax)/2,
      y = y + 0.015,
      label = stars
    ),
    inherit.aes = FALSE)
flower_presence_plot2

# ggsave(filename = "Output/Biomass/Figures_for_manuscript/Flower_presence_absence_NOR.png", 
#        plot = flower_presence_plot2,
#        width = 9, height = 8, units = "cm", dpi = 600)







# get max number of flowers -----------------------------------------------
max_flower_per_plant <- phenology_flower_count_nor |>
  group_by(site, treatment_site_temp, species, treat_competition, block_ID, unique_plant_ID) |>
  summarise(
    max_flower_number = max(value, na.rm = TRUE),
    .groups = "drop")



max_flower_per_plant <- 
  max_flower_per_plant |>
  filter(max_flower_number > 0)





# fit the model negative binomial  ------------------------------------
m_flower_number <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition + 
                              (1|species) + (1|block_ID),
                            data = max_flower_per_plant,
                            control = glmerControl(optimizer = "bobyqa"))

summary(m_flower_number)

car::Anova(m_flower_number)


tbl_regression(
  m_flower_number,
  exponentiate = TRUE
)

t1 <- tbl_regression(
  m_flower_number,
  exponentiate = TRUE) |>
  as_gt() |>
  tab_header(
    title = md("**Number of flowers NOR glm.nb**"))
t1
#gtsave(t1, "Output/Biomass/Number_flowers_model_summary_NOR.docx")




emm <- emmeans(m_flower_number,
               pairwise ~ treatment_site_temp * treat_competition)
emm


get_signi <- function(model) {
  
  signi_letters <- cld(emmeans(
    model,
    pairwise ~ treatment_site_temp * treat_competition),
    Letters = letters)
  
  print(signi_letters)
  return(signi_letters)
}

emm_sig <- get_signi(m_flower_number)






performance::model_performance(m_flower_number)

pred <- ggpredict(m_flower_number, c("treatment_site_temp", "treat_competition"))
plot(pred)


make_predictions <- function(model) {
  
  newdat <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_warm", "hi_ambi"),
    treat_competition = c("with", "without")) |> 
    as_tibble()
  
  pred <- predict(
    model,
    newdat = newdat,
    re.form = NA,
    type = "link", #  log scale
    se.fit = TRUE) |> 
    
    as_tibble() |>
    mutate(upper = exp(fit + 1.96 * se.fit),
           lower = exp(fit - 1.96 * se.fit), 
           fit = exp(fit)) |>
    bind_cols(newdat)
  
  return(pred)
}



flower_number_pred <- make_predictions(m_flower_number)
flower_number_pred


# rename treat info -------------------------------------------------------
# flower_number_pred$treatment_site_temp <- 
#   recode(flower_number_pred$treatment_site_temp,
#   "lo_ambi" = "low ambient",
#   "hi_ambi" = "high ambient",
#   "hi_warm" = "high warm")
# 
# max_flower_per_plant$treatment_site_temp <- 
#   recode(max_flower_per_plant$treatment_site_temp,
#   "lo_ambi" = "low ambient",
#   "hi_ambi" = "high ambient",
#   "hi_warm" = "high warm")



# this should be new script
# plot --------------------------------------------------
pd <- position_dodge(width = 0.9) 

flow_no <- ggplot(flower_number_pred, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,              
  shape = treatment_site_temp  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = .2,
                position = pd) +
  
  geom_jitter(
    data = max_flower_per_plant,
    aes(
      x = treatment_site_temp,
      y = max_flower_number,
      color = treat_competition
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9
    ),
    alpha = 0.05,
    size = 2,
    inherit.aes = FALSE
  )+
  
  scale_color_manual(values = c("#528B8B", "#CD950C")) +   
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "(Predicted) max flower number",
    title = "Effect of transplantation and warming on flower number",
    shape = "Treatment site × warming",
    color = "Biotic interactions"
  )+
  guides(shape = "none")
flow_no

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb.png", 
#       plot = flow_no, width = 12, height = 8, units = "in")



# Plot with raw as violin ------------------------------------------------

flow_no2 <- ggplot(flower_number_pred, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # model estimates
  geom_point(position = pd, size = 4, stroke = 1.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = .2,
    position = pd
  ) +
  
  # distribution of raw data
  geom_violin(
    data = max_flower_per_plant,
    aes(
      x = treatment_site_temp,
      y = max_flower_number,
      fill = treat_competition
    ),
    position = position_dodge(width = 0.9),
    alpha = 0.25,
    color = NA,
    trim = TRUE # cuts of at 0 because we dont have negative flower no
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C")) +
  scale_fill_manual(values = c("#528B8B", "#CD950C")) +
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 2
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "Log (max flower number)",
    title = "Effect of transplantation and warming on flower number",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  
  guides(shape = "none")+
  coord_transform(y = "log1p")+
  theme(legend.position = "bottom")

flow_no2

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb_violin.png", 
#       plot = flow_no2, width = 12, height = 8, units = "in")





# plot with significance letters ------------------------------------------
emm_sig <- emm_sig |>
  mutate(
    y = asymp.UCL + 10
  )
emm_sig


flow_no3 <- flow_no2 +
  geom_text(
    data = emm_sig,
    aes(
      x = treatment_site_temp,
      y = y,
      label = .group, 
      group = treat_competition),
    color = "grey43",
    position = pd,
    size = 5,
    show.legend = FALSE)
flow_no3







# plot with stars as significances ----------------------------------------

emm_contr <- emm$contrasts |> 
  as.data.frame()

emm_contr <- emm_contr |> mutate(stars = case_when(p.value < 0.001 ~ "***",
                                                   p.value < 0.01 ~ "**",
                                                   p.value < 0.05 ~ "*",
                                                   TRUE ~ ""))
emm_contr

emm_contr2 <- emm_contr |> 
  filter(contrast %in% c("lo_ambi with - lo_ambi without",
                         "hi_warm with - hi_warm without",
                         "hi_ambi with - hi_ambi without",
                         "lo_ambi with - hi_ambi with",
                         "lo_ambi without - hi_ambi without",
                         "hi_warm with - hi_ambi with",
                         "hi_warm without - hi_ambi without"))

pos_df <- tibble(
  group = c(
    "lo_ambi with",
    "lo_ambi without",
    "hi_warm with",
    "hi_warm without",
    "hi_ambi with",
    "hi_ambi without"
  ),
  x = c(
    0.775, 1.225,
    1.775, 2.225,
    2.775, 3.225
  )
)



brackets <- emm_contr2 |>
  separate(
    contrast,
    into = c("group1", "group2"),
    sep = " - "
  )
brackets <- brackets |>
  left_join(pos_df, by = c("group1" = "group")) |>
  rename(xmin = x) |>
  left_join(pos_df, by = c("group2" = "group")) |>
  rename(xmax = x)
brackets


brackets |>
  select(group1, group2, xmin, xmax, stars)

top_y <- max(flower_number_pred$upper)

brackets <- brackets |>
  mutate(
    y = top_y + c(18, 6, 27, 6, 6, 13, 22)
  )
brackets

theme_set(theme_bw(base_size = 20))

flow_no2_sig <- flow_no2 +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmax,
        y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmin,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmax, xend = xmax,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = brackets,
    aes(
      x = (xmin + xmax)/2,
      y = y + 0.5,
      label = stars
    ),
    inherit.aes = FALSE
  )
flow_no2_sig


# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb_violin_sig.png", 
#       plot = flow_no2_sig, width = 12, height = 9, units = "in")








# figure for manuscript ---------------------------------------------------
# with jitter instead of violins for raw data

flower_number_pred <- flower_number_pred |>
  mutate(
    point_fill = case_when(
      treatment_site_temp == "hi_warm" ~ "white",
      treat_competition == "with" ~ "#528B8B",
      treat_competition == "without" ~ "#CD950C"))

pd <- position_dodge(width = 0.7) 

theme_set(theme_bw(base_size = 8))


flow_no_4 <- ggplot(flower_number_pred, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp)) +
  
  # distribution of raw data
  geom_jitter(
    data = max_flower_per_plant,
    aes(
      x = treatment_site_temp,
      y = max_flower_number,
      color = treat_competition),
    position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.25),
    alpha = 0.14) +
  
  # model estimates
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    position = pd) +
  geom_point(position = pd, size = 2.5, stroke = 0.8,
             aes(fill = point_fill)) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C")) +
  scale_fill_manual(values = c("#528B8B", "#CD950C", "white")) +
  
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 24)) +
  scale_x_discrete(
    labels = c(
      "lo_ambi" = "low ambient",
      "hi_ambi" = "high ambient",
      "hi_warm" = "high warmed"))+
  labs(
    x = "Site temperature treatment",
    y = "Log (max flower number)",
    title = "Effect of transplantation and warming on flower number adjusted for biomass interactive",
    shape = "Treatment site × warming",
    color = "Biotic interactions") +
  
  guides(shape = "none",
         fill = "none")+
  coord_transform(y = "log1p")+
  theme(legend.position = "bottom",
        plot.title = element_blank())
flow_no_4


flow_no_4_sig <- flow_no_4 +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmax,
        y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmin,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmax, xend = xmax,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = brackets,
    aes(
      x = (xmin + xmax)/2,
      y = y + 0.5,
      label = stars
    ),
    inherit.aes = FALSE)
flow_no_4_sig

ggsave(filename = "Output/Biomass/Figures_for_manuscript/Flower_number_NOR.png", 
       plot = flow_no_4_sig,
       width = 9, height = 8, units = "cm", dpi = 600)


tab2 <- emm_contr2 |>
  gt() |>
  fmt_number(
    columns = c(estimate, SE),
    decimals = 2
  ) |>
  fmt_number(
    columns = p.value,
    decimals = 3
  ) |>
  cols_label(
    contrast = "Contrast",
    estimate = "Estimate",
    SE = "SE",
    p.value = "p",
    stars = ""
  ) |>
  cols_hide(c(df, z.ratio)) |>
  tab_header(
    title = md("**Number of flowers comparisons**")
  )
tab2

#gtsave(tab2, "Output/Biomass/Number_flowers_signif_NOR.docx")


