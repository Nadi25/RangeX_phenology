

# 03_temp_sens TMS4 species ------------------------------------------------------------

# Temperature sensitivity analysis -------------------------------------------------

# RangeX phenology effect of transplantation and warming on temperature sensitivity 

# Species level -----------------------------------------------------------


library(purrr)
library(ggrepel)

source("Temperature_sensitivity_analysis_TMS_2.R")

theme_set(theme_bw(base_size = 20))

# Filter correct onset dataset species --------------------------------------------------
# only flowering
# NOR
d_flower_cn <- filter_data_species(onset_all_temp, "cennig", "Flowering", "ambi")
d_flower_cc <- filter_data_species(onset_all_temp, "cyncri", "Flowering", "ambi")
d_flower_hm <- filter_data_species(onset_all_temp, "hypmac", "Flowering", "ambi")
d_flower_lv <- filter_data_species(onset_all_temp, "leuvul", "Flowering", "ambi")
d_flower_lm <- filter_data_species(onset_all_temp, "luzmul", "Flowering", "ambi")
d_flower_ps <- filter_data_species(onset_all_temp, "pimsax", "Flowering", "ambi")
d_flower_pl <- filter_data_species(onset_all_temp, "plalan", "Flowering", "ambi")
d_flower_sup <- filter_data_species(onset_all_temp, "sucpra", "Flowering", "ambi")
d_flower_tp <- filter_data_species(onset_all_temp, "tripra", "Flowering", "ambi")
d_flower_sd <- filter_data_species(onset_all_temp, "sildio", "Flowering", "ambi")

d_flower_bp <- filter_data_species(onset_all_temp, "brapin", "Flowering", "ambi")
d_flower_be <- filter_data_species(onset_all_temp, "broere", "Flowering", "ambi")
d_flower_cj <- filter_data_species(onset_all_temp, "cenjac", "Flowering", "ambi")
d_flower_dc <- filter_data_species(onset_all_temp, "daucar", "Flowering", "ambi")
d_flower_hp <- filter_data_species(onset_all_temp, "hypper", "Flowering", "ambi")
d_flower_ml <- filter_data_species(onset_all_temp, "medlup", "Flowering", "ambi")
d_flower_sv <- filter_data_species(onset_all_temp, "silvul", "Flowering", "ambi")
d_flower_pm <- filter_data_species(onset_all_temp, "plamed", "Flowering", "ambi")
d_flower_sp <- filter_data_species(onset_all_temp, "salpra", "Flowering", "ambi")
d_flower_sc <- filter_data_species(onset_all_temp, "scacol", "Flowering", "ambi")



# sensitivity models ------------------------------------------------------------------
flower_species <- list(
  cennig = d_flower_cn,
  cyncri = d_flower_cc,
  hypmac = d_flower_hm,
  leuvul = d_flower_lv,
  luzmul = d_flower_lm,
  pimsax = d_flower_ps,
  plalan = d_flower_pl,
  sucpra = d_flower_sup,
  tripra = d_flower_tp,
  sildio = d_flower_sd,
  brapin = d_flower_bp,
  broere = d_flower_be,
  cenjac = d_flower_cj,
  daucar = d_flower_dc,
  hypper = d_flower_hp,
  medlup = d_flower_ml,
  silvul = d_flower_sv,
  plamed = d_flower_pm,
  #salpra = d_flower_sp, # fails the loop to model
  scacol = d_flower_sc
)
flower_species


# loop through all species ------------------------------------------------
models_species <- map(flower_species, fit_model_sens_species)


flower_species$hypmac |> 
  ggplot(aes(Tmean, onset)) +
  geom_point() +
  geom_smooth(method = "lm")

flower_species$hypmac |>
  summarise(n = n())

flower_species$hypmac |>
  summarise(
    n_total = n(),
    n_complete = sum(complete.cases(Tmean, onset))
  )

flower_species$hypmac |>
  count(treat_competition, Tmean)

flower_species$cyncri |>
  count(treat_competition, Tmean)

flower_species$sildio |>
  count(treat_competition, Tmean)


imap_dfr(
  flower_species,
  ~ tibble(
    species = .y,
    n = nrow(.x)
  )
) |> 
  arrange(n)


# get model outputs -------------------------------------------------------
model_summaries <- map(models_species, summary)

model_summaries

model_summaries$cennig

model_summaries$hypmac



# get the actual sensitivity from coefficient -----------------------------
temp_sens <- map(models_species, get_temp_sens_coef)

temp_sens$cennig
temp_sens$hypmac


# make sensitivity into one dataframe -------------------------------------
temp_sens_all_species <- bind_rows(temp_sens, .id = "species")
temp_sens_all_species


# add region
temp_sens_all_species <- temp_sens_all_species |>
  mutate(
    region = ifelse(
      species %in% c(
        "cennig", "cyncri", "hypmac", "leuvul", "luzmul",
        "pimsax", "plalan", "sucpra", "tripra", "sildio"
      ),
      "Norway",
      "Switzerland"
    )
  )


# Plot --------------------------------------------------------------------
ts_species <- ggplot(
  temp_sens_all_species,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species)) +
  geom_point(size = 3.5) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ region) +
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Biotic interactions",
    title = "Temperature sensitivity by species hi vs lo ambi"
  ) +
  theme(legend.position = "none")
ts_species


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_lo_hi_species.png", 
#       plot = ts_species,
#       width = 10, height = 20, units = "in")





pd <- position_dodge(width = 0.2)

ts_species <- ggplot(
  temp_sens_all_species,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species
  )
) +
  geom_pointrange(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = pd,
    linewidth = 0.5
  ) +
  geom_line(
    position = pd,
    linewidth = 0.8
  ) +
  
  geom_text_repel(
    data = subset(
      temp_sens_all_species,
      treat_competition == "with"
    ),
    aes(label = species),
    segment.color = NA,
    show.legend = FALSE,
    max.overlaps = Inf
  )+
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ region) +
  theme(legend.position = "none")
ts_species


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_lo_hi_species.png", 
#       plot = ts_species,
#       width = 15, height = 20, units = "in")



# plot with label of species ----------------------------------------------
label_data <- temp_sens_all_species |> 
  group_by(region, species) |>
  arrange(
    is.na(Tmean.trend),          # valid values first
    treat_competition != "with"  # prefer with
  ) |>
  slice(1) |>
  ungroup()


label_data |>
  select(species, region, treat_competition, Tmean.trend)


ts_species_label <- ggplot(
  temp_sens_all_species,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species)) +
  
  geom_pointrange(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = pd,
    linewidth = 0.5) +
  geom_line(
    position = pd,
    linewidth = 0.8) +
  
  geom_text_repel(
    data = label_data,
    aes(label = species),
    size = 7,
    segment.color = NA,
    show.legend = FALSE,
    max.overlaps = Inf )+
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  facet_wrap(~ region) +
  theme(legend.position = "none")+
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Biotic interactions",
    title = "Temperature sensitivity flowering per species")
ts_species_label

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_lo_hi_flowering_species.png", 
#       plot = ts_species_label,
#       width = 18, height = 15, units = "in")




# Early vs late -----------------------------------------------------------

# get mean flowering onset per species ------------------------------------
onset_flower_general <- phenology2 |>
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage, functional_group) |>
  summarise(onset = min(jday), .groups = "drop") |>
  filter(is.finite(onset))
onset_flower_general


mean_onset_flower_species_region <- onset_flower_general |>
  group_by(region, species) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_flower_species_region



# define early and late flowering species ---------------------------------
# based on mean onset per species and region
species_flowering <- mean_onset_flower_species_region |>
  group_by(region) |>
  mutate(cutoff = median(mean_onset),
         flowering_time = if_else(mean_onset <= cutoff, "early", "late"))
species_flowering


table(species_flowering$region, species_flowering$flowering_time)


# join label data with early late type ------------------------------------
temp_sens_all_species <- temp_sens_all_species |> 
  left_join(species_flowering, by = c("region", "species")) |> 
  select(-c(cutoff, mean_onset))
temp_sens_all_species


label_data <- label_data |> 
  left_join(species_flowering, by = c("region", "species")) |> 
  select(-c(cutoff, mean_onset))

ts_species_label <- ggplot(
  temp_sens_all_species,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species,
    shape = flowering_time)) +
  
  geom_pointrange(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = pd,
    linewidth = 0.5,
    size = 1.5) +
  geom_line(
    position = pd,
    linewidth = 0.8) +
  
  geom_text_repel(
    data = label_data,
    aes(label = species),
    size = 7,
    segment.color = NA,
    show.legend = FALSE,
    max.overlaps = Inf )+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_shape_manual(values = c("early" = 1, "late" = 18))+
  facet_wrap(~ region) +
  guides(color = "none")+
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Biotic interactions",
    title = "Temperature sensitivity flowering per species hi vs lo",
    shape = "Flowering time")+
  theme(legend.position = "bottom")
ts_species_label


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_lo_hi_flowering_species.png", 
#       plot = ts_species_label,
#       width = 18, height = 15, units = "in")




# Ambient vs warmed  ------------------------------------------------------

# Filter correct onset dataset species --------------------------------------------------
# only flowering
# NOR
d_flower_cn <- filter_data_species_aw(onset_all_temp, "cennig", "Flowering", "hi")
d_flower_cc <- filter_data_species_aw(onset_all_temp, "cyncri", "Flowering", "hi")
d_flower_hm <- filter_data_species_aw(onset_all_temp, "hypmac", "Flowering", "hi")
d_flower_lv <- filter_data_species_aw(onset_all_temp, "leuvul", "Flowering", "hi")
d_flower_lm <- filter_data_species_aw(onset_all_temp, "luzmul", "Flowering", "hi")
d_flower_ps <- filter_data_species_aw(onset_all_temp, "pimsax", "Flowering", "hi")
d_flower_pl <- filter_data_species_aw(onset_all_temp, "plalan", "Flowering", "hi")
d_flower_sup <- filter_data_species_aw(onset_all_temp, "sucpra", "Flowering", "hi")
d_flower_tp <- filter_data_species_aw(onset_all_temp, "tripra", "Flowering", "hi")
d_flower_sd <- filter_data_species_aw(onset_all_temp, "sildio", "Flowering", "hi")

d_flower_bp <- filter_data_species_aw(onset_all_temp, "brapin", "Flowering", "hi")
d_flower_be <- filter_data_species_aw(onset_all_temp, "broere", "Flowering", "hi")
d_flower_cj <- filter_data_species_aw(onset_all_temp, "cenjac", "Flowering", "hi")
d_flower_dc <- filter_data_species_aw(onset_all_temp, "daucar", "Flowering", "hi")
d_flower_hp <- filter_data_species_aw(onset_all_temp, "hypper", "Flowering", "hi")
d_flower_ml <- filter_data_species_aw(onset_all_temp, "medlup", "Flowering", "hi")
d_flower_sv <- filter_data_species_aw(onset_all_temp, "silvul", "Flowering", "hi")
d_flower_pm <- filter_data_species_aw(onset_all_temp, "plamed", "Seeds", "hi")
d_flower_sp <- filter_data_species_aw(onset_all_temp, "salpra", "Seeds", "hi")
d_flower_sc <- filter_data_species_aw(onset_all_temp, "scacol", "Seeds", "hi")



# sensitivity models ------------------------------------------------------------------
flower_species <- list(
  cennig = d_flower_cn,
  cyncri = d_flower_cc,
  #hypmac = d_flower_hm,
  leuvul = d_flower_lv,
  luzmul = d_flower_lm,
  pimsax = d_flower_ps,
  plalan = d_flower_pl,
  sucpra = d_flower_sup,
  tripra = d_flower_tp,
  sildio = d_flower_sd,
  #brapin = d_flower_bp,
  broere = d_flower_be,
  cenjac = d_flower_cj,
  daucar = d_flower_dc,
  hypper = d_flower_hp,
  medlup = d_flower_ml,
  silvul = d_flower_sv,
  #plamed = d_flower_pm,
  #salpra = d_flower_sp, # fails the loop to model
  scacol = d_flower_sc
)
flower_species


# loop through all species ------------------------------------------------
models_species_aw <- map(flower_species, fit_model_sens_species)


imap_dfr(
  flower_species,
  ~ tibble(
    species = .y,
    n = nrow(.x)
  )
) |> 
  arrange(n)


# get outputs -------------------------------------------------------------
model_summaries_aw <- map(models_species_aw, summary)

model_summaries_aw$cennig

model_summaries_aw$hypmac



# get temp sens ----------------------------------------------------------
temp_sens_aw <- map(models_species_aw, get_temp_sens_coef)

temp_sens_aw$cennig
temp_sens_aw$hypmac


# make one dataframe ------------------------------------------------------
temp_sens_all_species_aw <- bind_rows(temp_sens_aw, .id = "species")
temp_sens_all_species_aw


temp_sens_all_species_aw <- temp_sens_all_species_aw |>
  mutate(
    region = ifelse(
      species %in% c(
        "cennig", "cyncri", "hypmac", "leuvul", "luzmul",
        "pimsax", "plalan", "sucpra", "tripra", "sildio"
      ),
      "Norway",
      "Switzerland"
    )
  )


pd <- position_dodge(width = 0.2)

# plot with label of species ----------------------------------------------
label_data_aw <- temp_sens_all_species_aw |> 
  group_by(region, species) |>
  arrange(
    is.na(Tmean.trend),          # valid values first
    treat_competition != "with"  # prefer with
  ) |>
  slice(1) |>
  ungroup()


label_data_aw |>
  select(species, region, treat_competition, Tmean.trend)

label_data_aw <- label_data_aw |> 
  left_join(species_flowering, by = c("region", "species")) |> 
  select(-c(cutoff, mean_onset))




# add early vs late flowering data ----------------------------------------
temp_sens_all_species_aw <- temp_sens_all_species_aw |>
  left_join(species_flowering, by = c("region", "species")) |>
  select(-c(cutoff, mean_onset))



# Plot ambi vs warm per species flowering ---------------------------------
ts_species_label_aw <- ggplot(
  temp_sens_all_species_aw,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species,
    shape = flowering_time
  )
) +
  geom_pointrange(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = pd,
    linewidth = 0.5
  ) +
  geom_line(
    position = pd,
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = label_data_aw,
    aes(label = species),
    size = 7,
    segment.color = NA,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_shape_manual(
    values = c(
      "early" = 1,
      "late" = 18
    )
  ) +
  facet_wrap(~region) +
  guides(color = "none") +
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Biotic interactions",
    title = "Temperature sensitivity flowering per species ambi vs warm",
    shape = "Flowering time"
  ) +
  theme(legend.position = "bottom")
ts_species_label_aw

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_ambi_warm_flowering_species.png", 
#        plot = ts_species_label_aw,
#        width = 18, height = 15, units = "in")



# Joined figure flowering -------------------------------------------------
temp_sens_all_species$comparison <- "low high"
temp_sens_all_species_aw$comparison <- "ambient warmed"

temp_sens_combined <- bind_rows(
  temp_sens_all_species,
  temp_sens_all_species_aw
)

label_data$comparison <- "low high"
label_data_aw$comparison <- "ambient warmed"

label_data_combined <- bind_rows(
  label_data,
  label_data_aw
)


temp_sens_flo_comb <- ggplot(
  temp_sens_combined,
  aes(
    x = treat_competition,
    y = Tmean.trend,
    color = species,
    group = species,
    shape = flowering_time
  )
) +
  geom_pointrange(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = pd,
    linewidth = 0.5
  ) +
  geom_line(
    position = pd,
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = label_data_combined,
    aes(label = species),
    size = 7,
    segment.color = NA,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_shape_manual(
    values = c(
      "early" = 1,
      "late" = 18
    )
  ) +
  facet_grid(region~ comparison, scales = "free") +
  guides(color = "none") +
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Biotic interactions",
    title = "Temperature sensitivity flowering per species hi vs lo and ambi vs warm",
    shape = "Flowering time"
  ) +
  theme(legend.position = "bottom")
temp_sens_flo_comb


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_lo_ambi_warm_flowering_species.png", 
#        plot = temp_sens_flo_comb,
#        width = 18, height = 18, units = "in")

