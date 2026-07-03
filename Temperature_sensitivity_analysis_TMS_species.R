

library(purrr)
library(ggrepel)

source("Temperature_sensitivity_analysis_TMS.R")

theme_set(theme_bw(base_size = 20))

# Filter correct onset dataset species --------------------------------------------------
# only flowering
# NOR
d_flower_cn <- filter_data_species(onset_all_gs, "cennig", "Flowering", "ambi")
d_flower_cc <- filter_data_species(onset_all_gs, "cyncri", "Flowering", "ambi")
d_flower_hm <- filter_data_species(onset_all_gs, "hypmac", "Flowering", "ambi")
d_flower_lv <- filter_data_species(onset_all_gs, "leuvul", "Flowering", "ambi")
d_flower_lm <- filter_data_species(onset_all_gs, "luzmul", "Flowering", "ambi")
d_flower_ps <- filter_data_species(onset_all_gs, "pimsax", "Flowering", "ambi")
d_flower_pl <- filter_data_species(onset_all_gs, "plalan", "Flowering", "ambi")
d_flower_sup <- filter_data_species(onset_all_gs, "sucpra", "Flowering", "ambi")
d_flower_tp <- filter_data_species(onset_all_gs, "tripra", "Flowering", "ambi")
d_flower_sd <- filter_data_species(onset_all_gs, "sildio", "Flowering", "ambi")

d_flower_bp <- filter_data_species(onset_all_gs, "brapin", "Flowering", "ambi")
d_flower_be <- filter_data_species(onset_all_gs, "broere", "Flowering", "ambi")
d_flower_cj <- filter_data_species(onset_all_gs, "cenjac", "Flowering", "ambi")
d_flower_dc <- filter_data_species(onset_all_gs, "daucar", "Flowering", "ambi")
d_flower_hp <- filter_data_species(onset_all_gs, "hypper", "Flowering", "ambi")
d_flower_ml <- filter_data_species(onset_all_gs, "medlup", "Flowering", "ambi")
d_flower_sv <- filter_data_species(onset_all_gs, "silvul", "Flowering", "ambi")
d_flower_pm <- filter_data_species(onset_all_gs, "plamed", "Flowering", "ambi")
d_flower_sp <- filter_data_species(onset_all_gs, "salpra", "Flowering", "ambi")
d_flower_sc <- filter_data_species(onset_all_gs, "scacol", "Flowering", "ambi")



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

models_species <- map(flower_species, fit_model_sens_species)


imap_dfr(
  flower_species,
  ~ tibble(
    species = .y,
    n = nrow(.x)
  )
) |> 
  arrange(n)


model_summaries <- map(models_species, summary)

model_summaries

model_summaries$cennig

model_summaries$hypmac



temp_sens <- map(models_species, get_temp_sens_coef)

temp_sens$cennig
temp_sens$hypmac


temp_sens_all_species <- bind_rows(temp_sens, .id = "species")
temp_sens_all_species


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


ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_lo_hi_species.png", 
      plot = ts_species,
      width = 15, height = 20, units = "in")



