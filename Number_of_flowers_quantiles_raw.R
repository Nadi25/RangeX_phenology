

# Number of flowers raw data quantiles ------------------------------------


# load library ------------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(ggplot2)

# set theme for plots ------------------------------------
theme_set(theme_bw())


# source the phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use phenology2 which has combined site and temperature treatment

# compare low site ambient with high site ambient = site effect



# factor treatment --------------------------------------------------------
phenology2$treat_competition <- factor(phenology2$treat_competition)


# filter only Norway ------------------------------------------------------
phenology_count_nor <- phenology2 |>
  filter(region == "Norway")


# exclude cennig and sildio -----------------------------------------------
# cennig is a calculation
# sildio is not correct counts because we also counted stems
# 
phenology_count_nor <- phenology_count_nor |> 
  filter(!species %in% c("cennig", "sildio"))



## don't use seeds_collected
phenology_count_nor <- phenology_count_nor |>
  filter(phenology_stage != "No_Seeds")



## calculate median of number of buds, flowers, infructescences
## calculate quantiles (fits to median) instead of sd to avoid having negative values for mean-/+ sd. 
## Negative ranges come when we have many low values and a few high ones
phenology_median_quant <- phenology_count_nor |> 
  group_by(species, date_measurement, site, treat_warming, treat_competition, treatment_site_temp, phenology_stage) |> 
  #group_by(species, date, treatment, Stage) |> 
  summarise(median = median(value),
            lo = quantile(value, probs = 0.1, na.rm = TRUE), 
            hi = quantile(value, probs = 0.9, na.rm = TRUE), 
            .groups = "drop")



# function to filter by species -------------------------------------------



## filter only plalan
phenology_clean_median_quant_plalan <- phenology_median_quant |> 
  filter(species == "plalan")


ggplot(phenology_clean_median_quant_plalan, aes(x = date_measurement)) +
  geom_line(aes(y = median, color = phenology_stage)) +
  geom_ribbon(aes(y = median, ymin = lo, ymax = hi, fill = phenology_stage), linewidth = 2, alpha = 0.5) +
  facet_grid(rows = vars(treatment_site_temp), cols = vars(treat_competition)) +  # Added `scales = "free"`
  labs(y = "Median +/- quantiles", x = "", title = "Plantago lanceolata") +
  scale_color_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
  scale_fill_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
  theme(legend.position = "top")




# Leuvul ------------------------------------------------------------------

phenology_clean_median_quant_leuvul <- phenology_median_quant |> 
  filter(species == "leuvul")


ggplot(phenology_clean_median_quant_leuvul, aes(x = date_measurement)) +
  geom_line(aes(y = median, color = phenology_stage)) +
  geom_ribbon(aes(y = median, ymin = lo, ymax = hi, fill = phenology_stage), linewidth = 2, alpha = 0.5) +
  facet_grid(rows = vars(treatment_site_temp), cols = vars(treat_competition)) +  # Added `scales = "free"`
  labs(y = "Median +/- quantiles", x = "", title = "Leucanthemum vulgare") +
  scale_color_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
  scale_fill_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
  theme(legend.position = "top")





filter_data_species_aw <- function(phenology_data, species_name) {
  
  phenology_median_quant_s <- phenology_data |> 
    filter(species == "species_name")
  
  p <- ggplot(phenology_median_quant_s, aes(x = date_measurement)) +
    geom_line(aes(y = median, color = phenology_stage)) +
    geom_ribbon(aes(y = median, ymin = lo, ymax = hi, fill = phenology_stage), linewidth = 2, alpha = 0.5) +
    facet_grid(rows = vars(treatment_site_temp), cols = vars(treat_competition)) +  # Added `scales = "free"`
    labs(y = "Median +/- quantiles", x = "", title = "Leucanthemum vulgare") +
    scale_color_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
    scale_fill_manual(values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")) +
    theme(legend.position = "top")
  
  
  return(phenology_median_quant_s)
  return(p)
}

plot_phenology_species <- function(data, species_name, species_title) {
  
  data |>
    filter(species == species_name) |>
    ggplot(aes(x = date_measurement)) +
    geom_line(aes(y = median, color = phenology_stage)) +
    geom_ribbon(
      aes(y = median, ymin = lo, ymax = hi, fill = phenology_stage),
      alpha = 0.5
    ) +
    facet_grid(
      rows = vars(treatment_site_temp),
      cols = vars(treat_competition)
    ) +
    labs(
      y = "Median +/- quantiles",
      x = "",
      title = species_title
    ) +
    scale_color_manual(
      values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")
    ) +
    scale_fill_manual(
      values = c("#00BB00FF", "#500050FF", "#FFBBFFFF")
    ) +
    theme(legend.position = "top")
}



p_cyncri <- plot_phenology_species(
  phenology_median_quant,
  "cyncri",
  "Cynosurus cristatus")
p_cyncri

p_hypmac <- plot_phenology_species(
  phenology_median_quant,
  "hypmac",
  "Hypericum maculatum")
p_hypmac

p_leuvul <- plot_phenology_species(
  phenology_median_quant,
  "leuvul",
  "Leucanthemum vulgare")
p_leuvul

p_luzmul <- plot_phenology_species(
  phenology_median_quant,
  "luzmul",
  "Luzula multiflora")
p_luzmul

p_pimsax <- plot_phenology_species(
  phenology_median_quant,
  "pimsax",
  "Pimpinella saxifraga")
p_pimsax

p_plalan <- plot_phenology_species(
  phenology_median_quant,
  "plalan",
  "Plantago lanceolata")
p_plalan

p_sucpra <- plot_phenology_species(
  phenology_median_quant,
  "sucpra",
  "Succisa pratensis")
p_plalan

p_tripra <- plot_phenology_species(
  phenology_median_quant,
  "tripra",
  "Trifolium pratense")
p_tripra












