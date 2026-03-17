
# BIOMASS 10 ---------------------------------------------------------------

# Adding predicted biomass species method to phenology ------------

## Data used: 
## Date:      03.03.26
## Author:    Nadine Arzt
## Purpose:   Add pred log biomass species method to phenology max flowers per plant


# colors ------------------------------------------------------------------
colors <- c(
  "#000000",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#8DD3C7",
  "#FB8072",
  "#80B1D3"
)

colors2 <- c(
  "#000000",  # black
  "#E69F00",  # orange
  "#D55E00",  # vermillion
  "#CC79A7",  # magenta
  "#009E73",  # green
  "#F0E442",  # yellow
  "#56B4E9",  # sky blue
  "#0072B2",  # blue
  "#999999",  # grey
  "#A6761D"   # brown
)

# load library ------------------------------------------------------------
library(ggeffects)


# source script that predicts the 23 biomass
source("Biomass_prediction_23_per_species.R")

theme_set(theme_bw(base_size = 22))


# use
bio_pred_23_species


# import pheno 23 data --------------------------------------------------------
# this is the latest clean version of phenology data
# what changed is the unique_plant_ID because I added the number in the end 
# a 1 since all plants survived without being replaced
phenol <- read.csv("Data/Clean/RangeX_clean_Phenology_2023_NOR.csv")


# import metadata NOR -----------------------------------------------------
meta_NOR <- read.csv("Data/RangeX_clean_MetadataFocal_NOR.csv")


# combine pheno with metadata ---------------------------------------------
phenol_23_NOR <- left_join(meta_NOR, phenol, by = c("unique_plant_ID", "species"))


# rename pheno stages to match regions ------------------------------------
# "number_infructescences" in NOR correpsonfs to "No_FloWithrd" in CHE?
phenol_23_NOR <- phenol_23_NOR |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "number_buds" = "No_Buds",
                                  "number_flowers" = "No_FloOpen",
                                  #"number_infructescences" = "No_FloWithrd",
                                  "seeds_collected" = "No_Seeds"))



bio_23_species <- bio_pred_23_species[, c("unique_plant_ID",
                                  "pred_log_biomass")]

# combine pheno with pred biomass -----------------------------------------
phenol_23_NOR_bio_species <- phenol_23_NOR |>
  left_join(bio_23_species,
            by = "unique_plant_ID") |> 
  rename(pred_log_biomass_species = pred_log_biomass)

## 

# plot pred biomass against max no of flowers -----------------------------
# get max number of flowers per plant
flowers_per_plant <- phenol_23_NOR_bio_species |>
  group_by(unique_plant_ID, species) |>
  summarise(flowers_max = max(value, na.rm = TRUE))



bio_flower_species <- phenol_23_NOR_bio_species |>
  left_join(flowers_per_plant, by = c("unique_plant_ID", "species"))


# vizualization
ggplot(bio_flower_species, aes(x = flowers_max, y = pred_log_biomass_species)) +
  geom_point(alpha = 0.5, color = "darkred") +
  facet_wrap(~ species) +
  labs(
    x = "Maximum number of flowers",
    y = "Predicted log(biomass)",
    title = "Predicted biomass vs. flower production"
  )

# filter out max flower number 0 ----------------------------------------------
bio_flower_species <- bio_flower_species |>
  filter(flowers_max > 0)

hist(bio_flower_species$flowers_max)
hist(log(bio_flower_species$flowers_max))


# log max flower -----------------------------------------------------------
qqnorm(bio_flower_species$flowers_max); qqline(bio_flower_species$flowers_max)
qqnorm(log(bio_flower_species$flowers_max)); qqline(log(bio_flower_species$flowers_max))

# log 
bio_flower_species <- bio_flower_species |>
  mutate(log_max_flower = log(flowers_max))

# fit model ---------------------------------------------------------------
# to see if the observed trend is significant
m_flowers_species <- lmerTest::lmer(pred_log_biomass_species ~ log_max_flower + (1|species) + (1|block_ID),
                            data = bio_flower_species)
summary(m_flowers_species)

# it is significant: <2e-16 ***

# plot per species
# plot biomass vs max flowers and the trend with a lm
ggplot(bio_flower_species, aes(x = log_max_flower, y = pred_log_biomass_species)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  #facet_wrap(~ species) +
  labs(x = "Log (max number of flowers)",
       y = "Predicted log(biomass)")

# plot with correct model  ------------------------------------------------
# and confidence intervals
pred_df_species <- ggpredict(m_flowers_species, terms = "log_max_flower")

g <- ggplot() +
  geom_point(data = bio_flower_species,
             aes(x = log_max_flower, y = pred_log_biomass_species),
             alpha = 0.4) +
  geom_line(data = pred_df_species,
            aes(x = x, y = predicted),
            color = "turquoise4", size = 1.2) +
  geom_ribbon(data = pred_df_species,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "turquoise4") +
  labs(x = "Log(max number of flowers)",
       y = "Predicted log(biomass)",
       title = "Species models")
g
# the points are the actual log data, the line is the model
# title = "Mixed-model prediction with 95% CI"

# ggsave(filename = "Output/Biomass/Log_pred_biomass_log_max_flowers_lmm_species_models.png", 
#        plot = g, width = 13, height = 9, units = "in")

h <- ggplot() +
  geom_point(data = bio_flower_species,
             aes(x = log_max_flower,
                 y = pred_log_biomass_species,
                 color = species,
                 shape = functional_group),
             alpha = 0.6) +
  geom_line(data = pred_df_species,
            aes(x = x, y = predicted),
            color = "grey7", size = 1.2) +
  geom_ribbon(data = pred_df_species,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "grey7") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16, 17, 15))+
  labs(x = "Log(max number of flowers)",
       y = "Predicted log(biomass)",
       title = "Species models",
       color = "Species",
       shape = "Functional group")
h

# ggsave(filename = "Output/Biomass/Log_pred_biomass_log_max_flowers_lmm_species_models_color.png", 
#        plot = h, width = 13, height = 9, units = "in")

