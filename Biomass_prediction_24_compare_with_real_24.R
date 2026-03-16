
# BIOMASS 4 ---------------------------------------------------------------

# RangeX biomass predictions 24 ------------

## Data used: 
## Date:      07.03.26
## Author:    Nadine Arzt
## Purpose:   Figure out which prediction method to use - one general model across species
##            or species specific models


# set theme and colors ----------------------------------------------------
theme_set(theme_bw(base_size = 22))

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


# source script with all the species traits models ------------------------
source("Biomass_traits_correlation_per_species.R")


# make list with all models -----------------------------------------------
species_models <- list(
  sucpra = m_sucpra,
  cennig = m_cennig,
  pimsax = m_pimsax,
  luzmul = m_luzmul,
  leuvul = m_leuvul,
  tripra = m_tripra,
  hypmac = m_hypmac,
  plalan = m_plalan,
  cyncri = m_cyncri,
  sildio = m_sildio
)



# predict biomass 24 with function per species ----------------------------
predict_species <- function(sp, model, data){
  
  df <- data |> filter(species == sp)
  
  # detect variables used in the model
  vars <- all.vars(formula(model))[-1]
  
  df <- df |> drop_na(all_of(vars))
  
  df$pred_log_biomass_species <- predict(
    model,
    newdata = df,
    re.form = NA
  )
  
  df
}



# perform prediction and have one df --------------------------------------
df_2024_pred_species <- purrr::map2_dfr(
  names(species_models),
  species_models,
  predict_species,
  data = analysis_data_24_log
)


# plot all species together -----------------------------------------------
f <- ggplot(df_2024_pred_species,
            aes(log_biomass,
                pred_log_biomass_species,
                color = species,
                shape = functional_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8))+
  labs(x = "log(biomass real 24)",
       y = "log(pred biomass 24 species model)",
       color = "Species",
       shape = "Functional group")
f


# plot separated per species ----------------------------------------------
g <- ggplot(df_2024_pred_species,
            aes(log_biomass,
                pred_log_biomass_species,
                color = species,
                shape = functional_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(16, 17, 15, 3, 8))+
  labs(x = "log(biomass real 24)",
       y = "log(pred biomass 24 species model)",
       color = "Species",
       shape = "Functional group")+
  theme(legend.position = "none")+
  facet_wrap(~ species)
g

ggsave(filename = "Output/Biomass/Log(biomass24)_log(pred_biomass24_species_seperated).png", 
       plot = g, width = 12, height = 9, units = "in")













