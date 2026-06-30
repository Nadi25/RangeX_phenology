


# 01_temp_sens ------------------------------------------------------------


# Functions for temperature sensitivity analysis --------------------------


# function to calculate mean onset per plot ------------------------------
# take the average onset for the three individuals per species in one plot
# first budding date per individual
# then average per plot

# data = phenology2
# stage_name = No_Buds, No_FloOpen, No_FloWithrd, No_Seeds
# onset_type = jday or GDD_cum

get_mean_onset <- function(data, stage_name, onset_type) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, treatment_site_temp, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID, functional_group
    ) |>
    summarise(
      onset = min(.data[[onset_type]]),
      .groups = "drop"
    )
}



# # Average across species, site and treat -----------------------------------------
# # this is not needed or at least it doesnt change from get_mean_onset
# # because unique_plot_ID is essentially treat_comp
# # but we still need an average across the plots
# get_plot_onset <- function(onset_data) {
#   onset_data |>
#     group_by(region, treatment_site_temp, treat_competition, species, block_ID, 
#              unique_plot_ID, functional_group) |>
#     summarise(
#       onset = mean(onset, na.rm = TRUE),
#       .groups = "drop"
#     )
# }




# filter dataset ----------------------------------------------------------
filter_data <- function(onset_data, region_name, stage_name, temp_name) {
  onset_region_stage <- onset_data |> 
    filter(region == region_name) |> 
    filter(stage == stage_name) |> 
    filter(treat_warming == temp_name)
  
  return(onset_region_stage)
}


# Model -------------------------------------------------------------------
# separately per region
# filter sensitivity dataset per region first


# lo ambi vs hi ambi ------------------------------------------------------

fit_model_sens_lo_hi <- function(onset_data) {
  
  model <- lmerTest::lmer(
    onset ~ treat_competition * Tmean + (1 | species) + (1 | block_ID),
    data = onset_data
  )
  
  return(model)
}





# predictions per stage ---------------------------------------------------

make_sens_predictions_lh <- function(temp_data, model) {
  
  newdata <- temp_data |>
    select(treat_competition, Tmean)|> 
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
#



# get the actual temperature sensitivity ----------------------------------
# get the coefficient and slope

get_temp_sens_coef <- function(model) {
  
  slopes <- as.data.frame(emtrends(model, ~ treat_competition, var = "Tmean"))
  
  return(slopes)
}

#



































