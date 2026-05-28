
# 01_temp_sens ------------------------------------------------------------


# Functions for temperature sensitivity analysis --------------------------


# function to calculate mean onset per plot ------------------------------
# take the average onset for the three individuals per species in one plot
# first budding date per individual
# then average per plot

# data = phenology3
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
      region, site, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID
    ) |>
    summarise(
      first_onset = min(.data[[onset_type]]),
      .groups = "drop"
    ) |>
    
    # mean onset across 3 individuals within plot
    group_by(
      region, site, treat_competition,
      species, block_ID, unique_plot_ID
    ) |>
    summarise(
      onset = mean(first_onset),
      .groups = "drop"
    )
}



# Average across species, site and treat -----------------------------------------
get_species_onset <- function(onset_data) {
  onset_data |>
    group_by(region, site, treat_competition, species, block_ID) |>
    summarise(
      onset = mean(onset, na.rm = TRUE),
      .groups = "drop"
    )
}



# Calculate temperature sensitivity ---------------------------------------
# function to get temp sens
get_temp_sens <- function(onset_mean_data, temperature_data) {
  onset_mean_data |>
    left_join(temperature_data, by = c("region","site")) |>
    group_by(region, species, treat_competition, block_ID) |>  # or group by site and species?
    pivot_wider(names_from = site,
                values_from = c(onset, Tmean)) |>
    mutate(
      temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
    )
}


# control plots sensitivity -----------------------------------------------
control_plot_temp_sens <- function(sensitivity_data) {
  ggplot(sensitivity_data, aes(x = treat_competition, y = temp_sens, color = species)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
    facet_wrap(~region)
}




# sensitivity model ------------------------------------------------------------------
# function for fitting the model per stage
# fit_sens_model <- function(sens_data) {
#   
#   model <- lmerTest::lmer(
#     temp_sens ~ region * treat_competition + (1 | species) + (1 | block_ID),
#     data = sens_data
#   )
#   
#   return(model)
# }

# separately per region
# filter sensitivity dataset per region first
fit_sens_model <- function(sens_data, region_name) {
  sens_region <- sens_data |> 
    filter(region == region_name)
  
  model <- lmerTest::lmer(
    temp_sens ~ treat_competition + (1 | species) + (1 | block_ID),
    data = sens_region
  )
  
  return(model)
}



# predict sensitivity for each stage with function ------------------------
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

make_sens_predictions <- function(model) {
  
  # 1. new data
  newdat <- expand.grid(
    treat_competition = c("with", "without"),
    temp_sens = 0
  )
  
  # 2. fixed-effect predictions
  newdat$temp_sens <- predict(
    model,
    newdata = newdat,
    re.form = NA
  )
  
  # 3. model matrix
  mm <- model.matrix(terms(model), newdat)
  
  # 4. fixed-effect variance
  pvar <- diag(mm %*% tcrossprod(vcov(model), mm))
  
  # 5. random effects
  var_species <- as.numeric(VarCorr(model)$species)
  var_block <- as.numeric(VarCorr(model)$block_ID)
  
  # 6. total variance
  tvar1 <- pvar + var_species + var_block + sigma(model)^2
  
  # 7. confidence intervals
  cmult <- 2
  
  newdat <- newdat |>
    mutate(
      plo = temp_sens - cmult * sqrt(pvar),
      phi = temp_sens + cmult * sqrt(pvar),
      tlo = temp_sens - cmult * sqrt(tvar1), # this is including the random effects
      thi = temp_sens + cmult * sqrt(tvar1)
      
    )
  
  return(newdat)
}



# but this does not include random effects
make_sens_predictions2 <- function(model) {
  
  pred <- ggpredict(
    model,
    terms = "treat_competition"
  )
  
  pred <- as.data.frame(pred)
  
  return(pred)
}


make_sens_predictions3 <- function(model) {
  
  newdat <- data.frame(treat_competition = c("with", "without"))
  
  pred <- predict(
    model,
    newdat = newdat,
    re.form = NA,
    se.fit = TRUE
  )
  
  newdat$fit <- pred$fit
  newdat$se <- pred$se.fit
  
  newdat$lo <- newdat$fit - 1.96 * newdat$se
  newdat$hi <- newdat$fit + 1.96 * newdat$se
  
  return(newdat)
}


# plot --------------------------------------------------------------------
plot_temp_sens_predictions <- function(sens_data_raw, sens_predictions) {
  
  ggplot() +
    
    # raw species values
    geom_jitter(
      data = sens_data_raw,
      aes(
        x = treat_competition,
        y = temp_sens,
        color = treat_competition
      ),
      width = 0.2,
      alpha = 0.2
    ) +
    
    # model predictions
    geom_point(
      data = sens_predictions,
      aes(
        x = treat_competition,
        y = temp_sens,
        color = treat_competition
      ),
      size = 3
    ) +
    
    geom_errorbar(
      data = sens_predictions,
      aes(
        x = treat_competition,
        ymin = plo,
        ymax = phi,
        color = treat_competition
      ),
      width = 0.12
    ) +
    
    facet_grid(region ~ stage) +
    
    scale_color_manual(values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"
    )) +
    
    labs(
      x = "Biotic interactions",
      y = expression("Temperature sensitivity (days/"*degree*"C)")
    ) +
    theme(
      legend.position = "none"
    )+
    geom_hline(yintercept=0, linetype = "dashed")
}







