


# Functions for onset -----------------------------------------------------
# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
conflict_prefer_all("lmerTest", quiet = TRUE)

# function to get first onset per individual ----------------------------------------------
get_onset <- function(data, stage_name, onset_type) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, treatment_site_temp, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID
    ) |>
    summarise(
      onset = min(.data[[onset_type]]),
      .groups = "drop"
    ) |> 
    # remove groups where budding never occurred
    filter(is.finite(onset))
}



# function of model per stage and region ---------------------------------------------
fit_onset_model <- function(onset_data, region_name) {
  region <- onset_data |> 
    filter(region == region_name)
  
  model <- lmer(
    onset ~ treatment_site_temp * treat_competition + (1 | species) + (1 | block_ID),
    data = region
  )
  
  return(model)
}

fit_onset_model_species <- function(onset_data, region_name) {
  region <- onset_data |> 
    filter(region == region_name)
  
  model <- lmer(
    onset ~ treatment_site_temp * treat_competition + 
      (treatment_site_temp | species) + (treat_competition | species) + (1 | block_ID),
    data = region
  )
  
  return(model)
}


# function for predictions ------------------------------------------------
make_onset_predictions <- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without")
  ) |>
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


make_onset_predictions_gdd <- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without")
  ) |>
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

make_onset_predictions_species3 <- function(model) {
  
  pred <- ggpredict(model,
                    terms = c("treatment_site_temp",
                              "treat_competition",
                              "species"),
                    type = "random")
  
  return(pred)
}


# Per species -------------------------------------------------------------
filter_species <- function(onset_data, species_name) {
  
  onset_data_species <- onset_data |>
    filter(species == species_name)
  
  return(onset_data_species)
}

fit_species_model <- function(onset_data_species) {
  
  lmer(onset ~ treatment_site_temp * treat_competition +
      (1 | block_ID),
    data = onset_data_species)
  
}

make_onset_predictions_species<- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without"),
    #species = levels(model.frame(model)$species),
    block_ID = NA
  ) |>
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

make_onset_predictions_species2 <- function(model) {
  
  pred <- ggpredict(model,
            terms = c("treatment_site_temp",
              "treat_competition"))
  
  return(pred)
}





make_onset_predictions_species_manual <- function(model) {
  
  # new data
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without"),
    species = levels(model.frame(model)$species)
  )
  
  # species-specific predictions
  newdata$fit <- predict(
    model,
    newdata = newdata,
    re.form = ~(1 | species)
  )
  
  # fixed-effect model matrix
  mm <- model.matrix(
    delete.response(terms(model)),
    newdata
  )
  
  # fixed-effect variance
  pvar <- diag(
    mm %*% tcrossprod(vcov(model), mm)
  )
  
  # random-effect variances
  var_species <- as.numeric(VarCorr(model)$species)
  var_block   <- as.numeric(VarCorr(model)$block_ID)
  
  # total variance
  tvar <- pvar + var_species + var_block + sigma(model)^2
  
  cmult <- 1.96
  
  newdata <- newdata |>
    mutate(
      se_fit = sqrt(pvar),
      
      lower_fixed = fit - cmult * sqrt(pvar),
      upper_fixed = fit + cmult * sqrt(pvar),
      
      lower = fit - cmult * sqrt(tvar),
      upper = fit + cmult * sqrt(tvar)
    ) |>
    tibble::as_tibble()
  
  return(newdata)
}




# test region differences -------------------------------------------------

# function of model per stage and region ---------------------------------------------
fit_onset_model_region <- function(onset_data) {
  
  model <- lmer(
    onset ~ region * treatment_site_temp * treat_competition + (1 | species) + (1 | block_ID),
    data = onset_data
  )
  
  return(model)
}


