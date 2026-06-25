


# Functions for onset -----------------------------------------------------


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
  
  model <- lmerTest::lmer(
    onset ~ treatment_site_temp + treat_competition + (1 | species) + (1 | block_ID),
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





# Per species -------------------------------------------------------------
fit_onset_model_species <- function(onset_data, region_name) {
  region <- onset_data |> 
    filter(region == region_name)
  
  model <- lmerTest::lmer(
    onset ~ treatment_site_temp * treat_competition * species + (1 | block_ID),
    data = region
  )
  
  return(model)
}


make_onset_predictions_species <- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without"),
    species = c("cennig", "cyncri", "hypmac", "leuvul", "luzmul", "pimsax",
                "plalan", "sildio", "sucpra", "tripra"),
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



make_onset_predictions_species_che <- function(model) {
  
  newdata <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without"),
    species = c("brapin", "broere", "cenjac", "daucar",
                "hypper", "medlup", "silvul", "plamed", "salpra", "scacol"),
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





















