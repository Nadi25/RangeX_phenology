
# BIOMASS 3 ---------------------------------------------------------------

source("Data_preparation_traits_biomass_NOR_24.R")


# library -----------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)

library(lmerTest)
conflicts_prefer(lmerTest::lmer)


# Fit individual models per species ---------------------------------------



# log transform the traits ------------------------------------------------
analysis_data_24_log <- biomass_traits_NOR_24 |>
  select(site, block_ID, treat_warming, treat_competition, functional_group, species, 
    log_biomass,
    height_vegetative,
    height_vegetative_str,
    height_reproductive,
    height_reproductive_str,
    no_stems,
    number_leaves
  ) |>
  mutate(
    log_height_vegetative       = log1p(height_vegetative),
    log_height_vegetative_str   = log1p(height_vegetative_str),
    log_height_reproductive     = log1p(height_reproductive),
    log_height_reproductive_str = log1p(height_reproductive_str),
    log_no_stems                = log1p(no_stems),
    log_number_leaves           = log1p(number_leaves)
  )


# split data by species ---------------------------------------------------
# use analysis_data because it has only rows that have the traits
data_species <- analysis_data_24_log |> 
  group_split(species)


# function with all model options -----------------------------------------
fit_models <- function(df){
  
  # vegetative height
  m_height_v <- lmerTest::lmer(
    log_biomass ~ height_vegetative + (1|block_ID),
    data = df)
  
  m_height_v_log <- lmerTest::lmer(
    log_biomass ~ log_height_vegetative + (1|block_ID),
    data = df)
  
  # stretched vegetative height
  m_height_vs <- lmerTest::lmer(
    log_biomass ~ height_vegetative_str + (1|block_ID),
    data = df)
  
  m_height_vs_log <- lmerTest::lmer(
    log_biomass ~ log_height_vegetative_str + (1|block_ID),
    data = df)
  
  # reproductive height
  m_height_r <- lmerTest::lmer(
    log_biomass ~ height_reproductive + (1|block_ID),
    data = df)
  
  m_height_r_log <- lmerTest::lmer(
    log_biomass ~ log_height_reproductive + (1|block_ID),
    data = df)
  
  # stretched reproductive height
  m_height_rs <- lmerTest::lmer(
    log_biomass ~ height_reproductive_str + (1|block_ID),
    data = df)
  
  m_height_rs_log <- lmerTest::lmer(
    log_biomass ~ log_height_reproductive_str + (1|block_ID),
    data = df)
  
  # stems
  m_stems <- lmerTest::lmer(
    log_biomass ~ no_stems + (1|block_ID),
    data = df)
  
  m_stems_log <- lmerTest::lmer(
    log_biomass ~ log_no_stems + (1|block_ID),
    data = df)
  
  # leaves
  m_leaves <- lmerTest::lmer(
    log_biomass ~ number_leaves + (1|block_ID),
    data = df)
  
  m_leaves_log <- lmerTest::lmer(
    log_biomass ~ log_number_leaves + (1|block_ID),
    data = df)
  
  # AIC comparison
  AIC(
    m_height_v,
    m_height_v_log,
    m_height_vs,
    m_height_vs_log,
    m_height_r,
    m_height_r_log,
    m_height_rs,
    m_height_rs_log,
    m_stems,
    m_stems_log,
    m_leaves,
    m_leaves_log
  )
}


# run the omdels for all species ------------------------------------------
results <- lapply(data_species, fit_models)
results



# attach species names ----------------------------------------------------
species_names <- analysis_data_24_log |> 
  distinct(species) |> 
  pull()

names(results) <- species_names


# inspect results for all species -----------------------------------------
results[["sucpra"]]
results[["cennig"]]
results[["pimsax"]]
results[["luzmul"]]
results[["leuvul"]]
results[["tripra"]]
results[["hypmac"]]
results[["plalan"]]
results[["cyncri"]]
results[["sildio"]]



# overview of the best model per species ----------------------------------
best_models <- do.call(rbind, lapply(names(results), function(sp){
  
  x <- results[[sp]]
  best <- x[which.min(x$AIC), ]
  
  data.frame(
    species = sp,
    model   = rownames(best),
    AIC     = best$AIC
  )
}))

best_models

# species           model       AIC
# 1   sucpra     m_stems_log 213.08173
# 2   cennig m_height_rs_log 119.49785
# 3   cyncri  m_height_r_log 165.59830
# 4   pimsax     m_stems_log 235.41082
# 5   luzmul     m_stems_log 249.14759
# 6   plalan m_height_rs_log 231.24656
# 7   leuvul     m_stems_log 308.15681
# 8   tripra  m_height_r_log  15.97995
# 9   hypmac     m_stems_log 129.95020
# 10  sildio    m_leaves_log 302.45545


## old results where na.omit was used in the log transform step
# species           model       AIC
# 1   sucpra     m_stems_log 212.85085
# 2   cennig m_height_rs_log 212.55040
# 3   pimsax    m_leaves_log 261.83000
# 4   luzmul    m_leaves_log 482.49088
# 5   leuvul     m_stems_log 355.26236
# 6   tripra m_height_rs_log 226.53698
# 7   hypmac    m_leaves_log 351.22768
# 8   plalan m_height_vs_log  15.15047
# 9   cyncri     m_stems_log 130.22105
# 10  sildio    m_leaves_log 345.35538


# this list the AIC for all models per species
# 0 is the best model
results_delta <- lapply(results, function(x){
  x$deltaAIC <- x$AIC - min(x$AIC)
  x[order(x$deltaAIC), ]
})
results_delta



# test if certain trait combinations are better than single traits --------
# use the two best single trait models for this

fit_top2_models <- function(species_name, data, results_delta){
  
  # filter dataset for the species
  df <- data |> filter(species == species_name)
  
  # get the top 2 traits (lowest deltaAIC)
  top2 <- results_delta[[species_name]] |> 
    slice_head(n = 2) |> 
    rownames()  # model names, e.g., "m_stems_log", "m_leaves_log"
  
  # map model names to actual variable names in the dataset
  trait_map <- c(
    m_height_v       = "height_vegetative",
    m_height_v_log   = "log_height_vegetative",
    m_height_vs      = "height_vegetative_str",
    m_height_vs_log  = "log_height_vegetative_str",
    m_height_r       = "height_reproductive",
    m_height_r_log   = "log_height_reproductive",
    m_height_rs      = "height_reproductive_str",
    m_height_rs_log  = "log_height_reproductive_str",
    m_stems          = "no_stems",
    m_stems_log      = "log_no_stems",
    m_leaves         = "number_leaves",
    m_leaves_log     = "log_number_leaves"
  )
  
  traits <- trait_map[top2]
  
  # store models in a list
  model_list <- list()
  
  # single best trait
  m_single <- lmer(as.formula(paste("log_biomass ~", traits[1], "+ (1|block_ID)")), data=df)
  model_list[["single"]] <- m_single
  
  # additive model (if 2 traits)
  if(length(traits) == 2){
    m_add <- lmer(as.formula(paste("log_biomass ~", paste(traits, collapse = " + "), "+ (1|block_ID)")), data=df)
    m_int <- lmer(as.formula(paste("log_biomass ~", paste(traits, collapse = " * "), "+ (1|block_ID)")), data=df)
    model_list[["additive"]] <- m_add
    model_list[["interaction"]] <- m_int
  }
  
  # return list of fitted models
  return(model_list)
}


top_models_species <- lapply(species_names, function(sp){
  fit_top2_models(sp, analysis_data_24_log, results_delta)
})

names(top_models_species) <- species_names



# sucpra ------------------------------------------------------------------
# single-trait model for sucpra
summary(top_models_species[["sucpra"]][["single"]])

# additive model for sucpra
summary(top_models_species[["sucpra"]][["additive"]])

# interaction model for sucpra
summary(top_models_species[["sucpra"]][["interaction"]])

# which model is best
AIC(
  top_models_species[["sucpra"]][["single"]],
  top_models_species[["sucpra"]][["additive"]],
  top_models_species[["sucpra"]][["interaction"]]
)

# top_models_species[["sucpra"]][["single"]]       4 129.95020
# top_models_species[["sucpra"]][["additive"]]     5  82.61082 ## best
# top_models_species[["sucpra"]][["interaction"]]  6  86.29368

# final model for sucpra --------------------------------------------------
# dataset for sucpra
df_sucpra <- analysis_data_24_log |> filter(species == "sucpra")

# final model
m_sucpra <- lmerTest::lmer(log_biomass ~ log_no_stems + log_number_leaves + (1 | block_ID),
                 data = df_sucpra)
summary(m_sucpra)

# cennig ------------------------------------------------------------------
# single-trait model for cennig
summary(top_models_species[["cennig"]][["single"]])

# additive model for cennig
summary(top_models_species[["cennig"]][["additive"]])

# interaction model for cennig
summary(top_models_species[["cennig"]][["interaction"]])

# which model is best
AIC(
  top_models_species[["cennig"]][["single"]],
  top_models_species[["cennig"]][["additive"]],
  top_models_species[["cennig"]][["interaction"]]
)

# top_models_species[["cennig"]][["single"]]       4 465.5881
# top_models_species[["cennig"]][["additive"]]     5 150.8483
# top_models_species[["cennig"]][["interaction"]]  6 142.0404 ## best

# final model for cennig --------------------------------------------------
# dataset for cennig
df_cennig <- analysis_data_24_log |> filter(species == "cennig")

# final model
m_cennig <- lmerTest::lmer(log_biomass ~ log_height_reproductive_str * log_no_stems + 
                             (1 | block_ID), data = df_cennig)
summary(m_cennig)


# pimsax ------------------------------------------------------------------
# single-trait model for pimsax
summary(top_models_species[["pimsax"]][["single"]])

# additive model for pimsax
summary(top_models_species[["pimsax"]][["additive"]])

# interaction model for pimsax
summary(top_models_species[["pimsax"]][["interaction"]])

# which model is best
AIC(
  top_models_species[["pimsax"]][["single"]],
  top_models_species[["pimsax"]][["additive"]],
  top_models_species[["pimsax"]][["interaction"]]
)

# top_models_species[["pimsax"]][["single"]]       4 290.6022
# top_models_species[["pimsax"]][["additive"]]     5 133.2130  ## best
# top_models_species[["pimsax"]][["interaction"]]  6 135.6647

# final model for pimsax --------------------------------------------------
# dataset for pimsax
df_pimsax <- analysis_data_24_log |> filter(species == "pimsax")

# final model
m_pimsax <- lmerTest::lmer(log_biomass ~ log_no_stems + log_number_leaves+ 
                             (1 | block_ID), data = df_pimsax)
summary(m_pimsax)


# luzmul ------------------------------------------------------------------
df_luzmul <- analysis_data_24_log |> filter(species == "luzmul")

# single-trait model
summary(top_models_species[["luzmul"]][["single"]])

# additive model
summary(top_models_species[["luzmul"]][["additive"]])

# interaction model
summary(top_models_species[["luzmul"]][["interaction"]])

# compare models
AIC(
  top_models_species[["luzmul"]][["single"]],
  top_models_species[["luzmul"]][["additive"]],
  top_models_species[["luzmul"]][["interaction"]]
)

# top_models_species[["luzmul"]][["single"]]       4 249.1476
# top_models_species[["luzmul"]][["additive"]]     5 144.3556  ## best
# top_models_species[["luzmul"]][["interaction"]]  6 145.2263

# final model for luzmul --------------------------------------------------
m_luzmul <- lmerTest::lmer(
  log_biomass ~ log_no_stems + log_number_leaves + (1 | block_ID),
  data = df_luzmul
)
summary(m_luzmul)

# leuvul ------------------------------------------------------------------
df_leuvul <- analysis_data_24_log |> filter(species == "leuvul")

summary(top_models_species[["leuvul"]][["single"]])
summary(top_models_species[["leuvul"]][["additive"]])
summary(top_models_species[["leuvul"]][["interaction"]])

AIC(
  top_models_species[["leuvul"]][["single"]],
  top_models_species[["leuvul"]][["additive"]],
  top_models_species[["leuvul"]][["interaction"]]
)

# top_models_species[["leuvul"]][["single"]]       4 235.41082
# top_models_species[["leuvul"]][["additive"]]     5  71.46198
# top_models_species[["leuvul"]][["interaction"]]  6  69.69727  ## best


# final model for leuvul --------------------------------------------------
m_leuvul <- lmerTest::lmer(
  log_biomass ~ log_no_stems * log_number_leaves + (1 | block_ID),
  data = df_leuvul
)
summary(m_leuvul)



# tripra ------------------------------------------------------------------
df_tripra <- analysis_data_24_log |> filter(species == "tripra")

summary(top_models_species[["tripra"]][["single"]])
summary(top_models_species[["tripra"]][["additive"]])
summary(top_models_species[["tripra"]][["interaction"]])

AIC(
  top_models_species[["tripra"]][["single"]],
  top_models_species[["tripra"]][["additive"]],
  top_models_species[["tripra"]][["interaction"]]
)

# top_models_species[["tripra"]][["single"]]       4 412.6296
# top_models_species[["tripra"]][["additive"]]     5 353.0691  ## best
# top_models_species[["tripra"]][["interaction"]]  6 354.6496

m_tripra2 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive + log_height_reproductive_str + (1 | block_ID),
  data = df_tripra
)
summary(m_tripra2)


# final model for tripra --------------------------------------------------
m_tripra <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str + log_no_stems + (1 | block_ID),
  data = df_tripra
)
summary(m_tripra)

AIC(m_tripra, m_tripra2)

# rep height str and not str are the two best
# but they are closely correlated so we take the next which is more independent 
## log no stems


# hypmac ------------------------------------------------------------------
df_hypmac <- analysis_data_24_log |> filter(species == "hypmac")

summary(top_models_species[["hypmac"]][["single"]])
summary(top_models_species[["hypmac"]][["additive"]])
summary(top_models_species[["hypmac"]][["interaction"]])

AIC(
  top_models_species[["hypmac"]][["single"]],
  top_models_species[["hypmac"]][["additive"]],
  top_models_species[["hypmac"]][["interaction"]]
)

# top_models_species[["hypmac"]][["single"]]       4 492.8488 ## best
# top_models_species[["hypmac"]][["additive"]]     5 500.2511
# top_models_species[["hypmac"]][["interaction"]]  6 502.7848

# but that is the same trait with and without log, so we take the next
# leaves log

m_hypmac <- lmerTest::lmer(
  log_biomass ~ log_no_stems + log_number_leaves + (1 | block_ID),
  data = df_hypmac
)
summary(m_hypmac)


m_hypmac2 <- lmerTest::lmer(
  log_biomass ~ log_no_stems * log_number_leaves + (1 | block_ID),
  data = df_hypmac
)
summary(m_hypmac2)


m_hypmac3 <- lmerTest::lmer(
  log_biomass ~ log_no_stems + (1 | block_ID),
  data = df_hypmac
)
summary(m_hypmac3)

m_hypmac4 <- lmerTest::lmer(
  log_biomass ~ log_number_leaves + (1 | block_ID),
  data = df_hypmac
)
summary(m_hypmac4)

AIC(m_hypmac, m_hypmac2, m_hypmac3, m_hypmac4)

# the best would be only log_number_leaves
# but if we want multi-trait model, we need two 

# final model for hypmac --------------------------------------------------
m_hypmac <- lmerTest::lmer(
  log_biomass ~ log_no_stems + log_number_leaves + (1 | block_ID),
  data = df_hypmac
)
summary(m_hypmac)



# plalan ------------------------------------------------------------------
df_plalan <- analysis_data_24_log |> filter(species == "plalan")

summary(top_models_species[["plalan"]][["single"]])
summary(top_models_species[["plalan"]][["additive"]]) # just with and without log
summary(top_models_species[["plalan"]][["interaction"]])

AIC(
  top_models_species[["plalan"]][["single"]],
  top_models_species[["plalan"]][["additive"]],
  top_models_species[["plalan"]][["interaction"]]
)

# top_models_species[["plalan"]][["single"]]       4 384.6810  ## best
# top_models_species[["plalan"]][["additive"]]     5 391.6417
# top_models_species[["plalan"]][["interaction"]]  6 390.0806

# best is height_rep_str
# second best is height_rep 

m_plalan <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str + log_height_reproductive + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan)

m_plalan2 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str * log_height_reproductive + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan2)

m_plalan3 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan3)

m_plalan4 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str + number_leaves + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan4)

m_plalan5 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str * number_leaves + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan5)

m_plalan6 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str * log_number_leaves + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan6)

AIC(m_plalan, m_plalan2, m_plalan3, m_plalan4, m_plalan5, m_plalan6)


# final model for plalan --------------------------------------------------
m_plalan <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive_str * log_number_leaves + (1 | block_ID),
  data = df_plalan
)
summary(m_plalan)



# cyncri ------------------------------------------------------------------
df_cyncri <- analysis_data_24_log |> filter(species == "cyncri")

summary(top_models_species[["cyncri"]][["single"]])
summary(top_models_species[["cyncri"]][["additive"]]) # same trait, with and without log
summary(top_models_species[["cyncri"]][["interaction"]])

AIC(
  top_models_species[["cyncri"]][["single"]],
  top_models_species[["cyncri"]][["additive"]],
  top_models_species[["cyncri"]][["interaction"]]
)

# top_models_species[["cyncri"]][["single"]]       4 127.3636  ## best
# top_models_species[["cyncri"]][["additive"]]     5 130.2517
# top_models_species[["cyncri"]][["interaction"]]  6 135.9990


# try with second best trait, height rep str log

m_cyncri3 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive + log_height_reproductive_str + (1 | block_ID),
  data = df_cyncri
)
summary(m_cyncri3)

m_cyncri2 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive * log_height_reproductive_str + (1 | block_ID),
  data = df_cyncri
)
summary(m_cyncri2)

m_cyncri <- lmerTest::lmer(
  log_biomass ~ log_no_stems * log_number_leaves + (1 | block_ID),
  data = df_cyncri
)
summary(m_cyncri)

m_cyncri4 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive + (1 | block_ID),
  data = df_cyncri
)
summary(m_cyncri4)


m_cyncri5 <- lmerTest::lmer(
  log_biomass ~ log_height_reproductive * log_number_leaves + (1 | block_ID),
  data = df_cyncri
)


AIC(m_cyncri3, m_cyncri2, m_cyncri, m_cyncri4, m_cyncri5)




# final model for cyncri --------------------------------------------------
m_cyncri <- lmerTest::lmer(
  log_biomass ~ log_no_stems * log_number_leaves + (1 | block_ID),
  data = df_cyncri
)
summary(m_cyncri)
# use the interactive model 
# first thought that the additive model is causing a split in the biomass

# sildio ------------------------------------------------------------------
df_sildio <- analysis_data_24_log |> filter(species == "sildio")

summary(top_models_species[["sildio"]][["single"]])
summary(top_models_species[["sildio"]][["additive"]])
summary(top_models_species[["sildio"]][["interaction"]])

AIC(
  top_models_species[["sildio"]][["single"]],
  top_models_species[["sildio"]][["additive"]],
  top_models_species[["sildio"]][["interaction"]]
)

# top_models_species[["sildio"]][["single"]]       4 169.55222
# top_models_species[["sildio"]][["additive"]]     4  13.28415  ## best
# top_models_species[["sildio"]][["interaction"]]  4  13.28415



cor(df_sildio[, c("log_number_leaves", "log_height_reproductive_str")], use = "pairwise.complete.obs")

m_sildio3 <- lmerTest::lmer(
  log_biomass ~ log_number_leaves + (1 | block_ID),
  data = df_sildio
)


m_sildio2 <- lmerTest::lmer(
  log_biomass ~ log_number_leaves + log_height_reproductive_str + (1 | block_ID),
  data = df_sildio
)


m_sildio <- lmerTest::lmer(
  log_biomass ~ log_number_leaves * log_height_reproductive_str + (1 | block_ID),
  data = df_sildio
)

m_sildio4 <- lmerTest::lmer(
  log_biomass ~ log_number_leaves + log_no_stems + (1 | block_ID),
  data = df_sildio
)

AIC(m_sildio, m_sildio2, m_sildio3, m_sildio4)

# final model for sildio --------------------------------------------------
m_sildio <- lmerTest::lmer(
  log_biomass ~ log_number_leaves + log_no_stems + (1 | block_ID),
  data = df_sildio
)
summary(m_sildio)

# before it looked better with log_number_leaves * log_height_reproductive_str
# in the plot but 
# this doesn't work now




