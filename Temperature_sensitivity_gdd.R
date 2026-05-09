



source("Temperature_sensitivity.R")





# sensitivity models gdd nor ------------------------------------------------------------------

# bud
m_sens_bud_gdd <- lmerTest::lmer(temp_sens ~ treat_competition + 
                               (1|species),
                             data = sens_bud_gdd_nor)
summary(m_sens_bud_gdd)
anova(m_sens_bud_gdd)

model_performance(m_sens_bud_gdd)

# flower
m_sens_flower_gdd <- lmerTest::lmer(temp_sens ~ treat_competition + 
                                  (1|species),
                                data = sens_flower_gdd_nor)
summary(m_sens_flower_gdd)
anova(m_sens_flower_gdd)

model_performance(m_sens_flower_gdd)

# fruit
m_sens_fruit_gdd <- lmerTest::lmer(temp_sens ~ treat_competition + 
                                 (1|species),
                               data = sens_fruit_gdd_nor)
summary(m_sens_fruit_gdd)
anova(m_sens_fruit_gdd)

model_performance(m_sens_fruit_gdd)

# seeds
m_sens_seed_gdd <- lmerTest::lmer(temp_sens ~ treat_competition + 
                                (1|species),
                              data = sens_seed_gdd_nor)
summary(m_sens_seed_gdd)
anova(m_sens_seed_gdd)

model_performance(m_sens_seed_gdd)




# make predictions for each stage -----------------------------------------

# bud ---------------------------------------------------------------------

# create new matrix for predicted data ------------------------------------
newdat_sens_bud_gdd <- expand.grid(
  treat_competition = c("with", "without"),
  temp_sens = 0
)
newdat_sens_bud_gdd



# predict  ----------------------------------------------------------------
newdat_sens_bud_gdd$temp_sens <- predict(
  m_sens_bud_gdd,
  newdata = newdat_sens_bud_gdd,
  re.form = NA 
)
newdat_sens_bud_gdd


# make model matrix -------------------------------------------------------
mm_sens_bud_gdd <- model.matrix(terms(m_sens_bud_gdd), newdat_sens_bud_gdd)
mm_sens_bud_gdd


pvar1_sens_bud_gdd <- diag(mm_sens_bud_gdd %*% tcrossprod(vcov(m_sens_bud_gdd), mm_sens_bud_gdd))
# tvar1 <- pvar1+VarCorr(m_onset_bud_cooling)$species[1]  ## must be adapted for more complex models


# 2. EXTRACT RANDOM EFFECT VARIANCES
# VarCorr returns variance-covariance matrices for each group
tvar1_sens_bud_gdd <- as.numeric(VarCorr(m_sens_bud_gdd)$species)
#var_block_sens_bud <- as.numeric(VarCorr(m_sens_bud)$block_ID)

# 3. CALCULATE TOTAL VARIANCE
# This is fixed-effect uncertainty + variance between sites + variance between blocks
# If you want the interval for a NEW observation (individual point), add sigma(fm1)^2 as well
#tvar1_sens_bud <- pvar1_sens_bud + var_species_sens_bud + var_block_sens_bud + sigma(m_sens_bud)^2

# 4. CALCULATE INTERVALS
cmult <- 2 # is roughly twice standard error


newdat_sens_bud_gdd <- data.frame(
  newdat_sens_bud_gdd
  , plo = newdat_sens_bud_gdd$temp_sens-cmult*sqrt(pvar1_sens_bud_gdd) # fixed effects only, confidence interval
  , phi = newdat_sens_bud_gdd$temp_sens+cmult*sqrt(pvar1_sens_bud_gdd)
  , tlo = newdat_sens_bud_gdd$temp_sens-cmult*sqrt(tvar1_sens_bud_gdd) # takes fixed and random effects, prediction intervall
  , thi = newdat_sens_bud_gdd$temp_sens+cmult*sqrt(tvar1_sens_bud_gdd)
)
newdat_sens_bud_gdd


# plot predicted onset ----------------------------------------------------
pd <- position_dodge(width = 0.4) 

b<- ggplot(newdat_sens_bud_gdd, aes(x=treat_competition, y= temp_sens, 
                                color=treat_competition)) +
  geom_point(position = pd)+
  geom_errorbar(aes(ymin= plo, ymax= phi), width=.2,
                position = pd)+
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  ))+
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (budding)",
    title = "Effect of transplantation on budding onset nor"
  )
print(b)


# flower ---------------------------------------------------------------------
newdat_sens_flower_gdd <- expand.grid(
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_flower_gdd$temp_sens <- predict(
  m_sens_flower_gdd,
  newdata = newdat_sens_flower_gdd,
  re.form = NA 
)
newdat_sens_flower_gdd

mm_sens_flower_gdd <- model.matrix(terms(m_sens_flower_gdd), newdat_sens_flower_gdd)

pvar1_sens_flower_gdd <- diag(mm_sens_flower_gdd %*% tcrossprod(vcov(m_sens_flower_gdd), mm_sens_flower_gdd))

tvar1_sens_flower_gdd <- as.numeric(VarCorr(m_sens_flower_gdd)$species)
#var_block_sens_flower <- as.numeric(VarCorr(m_sens_flower)$block_ID)

#tvar1_sens_flower <- pvar1_sens_flower + var_species_sens_flower + var_block_sens_flower + sigma(m_sens_flower)^2

cmult <- 2

newdat_sens_flower_gdd <- data.frame(
  newdat_sens_flower_gdd,
  plo = newdat_sens_flower_gdd$temp_sens - cmult * sqrt(pvar1_sens_flower_gdd),
  phi = newdat_sens_flower_gdd$temp_sens + cmult * sqrt(pvar1_sens_flower_gdd),
  tlo = newdat_sens_flower_gdd$temp_sens - cmult * sqrt(tvar1_sens_flower_gdd),
  thi = newdat_sens_flower_gdd$temp_sens + cmult * sqrt(tvar1_sens_flower_gdd)
)
newdat_sens_flower_gdd

pd <- position_dodge(width = 0.4)

f <- ggplot(newdat_sens_flower_gdd, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (flowering)",
    title = "Effect of transplantation on flowering onset nor"
  )
print(f)

# fruit ---------------------------------------------------------------------
newdat_sens_fruit_gdd <- expand.grid(
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_fruit_gdd$temp_sens <- predict(
  m_sens_fruit_gdd,
  newdata = newdat_sens_fruit_gdd,
  re.form = NA 
)

mm_sens_fruit_gdd <- model.matrix(terms(m_sens_fruit_gdd), newdat_sens_fruit_gdd)

pvar1_sens_fruit_gdd <- diag(mm_sens_fruit_gdd %*% tcrossprod(vcov(m_sens_fruit_gdd), mm_sens_fruit_gdd))

tvar1_sens_fruit_gdd <- as.numeric(VarCorr(m_sens_fruit_gdd)$species)
#var_block_sens_fruit <- as.numeric(VarCorr(m_sens_fruit)$block_ID)

#tvar1_sens_fruit <- pvar1_sens_fruit + var_species_sens_fruit + var_block_sens_fruit + sigma(m_sens_fruit)^2

cmult <- 2

newdat_sens_fruit_gdd <- data.frame(
  newdat_sens_fruit_gdd,
  plo = newdat_sens_fruit_gdd$temp_sens - cmult * sqrt(pvar1_sens_fruit_gdd),
  phi = newdat_sens_fruit_gdd$temp_sens + cmult * sqrt(pvar1_sens_fruit_gdd),
  tlo = newdat_sens_fruit_gdd$temp_sens - cmult * sqrt(tvar1_sens_fruit_gdd),
  thi = newdat_sens_fruit_gdd$temp_sens + cmult * sqrt(tvar1_sens_fruit_gdd)
)

pd <- position_dodge(width = 0.4)

fr <- ggplot(newdat_sens_fruit_gdd, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (fruiting)",
    title = "Effect of transplantation on fruiting onset nor"
  )
print(fr)


# seed ---------------------------------------------------------------------
newdat_sens_seed_gdd <- expand.grid(
  treat_competition = c("with", "without"),
  temp_sens = 0
)

newdat_sens_seed_gdd$temp_sens <- predict(
  m_sens_seed_gdd,
  newdata = newdat_sens_seed_gdd,
  re.form = NA 
)
newdat_sens_seed_gdd

mm_sens_seed_gdd <- model.matrix(terms(m_sens_seed_gdd), newdat_sens_seed_gdd)

pvar1_sens_seed_gdd <- diag(mm_sens_seed_gdd %*% tcrossprod(vcov(m_sens_seed_gdd), mm_sens_seed_gdd))

tvar1_sens_seed_gdd <- as.numeric(VarCorr(m_sens_seed_gdd)$species)
#var_block_sens_seed <- as.numeric(VarCorr(m_sens_seed)$block_ID)

#tvar1_sens_seed <- pvar1_sens_seed + var_species_sens_seed + var_block_sens_seed + sigma(m_sens_seed)^2

cmult <- 2

newdat_sens_seed_gdd<- data.frame(
  newdat_sens_seed_gdd,
  plo = newdat_sens_seed_gdd$temp_sens - cmult * sqrt(pvar1_sens_seed_gdd),
  phi = newdat_sens_seed_gdd$temp_sens + cmult * sqrt(pvar1_sens_seed_gdd),
  tlo = newdat_sens_seed_gdd$temp_sens - cmult * sqrt(tvar1_sens_seed_gdd),
  thi = newdat_sens_seed_gdd$temp_sens + cmult * sqrt(tvar1_sens_seed_gdd)
)
newdat_sens_seed_gdd


pd <- position_dodge(width = 0.4)

s <- ggplot(newdat_sens_seed_gdd, aes(x = treat_competition, y = temp_sens, color = treat_competition)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = plo, ymax = phi), width = .2, position = pd) +
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (seed set)",
    title = "Effect of transplantation on seed onset nor"
  )
print(s)



# combine all stages into one plot ----------------------------------------

stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636",
  "Seeds" = "grey50"
)

plot_df_sens_bud_gdd  <- newdat_sens_bud_gdd  |> mutate(stage = "Budding")

plot_df_sens_flower_gdd  <- newdat_sens_flower_gdd  |> mutate(stage = "Flowering")

plot_df_sens_fruit_gdd  <- newdat_sens_fruit_gdd  |> mutate(stage = "Fruiting")

plot_df_sens_seed_gdd  <- newdat_sens_seed_gdd  |> mutate(stage = "Seeds")

plot_df_sens_all_gdd  <- bind_rows(
  plot_df_sens_bud_gdd ,
  plot_df_sens_flower_gdd ,
  plot_df_sens_fruit_gdd ,
  plot_df_sens_seed_gdd 
)
plot_df_sens_all_gdd 




# plot --------------------------------------------------------------------
pd <- position_dodge(width = 0.6) 

b_f_fr_s <- ggplot(plot_df_sens_all_gdd, aes(
  x = treat_competition,
  y = temp_sens,
  color = stage
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = plo, ymax = phi),
                width = .2,
                position = pd) +
  
  scale_color_manual(values = stage_colors) +
  labs(
    x = "Biotic interactions",
    y = "Predicted temperature sensitivity (GDD per °C)",
    title = "Effect of transplantation on onset nor",
    color = "Phenological stage"
  )
print(b_f_fr_s)



# ggsave(filename = "Output/Onset/Temperature_sensitivity_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#        plot = b_f_fr_s2,
#        width = 15, height = 10, units = "in")























