

# 02_temp_sens TMS4 ------------------------------------------------------------

# Temperature sensitivity analysis -------------------------------------------------


# RangeX phenology effect of transplantation on temperature sensitivity ------------

## Data used: RangeX_clean_phenology_2023_NOR.csv
##            
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
##            
## Date:      07.07.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on temperature sensitivity

## 
# using mean temperatures per stage
# season start - last budding, 
# season start - last flowering

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lme4)
library(performance)
library(see)
library(emmeans)
library(ggeffects)
library(sjPlot)
library(lattice)
library(gt)

theme_set(theme_bw(base_size = 20))


# source script with functions --------------------------------------------
source("Functions_temperature_sensitivity.R")


# source script with final tms data ---------------------------------------
# this one sources the NOR and CHE preparation scripts
source("DOY_GDD_TMS4_NOR_CHE.R")


# add region and combine the datasets -------------------------------------

# > season_start_nor_tms
# # A tibble: 6 × 4
# site  treat_warming treat_competition season_start
# <chr> <chr>         <chr>             <date>      
# 1 hi    ambi          bare              2023-04-12  
# 2 hi    ambi          vege              2023-04-12  
# 3 hi    warm          bare              2023-04-11  
# 4 hi    warm          vege              2023-04-11  
# 5 lo    ambi          bare              2023-04-11  
# 6 lo    ambi          vege              2023-04-11 
season_start_nor_tms$region <- "Norway"

# season_start_che_tms
# # A tibble: 6 × 4
# site  treat_warming treat_competition season_start
# <chr> <chr>         <chr>             <date>      
# 1 hi    ambi          bare              2022-04-16  
# 2 hi    ambi          vege              2022-04-14  
# 3 hi    warm          bare              2022-04-16  
# 4 hi    warm          vege              2022-04-13  
# 5 lo    ambi          bare              2022-03-27  
# 6 lo    ambi          vege              2022-03-27  
season_start_che_tms$region <- "Switzerland"



# combine nor and che tms data --------------------------------------------
# and change names of treat_competition to match
# this is the season start dataset for
season_start_all <- bind_rows(
  season_start_nor_tms,
  season_start_che_tms) |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))


tms_final_nor_che <- tms_final_nor_che |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))




# season_start_all <- season_start_all |> 
#   mutate(
#     temp_start = case_when(
#       region == "Switzerland" ~ as.Date("2022-04-28"),
#       #region == "Switzerland" & site == "lo" ~ season_start,
#       TRUE ~ season_start))






# source clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)


# rename infructescence stage of NOR to FlowWithrd ------------------------
# this will be the fruiting stage
# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology <- phenology |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "No_Infructescences" = "No_FloWithrd"))

# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology4 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional 



phenology4 <- phenology4 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))

phenology4$treatment_site_temp <- paste(phenology4$site, phenology4$treat_warming, sep = "_")




# Get first onset per individual ------------------------------------------
onset_ind <- phenology4 |>
  filter(value > 0) |>
  group_by(
    region,
    site,
    treat_competition,
    treat_warming,
    phenology_stage,
    unique_plant_ID
  ) |>
  summarise(
    onset = min(date_measurement),
    .groups = "drop"
  )


# get mean onset per treatment --------------------------------------------
onset_treat <- onset_ind |>
  group_by(
    region, site,
    treat_competition, treat_warming,
    phenology_stage
  ) |>
  summarise(
    start = min(onset),
    end = max(onset),
    .groups = "drop"
  )


# Get all onsets together -------------------------------------------------
stage_windows <- onset_treat |>
  pivot_wider(
    names_from = phenology_stage,
    values_from = c(start, end)
  )
stage_windows

stage_windows <- stage_windows |>
  left_join(
    season_start_all,
    by = c(
      "region",
      "site",
      "treat_warming",
      "treat_competition"
    )
  )

stage_windows <- stage_windows |>
  mutate(
    bud_start    = season_start,
    bud_end      = end_No_Buds,
    
    flower_start = season_start,
    flower_end   = end_No_FloOpen,
    
    fruit_start  = season_start,
    fruit_end    = end_No_FloWithrd,
    
    seed_start   = season_start,
    seed_end     = end_No_Seeds
  )
stage_windows

#write.csv(stage_windows, "Data/Stage_windows_NOR_CHE.csv")

# Function to calculate mean temp per stage -------------------------------
calc_mean_temp <- function(data, windows, start_col, end_col, label) {
  
  data |>
    inner_join(
      windows,
      by = c("region", "site", "treat_warming", "treat_competition")
    ) |>
    filter(
      date >= .data[[start_col]],
      date <= .data[[end_col]]
    ) |>
    group_by(region, site, treat_warming, treat_competition) |>
    summarise(
      mean_temp = mean(temp_mean, na.rm = TRUE),
      stage = label,
      .groups = "drop"
    )
}


# apply the function per stage -------------------------------------------
# bud
temp_bud <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "bud_start",
  "bud_end",
  "budding"
)

temp_bud$treatment_site_temp <- paste(temp_bud$site, temp_bud$treat_warming, sep = "_")


# flower
temp_flower <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "flower_start",
  "flower_end",
  "flowering"
)

temp_flower$treatment_site_temp <- paste(temp_flower$site, temp_flower$treat_warming, sep = "_")

# fruit
temp_fruit <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "fruit_start",
  "fruit_end",
  "fruiting"
)

temp_fruit$treatment_site_temp <- paste(temp_fruit$site, temp_fruit$treat_warming, sep = "_")

# seed
temp_seed <- calc_mean_temp(
  tms_final_nor_che,
  stage_windows,
  "seed_start",
  "seed_end",
  "seeds"
)

temp_seed$treatment_site_temp <- paste(temp_seed$site, temp_seed$treat_warming, sep = "_")


# mean temp per stage and treatment - combine all stages ------------------------------------------------------
temp_all <- bind_rows(temp_bud, temp_flower, temp_fruit, temp_seed)
temp_all

# add combined treatment
temp_all$treatment_site_temp <- paste(temp_all$site, temp_all$treat_warming, sep = "_")

# plot
ggplot(temp_all,
       aes(x = stage,
           y = mean_temp,
           color = treat_warming,
           group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(region + site ~ treat_competition) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature during stage",
    color = "Warming"
  )


temp_all_ambi <- temp_all |> 
  filter(treat_warming == "ambi")

ggplot(temp_all_ambi,
       aes(x = stage,
           y = mean_temp,
           color = site,
           group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(region ~ treat_competition) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature during stage",
    color = "Warming"
  )

temp_all_hi <- temp_all |> 
  filter(site == "hi")

ggplot(temp_all_hi,
       aes(x = stage,
           y = mean_temp,
           color = treat_warming,
           group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(region ~ treat_competition) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature during stage",
    color = "Warming"
  )

# hi site: bud: ambi is warmer than warmed


# Plot mean temp per stage ------------------------------------------------
theme_set(theme_bw(base_size = 14))
temp <- ggplot(temp_all,
       aes(x = stage,
           y = mean_temp,
           color = treatment_site_temp,
           shape = treatment_site_temp,
           group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(region ~ treat_competition,
             labeller = labeller(treat_competition = c("with" = "With biotic interactions",
                                   "without" = "Without biotic interactions"))) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature °C",
    shape = "Treatment site warming",
    title = "Mean temperature per stage (season start - stage end)")+
  scale_color_manual(values = c("hi_ambi" = "turquoise4", "hi_warm" = "darkred", "lo_ambi"= "grey34"))+
  scale_shape_manual(values = c("lo_ambi" = 16, "hi_ambi" = 17,  "hi_warm" = 2))+
  guides(color = "none")+
  theme(legend.position = "bottom")
temp

# ggsave(filename = "Output/Sensitivity/Temperature_mean_per_stage_NOR_CHE.png", 
#       plot = temp,
#       width = 8, height = 6, units = "in")


define_colors <- c(
  "hi_ambi.with" = "turquoise4",
  "hi_ambi.without" = lighten("turquoise4", 0.4),
  
  "hi_warm.with" = "darkred",
  "hi_warm.without" = lighten("darkred", 0.4),
  
  "lo_ambi.with" = "grey34",
  "lo_ambi.without" = lighten("grey34", 0.4)
)


temp <- ggplot(
  temp_all,
  aes(
    x = stage,
    y = mean_temp,
    color = interaction(treatment_site_temp, treat_competition),
    shape = treatment_site_temp,
    group = interaction(site, treat_warming, treat_competition))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_grid(treat_competition
     ~ region,
    labeller = labeller(
      treat_competition = c(
        "with" = "With biotic interactions",
        "without" = "Without biotic interactions"))) +
  labs(
    x = "Phenological stage",
    y = "Mean temperature °C",
    shape = "Treatment site warming",
    title = "Mean temperature per stage (season start - stage end) 12h") +
  scale_color_manual(values = define_colors) +
  scale_shape_manual(
    values = c(
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 2)) +
  theme(legend.position = "bottom")+
  guides(color = "none")
temp

# ggsave(filename = "Output/Sensitivity/Temperature_mean_per_stage_NOR_CHE_12h.png", 
#       plot = temp,
#       width = 9, height = 7, units = "in")






temp_all |>
  select(
    stage,
    region,
    treatment_site_temp,
    treat_competition,
    mean_temp
  ) |>
  arrange(stage, region, treatment_site_temp, treat_competition) |>
  pivot_wider(names_from = region,
              values_from = mean_temp) |> 
  gt(groupname_col = "stage") |> 
  cols_label(treatment_site_temp = "Site temp treatment",
             treat_competition = "Biotic interactions")




pivot_wider(
  names_from = region,
  values_from = c(estimate_ci, letters)
) |>
  gt(rowname_col = NULL,
     groupname_col = "stage") |>
  tab_spanner(
    label = "Norway",
    columns = c(
      estimate_ci_Norway,
      letters_Norway
    )
  ) |>
  tab_spanner(
    label = "Switzerland",
    columns = c(
      estimate_ci_Switzerland,
      letters_Switzerland
    )
  )





temp_tab <- temp_all |>
  select(
    stage,
    region,
    site,
    treatment_site_temp,
    treat_competition,
    mean_temp
  ) |>
  arrange(stage, region, treatment_site_temp, treat_competition) |>
  gt(
    groupname_col = "stage"
  ) |>
  cols_label(
    region = "Region",
    site = "Site",
    treatment_site_temp = "Site temp treatment",
    treat_competition = "Biotic interactions",
    mean_temp = "Mean temperature (°C)"
  ) |>
  fmt_number(
    columns = mean_temp,
    decimals = 1
  ) |>
  cols_align(
    align = "center",
    columns = c(
      site,
      treatment_site_temp,
      treat_competition,
      mean_temp
    )
  ) |>
  tab_header(
    title = md("**Mean temperatures by treatment**"),
    subtitle = "Average temperatures for each region, stage and treatment combination"
  ) |>
  opt_table_font(
    font = list(
      google_font("Source Sans Pro"),
      default_fonts()
    )
  ) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(18),
    row_group.font.weight = "bold",
    row_group.background.color = "grey95",
    data_row.padding = px(5)
  )

temp_tab


# calculate mean onset per species and individual for all stages ----------------
onset_bud    <- get_mean_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_mean_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_mean_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_mean_onset(phenology2, "No_Seeds", "jday")



# Join phenology with temperature data ------------------------------------
onset_bud_temp <- onset_bud |>
  left_join(temp_bud,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_flower_temp <- onset_flower |>
  left_join(temp_flower,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_fruit_temp <- onset_fruit |>
  left_join(temp_fruit,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_seed_temp <- onset_seed |>
  left_join(temp_seed,
            by = c("region", "treatment_site_temp", "treat_competition"))


# combine from all stages -------------------------------------------
onset_all_temp <- bind_rows(
  onset_bud_temp   |> mutate(stage = "Budding"),
  onset_flower_temp|> mutate(stage = "Flowering"),
  onset_fruit_temp |> mutate(stage = "Fruiting"),
  onset_seed_temp  |> mutate(stage = "Seeds")
)
onset_all_temp

onset_all_temp <- onset_all_temp |> 
  rename(Tmean = mean_temp)


# quick control plot ------------------------------------------------------
# plot the raw sens data 
ggplot(onset_all_temp,
       aes(x = stage,
           y = onset,
           color = treat_competition)) +
  
  # individual species
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.4
    ),
    alpha = 0.3
  ) +
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.4),
    size = 4
  ) +
  
  labs(
    x = "Phenological stage",
    y = expression("Onset"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)




# Low vs high site --------------------------------------------------------


# Filter correct onset dataset low vs hi --------------------------------------------------
# per stage
# only ambi

# NOR
d_bud_lh_nor <- filter_data_ambi(onset_all_temp, "Norway", "Budding", "ambi")
d_flower_lh_nor <- filter_data_ambi(onset_all_temp, "Norway", "Flowering", "ambi")
d_fruit_lh_nor <- filter_data_ambi(onset_all_temp, "Norway", "Fruiting", "ambi")
d_seed_lh_nor <- filter_data_ambi(onset_all_temp, "Norway", "Seeds", "ambi")


# CHE
d_bud_lh_che <- filter_data_ambi(onset_all_temp, "Switzerland", "Budding", "ambi")
d_flower_lh_che <- filter_data_ambi(onset_all_temp, "Switzerland", "Flowering", "ambi")
d_fruit_lh_che <- filter_data_ambi(onset_all_temp, "Switzerland", "Fruiting", "ambi")
d_seed_lh_che <- filter_data_ambi(onset_all_temp, "Switzerland", "Seeds", "ambi")



# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_lh_nor    <- fit_model_sens(d_bud_lh_nor)
m_sens_flower_gs_lh_nor <- fit_model_sens(d_flower_lh_nor)
m_sens_fruit_gs_lh_nor  <- fit_model_sens(d_fruit_lh_nor)
m_sens_seed_gs_lh_nor   <- fit_model_sens(d_seed_lh_nor)


summary(m_sens_bud_gs_lh_nor)
summary(m_sens_flower_gs_lh_nor)
summary(m_sens_fruit_gs_lh_nor)
summary(m_sens_seed_gs_lh_nor)

summary(m_sens_bud_gs_lh_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_lh_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_lh_che    <- fit_model_sens(d_bud_lh_che)
m_sens_flower_gs_lh_che <- fit_model_sens(d_flower_lh_che)
m_sens_fruit_gs_lh_che  <- fit_model_sens(d_fruit_lh_che)
m_sens_seed_gs_lh_che   <- fit_model_sens(d_seed_lh_che)


summary(m_sens_bud_gs_lh_che)
summary(m_sens_flower_gs_lh_che)
summary(m_sens_fruit_gs_lh_che)
summary(m_sens_seed_gs_lh_che)

summary(m_sens_bud_gs_lh_che)$coefficients["Tmean", ]



# filter temp data --------------------------------------------------------
# this is to get the newdataframe for the predictions
temp_lh_nor_b <- temp_bud |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor_b

temp_lh_nor_f <- temp_flower |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor_f

temp_lh_nor_fr <- temp_fruit |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor_fr

temp_lh_nor_s <- temp_seed |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor_s


temp_lh_che_b <- temp_bud |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che_b

temp_lh_che_f <- temp_flower |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che_f

temp_lh_che_fr <- temp_fruit |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che_fr

temp_lh_che_s <- temp_seed |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che_s




# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_lh_nor <- get_temp_sens_coef(m_sens_bud_gs_lh_nor)
ts_flower_lh_nor <- get_temp_sens_coef(m_sens_flower_gs_lh_nor)
ts_fruit_lh_nor <- get_temp_sens_coef(m_sens_fruit_gs_lh_nor)
ts_seed_lh_nor <- get_temp_sens_coef(m_sens_seed_gs_lh_nor)


# CHE ---------------------------------------------------------------------

ts_bud_lh_che <- get_temp_sens_coef(m_sens_bud_gs_lh_che)
ts_flower_lh_che <- get_temp_sens_coef(m_sens_flower_gs_lh_che)
ts_fruit_lh_che <- get_temp_sens_coef(m_sens_fruit_gs_lh_che)
ts_seed_lh_che <- get_temp_sens_coef(m_sens_seed_gs_lh_che)



# combine sens from all stages -------------------------------------------
sens_all_lh <- bind_rows(
  ts_bud_lh_nor$slopes    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_lh_nor$slopes |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_lh_nor$slopes  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_lh_nor$slopes   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_lh_che$slopes    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_lh_che$slopes |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_lh_che$slopes  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_lh_che$slopes   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_lh

sens_all_lh <- sens_all_lh |>
  mutate(
    slope_stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""))

sens_all_lh

sig_all_lh <- bind_rows(
  ts_bud_lh_nor$pairs    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_lh_nor$pairs |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_lh_nor$pairs  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_lh_nor$pairs   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_lh_che$pairs    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_lh_che$pairs |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_lh_che$pairs  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_lh_che$pairs   |> mutate(stage = "Seeds", region = "Switzerland")
)
sig_all_lh

sig_all_lh <- sig_all_lh |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"))
sig_all_lh


# Plot sensitivity --------------------------------------------------------
ts <- ggplot(sens_all_lh, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity low vs high ambient")+
  facet_grid(region ~ stage, scales = "free")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_low_high_NOR_CHE.png", 
#       plot = ts,
#        width = 15, height = 10, units = "in")




# High ambi vs warm -------------------------------------------------------


# Filter correct onset dataset hi ambi vs warm --------------------------------------------------
# per stage
# only high

# NOR
d_bud_aw_nor <- filter_data_hi(onset_all_temp, "Norway", "Budding", "hi")
d_flower_aw_nor <- filter_data_hi(onset_all_temp, "Norway", "Flowering", "hi")
d_fruit_aw_nor <- filter_data_hi(onset_all_temp, "Norway", "Fruiting", "hi")
d_seed_aw_nor <- filter_data_hi(onset_all_temp, "Norway", "Seeds", "hi")


# CHE
d_bud_aw_che <- filter_data_hi(onset_all_temp, "Switzerland", "Budding", "hi")
d_flower_aw_che <- filter_data_hi(onset_all_temp, "Switzerland", "Flowering", "hi")
d_fruit_aw_che <- filter_data_hi(onset_all_temp, "Switzerland", "Fruiting", "hi")
d_seed_aw_che <- filter_data_hi(onset_all_temp, "Switzerland", "Seeds", "hi")



# sensitivity models ambi warm ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_aw_nor    <- fit_model_sens(d_bud_aw_nor)
m_sens_flower_gs_aw_nor <- fit_model_sens(d_flower_aw_nor)
m_sens_fruit_gs_aw_nor  <- fit_model_sens(d_fruit_aw_nor)
m_sens_seed_gs_aw_nor   <- fit_model_sens(d_seed_aw_nor)


summary(m_sens_bud_gs_aw_nor)
summary(m_sens_flower_gs_aw_nor)
summary(m_sens_fruit_gs_aw_nor)
summary(m_sens_seed_gs_aw_nor)

summary(m_sens_bud_gs_aw_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_aw_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_aw_che    <- fit_model_sens(d_bud_aw_che)
m_sens_flower_gs_aw_che <- fit_model_sens(d_flower_aw_che)
m_sens_fruit_gs_aw_che  <- fit_model_sens(d_fruit_aw_che)
m_sens_seed_gs_aw_che   <- fit_model_sens(d_seed_aw_che)


summary(m_sens_bud_gs_aw_che)
summary(m_sens_flower_gs_aw_che)
summary(m_sens_fruit_gs_aw_che)
summary(m_sens_seed_gs_aw_che)

summary(m_sens_bud_gs_aw_che)$coefficients["Tmean", ]




# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_aw_nor <- get_temp_sens_coef(m_sens_bud_gs_aw_nor)
ts_flower_aw_nor <- get_temp_sens_coef(m_sens_flower_gs_aw_nor)
ts_fruit_aw_nor <- get_temp_sens_coef(m_sens_fruit_gs_aw_nor)
ts_seed_aw_nor <- get_temp_sens_coef(m_sens_seed_gs_aw_nor)


# CHE ---------------------------------------------------------------------

ts_bud_aw_che <- get_temp_sens_coef(m_sens_bud_gs_aw_che)
ts_flower_aw_che <- get_temp_sens_coef(m_sens_flower_gs_aw_che)
ts_fruit_aw_che <- get_temp_sens_coef(m_sens_fruit_gs_aw_che)
ts_seed_aw_che <- get_temp_sens_coef(m_sens_seed_gs_aw_che)



# combine sens from all stages -------------------------------------------
sens_all_aw <- bind_rows(
  ts_bud_aw_nor$slopes    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_aw_nor$slopes |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_aw_nor$slopes  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_aw_nor$slopes   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_aw_che$slopes    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_aw_che$slopes |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_aw_che$slopes  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_aw_che$slopes   |> mutate(stage = "Seeds", region = "Switzerland")
)

sens_all_aw <- sens_all_aw |>
  mutate(
    slope_stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""))
sens_all_aw

sig_all_aw <- bind_rows(
  ts_bud_aw_nor$pairs    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_aw_nor$pairs |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_aw_nor$pairs  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_aw_nor$pairs   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_aw_che$pairs    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_aw_che$pairs |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_aw_che$pairs  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_aw_che$pairs   |> mutate(stage = "Seeds", region = "Switzerland")
)

sig_all_aw <- sig_all_aw |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"))
sig_all_aw




# Plot sensitivity --------------------------------------------------------
ts_aw <- ggplot(sens_all_aw, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity high ambient vs warmed")+
  facet_grid(region ~ stage, scales = "free")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_aw

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_ambient_warmed_NOR_CHE.png", 
#       plot = ts_aw,
#        width = 15, height = 10, units = "in")




# joined figure for hi vs lo and ambi vs warm -----------------------------

sens_all_lh$type <- "low-high"
sens_all_aw$type <- "warmed-ambient"

sens_combined <- bind_rows(sens_all_lh, sens_all_aw)

pd <- position_dodge(width = 0.4)

ts_hl_aw <- ggplot(
  sens_combined,
  aes(
    x = type,
    y = Tmean.trend,
    color = treat_competition,
    shape = type
  )
) +
  geom_point(
    size = 3.5,
    position = pd
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.1,
    position = pd
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(region ~ stage, scales = "free") +
  labs(
    y = "Temperature sensitivity (days / °C)",
    x = "Temperature shift",
    title = "Temperature sensitivity 12h daily mean per stage",
    color = "Biotic interactions") +
  theme(legend.position = "right") +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"))+
  guides(shape = "none")+
  theme(legend.position = "bottom")
ts_hl_aw


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_lo_ambi_warm_NOR_CHE.png", 
#       plot = ts_hl_aw,
#       width = 15, height = 10, units = "in")


# add significance stars --------------------------------------------------

sig_all_lh$type <- "low-high"
sig_all_aw$type <- "warmed-ambient"

sig_combined <- bind_rows(sig_all_lh, sig_all_aw)
sig_combined


brackets <- sens_combined |>
  group_by(region, stage, type) |>
  summarise(
    y_bracket = max(upper.CL) + 4,
    .groups = "drop"
  ) |>
  left_join(
    sig_combined |>
      select(region, stage, type, stars),
    by = c("region", "stage", "type")
  ) |>
  mutate(
    x = ifelse(type == "low-high", 1, 2),
    xmin = x - 0.2,
    xmax = x + 0.2
  )
brackets

ts_hl_aw_sig <- ts_hl_aw +
  geom_text(
    aes(
      label = slope_stars,
      y = upper.CL + 0.4
    ),
    position = pd,
    show.legend = FALSE
  ) +
  
  # horizontal line
  geom_segment(
    data = brackets,
    aes(
      x = xmin,
      xend = xmax,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE
  ) +
  
  # left tick
  geom_segment(
    data = brackets,
    aes(
      x = xmin,
      xend = xmin,
      y = y_bracket,
      yend = y_bracket - 0.3
    ),
    inherit.aes = FALSE
  ) +
  
  # right tick
  geom_segment(
    data = brackets,
    aes(
      x = xmax,
      xend = xmax,
      y = y_bracket,
      yend = y_bracket - 0.3
    ),
    inherit.aes = FALSE
  ) +
  
  # significance text
  geom_text(
    data = brackets,
    aes(
      x = x,
      y = y_bracket + 4,
      label = stars
    ),
    inherit.aes = FALSE,
    size = 5
  )
ts_hl_aw_sig


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_lo_ambi_warm_signif_NOR_CHE_12h.png", 
#       plot = ts_hl_aw_sig,
#       width = 15, height = 10, units = "in")





# make significance table -------------------------------------------------

tab_slopes <- sens_combined |> 
  select(region, stage, type, treat_competition, Tmean.trend, 
         lower.CL, upper.CL, p.value, slope_stars) |> 
  rename(slope = Tmean.trend,
         p_slope = p.value)
tab_slopes

tab_comp <- sig_combined |>
  select(region, stage, type, estimate, p.value, stars) |>
  rename(slope_difference = estimate,
    p_difference = p.value,
    difference_stars = stars)
tab_comp


tab_all <- tab_slopes |>
  left_join(
    tab_comp,
    by = c("region", "stage", "type"))
tab_all




gt(tab_all) |>
  fmt_number(
    columns = c(
      slope,
      lower.CL,
      upper.CL,
      slope_difference
    ),
    decimals = 2
  ) |>
  fmt_scientific(
    columns = c(
      p_slope,
      p_difference
    ),
    decimals = 2
  ) |>
  cols_label(
    region = "Region",
    stage = "Stage",
    type = "Temperature shift",
    treat_competition = "Biotic interactions",
    slope = "Temp. sensitivity",
    lower.CL = "Lower CI",
    upper.CL = "Upper CI",
    p_slope = "P(slope)",
    slope_stars = "Slope sig.",
    p_difference = "P(with vs without)",
    difference_stars = "Difference sig.")







# GDD ---------------------------------------------------------------------

source("GDD_Transplantation_onset_predictions_NOR_CHE.R")


# calculate first onset per species and plot for all stages ----------------
onset_bud_gdd    <- get_onset(phenology_gdd_nor_che, "No_Buds", "GDD_cum")
onset_flower_gdd <- get_onset(phenology_gdd_nor_che, "No_FloOpen", "GDD_cum")
onset_fruit_gdd  <- get_onset(phenology_gdd_nor_che, "No_FloWithrd", "GDD_cum")
onset_seed_gdd   <- get_onset(phenology_gdd_nor_che, "No_Seeds", "GDD_cum")



# Join phenology with temperature data ------------------------------------
onset_bud_temp_gdd <- onset_bud_gdd |>
  left_join(temp_bud,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_flower_temp_gdd <- onset_flower_gdd |>
  left_join(temp_flower,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_fruit_temp_gdd <- onset_fruit_gdd |>
  left_join(temp_fruit,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_seed_temp_gdd <- onset_seed_gdd |>
  left_join(temp_seed,
            by = c("region", "treatment_site_temp", "treat_competition"))


# combine from all stages -------------------------------------------
onset_all_gdd <- bind_rows(
  onset_bud_temp_gdd   |> mutate(stage = "Budding"),
  onset_flower_temp_gdd|> mutate(stage = "Flowering"),
  onset_fruit_temp_gdd |> mutate(stage = "Fruiting"),
  onset_seed_temp_gdd |> mutate(stage = "Seeds")
)
onset_all_gdd

onset_all_gdd <- onset_all_gdd |> 
  rename(Tmean = mean_temp)


# Filter correct onset dataset low vs hi --------------------------------------------------
# per stage
# only ambi

# NOR
d_bud_gdd_lh_nor <- filter_data_ambi(onset_all_gdd, "Norway", "Budding", "ambi")
d_flower_gdd_lh_nor <- filter_data_ambi(onset_all_gdd, "Norway", "Flowering", "ambi")
d_fruit_gdd_lh_nor <- filter_data_ambi(onset_all_gdd, "Norway", "Fruiting", "ambi")
d_seed_gdd_lh_nor <- filter_data_ambi(onset_all_gdd, "Norway", "Seeds", "ambi")


# CHE
d_bud_gdd_lh_che <- filter_data_ambi(onset_all_gdd, "Switzerland", "Budding", "ambi")
d_flower_gdd_lh_che <- filter_data_ambi(onset_all_gdd, "Switzerland", "Flowering", "ambi")
d_fruit_gdd_lh_che <- filter_data_ambi(onset_all_gdd, "Switzerland", "Fruiting", "ambi")
d_seed_gdd_lh_che <- filter_data_ambi(onset_all_gdd, "Switzerland", "Seeds", "ambi")




# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gdd_lh_nor    <- fit_model_sens(d_bud_gdd_lh_nor)
m_sens_flower_gdd_lh_nor <- fit_model_sens(d_flower_gdd_lh_nor)
m_sens_fruit_gdd_lh_nor  <- fit_model_sens(d_fruit_gdd_lh_nor)
m_sens_seed_gdd_lh_nor   <- fit_model_sens(d_seed_gdd_lh_nor)


summary(m_sens_bud_gdd_lh_nor)
summary(m_sens_flower_gdd_lh_nor)
summary(m_sens_fruit_gdd_lh_nor)
summary(m_sens_seed_gdd_lh_nor)

summary(m_sens_bud_gdd_lh_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gdd_lh_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gdd_lh_che    <- fit_model_sens(d_bud_gdd_lh_che)
m_sens_flower_gdd_lh_che <- fit_model_sens(d_flower_gdd_lh_che)
m_sens_fruit_gdd_lh_che  <- fit_model_sens(d_fruit_gdd_lh_che)
m_sens_seed_gdd_lh_che   <- fit_model_sens(d_seed_gdd_lh_che)


summary(m_sens_bud_gdd_lh_che)
summary(m_sens_flower_gdd_lh_che)
summary(m_sens_fruit_gdd_lh_che)
summary(m_sens_seed_gdd_lh_che)

summary(m_sens_bud_gdd_lh_che)$coefficients["Tmean", ]





# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_gdd_lh_nor <- get_temp_sens_coef(m_sens_bud_gdd_lh_nor)
ts_flower_gdd_lh_nor <- get_temp_sens_coef(m_sens_flower_gdd_lh_nor)
ts_fruit_gdd_lh_nor <- get_temp_sens_coef(m_sens_fruit_gdd_lh_nor)
ts_seed_gdd_lh_nor <- get_temp_sens_coef(m_sens_seed_gdd_lh_nor)


# CHE ---------------------------------------------------------------------

ts_bud_gdd_lh_che <- get_temp_sens_coef(m_sens_bud_gdd_lh_che)
ts_flower_gdd_lh_che <- get_temp_sens_coef(m_sens_flower_gdd_lh_che)
ts_fruit_gdd_lh_che <- get_temp_sens_coef(m_sens_fruit_gdd_lh_che)
ts_seed_gdd_lh_che <- get_temp_sens_coef(m_sens_seed_gdd_lh_che)



# combine sens from all stages -------------------------------------------
sens_all_gdd_lh <- bind_rows(
  ts_bud_gdd_lh_nor$slopes    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_lh_nor$slopes |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_lh_nor$slopes  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_lh_nor$slopes   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_lh_che$slopes    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_lh_che$slopes |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_lh_che$slopes  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_lh_che$slopes   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_gdd_lh


sens_all_gdd_lh <- sens_all_gdd_lh |>
  mutate(
    slope_stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""))
sens_all_gdd_lh

sig_all_gdd_lh <- bind_rows(
  ts_bud_gdd_lh_nor$pairs    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_lh_nor$pairs |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_lh_nor$pairs  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_lh_nor$pairs   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_lh_che$pairs    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_lh_che$pairs |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_lh_che$pairs  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_lh_che$pairs   |> mutate(stage = "Seeds", region = "Switzerland")
)

sig_all_gdd_lh <- sig_all_gdd_lh |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"))
sig_all_gdd_lh


# Plot sensitivity --------------------------------------------------------
ts_gdd_lh <- ggplot(sens_all_gdd_lh, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (GDD2 / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity low vs high ambient")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_gdd_lh

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_GDD_low_high_NOR_CHE.png", 
#       plot = ts_gdd_lh,
#        width = 15, height = 10, units = "in")




# GDD hi ambi vs warm -----------------------------------------------------


# Filter correct onset dataset hi ambi vs warm --------------------------------------------------
# per stage
# only hi

# NOR
d_bud_gdd_aw_nor <- filter_data_hi(onset_all_gdd, "Norway", "Budding", "hi")
d_flower_gdd_aw_nor <- filter_data_hi(onset_all_gdd, "Norway", "Flowering", "hi")
d_fruit_gdd_aw_nor <- filter_data_hi(onset_all_gdd, "Norway", "Fruiting", "hi")
d_seed_gdd_aw_nor <- filter_data_hi(onset_all_gdd, "Norway", "Seeds", "hi")


# CHE
d_bud_gdd_aw_che <- filter_data_hi(onset_all_gdd, "Switzerland", "Budding", "hi")
d_flower_gdd_aw_che <- filter_data_hi(onset_all_gdd, "Switzerland", "Flowering", "hi")
d_fruit_gdd_aw_che <- filter_data_hi(onset_all_gdd, "Switzerland", "Fruiting", "hi")
d_seed_gdd_aw_che <- filter_data_hi(onset_all_gdd, "Switzerland", "Seeds", "hi")




# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gdd_aw_nor    <- fit_model_sens(d_bud_gdd_aw_nor)
m_sens_flower_gdd_aw_nor <- fit_model_sens(d_flower_gdd_aw_nor)
m_sens_fruit_gdd_aw_nor  <- fit_model_sens(d_fruit_gdd_aw_nor)
m_sens_seed_gdd_aw_nor   <- fit_model_sens(d_seed_gdd_aw_nor)


summary(m_sens_bud_gdd_aw_nor)
summary(m_sens_flower_gdd_aw_nor)
summary(m_sens_fruit_gdd_aw_nor)
summary(m_sens_seed_gdd_aw_nor)

summary(m_sens_bud_gdd_aw_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gdd_aw_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gdd_aw_che    <- fit_model_sens(d_bud_gdd_aw_che)
m_sens_flower_gdd_aw_che <- fit_model_sens(d_flower_gdd_aw_che)
m_sens_fruit_gdd_aw_che  <- fit_model_sens(d_fruit_gdd_aw_che)
m_sens_seed_gdd_aw_che   <- fit_model_sens(d_seed_gdd_aw_che)


summary(m_sens_bud_gdd_aw_che)
summary(m_sens_flower_gdd_aw_che)
summary(m_sens_fruit_gdd_aw_che)
summary(m_sens_seed_gdd_aw_che)

summary(m_sens_bud_gdd_aw_che)$coefficients["Tmean", ]



# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_gdd_aw_nor <- get_temp_sens_coef(m_sens_bud_gdd_aw_nor)
ts_flower_gdd_aw_nor <- get_temp_sens_coef(m_sens_flower_gdd_aw_nor)
ts_fruit_gdd_aw_nor <- get_temp_sens_coef(m_sens_fruit_gdd_aw_nor)
ts_seed_gdd_aw_nor <- get_temp_sens_coef(m_sens_seed_gdd_aw_nor)


# CHE ---------------------------------------------------------------------

ts_bud_gdd_aw_che <- get_temp_sens_coef(m_sens_bud_gdd_aw_che)
ts_flower_gdd_aw_che <- get_temp_sens_coef(m_sens_flower_gdd_aw_che)
ts_fruit_gdd_aw_che <- get_temp_sens_coef(m_sens_fruit_gdd_aw_che)
ts_seed_gdd_aw_che <- get_temp_sens_coef(m_sens_seed_gdd_aw_che)



# combine sens from all stages -------------------------------------------
sens_all_gdd_aw <- bind_rows(
  ts_bud_gdd_aw_nor$slopes    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_aw_nor$slopes |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_aw_nor$slopes  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_aw_nor$slopes   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_aw_che$slopes    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_aw_che$slopes |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_aw_che$slopes  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_aw_che$slopes   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_gdd_aw


sens_all_gdd_aw <- sens_all_gdd_aw |>
  mutate(
    slope_stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""))
sens_all_gdd_aw

sig_all_gdd_aw <- bind_rows(
  ts_bud_gdd_aw_nor$pairs    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_aw_nor$pairs |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_aw_nor$pairs  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_aw_nor$pairs   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_aw_che$pairs    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_aw_che$pairs |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_aw_che$pairs  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_aw_che$pairs   |> mutate(stage = "Seeds", region = "Switzerland")
)

sig_all_gdd_aw <- sig_all_gdd_aw |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"))
sig_all_gdd_aw




# Plot sensitivity --------------------------------------------------------
ts_gdd_aw <- ggplot(sens_all_gdd_aw, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (GDD5 / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity high ambient vs warmed")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_gdd_aw

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_GDD_ambi_warm_NOR_CHE.png", 
#       plot = ts_gdd_aw,
#       width = 15, height = 10, units = "in")






# joined figure for GDD hi vs lo and ambi vs warm -----------------------------

sens_all_gdd_lh$type <- "low-high"
sens_all_gdd_aw$type <- "warmed-ambient"

sens_combined_gdd <- bind_rows(sens_all_gdd_lh, sens_all_gdd_aw)

pd <- position_dodge(width = 0.4)

ts_hl_aw_gdd <- ggplot(
  sens_combined_gdd,
  aes(
    x = type,
    y = Tmean.trend,
    color = treat_competition,
    shape = type
  )
) +
  geom_point(
    size = 3.5,
    position = pd
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.1,
    position = pd
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(region ~ stage, scales = "free") +
  labs(
    y = "Temperature sensitivity (GDD5 / °C)",
    x = "Temperature shift",
    title = "Temperature sensitivity 12h daily mean per stage 5°C base",
    color = "Biotic interactions") +
  theme(legend.position = "right") +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"))+
  guides(shape = "none")+
  theme(legend.position = "bottom")
ts_hl_aw_gdd


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_lo_ambi_warm_GDD_NOR_CHE.png", 
#       plot = ts_hl_aw_gdd,
#       width = 15, height = 10, units = "in")




# add significance stars --------------------------------------------------

sig_all_gdd_lh$type <- "low-high"
sig_all_gdd_aw$type <- "warmed-ambient"

sig_combined_gdd <- bind_rows(sig_all_gdd_lh, sig_all_gdd_aw)
sig_combined_gdd


brackets_gdd <- sens_combined_gdd |>
  group_by(region, stage, type) |>
  summarise(
    y_bracket = max(upper.CL) + 50,
    .groups = "drop"
  ) |>
  left_join(
    sig_combined_gdd |>
      select(region, stage, type, stars),
    by = c("region", "stage", "type")
  ) |>
  mutate(
    x = ifelse(type == "low-high", 1, 2),
    xmin = x - 0.2,
    xmax = x + 0.2
  )
brackets_gdd

ts_hl_aw_sig_gdd <- ts_hl_aw_gdd +
  geom_text(
    aes(
      label = slope_stars,
      y = upper.CL + 0.4
    ),
    position = pd,
    show.legend = FALSE
  ) +
  
  # horizontal line
  geom_segment(
    data = brackets_gdd,
    aes(
      x = xmin,
      xend = xmax,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE
  ) +
  
  # left tick
  geom_segment(
    data = brackets_gdd,
    aes(
      x = xmin,
      xend = xmin,
      y = y_bracket,
      yend = y_bracket - 0.3
    ),
    inherit.aes = FALSE
  ) +
  
  # right tick
  geom_segment(
    data = brackets_gdd,
    aes(
      x = xmax,
      xend = xmax,
      y = y_bracket,
      yend = y_bracket - 0.3
    ),
    inherit.aes = FALSE
  ) +
  
  # significance text
  geom_text(
    data = brackets_gdd,
    aes(
      x = x,
      y = y_bracket + 30,
      label = stars
    ),
    inherit.aes = FALSE,
    size = 5
  )
ts_hl_aw_sig_gdd


# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_lo_ambi_warm_GDD5_NOR_CHE_12h.png", 
#       plot = ts_hl_aw_sig_gdd,
#       width = 15, height = 10, units = "in")








