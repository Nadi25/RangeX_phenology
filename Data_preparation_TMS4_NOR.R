




# Tomst loggers -----------------------------------------------------------


# read in tms data 23 -----------------------------------------------------
tms_nor <- read.csv("Data/Clean/RangeX_clean_EnvTMS4_2023_NOR.csv")



# import metadata ---------------------------------------------------------
meta_p_nor <- read.csv("Data/RangeX_clean_MetadataPlot_NOR.csv")



# combine metadata with tms data ------------------------------------------
env_nor <- tms_nor |> 
  left_join(meta_p_nor, by = "unique_plot_ID")


# Date as date ------------------------------------------------------------
env_nor <- env_nor |> 
  mutate(date_time = parse_date_time(date_time, orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d")))

# check if CHE has NAs in date_time ------------------------------------------------
date_na <- env_nor |> 
  filter(is.na(date_time))
date_na

# more than expected

# plot_ID_original as chr -------------------------------------------------
env_nor <- env_nor |> mutate(plot_ID_original = as.character(plot_ID_original))

# block_ID_original as chr -------------------------------------------------
env_nor <- env_nor |> mutate(block_ID_original  = as.character(block_ID_original ))


# block_ID as chr -------------------------------------------------
env_nor <- env_nor |> mutate(block_ID  = as.character(block_ID ))


# data availability 2023 ---------------------------------------------------
env_nor_23 <- env_nor |>
  mutate(
    date = as.Date(date_time),
    year = lubridate::year(date_time)
  ) |>
  group_by(site, unique_plot_ID, treat_warming, treat_competition) |>
  summarise(
    start_date = min(date, na.rm = TRUE),
    end_date   = max(date, na.rm = TRUE),
    n_records  = n(),
    .groups = "drop"
  ) |>
  arrange(site, start_date)
env_nor_23

env_nor_23 |>
  count(site, treat_warming, treat_competition)




env_nor_daily <- env_nor |>
  filter(!is.na(date_time)) |>
  mutate(date = as.Date(date_time)) |>
  group_by(
    region, site, added_focals,
    unique_plot_ID,
    treat_warming,
    treat_competition,
    date
  ) |>
  summarise(
    T1_mean = mean(TMS_T1, na.rm = TRUE),
    T1_min  = min(TMS_T1, na.rm = TRUE),
    T1_max  = max(TMS_T1, na.rm = TRUE),
    
    T2_mean = mean(TMS_T2, na.rm = TRUE),
    T2_min  = min(TMS_T2, na.rm = TRUE),
    T2_max  = max(TMS_T2, na.rm = TRUE),
    
    T3_mean = mean(TMS_T3, na.rm = TRUE),
    T3_min  = min(TMS_T3, na.rm = TRUE),
    T3_max  = max(TMS_T3, na.rm = TRUE),
    
    .groups = "drop"
  )



env_nor_daily_long <- env_nor_daily |>
  pivot_longer(
    cols = c(T1_mean, T2_mean, T3_mean),
    names_to = "sensor",
    values_to = "temp"
  )


ggplot(
  env_nor_daily_long,
  aes(x = date, y = temp,
      color = treat_warming)
) +
  geom_line(alpha = 0.7) +
  facet_grid(sensor ~ site)




env_plot <- env_nor_daily |>
  group_by(
    date,
    site,
    treat_warming
  ) |>
  summarise(
    temp = mean(T1_mean, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  env_plot,
  aes(date, temp,
      color = treat_warming)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~site)


env_plot <- env_nor_daily |>
  group_by(
    date,
    site,
    treat_warming
  ) |>
  summarise(
    mean_temp = mean(T1_mean),
    min_temp = mean(T1_min),
    max_temp = mean(T1_max),
    .groups = "drop"
  )

ggplot(
  env_plot,
  aes(date, mean_temp,
      color = treat_warming)
) +
  geom_ribbon(
    aes(
      ymin = min_temp,
      ymax = max_temp,
      fill = treat_warming
    ),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(~site)



