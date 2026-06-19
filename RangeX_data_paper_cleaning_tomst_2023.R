
# Climate data TOMST loggers NOR 2023 --------------------------------------------

## Data used: Data/Data_tomst_loggers/tomst_2023/,
##            tomst_plot_codes_2023.csv,
##            RangeX_metadata_plot_NOR.csv
##            Sunrise_sundown_Voss_2023.csv
## Date:      03.03.2025
## Author:    Nadine Arzt
## Purpose:   Clean TOMST logger data 2023

# load library ------------------------------------------------------------
library(tidyverse)
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(janitor)

# default theme
theme_set(theme_bw())

# comments --------------------------------------------------------------
# https://tomst.com/web/en/systems/tms/software/
# temp1 is in soil, temp2 at 0cm and temp3 20cm above ground

# 2023: high out: "08.06.2023" and "20.06.2023" - "23.10.2023"
# 2023: low out: "12.05.2023" and "19.06.2023" - "24.10.2023"
# so take the later one each?

# deleted tomst 94201723 and 94217314 because of impossible values
# should I delete 94217320? it has a drop in early August: 2023-08-19 11:15:00

# calculate max, min and daily amplitude
# calculate rolling average with rollmean()

# calculate soil moisture: https://github.com/audhalbritter/Three-D/blob/master/R/functions/soilmoisture_correction.R


# source script with functions --------------------------------------------
source("RangeX_data_paper_functions.R")

# import data 2023 --------------------------------------------------------
# List all files in the 'Data_tomst_loggers' folder that start with 'data'
tomst_23 <- list.files(path = "Data/Raw/tomst_2023/", pattern = "^data_\\d+.*\\.csv$", full.names = TRUE)

# test to see structure of files
test_file <- read_delim(tomst_23[1], delim = ";", skip = 1)
head(test_file)
# has , 

# define colnames as header
column_names <- c("number", "date_time", "Column1", "TMS_T1", "TMS_T2", "TMS_T3", "Soilmoisture_raw", "Column6", "Column7")

# define coltypes to have the values correct later
column_types <- cols(
  number = col_double(),
  date_time = col_character(),   # Read as character first, convert to datetime later
  Column1 = col_double(),
  TMS_T1 = col_character(),  # Read as character to handle commas, convert later
  TMS_T2 = col_character(),
  TMS_T3 = col_character(),
  Soilmoisture_raw = col_double(),
  Column6 = col_double(),
  Column7 = col_double()
)

# get one dataframe with data from all files using the list of files (tomst_23) with a loop
tomst_data_23 <- map(tomst_23, read_tomst_file_23) |> 
  list_rbind()
head(tomst_data_23)

# get plot codes 23 ------------------------------------------------------
plot_codes_23 <- read.csv2("Data/Raw/tomst_plot_codes_2023.csv", skip = 1)
plot_codes_23

# split dataset into low and high 
# high
plot_high <- plot_codes_23 |> 
  select(block:comment)
# low
plot_low <- plot_codes_23 |> 
  select(block.1:comment.1)

plot_low <- plot_low |> 
  filter(if_any(everything(), ~ !is.na(.) & . != "")) |> 
  rename(block = block.1,
        treat = treat.1,
        tomst = tomst.1,
        date_out = date_out.1,
        comment = comment.1)

# combine again under each other
plot_codes_clean <- bind_rows(hi = plot_high, lo = plot_low, .id = "site")
head(plot_codes_clean)
str(plot_codes_clean) # tomst = chr 

# tomst into character
plot_codes_clean <- plot_codes_clean |> 
  mutate(tomst = as.character(tomst))

# combine tomst data with plot labels -------------------------------------
tomst_23_raw <- left_join(plot_codes_clean, tomst_data_23, by = "tomst")
head(tomst_23_raw)


# Add treat_warming and treat_competition columns based on treat ----------
# why does it not work with .default ~ NA
tomst_23_raw <- tomst_23_raw |> 
  mutate(
    treat_warming = case_when(
      site == "hi" & treat %in% c("A", "C", "E") ~ "warm",
      site == "hi" & treat %in% c("B", "D", "F") ~ "ambi",
      site == "lo" & treat %in% c("A", "B") ~ "ambi",  
      TRUE ~ NA  # Assign NA to unexpected cases
    ),
    treat_competition = case_when(
      site == "hi" & treat %in% c("A", "B") ~ "vege",
      site == "hi" & treat %in% c("C", "D") ~ "control",
      site == "hi" & treat %in% c("E", "F") ~ "bare",
      site == "lo" & treat == "A" ~ "vege",  
      site == "lo" & treat == "B" ~ "bare",
      TRUE ~ NA
    )
  )

# delete empty columns -----------------------------------------------------
tomst_23_raw <- tomst_23_raw |> 
  select(where(~ !all(is.na(.))))

# delete tomst 94201723 because of impossible values --------------------------------------------------
tomst_23_raw <- tomst_23_raw |> 
  filter(tomst != "94201723")

# delete 94217314 ---------------------------------------------------------
tomst_23_raw <- tomst_23_raw |> 
  filter(tomst != "94217314")

# calculate soil moisture -----------------------------------------------------------
# function to calculate soil moisture from raw values in function script
# RangeX_data_paper_functions.R
# apply the soil moisture function here
tomst_23_raw <- tomst_23_raw |> 
  mutate(VWC = calc_soil_moist(rawsoilmoist = Soilmoisture_raw, 
                                    soil_temp = TMS_T1, 
                                    soilclass ="silt_loam"))

# rename raw soil moisture into TMS_moist ---------------------------------
# these are the raw moisture count values
tomst_23_raw <- tomst_23_raw |> 
  rename(TMS_moist = Soilmoisture_raw)


# combined treatment column -----------------------------------------------
tomst_23_raw <- tomst_23_raw |> 
  mutate(treat_combined = paste(site, treat_warming, treat_competition, sep = "_"))


unique(tomst_23_raw$date_out)



tomst_23_raw2 <- tomst_23_raw |>
  mutate(date_out = as.Date(date_out, format = "%d.%m.%Y"))

logger_time <- tomst_23_raw2 |>
  group_by(tomst, site, treat_warming, treat_competition) |>
  summarise(
    deployed = min(date_out, na.rm = TRUE),
    .groups = "drop"
  )
logger_time


tomst_joined <- tomst_23_raw2 |>
  left_join(logger_time, by = c("tomst", "site", "treat_warming", "treat_competition"))


tomst_filtered <- tomst_joined |>
  filter(
    date_time >= (deployed + 1),
    date_time <= as.Date("2023-10-23") # were collected on 24.10
  )






# import metadata -------------------------------------------------------
metadata <- read.csv("Data/RangeX_clean_MetadataPlot_NOR.csv")


# fix col names --------------------------------------------------------
names(metadata)
names(tomst_filtered)

tomst_filtered <- tomst_filtered |> 
  rename("block_ID_original" = "block",
         "plot_ID_original" = "treat")

# to match plot_ID_original
tomst_filtered <- tomst_filtered |> 
  mutate(plot_ID_original = case_when(
    plot_ID_original == "A" ~ "a",
    plot_ID_original == "B" ~ "b",
    plot_ID_original == "C" ~ "c",
    plot_ID_original == "D" ~ "d",
    plot_ID_original == "E" ~ "e",
    plot_ID_original == "F" ~ "f"
  ))


# merge metadata with tomst_23_raw_filtered -------------------------------
metadata <- metadata |> 
  mutate(across(c(block_ID_original, plot_ID_original), as.character))

tomst_filtered <- tomst_filtered |> 
  mutate(across(c(block_ID_original, plot_ID_original), as.character))

tomst_23_clean<- left_join(metadata, tomst_filtered, 
                         by = c( "site", "block_ID_original",
                                 "plot_ID_original",
                                 "treat_warming", "treat_competition"))


# select only columns needed for clean data on OSF ------------------------
rx_tomst_23_clean <- tomst_23_clean |> 
  select(site,treat_competition, treat_warming, tomst,
         unique_plot_ID, date_time, 
         TMS_T1, TMS_T2, TMS_T3, 
         TMS_moist, VWC)


# save clean data ---------------------------------------------------------
#write.csv(rx_tomst_23_clean, "Data/Data_tomst_loggers/CleanEnvTMS4/RangeX_clean_EnvTMS4_2023_NOR.csv", row.names = FALSE)

# when does which logger start
logger_start_dates <- tomst_23_clean |>
  group_by(site,treat_competition, treat_warming, tomst) |>
  summarise(
    first_record = min(date_time, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(site, treat_competition, treat_warming, first_record)

logger_start_dates














