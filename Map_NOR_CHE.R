

# Map of sites ------------------------------------------------------------

# Load required libraries
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)

theme_set(theme_bw(base_size = 5))

world <- ne_countries(scale = "medium", returnclass = "sf")

locations_nor_che <- data.frame(
  site = c("Norway", "Switzerland"),
  lat = c(60.5403, 46.888),
  lon = c(6.5676, 9.489)
)

locations_sf <- st_as_sf(
  locations_nor_che,
  coords = c("lon", "lat"),
  crs = 4326
)

map <- ggplot() +
  geom_sf(data = world,
          fill = "grey85",
          color = "white",
          linewidth = 0.2) +
  geom_sf(data = locations_sf,
          aes(color = site),
          size = 1) +
  # geom_label_repel(
  #   data = locations_nor_che,
  #   aes(lon, lat, label = site, color = site),
  #   size = 1.25, label.size = 0.01
  # ) +
  geom_text_repel(
    data = locations_nor_che,
    aes(lon, lat, label = site),
    size = 2.2,
    color = "aquamarine4",
    box.padding = 0.1,
    point.padding = 0.1,
    segment.color = "grey60",
    min.segment.length = 0
  )+

  coord_sf(
    xlim = c(2, 14),
    ylim = c(45, 62),
    expand = FALSE
  )+
  scale_color_manual(
    values = c(
      "Norway" = "aquamarine4",
      "Switzerland" = "aquamarine4"
    )
  ) +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  )
map

# ggsave(filename = "Output/Map/Map_NOR_CHE2.png",
#        plot = map,
#   width = 5,
#   height = 5,
#   units = "cm",
#   dpi = 600)

# ggsave(filename = "Output/Temperature/Figures_for_manuscript/DOY_GDD5_TMS4_NOR_CHE_12h.png", 
#       plot = doy_gdd3, width = 18, height = 12, units = "cm", dpi = 600)


