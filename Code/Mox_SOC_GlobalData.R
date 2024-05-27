## Global SOC - Alox ##
## 2024-05-23 ##
## Sophie von Fromm ##

library(tidyverse)
library(raster)
library(sf)
library(tmap)

#### Load merged dataset and add global data products

global_data <- read_csv("./Data/Database_all_merged_2024-05-27.csv")

skimr::skim_without_charts(global_data)

#### Add Holdridge Life Zones (HLZ)

hlz_shp <- sf::read_sf("./Data/HLZ/holdrid/holdrid.shp")

# plot(hlz_shp[,8])

st_crs(hlz_shp) <- 4326

data_sf <- sf::st_as_sf(global_data, coords = c("Longitude", "Latitude"), crs = 4326)

hlz_sf <- sf::st_join(data_sf, hlz_shp)

summary(hlz_sf)

### A few entries get duplicated: 38639-38617 = 22; because they are at the border of two zones
hlz_dup <- hlz_sf %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) 

view(hlz_dup)

# Have a look at duplicated entries
tmap_mode("view")
tm_shape(hlz_shp) +
  tm_polygons(col = "DESC") +
  tm_shape(hlz_sf %>%
             filter(ID == "Sukoi_14462"|
                      ID == "Schrumpf_2022_N55_N55"|
                      ID == "Schrumpf_2022_N68_N68"|
                      ID == "Schrumpf_2022_N73_N73")) +
  tm_dots()

# Convert HLZ spatial dataset back into normal datset with long/lat as columns
hlz_data <- st_set_geometry(hlz_sf, NULL)

hlz_all <- hlz_data %>% 
  left_join(global_data)

skimr::skim_without_charts(hlz_all)

## Remove duplicates that where created during merging
# Manually remove duplicates; only keep one entry:
# if 3 duplicates chose HLZ that is two times present
# if 2 duplicates make decision based on map

hlz_dup <- hlz_all %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) 

hlz_all_wo_dup <- hlz_all %>% 
  anti_join(hlz_dup)

Sukoi_14462 <- hlz_all %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) %>% 
  filter(ID == "Sukoi_14462" & DESC == "Polar rain tundra" & HOLDRIG_ID == 365) 

Schrumpf_2022_N55_N55 <- hlz_all %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) %>% 
  filter(ID == "Schrumpf_2022_N55_N55" & DESC == "Boreal moist forest") 

Schrumpf_2022_N68_N68 <- hlz_all %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) %>% 
  filter(ID == "Schrumpf_2022_N68_N68" & DESC == "Polar wet tundra")

Schrumpf_2022_N73_N73 <- hlz_all %>% 
  janitor::get_dupes(ID, hzn_top, hzn_bot, hzn_mid, carbon, al_ox, fe_ox) %>% 
  filter(ID == "Schrumpf_2022_N73_N73" & DESC == "Polar desert")

hlz_all_clean <- hlz_all_wo_dup %>% 
  full_join(Sukoi_14462) %>% 
  full_join(Schrumpf_2022_N55_N55) %>% 
  full_join(Schrumpf_2022_N68_N68) %>% 
  full_join(Schrumpf_2022_N73_N73) %>% 
  dplyr::select(-dupe_count)

## Add WorldClim data
# Need to download WorldClim data first

MAP_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP_raster <- raster(MAP_dir)
MAP <- raster::extract(MAP_raster, cbind(hlz_all_clean$Longitude,
                                         hlz_all_clean$Latitude))

MAT_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT_raster <- raster(MAT_dir)
MAT <- raster::extract(MAT_raster, cbind(hlz_all_clean$Longitude,
                                         hlz_all_clean$Latitude))

hlz_all_climate <- cbind(hlz_all_clean, MAP, MAT) %>% 
  tibble()

summary(hlz_all_climate)

hlz_all_climate %>% 
  filter(is.na(DESC)) %>% 
  summarise(n_id = n_distinct(ID, Longitude, Latitude))

hlz_all_climate %>%
  filter(is.na(DESC)) %>%
  count(ID, Origin) %>%
  view()

## Manually assign missing MAP/MAT
# Use MAP/MAT reported in study
hlz_all_climate$MAP <- replace(hlz_all_climate$MAP,
                               which(hlz_all_climate$Latitude == -66.283333 &
                                       is.na(hlz_all_climate$MAP)), 180)

hlz_all_climate$MAT <- replace(hlz_all_climate$MAT,
                               which(hlz_all_climate$Latitude == -66.283333 &
                                       is.na(hlz_all_climate$MAT)), -9.3)

hlz_all_climate$MAP <- replace(hlz_all_climate$MAP,
                               which(hlz_all_climate$ID == "Sasalaguan_10980" &
                                       is.na(hlz_all_climate$MAP)), 768)

hlz_all_climate$MAT <- replace(hlz_all_climate$MAT,
                               which(hlz_all_climate$ID == "Sasalaguan_10980" &
                                       is.na(hlz_all_climate$MAT)), 27.3)

hlz_all_climate$MAP <- replace(hlz_all_climate$MAP,
                               which(hlz_all_climate$ID == "(unnamed)_12506" &
                                       is.na(hlz_all_climate$MAP)), 307)

hlz_all_climate$MAT <- replace(hlz_all_climate$MAT,
                               which(hlz_all_climate$ID == "(unnamed)_12506" &
                                       is.na(hlz_all_climate$MAT)), 26.8)

hlz_all_climate$MAP <- replace(hlz_all_climate$MAP,
                               which(grepl("CB", hlz_all_climate$ID) &
                                       is.na(hlz_all_climate$MAP)), 3200)

hlz_all_climate$MAT <- replace(hlz_all_climate$MAT,
                               which(grepl("CB", hlz_all_climate$ID) &
                                       is.na(hlz_all_climate$MAT)), 9.5)


## Manually assign missing HLZ (based on MAP and MAT)
hlz_all_fill <- hlz_all_climate %>%  
  mutate(DESC_filled = case_when(
    is.na(DESC) & MAT < 1.5 ~ "Polar desert",
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP < 125 ~ "Polar dry tundra",
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 125 & MAP < 250 ~ "Polar moist tundra",
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 250 & MAP < 500 ~ "Polar wet tundra",
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 500 & MAP < 1000 ~ "Polar rain tundra",
    is.na(DESC) & MAT >= 3 & MAT < 6 & MAP < 125 ~ "Boreal desert",
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 125 & MAP < 250 ~ "Boreal dry bush",
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 250 & MAP < 500 ~ "Boreal moist forest",
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 500 & MAP < 1000 ~ "Boreal wet forest",
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 1000 ~ "Boreal rain forest",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP < 125 ~ "Cool temperate desert",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 125 & MAP < 250 ~ "Cool temperate desert bush",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 250 & MAP < 500 ~ "Cool temperate steppe",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 250 & MAP < 500 ~ "Cool temperate steppe",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 500 & MAP < 1000 ~ "Cool temperate moist forest",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 1000 & MAP < 2000 ~ "Cool temperate wet forest",
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 2000 ~ "Cool temperate rain forest",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP < 125 ~ "Warm temperate desert",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 125 & MAP < 250 ~ "Warm temperate desert bush",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 250 & MAP < 500 ~ "Warm temperate thorn steppe",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 500 & MAP < 1000 ~ "Warm temperate dry forest",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 1000 & MAP < 2000 ~ "Warm temperate moist forest",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 2000 & MAP < 4000 ~ "Warm temperate wet forest",
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 4000 ~ "Warm temperate rain forest",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP < 125 ~ "Subtropical desert",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 125 & MAP < 250 ~ "Subtropical desert bush",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 250 & MAP < 500 ~ "Subtropical thorn steppe",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 500 & MAP < 1000 ~ "Subtropical dry forest",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 1000 & MAP < 2000 ~ "Subtropical moist forest",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 2000 & MAP < 4000 ~ "Subtropical wet forest",
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 4000 ~ "Subtropical rain forest",
    is.na(DESC) & MAT >= 24 & MAP < 125 ~ "Tropical desert",
    is.na(DESC) & MAT >= 24 & MAP >= 125 & MAP < 250 ~ "Tropical desert bush",
    is.na(DESC) & MAT >= 24 & MAP >= 250 & MAP < 500 ~ "Tropical thorn steppe",
    is.na(DESC) & MAT >= 24 & MAP >= 500 & MAP < 1000 ~ "Tropical very dry forest",
    is.na(DESC) & MAT >= 24 & MAP >= 1000 & MAP < 2000 ~ "Tropical dry forest",
    is.na(DESC) & MAT >= 24 & MAP >= 2000 & MAP < 4000 ~ "Tropical moist forest",
    is.na(DESC) & MAT >= 24 & MAP >= 4000 & MAP < 8000 ~ "Tropical wet forest",
    is.na(DESC) & MAT >= 24 & MAP >= 8000 ~ "Tropical rain forest",
    #Manually assign missing HLZ (e.g. too close to water body)
    is.na(DESC) & grepl("Sogi", ID) ~ "Tropical moist forest",
    is.na(DESC) & grepl("Susannaberg", ID) ~ "Tropical dry forest",
    is.na(DESC) & grepl("Not named", ID) ~ "Tropical dry forest",
    is.na(DESC) & grepl("Bullards", ID) ~ "Cool temperate rain forest",
    is.na(DESC) & grepl("Haro_29452", ID) ~ "Cool temperate moist forest",
    is.na(DESC) & grepl("Hoypus", ID) ~ "Cool temperate moist forest",
    is.na(DESC) & grepl("Sarkar_10662", ID) ~ "Cool temperate rain forest",
    is.na(DESC) & grepl("SND_21395", ID) ~ "Tropical dry forest",
    is.na(DESC) & grepl("Peleliu", ID) ~ "Tropical moist forest",
    is.na(DESC) & grepl("Ngedebus", ID) ~ "Tropical moist forest",
    TRUE ~ DESC
  )) %>% 
  mutate(ZONE_filled = case_when(
    is.na(DESC) & MAT < 1.5 ~ 2,
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP < 125 ~ 2,
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 125 & MAP < 250 ~ 4,
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 250 & MAP < 500 ~ 5,
    is.na(DESC) & MAT >= 1.5 & MAT < 3 & MAP >= 500 & MAP < 1000 ~ 6,
    is.na(DESC) & MAT >= 3 & MAT < 6 & MAP < 125 ~ 7,
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 125 & MAP < 250 ~ 8,
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 250 & MAP < 500 ~ 9,
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 500 & MAP < 1000 ~ 10,
    is.na(DESC) & MAT >= 3 & MAT < 6  & MAP >= 1000 ~ 11,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP < 125 ~ 12,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 125 & MAP < 250 ~ 13,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 250 & MAP < 500 ~ 14,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 500 & MAP < 1000 ~ 15,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 1000 & MAP < 2000 ~ 16,
    is.na(DESC) & MAT >= 6 & MAT < 12 & MAP >= 2000 ~ 17,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP < 125 ~ 18,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 125 & MAP < 250 ~ 19,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 250 & MAP < 500 ~ 20,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 500 & MAP < 1000 ~ 21,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 1000 & MAP < 2000 ~ 22,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 2000 & MAP < 4000 ~ 23,
    is.na(DESC) & MAT >= 12 & MAT < 18 & MAP >= 4000 ~ 24,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP < 125 ~ 25,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 125 & MAP < 250 ~ 26,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 250 & MAP < 500 ~ 27,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 500 & MAP < 1000 ~ 28,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 1000 & MAP < 2000 ~ 29,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 2000 & MAP < 4000 ~ 30,
    is.na(DESC) & MAT >= 18 & MAT < 24 & MAP >= 4000 ~ 31,
    is.na(DESC) & MAT >= 24 & MAP < 125 ~ 32,
    is.na(DESC) & MAT >= 24 & MAP >= 125 & MAP < 250 ~ 33,
    is.na(DESC) & MAT >= 24 & MAP >= 250 & MAP < 500 ~ 34,
    is.na(DESC) & MAT >= 24 & MAP >= 500 & MAP < 1000 ~ 35,
    is.na(DESC) & MAT >= 24 & MAP >= 1000 & MAP < 2000 ~ 36,
    is.na(DESC) & MAT >= 24 & MAP >= 2000 & MAP < 4000 ~ 37,
    is.na(DESC) & MAT >= 24 & MAP >= 4000 & MAP < 8000 ~ 38,
    is.na(DESC) & MAT >= 24 & MAP >= 8000 ~ 39,
    #Manually assign missing HLZ (e.g. too close to water body)
    is.na(DESC) & grepl("Sogi", ID) ~ 37,
    is.na(DESC) & grepl("Susannaberg", ID) ~ 36,
    is.na(DESC) & grepl("Not named", ID) ~ 36,
    is.na(DESC) & grepl("Bullards", ID) ~ 17,
    is.na(DESC) & grepl("Haro_29452", ID) ~ 15,
    is.na(DESC) & grepl("Hoypus", ID) ~ 15,
    is.na(DESC) & grepl("Sarkar_10662", ID) ~ 17,
    is.na(DESC) & grepl("SND_21395", ID) ~ 36,
    is.na(DESC) & grepl("Peleliu", ID) ~ 37,
    is.na(DESC) & grepl("Ngedebus", ID) ~ 37,
    TRUE ~ ZONE
  ))

hlz_all_fill %>% 
  filter(is.na(DESC_filled)) %>%
  summarise(n_id = n_distinct(ID, Longitude, Latitude))
  
skimr::skim_without_charts(hlz_all_fill)

## Final dataset
hlz_final <- hlz_all_fill %>% 
  dplyr::select(-c(AREA, PERIMETER, HOLDRIG_, CASE_, FREQUENCY, SYMBOL))

write_csv(x = hlz_final, file = paste0("./Data/Database_HLZ_", 
                                       Sys.Date(), ".csv"))

