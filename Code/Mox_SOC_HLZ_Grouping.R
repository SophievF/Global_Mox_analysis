## Global SOC - Mox ##
## 2024-03-01 ##
## Sophie von Fromm ##
## Grouping of HLZ ##

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(sf)
library(igraph)

### NO NEED TO RUN CODE, FINAL DATASET HAS BEEN SAVED ####
# The code below allows to group the HLZ zones as proposed in the paper

### Units - carbon in wt-%; oxalate-extractable metals in g/kg 

## Load modified database
all_data <- read_csv("./Data/Database_HLZ_2024-08-14.csv")

all_data %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

## HLZ grouping based on temperature and moisture
all_data_grp <- all_data %>% 
  # drop_na(DESC_filled) %>% 
   mutate(HLZ_temp = case_when(
    grepl("Polar", DESC_filled) ~ "(sub-)polar",
    grepl("Boreal", DESC_filled) ~ "boreal",
    grepl("Tropical", DESC_filled) ~ "tropical",
    grepl("Subtropical", DESC_filled) ~ "subtropical",
    grepl("Cool temperate", DESC_filled) ~ "cool temperate",
    grepl("Warm temperate", DESC_filled) ~ "warm temperate",
  )) %>% 
  mutate(HLZ_temp1 = case_when(
    grepl("Polar desert", DESC_filled) ~ "polar",
    grepl("Polar", DESC_filled) ~ "subpolar",
    TRUE ~ HLZ_temp
    )) %>% 
  mutate(HLZ_moist_full = case_when(
    DESC_filled == "Tropical desert" ~ "superarid",
    DESC_filled == "Tropical desert bush" ~ "perarid",
    DESC_filled == "Subtropical desert" ~ "perarid",
    DESC_filled == "Warm temperate desert" ~ "perarid",
    DESC_filled == "Tropical thorn steppe" ~ "arid",
    DESC_filled == "Subtropical desert bush" ~ "arid",
    DESC_filled == "Warm temperate desert bush" ~ "arid",
    DESC_filled == "Cool temperate desert" ~ "arid",
    DESC_filled == "Tropical very dry forest" ~ "semiarid",
    DESC_filled == "Subtropical thorn steppe" ~ "semiarid",
    DESC_filled == "Warm temperate thorn steppe" ~ "semiarid",
    DESC_filled == "Cool temperate desert bush" ~ "semiarid",
    DESC_filled == "Boreal desert" ~ "semiarid",
    DESC_filled == "Tropical dry forest" ~ "subhumid",
    DESC_filled == "Subtropical dry forest" ~ "subhumid",
    DESC_filled == "Warm temperate dry forest" ~ "subhumid",
    DESC_filled == "Cool temperate steppe" ~ "subhumid",
    DESC_filled == "Boreal dry bush" ~ "subhumid",
    DESC_filled == "Polar dry tundra" ~ "subhumid",
    grepl("moist", DESC_filled) ~ "humid",
    DESC_filled == "Polar desert" ~ "humid",
    grepl("wet", DESC_filled) ~ "perhumid",
    grepl("rain", DESC_filled) ~ "superhumid",
  )) %>% 
  mutate(HLZ_moist = case_when(
    HLZ_moist_full == "superarid" ~ "arid",
    HLZ_moist_full == "perarid" ~ "arid",
    HLZ_moist_full == "arid" ~ "arid",
    HLZ_moist_full == "semiarid" ~ "semiarid",
    HLZ_moist_full == "subhumid" ~ "subhumid",
    HLZ_moist_full == "humid" ~ "humid",
    HLZ_moist_full == "perhumid" ~ "perhumid",
    HLZ_moist_full == "superhumid" ~ "perhumid"
  )) %>% 
  unite("HLZ_new", HLZ_temp1, HLZ_moist_full,remove = FALSE, sep = " ") %>% 
  # Group soils based on age + loess
  mutate(soil_age = case_when(
    grepl("Loess", sediment_class) ~ "young",
    lgm_vegetation_type == "Ice sheet and other permanent ice" ~ "young",
    lgm_vegetation_type == "Polar and alpine desert" ~ "young",
    # grepl("tundra", lgm_vegetation_type) ~ "young",
    # lgm_vegetation_type == "Tundra" ~ "young",
    TRUE ~ "old"
  ))

# Save final file
write_csv(x = all_data_grp, file = paste0("./Data/Database_HLZ_grp_", 
                                          Sys.Date(), ".csv"))

all_data_grp %>% 
  count(DESC_filled, HLZ_temp, HLZ_temp1, HLZ_moist_full, HLZ_moist, HLZ_new) %>% 
  view()

all_data_grp_sum <- all_data_grp %>% 
  group_by(ZONE_filled, DESC_filled, HLZ_new, HLZ_moist, HLZ_temp) %>% 
  summarise(n_profiles = n_distinct(ID))

view(all_data_grp_sum)

write_csv(x = all_data_grp_sum, file = paste0("./Output/all_data_HLZ_grp_sum_", 
                                              Sys.Date(), ".csv"))

#Subtropical desert, tropical desert and tropical desert bush are missing

## Check sequential extraction from Doetterl_2021
# This is the only study that we are aware of that used a sequential extraction, instead of a parallel one
# Data looks okay; no need to remove/modify
all_data %>% 
  filter(grepl("Doetterl_2021", ID)) %>% 
  count(DESC_filled)

all_data %>% 
  filter(DESC_filled == "Warm temperate moist forest"|
           DESC_filled == "Subtropical moist forest") %>%
  ggplot(aes(x = al_ox, y = carbon, color = grepl("Doetterl_2021", ID),
             shape = grepl("Doetterl_2021", ID),
             alpha = grepl("Doetterl_2021", ID))) +
  geom_point(size = 4) +
  facet_wrap(~DESC_filled) +
  theme_bw(base_size = 16) +
  scale_x_continuous(expression(paste("Al"[ox], " [wt-%]")),
                     labels = c(0.001,0.01,0.1,1,10), expand = c(0,0),
                     breaks = c(0.001,0.01,0.1,1,10),
                     trans = "log10", limits = c(0.001,20)) +
  scale_y_continuous("Soil organic carbon [wt-%]",
                     labels = c(0.01,0.1,1,10), expand = c(0,0),
                     breaks = c(0.01,0.1,1,10),
                     trans = "log10", limits = c(0.01,30))

