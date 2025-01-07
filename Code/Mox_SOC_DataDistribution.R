## Global SOC - Alox ##
## 2024-05-27 ##
## Sophie von Fromm ##

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(broom)
library(rstatix)
library(ggpattern)

##Units - carbon in wt-%; oxalate-extractable metals in g/kg

## Load final database
all_data <- read_csv("./Data/Database_HLZ_grp_2024-09-09.csv") %>% 
  mutate(al_fe_ox = al_ox + (1/2*fe_ox))

names(all_data)

# Count number of total samples and unique profiles
all_data %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

# Create ordered factors for the different HLZ groups
all_data$HLZ_temp <- factor(all_data$HLZ_temp, 
                            levels = unique(all_data$HLZ_temp[order(all_data$ZONE_filled)]), 
                            ordered = TRUE)

all_data$HLZ_temp1 <- factor(all_data$HLZ_temp1, 
                             levels = unique(all_data$HLZ_temp1[order(all_data$ZONE_filled)]), 
                             ordered = TRUE)

all_data$HLZ_moist <- factor(all_data$HLZ_moist, 
                             levels = c("arid", "semiarid", "subhumid", "humid", "perhumid"), 
                             ordered = TRUE)

all_data$HLZ_new <- factor(all_data$HLZ_new, 
                           levels = unique(all_data$HLZ_new[order(all_data$ZONE_filled)]), 
                           ordered = TRUE)

# Have a quick look at the entire data
skimr::skim_without_charts(all_data)

# Look for HLZ that have less than 10 profiles
all_data %>% 
  group_by(HLZ_new) %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID)) %>% 
  filter(n_profiles < 11)

## Define color schemes for plotting
polar <- brewer.pal(5, "Greys")
boreal <- brewer.pal(5, "Blues")
cool_temp <- brewer.pal(6, "Purples")
warm_temp <- brewer.pal(7, "Greens")
subtrop <- brewer.pal(6, "Oranges")
trop <- brewer.pal(5, "RdPu")

hlz_color <- c(polar, boreal, cool_temp, warm_temp, subtrop, trop)

hlz_color_temp <- c("#636363", "#3182BD", "#9E9AC8", "#74C476", "#FD8D3C", "#C51B8A")

hlz_color_moist <- c("#e5e5b7", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")

## Check depth distribution of entire dataset
plotly::ggplotly(
  all_data %>% 
    drop_na(al_ox, fe_ox) %>% 
    ggplot() +
    stat_ecdf(aes(x = hzn_bot), color = "red") +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_hline(yintercept = 0.75, linetype = "dashed") +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1050,100))
)

# based on hzn_bot: 95% of the data <= 200 cm
all_data %>% 
  filter(hzn_bot <= 200) %>% 
  skimr::skim_without_charts(hzn_top, hzn_mid, hzn_bot)

# count number of profiles and samples
all_data %>% 
  filter(hzn_bot <= 200) %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

# filter data hzn_bot <= 200 cm and create depth bins
all_data_depth <- all_data %>% 
  filter(hzn_bot <= 200) %>% 
  mutate(depth_cat = cut(hzn_mid, breaks = c(0, 20, 50, 100, 200), 
                         labels = c("0-20cm","20-50cm","50-100cm", "100-200cm"),
                         ordered = TRUE))

# Count number of total samples and unique profiles
all_data_depth %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

# Count number of total samples by HLZ
all_data_depth %>% 
  group_by(ZONE_filled, DESC_filled) %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID)) %>% view()

# Look for HLZ that have less than 10 profiles
all_data_depth %>% 
  group_by(HLZ_new) %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID)) %>% 
  filter(n_profiles < 11)

## Check data again
skimr::skim_without_charts(all_data_depth)

## Check Alox distribution
plotly::ggplotly(
  all_data_depth %>% 
    drop_na(al_ox) %>% 
    ggplot() +
    stat_ecdf(aes(x = al_ox), color = "red") +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_hline(yintercept = 0.75, linetype = "dashed") +
    theme_bw(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1050,100))
)

## Map sampling locations
world <- map_data("world")

global_map <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgrey", linewidth = 0.01)  +
  geom_point(data = all_data,
             aes(x = Longitude, y = Latitude, color = HLZ_new),
             size = 1.5, shape = 21) +
  theme_bw(base_size = 14) +
  theme(rect = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        plot.margin = margin(t = 5, b = 0, l = 5, r = 5),
        legend.position = "none") +
  coord_sf() +
  scale_color_manual("HLZ zones", values = hlz_color) +
  scale_x_continuous("", labels = c("100°W", "0", "100°E"), 
                     breaks = c(-100,0,100), expand = c(0,0)) +
  scale_y_continuous("",labels = c("50°S", "0", "50°N"), 
                     breaks = c(-50,0,50), expand = c(0,0))

ggsave(file = paste0("./Output/Figure1a_", Sys.Date(), ".jpeg"),
       width = 11, height = 6)

### Holdridge Life zones - distribution

## Relative data distribution compared to relative distribution of HLZ
# Load HLZ area (based on Jungkunst et al. 2021)
HLZ_area <- read_csv("./Data/HLZ_area.csv") %>% 
  #remove ice for which we don't have data anyway
  filter(DESC_filled != "Ice")

# Calculate total area and area coverd by our profiles for each HLZ
rel_dis_HLZ <- all_data_depth  %>%
  distinct(ID, .keep_all = TRUE) %>% 
  dplyr::select(ID, DESC_filled, HLZ_moist, HLZ_temp, HLZ_temp1, HLZ_new) %>% 
  dplyr::count(DESC_filled, HLZ_temp, HLZ_moist, HLZ_temp1, HLZ_new) %>% 
  right_join(HLZ_area) %>% 
  replace_na(list(n = 0)) %>% 
  mutate(n_sum = sum(n),
         area_sum = sum(area_km2_10_6)) %>% 
  rowwise() %>% 
  mutate(synthesis = n*100/n_sum,
         global = area_km2_10_6*100/area_sum,
         HLZ_temp = factor(HLZ_temp, 
                           levels = unique(HLZ_area$HLZ_temp[order(HLZ_area$ZONE)]), 
                           ordered = TRUE))

# Calculate relative distribution based on moisture groups
rel_dis_HLZ %>% 
  group_by(HLZ_moist) %>% 
  reframe(synthesis = sum(synthesis),
          global = sum(global)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(synthesis_rel = synthesis*100/global)

# Calculate relative distribution based on temperature groups
rel_dis_HLZ %>% 
  group_by(HLZ_temp) %>% 
  reframe(synthesis = sum(synthesis),
          global = sum(global)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(synthesis_rel = synthesis*100/global)

# Plot relative distribution based on global coverage and our profiles
p_dis_HLZ <- rel_dis_HLZ %>% 
  #combine temperature zones from the same moisture group and calculate their shared distribution
  unite("HLZ_moist_temp", c(HLZ_temp1, HLZ_moist), sep = " ", remove = FALSE) %>% 
  group_by(HLZ_moist_temp) %>% 
  reframe(HLZ_moist = HLZ_moist,
          HLZ_temp = HLZ_temp,
          synthesis = sum(synthesis),
          global = sum(global)) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(ZONE = case_when(
    grepl("boreal", HLZ_temp) ~ 2,
    grepl("cool temperate", HLZ_temp) ~ 3,
    grepl("warm temperate", HLZ_temp) ~ 4,
    grepl("subtropical", HLZ_temp) ~ 5,
    grepl("tropical", HLZ_temp) ~ 6,
    TRUE ~ 1
  )) %>%
  pivot_longer(synthesis:global, values_to = "percentage", names_to = "study") %>% 
  mutate(HLZ_moist = factor(HLZ_moist,
                           levels = c("arid", "semiarid", "subhumid",
                                      "humid", "perhumid"),
                           ordered = TRUE)) %>% 
  ggplot() +
  geom_bar(aes(x = reorder(HLZ_moist_temp, ZONE), y = percentage,
               fill = HLZ_temp, alpha = study),
           color = "darkgrey", stat = "identity", position = "dodge") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 5, b = 0, l = 25, r = 5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_text(color = "gray")) +
  facet_grid(cols = vars(HLZ_moist), scales = "free_x", space = "free") +
  scale_x_discrete("") +
  scale_y_continuous("Relative distribution [%]", trans = "log1p",
                     expand = c(0,0), limits = c(0,35)) +
  scale_alpha_manual("Data", values = c(0.5,1)) +
  scale_fill_manual( values = hlz_color_temp) +
  guides(fill = "none")

# Change background color of facets to match moisture groups
g <- ggplot_gtable(ggplot_build(p_dis_HLZ))
stript <- which(grepl('strip-t', g$layout$name))
fills <- hlz_color_moist
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

ggsave(g, file = paste0("./Output/Figure1b_", Sys.Date(), ".jpeg"), 
       width = 12, height = 6)

### Data distribution based on soil age
all_data_depth %>% 
  group_by(HLZ_new, soil_age) %>% 
  reframe(n = n()) %>% 
  group_by(HLZ_new) %>% 
  mutate(rel = n * 100 / sum(n)) %>% view()
  ggplot() +
  geom_bar(aes(x = HLZ_new, y = rel, fill = soil_age), stat = "identity") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_y_continuous("Relative distribution", expand = c(0,0))

#### Data distribution - SOC, Alox, Feox and Mox across HLZ
### All HLZ
plot_box <- function(dataset, y){
  dataset %>% 
    ggplot(aes(x = reorder(HLZ_new, ZONE_filled), 
               fill = reorder(HLZ_new, ZONE_filled))) +
    geom_boxplot(aes(y = .data[[y]])) +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_manual("HLZ zones", values = hlz_color) +
    scale_x_discrete("") 
}

## Boxplots
# SOC
p_soc <- plot_box(all_data_depth %>% 
                 drop_na(carbon), 
               y = "carbon") +
  theme(plot.margin = margin(t = 0, b = 0, l = 13, r = 5)) +
  scale_y_continuous("SOC [wt-%]",
                     expand = c(0,0), trans = "log1p", limits = c(0,21))

# Mox
p_mox <- plot_box(all_data_depth %>% 
                 drop_na(al_fe_ox), 
               y = "al_fe_ox") +
  scale_y_continuous(expression(paste("M"[ox], " [g/kg]")),
                     expand = c(0,0), trans = "log1p", limits = c(0,176))

# Plot Mox and SOC together
p_mox.1 <- p_mox +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(t = 5, b = 0, l = 5, r = 5))

ggarrange(p_mox.1, p_soc, nrow = 2, heights = c(1,1.35),
          labels = c("a)", "b)"), vjust = c(1.5, 1))

ggsave(file = paste0("./Output/Figure2_", Sys.Date(), ".jpeg"), 
       width = 13, height = 9)

# Alox
p_alox <- plot_box(all_data_depth %>% 
           drop_na(al_ox), 
         y = "al_ox") +
  scale_y_continuous(expression(paste("Al"[ox], " [g/kg]")),
                     expand = c(0,0), trans = "log1p", limits = c(0,176))

# Feox
p_feox <- plot_box(all_data_depth %>% 
                     drop_na(fe_ox), 
                   y = "fe_ox")+
  theme(plot.margin = margin(t = 0, b = 0, l = 5, r = 5)) +
  scale_y_continuous(expression(paste("Fe"[ox], " [g/kg]")),
                     expand = c(0,0), trans = "log1p", limits = c(0,176))

# Plot Alox and Feox together
p_alox.1 <- p_alox +
  theme(axis.text.x = element_blank(),
        plot.margin = margin(t = 5, b = 0, l = 5, r = 5))

ggarrange(p_alox.1, p_feox, nrow = 2, heights = c(1,1.35),
          labels = c("a)", "b)"), vjust = c(1.5, 1))

ggsave(file = paste0("./Output/FigureA3_", Sys.Date(), ".jpeg"), 
       width = 13, height = 9)

## Statistical test: Kruskal-Wallis & Dunn test
# SOC
all_data_depth %>% 
  kruskal_test(carbon ~ HLZ_new)

dunn_test_all_SOC <- all_data_depth %>%
  dunn_test(carbon ~ HLZ_new, p.adjust.method = "bonferroni")

dunn_test_all_SOC %>% 
  count(p.adj.signif)

dunn_test_all_SOC$group1 <- factor(dunn_test_all_SOC$group1, 
                                   levels = levels(all_data_depth$HLZ_new), 
                                   ordered = TRUE)

dunn_test_all_SOC$group2 <- factor(dunn_test_all_SOC$group2, 
                                   levels = levels(all_data_depth$HLZ_new), 
                                   ordered = TRUE)

dunn_test_all_SOC %>% 
  ggplot(aes(x = group1, y = group2, 
             fill = p.adj.signif)) +
  geom_tile(color = "black") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete("", expand = c(0,0)) +
  scale_y_discrete("", expand = c(0,0)) +
  scale_fill_manual(values = c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "white")) +
  geom_rect(xmin = 0.5, xmax = 5.5, ymin = 0.5, ymax = 4.5, fill = NA, color = hlz_color_temp[1],
            linewidth = 2) +
  geom_rect(xmin = 5.5, xmax = 10.5, ymin = 4.5, ymax = 9.5, fill = NA, color = hlz_color_temp[2],
            linewidth = 2) +
  geom_rect(xmin = 10.5, xmax = 16.5, ymin = 9.5, ymax = 15.5, fill = NA, color = hlz_color_temp[3],
            linewidth = 2) +
  geom_rect(xmin = 16.5, xmax = 23.5, ymin = 15.5, ymax = 22.5, fill = NA, color = hlz_color_temp[4],
            linewidth = 2) +
  geom_rect(xmin = 23.5, xmax = 29.5, ymin = 22.5, ymax = 28.5, fill = NA, color = hlz_color_temp[5],
            linewidth = 2) +
  geom_rect(xmin = 29.5, xmax = 33.5, ymin = 28.5, ymax = 33.5, fill = NA, color = hlz_color_temp[6],
            linewidth = 2) +
  geom_text(x = 4, y = 1.5, label = "(sub-)polar", color = hlz_color_temp[1]) +
  geom_text(x = 9.5, y = 5.5, label = "boreal", color = hlz_color_temp[2]) +
  geom_text(x = 14.5, y = 10.5, label = "cool temperate", color = hlz_color_temp[3]) +
  geom_text(x = 21.5, y = 16.5, label = "warm temperate", color = hlz_color_temp[4]) +
  geom_text(x = 28, y = 23.5, label = "subtropical", color = hlz_color_temp[5]) +
  geom_text(x = 32.5, y = 29.5, label = "tropical", color = hlz_color_temp[6])
ggsave(file = paste0("./Output/FigureA2_", Sys.Date(), ".jpeg"), 
       width = 14, height = 8)

# Mox
all_data_depth %>% 
  kruskal_test(al_fe_ox ~ HLZ_new)

# Dunn post-hoc test
dunn_test_all <- all_data_depth %>%
  dunn_test(al_fe_ox ~ HLZ_new, p.adjust.method = "bonferroni") 

dunn_test_all %>% 
  count(p.adj.signif)

dunn_test_all %>% 
  mutate(p_adj = p.adj/p) %>% 
  view()

dunn_test_all$group1 <- factor(dunn_test_all$group1, 
                               levels = levels(all_data_depth$HLZ_new), 
                               ordered = TRUE)

dunn_test_all$group2 <- factor(dunn_test_all$group2, 
                               levels = levels(all_data_depth$HLZ_new), 
                               ordered = TRUE)

dunn_test_all %>% 
  ggplot(aes(x = group1, y = group2, 
             fill = p.adj.signif)) +
  geom_tile(color = "black") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete("", expand = c(0,0)) +
  scale_y_discrete("", expand = c(0,0)) +
  scale_fill_manual(values = c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "white")) +
  geom_rect(xmin = 0.5, xmax = 5.5, ymin = 0.5, ymax = 4.5, fill = NA, color = hlz_color_temp[1],
            linewidth = 2) +
  geom_rect(xmin = 5.5, xmax = 10.5, ymin = 4.5, ymax = 9.5, fill = NA, color = hlz_color_temp[2],
            linewidth = 2) +
  geom_rect(xmin = 10.5, xmax = 16.5, ymin = 9.5, ymax = 15.5, fill = NA, color = hlz_color_temp[3],
            linewidth = 2) +
  geom_rect(xmin = 16.5, xmax = 23.5, ymin = 15.5, ymax = 22.5, fill = NA, color = hlz_color_temp[4],
            linewidth = 2) +
  geom_rect(xmin = 23.5, xmax = 29.5, ymin = 22.5, ymax = 28.5, fill = NA, color = hlz_color_temp[5],
            linewidth = 2) +
  geom_rect(xmin = 29.5, xmax = 33.5, ymin = 28.5, ymax = 33.5, fill = NA, color = hlz_color_temp[6],
            linewidth = 2) +
  geom_text(x = 4, y = 1.5, label = "(sub-)polar", color = hlz_color_temp[1]) +
  geom_text(x = 9.5, y = 5.5, label = "boreal", color = hlz_color_temp[2]) +
  geom_text(x = 14.5, y = 10.5, label = "cool temperate", color = hlz_color_temp[3]) +
  geom_text(x = 21.5, y = 16.5, label = "warm temperate", color = hlz_color_temp[4]) +
  geom_text(x = 28, y = 23.5, label = "subtropical", color = hlz_color_temp[5]) +
  geom_text(x = 32.5, y = 29.5, label = "tropical", color = hlz_color_temp[6])
ggsave(file = paste0("./Output/FigureA1_", Sys.Date(), ".jpeg"), 
       width = 14, height = 8)

### By moisture
plot_box_moist <- function(dataset, y){
  dataset %>% 
    ggplot(aes(x = HLZ_moist, fill = HLZ_moist)) +
    geom_boxplot(aes(y = .data[[y]])) +
    facet_wrap(~depth_cat) +
    theme_bw(base_size = 18) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_manual("HLZ zones", values = hlz_color_moist) +
    scale_x_discrete("") 
}

# SOC
plot_box_moist(all_data_depth %>% 
                 drop_na(carbon), 
               y = "carbon")  +
  scale_y_continuous("SOC [wt-%]",
                     expand = c(0,0), trans = "log1p", limits = c(0,21))
ggsave(file = paste0("./Output/FigureA5_", Sys.Date(), ".jpeg"), 
       width = 13, height = 7)

all_data_depth %>% 
  group_by(depth_cat) %>% 
  kruskal_test(carbon ~ HLZ_moist)

all_data_depth %>%
  group_by(depth_cat) %>% 
  dunn_test(carbon ~ HLZ_moist, p.adjust.method = "bonferroni") %>% 
  view()

all_data_depth %>%
  group_by(HLZ_moist) %>% 
  dunn_test(carbon ~ depth_cat, p.adjust.method = "bonferroni") %>% 
  view()

# Mox
plot_box_moist(dataset = all_data_depth %>% 
                 drop_na(al_fe_ox), y = "al_fe_ox") +
  scale_y_continuous(expression(paste("M"[ox], " [g/kg]")),
                     expand = c(0,0), trans = "log1p", limits = c(0,176))
ggsave(file = paste0("./Output/FigureA4_", Sys.Date(), ".jpeg"),
       width = 13, height = 7)

all_data_depth %>% 
  group_by(depth_cat) %>% 
  kruskal_test(al_fe_ox ~ HLZ_moist)

all_data_depth %>%
  group_by(depth_cat) %>% 
  dunn_test(al_fe_ox ~ HLZ_moist, p.adjust.method = "bonferroni") %>% 
  view()

all_data_depth %>%
  group_by(HLZ_moist) %>% 
  dunn_test(al_fe_ox ~ depth_cat, p.adjust.method = "bonferroni") %>% 
  view()

### By temperature
plot_box_temp <- function(dataset, y){
  dataset %>% 
    ggplot(aes(x = HLZ_temp, fill = HLZ_temp)) +
    geom_boxplot(aes(y = .data[[y]])) +
    facet_wrap(~depth_cat) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_manual("HLZ zones", values = hlz_color_temp) +
    scale_x_discrete("") 
}

# SOC
plot_box_temp(all_data_depth %>% 
                 drop_na(carbon), 
              y = "carbon") +
  scale_y_continuous("SOC [wt-%]",
                     expand = c(0,0), trans = "log1p", limits = c(0,21))
ggsave(file = paste0("./Output/FigureA7_", Sys.Date(), ".jpeg"), 
       width = 13, height = 7)

# SOC
all_data_depth %>% 
  group_by(depth_cat) %>% 
  kruskal_test(carbon ~ HLZ_temp)

all_data_depth %>%
  group_by(depth_cat) %>% 
  dunn_test(carbon ~ HLZ_temp, p.adjust.method = "bonferroni") %>% 
  view()

# Mox
plot_box_temp(dataset = all_data_depth %>% 
                drop_na(al_fe_ox), y = "al_fe_ox") +
  scale_y_continuous(expression(paste("M"[ox], " [g/kg]")),
                     expand = c(0,0), trans = "log1p", limits = c(0,176))
ggsave(file = paste0("./Output/FigureA6_", Sys.Date(), ".jpeg"), 
       width = 13, height = 7)

all_data_depth %>% 
  group_by(depth_cat) %>% 
  kruskal_test(al_fe_ox ~ HLZ_temp)

all_data_depth %>%
  group_by(depth_cat) %>% 
  dunn_test(al_fe_ox ~ HLZ_temp, p.adjust.method = "bonferroni") %>% 
  view()

