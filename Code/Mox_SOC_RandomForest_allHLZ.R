## Global SOC - Alox ##
## 2024-05-28 ##
## Sophie von Fromm ##

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(ranger)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(iml)
library(mlr3spatial)
library(mlr3spatiotempcv) 
library(iml)

#### Random Forest analysis ####

## Load and prepare data
all_data <- read_csv("./Data/Database_HLZ_grp_2024-09-09.csv") %>% 
  mutate(al_fe_ox = al_ox + (1/2*fe_ox)) %>% 
  filter(hzn_bot <= 200) %>% 
  mutate(depth_cat = cut(hzn_mid, breaks = c(0, 20, 50, 100, 200), 
                         labels = c("0-20cm","20-50cm","50-100cm", "100-200cm"),
                         ordered = TRUE))
names(all_data)

all_data %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

all_data$DESC_filled <- factor(all_data$DESC_filled, 
                               levels = unique(all_data$DESC_filled[order(all_data$ZONE_filled)]), 
                               ordered = TRUE)

all_data$HLZ_temp <- factor(all_data$HLZ_temp, 
                            levels = unique(all_data$HLZ_temp[order(all_data$ZONE_filled)]), 
                            ordered = TRUE)

all_data$HLZ_moist <- factor(all_data$HLZ_moist, 
                             levels = c("arid", "semiarid", "subhumid", "humid", "perhumid"), 
                             ordered = TRUE)

## Color Scheme
polar <- brewer.pal(5, "Greys")
boreal <- brewer.pal(5, "Blues")
cool_temp <- brewer.pal(6, "Purples")
warm_temp <- brewer.pal(7, "Greens")
subtrop <- brewer.pal(6, "Oranges")
trop <- brewer.pal(5, "RdPu")

hlz_color <- c(polar, boreal, cool_temp, warm_temp, subtrop, trop)

hlz_color_temp <- c("#636363", "#3182BD", "#9E9AC8", "#74C476", "#FD8D3C", "#C51B8A")

hlz_color_moist <- c("#e5e5b7", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")

### Random forest analysis
set.seed(42)
data_sf <- sf::st_as_sf(all_data, coords = c("Longitude", "Latitude"), crs = 4326)

data_rf <- data_sf %>% 
  drop_na(al_ox, fe_ox) %>% 
  dplyr::select(ID, carbon, depth_cat, al_ox, fe_ox, HLZ_new)

## Set-up random forest
task_rf_data <- as_task_regr_st(x = data_rf, target = "carbon")

lrn_rf_data <- lrn("regr.ranger", importance = "permutation",
                   num.trees = 1000)

# Add id as group for CV (same id kept together)
task_rf_data$set_col_roles("ID", roles = "group")
print(task_rf_data)

resampling_data <- rsmp("cv", folds = 10)
resampling_data$instantiate(task_rf_data)

## Train model & check performance
rr_data <- mlr3::resample(task = task_rf_data, learner = lrn_rf_data,
                          resampling = resampling_data, store_models = TRUE)

# rr_data <- mlr3::resample(task = task_rf_data, learner = lrn_rf_data, 
#                           resampling = resampling_data)

rr_data$aggregate(measures = msrs(c("regr.rmse", "regr.mse", "regr.rsq")))
rr_score_data_ind <- rr_data$score(measures = msrs(c("regr.rmse", "regr.mse", 
                                                     "regr.rsq")))

rr_data_measure <- data.frame(regr_rmse = rr_score_data_ind$regr.rmse,
                              regr_mse = rr_score_data_ind$regr.mse,
                              regr_rsq = rr_score_data_ind$regr.rsq)

skimr::skim_without_charts(rr_data_measure)

rr_pred_data <- rr_data$prediction(predict_sets = "test")
rr_data_prediction <- data.frame(truth = rr_pred_data$truth,
                                 response = rr_pred_data$response)

vi_data_all <- lapply(rr_data$learners, function(x) x$model$variable.importance) 

vi_data_sum <- vi_data_all %>%
  plyr::ldply() %>%
  mutate(row_sum = rowSums(pick(where(is.numeric)))) %>%
  mutate(HLZ_new = HLZ_new/row_sum*100,
         al_ox = al_ox/row_sum*100,
         fe_ox = fe_ox/row_sum*100,
         depth_cat = depth_cat/row_sum*100) %>%
  dplyr::select(-row_sum) %>%
  dplyr::summarise(across(everything(), list(median = median, mad = mad),
                          .names = "{.col}.{.fn}")) %>%
  pivot_longer(everything(), names_to = "names", values_to = "values") %>%
  separate_wider_delim(names, ".", names = c("predictor", "measure")) %>%
  pivot_wider(names_from = measure, values_from = values)

vi_data_sum  %>%
  ggplot() +
  geom_bar(aes(x = reorder(predictor, -median), y = median, fill = predictor),
           stat = "identity") +
  geom_errorbar(aes(ymin = median-mad, ymax = median+mad,
                    x = reorder(predictor, -median)), width = 0.2) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.x = element_blank()) +
  scale_y_continuous("Relative importance [%]", expand = c(0,0),
                     limits = c(0,40)) +
  scale_x_discrete(labels = c("Alox", "Depth", "Feox", "HLZ all"))
ggsave(file = paste0("./Output/FigureA14_",
                     Sys.Date(), ".jpeg"), width = 8, height = 6)

rr_data_prediction %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(shape = 21) +
  geom_rug() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC content [wt-%]", limits = c(0,20)) +
  scale_x_continuous("Predicted SOC content [wt-%]", 
                     limits = c(0,20)) +
  geom_text(aes(x = 16, y = 9), 
            label = paste("R² = ", round(median(rr_data_measure$regr_rsq), 2), 
                          " ± ", round(mad(rr_data_measure$regr_rsq), 2))) +
  geom_text(aes(x = 16, y = 8), 
            label = paste("MSE = ", round(median(rr_data_measure$regr_mse), 2), 
                          " ± ", round(mad(rr_data_measure$regr_mse), 2))) +
  geom_text(aes(x = 16.5, y = 7), 
            label = paste("RMSE = ", round(median(rr_data_measure$regr_rmse), 2), 
                          " ± ", round(mad(rr_data_measure$regr_rmse), 2), "wt-%"))

### Model interpretation (iml)
## HLZ moist
rf_data_pdp_m <- all_data %>%
  drop_na(al_ox, fe_ox) %>%
  arrange(ID) %>%
  dplyr::select(ID, carbon, depth_cat, al_ox, fe_ox, HLZ_new)

task_rf_data_pdp_m <- as_task_regr(x = rf_data_pdp_m %>%
                                     dplyr::select(-ID),
                                   target = "carbon")

lrn_rf_data_pdp_m <- lrn("regr.ranger", importance = "permutation",
                         num.trees = 1000) 

lrn_rf_data_pdp_m$train(task_rf_data_pdp_m)

lrn_rf_data_pdp_m$importance()

lrn_rf_data_pdp_m$model

model_data_m <- Predictor$new(lrn_rf_data_pdp_m, data = rf_data_pdp_m %>% 
                                dplyr::select(-ID, -carbon))

pdp_all <- FeatureEffects$new(model_data_m, method = "pdp")

#Extracting data
fun_df <- function(list_name, feat_name){
  df <- data.frame(response = list_name[[feat_name]][".value"],
                   predictor = list_name[[feat_name]][".borders"],
                   variable = feat_name)
}

names(rf_data_pdp_m)
feat_cont <- names(rf_data_pdp_m)[4:5]

df_pdp_cont <- map_df(feat_cont, ~fun_df(list_name = pdp_all$results, 
                                         feat_name = .x))

rf_data_pdp_m %>% 
  dplyr::select(al_ox, fe_ox, carbon) %>%
  pivot_longer(cols = al_ox:fe_ox,
               names_to = "variable", values_to = ".borders") %>% 
  mutate(variable = case_when(
    variable == "al_ox" ~ "oxalate Al",
    variable == "fe_ox" ~ "oxalate Fe"
  )) %>% 
  dplyr::rename(.value = carbon)

df_pdp_cont %>% 
  mutate(variable = case_when(
    variable == "al_ox" ~ "oxalate Al",
    variable == "fe_ox" ~ "oxalate Fe"
  )) %>% 
  ggplot(aes(y = .value, x = .borders)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  geom_rug(data = rf_data_pdp_m %>%
             dplyr::select(al_ox, fe_ox, carbon) %>%
             pivot_longer(cols = al_ox:fe_ox,
                          names_to = "variable", values_to = ".borders") %>%
             mutate(variable = case_when(
               variable == "al_ox" ~ "oxalate Al",
               variable == "fe_ox" ~ "oxalate Fe",
             )) %>% 
             dplyr::rename(.value = carbon), sides = "b") +
  facet_wrap(~variable, scales = "free_x") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.spacing.x = unit(1, "line")) +
  scale_x_continuous("Range of predictor", expand = c(0,0)) +
  scale_y_continuous("Predicted SOC content [wt-%]", expand = c(0,0), limits = c(0,8))
ggsave(file = paste0("./Output/FigureA13_",
                     Sys.Date(), ".jpeg"), width = 11, height = 7)

