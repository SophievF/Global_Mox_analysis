## Global SOC - Alox ##
## 2024-05-28 ##
## Sophie von Fromm ##

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(nlme)
library(MuMIn)
library(emmeans)
library(sp)
library(gstat)

#### Linear mixed-effects models: moisture groups

## Load and prepare data
all_data <- read_csv("./Data/Database_HLZ_grp_2024-09-09.csv") %>% 
  mutate(al_fe_ox = al_ox + (1/2*fe_ox)) %>% 
  filter(hzn_bot <= 200) %>% 
  mutate(depth_cat = cut(hzn_mid, breaks = c(0, 20, 50, 100, 200), 
                         labels = c("0-20cm","20-50cm","50-100cm", "100-200cm"),
                         ordered = TRUE))

names(all_data)

all_data$DESC_filled <- factor(all_data$DESC_filled, 
                               levels = unique(all_data$DESC_filled[order(all_data$ZONE_filled)]), 
                               ordered = TRUE)

all_data$HLZ_temp <- factor(all_data$HLZ_temp, 
                            levels = unique(all_data$HLZ_temp[order(all_data$ZONE_filled)]), 
                            ordered = TRUE)

all_data$HLZ_moist <- factor(all_data$HLZ_moist, 
                             levels = c("arid", "semiarid", "subhumid", "humid", "perhumid"), 
                             ordered = TRUE)

all_data %>% 
  drop_na(al_fe_ox) %>%
  group_by(HLZ_moist, depth_cat) %>% 
  summarise(n = n(),
            n_profiles = n_distinct(ID))

all_data %>% 
  drop_na(al_fe_ox) %>%
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = al_ox + 0.005,
                carbon = carbon + 0.005,
                fe_ox = fe_ox + 0.005,
                al_fe_ox = al_fe_ox + 0.005) %>% 
  group_by(HLZ_moist, depth_cat) %>% 
  skimr::skim_without_charts(al_fe_ox, carbon)

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

#### Linear mixed-effects models
### Prepare data
all_data_lm_moist <- all_data %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005)) 

# skimr::skim_without_charts(all_data_lm_moist)

all_data_lm_moist %>% 
  group_by(depth_cat, HLZ_moist) %>% 
  ggplot(aes(y = al_fe_ox, x = HLZ_moist)) +
  geom_boxplot() +
  facet_wrap(~depth_cat)

### Fit model
## Alox + 1/2 Feox
# Final model
lmm_al_fe_depth_moist <- lme(carbon ~ al_fe_ox*HLZ_moist + al_fe_ox*depth_cat + al_fe_ox*HLZ_moist*depth_cat, 
                             random = ~hzn_mid|ID,
                             data = all_data_lm_moist)

summary(lmm_al_fe_depth_moist)
r.squaredGLMM(lmm_al_fe_depth_moist)
anova(lmm_al_fe_depth_moist)
sqrt(mean(residuals(lmm_al_fe_depth_moist)^2))

LMM_moist_mox_df <- anova(lmm_al_fe_depth_moist) %>% 
  as.data.frame()

LMM_moist_mox_df$predictor <- rownames(LMM_moist_mox_df)

LMM_moist_mox_df

write_csv(x = LMM_moist_mox_df, file = paste0("./Output/TableA2_", 
                                              Sys.Date(), ".csv"))

# Graphical exploration of residuals
plot(lmm_al_fe_depth_moist)

F_Final <- fitted(lmm_al_fe_depth_moist)
R_Final <- residuals(lmm_al_fe_depth_moist) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ all_data_lm_moist$HLZ_moist)
boxplot(R_Final ~ all_data_lm_moist$depth_cat)
plot(R_Final ~ all_data_lm_moist$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_al_fe_depth_moist)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * all_data_lm_moist$al_fe_ox

{scatter.smooth(all_data_lm_moist$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*all_data_lm_moist$al_fe_ox ~ all_data_lm_moist$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous

# Check for spatial correlation
resid_norm <- resid(lmm_al_fe_depth_moist, type = "normalized")
data_moist <- data.frame(resid_norm, all_data_lm_moist$Longitude, 
                         all_data_lm_moist$Latitude)
coordinates(data_moist) <- c("all_data_lm_moist.Longitude",
                                    "all_data_lm_moist.Latitude")
bubble(data_moist, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario1 <- variogram(resid_norm ~ 1, data_moist)
plot(Vario1)

#--> Looks all good, no increase in semivariacen with distance

## Alox
# Random slope for depth
lmm_al_depth_moist <- lme(carbon ~ al_ox*HLZ_moist + al_ox*depth_cat + al_ox*HLZ_moist*depth_cat, 
                          random = ~hzn_mid|ID,
                          data = all_data_lm_moist)

summary(lmm_al_depth_moist)
r.squaredGLMM(lmm_al_depth_moist)
anova(lmm_al_depth_moist)
sqrt(mean(residuals(lmm_al_depth_moist)^2))

LMM_moist_alox_df <- anova(lmm_al_depth_moist) %>% 
  as.data.frame()

LMM_moist_alox_df$predictor <- rownames(LMM_moist_alox_df)

LMM_moist_alox_df

write_csv(x = LMM_moist_alox_df, file = paste0("./Output/TableA3_", 
                                              Sys.Date(), ".csv"))

# Graphical exploration of residuals
plot(lmm_al_depth_moist)

F_Final <- fitted(lmm_al_depth_moist)
R_Final <- residuals(lmm_al_depth_moist) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ all_data_lm_moist$HLZ_moist)
boxplot(R_Final ~ all_data_lm_moist$depth_cat)
plot(R_Final ~ all_data_lm_moist$al_ox)

# Partial residual plots
mox_c <- fixef(lmm_al_depth_moist)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * all_data_lm_moist$al_ox

{scatter.smooth(all_data_lm_moist$al_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*all_data_lm_moist$al_ox ~ all_data_lm_moist$al_ox), col = "red")}

# Check for spatial correlation
resid_norm <- resid(lmm_al_depth_moist, type = "normalized")
data_moist <- data.frame(resid_norm, all_data_lm_moist$Longitude, 
                         all_data_lm_moist$Latitude)
coordinates(data_moist) <- c("all_data_lm_moist.Longitude",
                             "all_data_lm_moist.Latitude")
bubble(data_moist, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario2 <- variogram(resid_norm ~ 1, data_moist)
plot(Vario2)

## Feox
# Random slope for depth
lmm_fe_depth_moist <- lme(carbon ~ fe_ox*HLZ_moist + fe_ox*depth_cat + fe_ox*HLZ_moist*depth_cat, 
                          random = ~hzn_mid|ID, 
                          data = all_data_lm_moist)

summary(lmm_fe_depth_moist)
anova(lmm_fe_depth_moist)
r.squaredGLMM(lmm_fe_depth_moist)
sqrt(mean(residuals(lmm_fe_depth_moist)^2))

LMM_moist_feox_df <- anova(lmm_fe_depth_moist) %>% 
  as.data.frame()

LMM_moist_feox_df$predictor <- rownames(LMM_moist_feox_df)

LMM_moist_feox_df

write_csv(x = LMM_moist_feox_df, file = paste0("./Output/TableA4_", 
                                               Sys.Date(), ".csv"))

# Graphical exploration of residuals
plot(lmm_fe_depth_moist)

F_Final <- fitted(lmm_fe_depth_moist)
R_Final <- residuals(lmm_fe_depth_moist, type = "response", scaled = TRUE) #type="response" 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ all_data_lm_moist$HLZ_moist)
boxplot(R_Final ~ all_data_lm_moist$depth_cat)
plot(R_Final ~ all_data_lm_moist$fe_ox)

mox_c <- fixef(lmm_fe_depth_moist)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * all_data_lm_moist$fe_ox

{scatter.smooth(all_data_lm_moist$fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*all_data_lm_moist$fe_ox ~ all_data_lm_moist$fe_ox), col = "red")}

# Check for spatial correlation
resid_norm <- resid(lmm_fe_depth_moist, type = "normalized")
data_moist <- data.frame(resid_norm, all_data_lm_moist$Longitude, 
                         all_data_lm_moist$Latitude)
coordinates(data_moist) <- c("all_data_lm_moist.Longitude",
                             "all_data_lm_moist.Latitude")
bubble(data_moist, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario3 <- variogram(resid_norm ~ 1, data_moist)
plot(Vario3)

## Compare all three models (based on method = Maximum Likelihood)
# lmm_al_fe_depth_moist_ML <- lme(carbon ~ al_fe_ox*HLZ_moist + al_fe_ox*depth_cat + al_fe_ox*HLZ_moist*depth_cat,
#                                 random = ~hzn_mid|ID, method = "ML",
#                                 data = all_data_lm_moist)
# lmm_al_depth_moist_ML <- lme(carbon ~ al_ox*HLZ_moist + al_ox*depth_cat + al_ox*HLZ_moist*depth_cat,
#                              random = ~hzn_mid|ID, method = "ML",
#                              data = all_data_lm_moist)
# lmm_fe_depth_moist_ML <- lme(carbon ~ fe_ox*HLZ_moist + fe_ox*depth_cat + fe_ox*HLZ_moist*depth_cat,
#                              random = ~hzn_mid|ID, method = "ML",
#                              data = all_data_lm_moist)
# anova(lmm_al_fe_depth_moist_ML, lmm_al_depth_moist_ML, lmm_fe_depth_moist_ML)

### Post-Hoc analysis
## Alox + 1/2 Feox
summary(lmm_al_fe_depth_moist)
anova(lmm_al_fe_depth_moist)
r.squaredGLMM(lmm_al_fe_depth_moist)

em_al_fe_moist <- emtrends(lmm_al_fe_depth_moist, pairwise ~ HLZ_moist | depth_cat, 
                           var = "al_fe_ox")

em_al_fe_moist$emtrends
#trend - slope direction
em_al_fe_moist$contrasts
#significant differences between slopes

as.data.frame(em_al_fe_moist$contrasts) %>% 
  filter(p.value <= 0.05) 

emip_al_fe_moist <- emmip(lmm_al_fe_depth_moist, HLZ_moist ~ al_fe_ox | depth_cat, cov.reduce = range, 
                          plotit = FALSE, CIs = TRUE)

emip_al_fe_moist$HLZ_moist <- as.ordered(emip_al_fe_moist$HLZ_moist)

emip_al_fe_data_moist <- emip_al_fe_moist %>% 
  as.data.frame() %>% 
  left_join(all_data %>% 
              distinct(HLZ_moist), 
            multiple = "all")

emip_al_fe_data_moist$HLZ_moist <- factor(emip_al_fe_data_moist$HLZ_moist, 
                                          levels = c("arid", "semiarid", "subhumid", 
                                                     "humid", "perhumid"),  
                                          ordered = TRUE)

emip_al_fe_moist %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = HLZ_moist), alpha = 0.25) +
  geom_line(aes(color = HLZ_moist), linewidth = 1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_blank()) +
  facet_wrap(~depth_cat) +
  scale_color_manual("Moisture group", values = hlz_color_moist) +
  scale_fill_manual("Moisture group", values = hlz_color_moist) +
  scale_x_continuous(expression(paste("Log-scaled concentration in M"[ox])), 
                     expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC content", 
                     expand = c(0,0), limits = c(-5,3))
ggsave(file = paste0("./Output/Figure3_", 
                     Sys.Date(), ".jpeg"), width = 14, height = 7)

emip_al_fe_moist %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_point(data = all_data_lm_moist, aes(x = al_fe_ox, y = carbon, 
                                           color = HLZ_moist),
             shape = 21) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = HLZ_moist), alpha = 0.25) +
  geom_line(aes(color = HLZ_moist), linewidth = 1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_blank(),
        legend.position = "none") +
  facet_wrap(vars(depth_cat, HLZ_moist)) +
  scale_color_manual("Moisture group", values = hlz_color_moist) +
  scale_fill_manual("Moisture group", values = hlz_color_moist) +
  scale_x_continuous(expression(paste("Log-scaled concentration in M"[ox])),
                     breaks = seq(-5,5,2.5)) +
  scale_y_continuous("Predicted log-scaled SOC content",
                     breaks = seq(-5,2.5,2.5))
ggsave(file = paste0("./Output/FigureA10_", 
                     Sys.Date(), ".jpeg"), width = 12, height = 10)

# Prepare data for summary table
all_data_mox_summary <- all_data %>% 
  drop_na(al_fe_ox) %>%
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = al_ox + 0.005,
                carbon = carbon + 0.005,
                fe_ox = fe_ox + 0.005,
                al_fe_ox = al_fe_ox + 0.005) %>% 
  group_by(HLZ_moist, depth_cat) %>% 
  skimr::skim_without_charts(al_fe_ox, carbon) %>% 
  as.data.frame() %>% 
  dplyr::select(skim_variable, HLZ_moist, depth_cat, numeric.p0, numeric.p50, numeric.p100) %>% 
  pivot_wider(names_from = skim_variable, values_from = c(numeric.p0, numeric.p50, numeric.p100))

mox_n_samples <- all_data %>% 
  drop_na(al_fe_ox) %>%
  group_by(HLZ_moist, depth_cat) %>% 
  summarise(n = n())

all_data_mox_summary$HLZ_moist <- factor(all_data_mox_summary$HLZ_moist, 
                                         levels = c("arid", "semiarid", "subhumid", 
                                                    "humid", "perhumid"),  
                                         ordered = TRUE)

em_al_fe_moist_slopes <- em_al_fe_moist$emtrends %>% 
  as.data.frame()

em_al_fe_moist_slopes$HLZ_moist <- factor(em_al_fe_moist_slopes$HLZ_moist, 
                                          levels = c("arid", "semiarid", "subhumid", 
                                                     "humid", "perhumid"),  
                                          ordered = TRUE)

em_al_fe_moist_slopes$depth_cat <- factor(em_al_fe_moist_slopes$depth_cat,
                                          levels = c("0-20cm","20-50cm","50-100cm", "100-200cm"),
                                          ordered = TRUE)

posthoc_mox_summary <- em_al_fe_moist_slopes %>% 
  left_join(all_data_mox_summary) %>% 
  left_join(mox_n_samples) %>% 
  mutate_if(is.numeric, round, 2)

write_csv(x = posthoc_mox_summary, file = paste0("./Output/Table1_", 
                                              Sys.Date(), ".csv"))

## Alox
summary(lmm_al_depth_moist)
anova(lmm_al_depth_moist)
r.squaredGLMM(lmm_al_depth_moist)

em_al_moist <- emtrends(lmm_al_depth_moist, pairwise ~ HLZ_moist | depth_cat, var = "al_ox")

em_al_moist$emtrends
#trend - slope direction
em_al_moist$contrasts
#significant differences between slopes

emip_al_moist <- emmip(lmm_al_depth_moist, HLZ_moist ~ al_ox | depth_cat, cov.reduce = range, 
                       plotit = FALSE, CIs = TRUE)

emip_al_moist$HLZ_moist <- as.ordered(emip_al_moist$HLZ_moist)

emip_al_data_moist <- emip_al_moist %>% 
  as.data.frame() %>% 
  left_join(all_data %>% 
              distinct(HLZ_moist), 
              multiple = "all")

emip_al_data_moist$HLZ_moist <- factor(emip_al_data_moist$HLZ_moist, 
                                 levels = c("arid", "semiarid", "subhumid", "humid", "perhumid"),  
                                 ordered = TRUE)

emip_al_moist %>%   
  ggplot(aes(x = al_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = HLZ_moist), alpha = 0.25) +
  geom_line(aes(color = HLZ_moist), linewidth = 1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_blank()) +
  facet_wrap(~depth_cat) +
  scale_color_manual("Moisture group", values = hlz_color_moist) +
  scale_fill_manual("Moisture group", values = hlz_color_moist) +
  scale_x_continuous(expression(paste("Log-scaled concentration in Al"[ox])), 
                     expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC content", 
                     expand = c(0,0), limits = c(-5,3))
ggsave(file = paste0("./Output/FigureA11_", 
                     Sys.Date(), ".jpeg"), width = 14, height = 7)

## Feox
summary(lmm_fe_depth_moist)
anova(lmm_fe_depth_moist)
r.squaredGLMM(lmm_fe_depth_moist)

em_fe_moist <- emtrends(lmm_fe_depth_moist, pairwise ~ HLZ_moist | depth_cat, var = "fe_ox")

em_fe_moist$emtrends
#trend - slope direction
em_fe_moist$contrasts
#significant differences between slopes

emip_fe_moist <- emmip(lmm_fe_depth_moist, HLZ_moist ~ fe_ox | depth_cat, cov.reduce = range, 
                       plotit = FALSE, CIs = TRUE)

emip_fe_moist$HLZ_moist <- as.ordered(emip_fe_moist$HLZ_moist)

emip_fe_data_moist <- emip_fe_moist %>% 
  as.data.frame() %>% 
  left_join(all_data %>% 
              distinct(HLZ_moist), 
            multiple = "all")

emip_fe_data_moist$HLZ_moist <- factor(emip_fe_data_moist$HLZ_moist, 
                                       levels = c("arid", "semiarid", "subhumid", "humid", "perhumid"),, 
                                       ordered = TRUE)

emip_fe_moist %>%   
  ggplot(aes(x = fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = HLZ_moist), alpha = 0.25) +
  geom_line(aes(color = HLZ_moist), linewidth = 1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_blank()) +
  facet_wrap(~depth_cat) +
  scale_color_manual("Moisture group", values = hlz_color_moist) +
  scale_fill_manual("Moisture group", values = hlz_color_moist) +
  scale_x_continuous(expression(paste("Log-scaled concentration in Fe"[ox])), 
                     expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC content", 
                     expand = c(0,0), limits = c(-5,3))
ggsave(file = paste0("./Output/FigureA12_", 
                     Sys.Date(), ".jpeg"), width = 14, height = 7)
