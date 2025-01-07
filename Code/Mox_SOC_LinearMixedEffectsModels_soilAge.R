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

#### Linear mixed-effects models: soil age ####

## Only focus on cool temperate where we have enough data and young vs old soil age

## Load and prepare data
all_data <- read_csv("./Data/Database_HLZ_grp_2024-09-09.csv") %>% 
  mutate(al_fe_ox = al_ox + (1/2*fe_ox)) %>% 
  filter(hzn_bot <= 200) %>% 
  mutate(depth_cat = cut(hzn_mid, breaks = c(0, 20, 50, 100, 200), 
                         labels = c("0-20cm","20-50cm","50-100cm", "100-200cm"),
                         ordered = TRUE)) 

all_data %>% 
  count(soil_age)

all_data %>% 
  count(HLZ_moist)

cool_data <- all_data %>% 
  filter(HLZ_temp == "cool temperate") %>% 
  #drop arid and semiarid (n = 4 for young), for which we have not enough data
  filter(HLZ_moist != "arid",
         HLZ_moist != "semiarid")

cool_data %>% 
  group_by(HLZ_moist) %>% 
  count(soil_age)

#### Linear mixed-effects models
### Fit model - cool temperate all
cool_data_lm <- cool_data %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005))
## Mox
lmm_cool <- lme(carbon ~ al_fe_ox*soil_age + al_fe_ox*depth_cat + 
                  al_fe_ox*soil_age*depth_cat + al_fe_ox*HLZ_moist + al_fe_ox*depth_cat + 
                  al_fe_ox*HLZ_moist*depth_cat, 
                random = ~hzn_mid|ID, 
                data = cool_data_lm)

summary(lmm_cool)
r.squaredGLMM(lmm_cool)
anova(lmm_cool)
sqrt(mean(residuals(lmm_cool)^2))

lmm_cool_hlz <- lme(carbon ~ al_fe_ox*HLZ_moist + al_fe_ox*depth_cat + 
                      al_fe_ox*HLZ_moist*depth_cat, 
                    random = ~hzn_mid|ID, 
                    data = cool_data_lm)

summary(lmm_cool_hlz)
r.squaredGLMM(lmm_cool_hlz)
anova(lmm_cool_hlz)
sqrt(mean(residuals(lmm_cool)^2))

#### Linear mixed-effects models
### Fit model - humid cool
cool_data_lm_humid <- cool_data %>% 
  filter(HLZ_moist == "humid") %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005))
## Mox
lmm_cool_humid <- lme(carbon ~ al_fe_ox*soil_age + al_fe_ox*depth_cat + 
                        al_fe_ox*soil_age*depth_cat, 
                      random = ~hzn_mid|ID, 
                      data = cool_data_lm_humid)

summary(lmm_cool_humid)
r.squaredGLMM(lmm_cool_humid)
anova(lmm_cool_humid)
sqrt(mean(residuals(lmm_cool_humid)^2))

# Graphical exploration of residuals
plot(lmm_cool_humid)

F_Final <- fitted(lmm_cool_humid)
R_Final <- residuals(lmm_cool_humid) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ cool_data_lm_humid$soil_age)
boxplot(R_Final ~ cool_data_lm_humid$depth_cat)
plot(R_Final ~ cool_data_lm_humid$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_cool_humid)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * cool_data_lm_humid$al_fe_ox

{scatter.smooth(cool_data_lm_humid$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*cool_data_lm_humid$al_fe_ox ~ cool_data_lm_humid$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous

# Check for spatial correlation
resid_norm <- resid(lmm_cool_humid, type = "normalized")
data_humid <- data.frame(resid_norm, cool_data_lm_humid$Longitude, 
                         cool_data_lm_humid$Latitude)
coordinates(data_humid) <- c("cool_data_lm_humid.Longitude",
                             "cool_data_lm_humid.Latitude")
bubble(data_humid, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario1 <- variogram(resid_norm ~ 1, data_humid)
plot(Vario1)

### Post-Hoc analysis
## Mox
em_al_fe_humid <- emtrends(lmm_cool_humid, pairwise ~ soil_age | depth_cat, 
                           var = "al_fe_ox")

em_al_fe_humid$emtrends
#trend - slope direction
em_al_fe_humid$contrasts
#significant differences between slopes

emip_al_fe_humid <- emmip(lmm_cool_humid, soil_age ~ al_fe_ox | depth_cat, 
                          cov.reduce = range, 
                          plotit = FALSE, CIs = TRUE)

em_al_fe_humid_data <- emip_al_fe_humid %>% 
  as.data.frame() %>% 
  left_join(cool_data_lm_humid %>% 
              droplevels() %>% 
              distinct(soil_age), 
            multiple = "all")

em_al_fe_humid_data %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = soil_age), alpha = 0.25) +
  geom_line(aes(color = soil_age), linewidth = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(vars(depth_cat)) +
  scale_color_discrete("Soil age") +
  scale_fill_discrete("Soil age") +
  scale_x_continuous(expression(paste("Log-scaled concentration in M"[ox])), expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC", expand = c(0,0), limits = c(-7,3))

### Fit model - subhumid cool
cool_data_lm_subhumid <- cool_data %>% 
  filter(HLZ_moist == "subhumid") %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005))
## Mox
lmm_cool_subhumid <- lme(carbon ~ al_fe_ox*soil_age + al_fe_ox*depth_cat + 
                           al_fe_ox*soil_age*depth_cat, 
                         random = ~hzn_mid|ID, 
                         data = cool_data_lm_subhumid)

summary(lmm_cool_subhumid)
r.squaredGLMM(lmm_cool_subhumid)
anova(lmm_cool_subhumid)
sqrt(mean(residuals(lmm_cool_subhumid)^2))

# Graphical exploration of residuals
plot(lmm_cool_subhumid)

F_Final <- fitted(lmm_cool_subhumid)
R_Final <- residuals(lmm_cool_subhumid) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ cool_data_lm_subhumid$soil_age)
boxplot(R_Final ~ cool_data_lm_subhumid$depth_cat)
plot(R_Final ~ cool_data_lm_subhumid$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_cool_subhumid)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * cool_data_lm_subhumid$al_fe_ox

{scatter.smooth(cool_data_lm_subhumid$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*cool_data_lm_subhumid$al_fe_ox ~ cool_data_lm_subhumid$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous

# Check for spatial correlation
resid_norm <- resid(lmm_cool_subhumid, type = "normalized")
data_subhumid <- data.frame(resid_norm, cool_data_lm_subhumid$Longitude, 
                         cool_data_lm_subhumid$Latitude)
coordinates(data_subhumid) <- c("cool_data_lm_subhumid.Longitude",
                             "cool_data_lm_subhumid.Latitude")

bubble(data_subhumid, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario1 <- variogram(resid_norm ~ 1, data_subhumid)
plot(Vario1)

#maybe not ideal?

### Post-Hoc analysis
## Mox
em_al_fe_subhumid <- emtrends(lmm_cool_subhumid, pairwise ~ soil_age | depth_cat, 
                           var = "al_fe_ox")

em_al_fe_subhumid$emtrends
#trend - slope direction
em_al_fe_subhumid$contrasts
#significant differences between slopes

emip_al_fe_subhumid <- emmip(lmm_cool_subhumid, soil_age ~ al_fe_ox | depth_cat, 
                          cov.reduce = range, 
                          plotit = FALSE, CIs = TRUE)

em_al_fe_subhumid_data <- emip_al_fe_subhumid %>% 
  as.data.frame() %>% 
  left_join(cool_data_lm_subhumid %>% 
              droplevels() %>% 
              distinct(soil_age), 
            multiple = "all")

em_al_fe_subhumid_data %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = soil_age), alpha = 0.25) +
  geom_line(aes(color = soil_age), linewidth = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(vars(depth_cat)) +
  scale_color_discrete("Soil age") +
  scale_fill_discrete("Soil age") +
  scale_x_continuous(expression(paste("Log-scaled concentration in M"[ox])), expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC", expand = c(0,0), limits = c(-7,3))

### Fit model - perhumid cool
cool_data_lm_perhumid <- cool_data %>% 
  filter(HLZ_moist == "perhumid") %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005))
## Mox
lmm_cool_perhumid <- lme(carbon ~ al_fe_ox*soil_age + al_fe_ox*depth_cat + 
                           al_fe_ox*soil_age*depth_cat, 
                         random = ~hzn_mid|ID, 
                         data = cool_data_lm_perhumid)

summary(lmm_cool_perhumid)
r.squaredGLMM(lmm_cool_perhumid)
anova(lmm_cool_perhumid)
sqrt(mean(residuals(lmm_cool_perhumid)^2))

# Graphical exploration of residuals
plot(lmm_cool_perhumid)

F_Final <- fitted(lmm_cool_perhumid)
R_Final <- residuals(lmm_cool_perhumid) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ cool_data_lm_perhumid$soil_age)
boxplot(R_Final ~ cool_data_lm_perhumid$depth_cat)
plot(R_Final ~ cool_data_lm_perhumid$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_cool_perhumid)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * cool_data_lm_perhumid$al_fe_ox

{scatter.smooth(cool_data_lm_perhumid$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*cool_data_lm_perhumid$al_fe_ox ~ cool_data_lm_perhumid$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous

# Check for spatial correlation
resid_norm <- resid(lmm_cool_perhumid, type = "normalized")
data_perhumid <- data.frame(resid_norm, cool_data_lm_perhumid$Longitude, 
                            cool_data_lm_perhumid$Latitude)
coordinates(data_perhumid) <- c("cool_data_lm_perhumid.Longitude",
                                "cool_data_lm_perhumid.Latitude")

bubble(data_perhumid, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario1 <- variogram(resid_norm ~ 1, data_perhumid)
plot(Vario1)

### Post-Hoc analysis
## Mox
em_al_fe_perhumid <- emtrends(lmm_cool_perhumid, pairwise ~ soil_age | depth_cat, 
                              var = "al_fe_ox")

em_al_fe_perhumid$emtrends
#trend - slope direction
em_al_fe_perhumid$contrasts
#significant differences between slopes

emip_al_fe_perhumid <- emmip(lmm_cool_perhumid, soil_age ~ al_fe_ox | depth_cat, 
                             cov.reduce = range, 
                             plotit = FALSE, CIs = TRUE)

em_al_fe_perhumid_data <- emip_al_fe_perhumid %>% 
  as.data.frame() %>% 
  left_join(cool_data_lm_perhumid %>% 
              droplevels() %>% 
              distinct(soil_age), 
            multiple = "all")

em_al_fe_perhumid_data %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = soil_age), alpha = 0.25) +
  geom_line(aes(color = soil_age), linewidth = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(vars(depth_cat)) +
  scale_color_discrete("Soil age") +
  scale_fill_discrete("Soil age") +
  scale_x_continuous(expression(paste("Log-scaled concentration in M"[ox])), expand = c(0,0)) +
  scale_y_continuous("Predicted log-scaled SOC", expand = c(0,0), limits = c(-7,3))

