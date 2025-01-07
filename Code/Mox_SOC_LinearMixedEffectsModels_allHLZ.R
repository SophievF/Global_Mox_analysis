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

#### Linear mixed-effects models: all HLZ ####

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

all_data$HLZ_new <- factor(all_data$HLZ_new, 
                           levels = unique(all_data$HLZ_new[order(all_data$ZONE_filled)]), 
                           ordered = TRUE)

all_data$HLZ_temp <- factor(all_data$HLZ_temp, 
                            levels = unique(all_data$HLZ_temp[
                              order(all_data$ZONE_filled)]), 
                            ordered = TRUE)

all_data$HLZ_moist <- factor(all_data$HLZ_moist, 
                             levels = c("arid", "semiarid", "subhumid", "humid", 
                                        "perhumid"), 
                             ordered = TRUE)

#### Linear mixed-effects models
### Possible auto-correlation
all_data %>% 
  drop_na(al_ox) %>% 
  ggplot(aes(x = al_ox, y = hzn_mid)) +
  geom_point()

cor.test(all_data$al_ox, all_data$hzn_mid, method = "spearman") 

all_data %>% 
  drop_na(al_ox) %>% 
  ggplot(aes(x = fe_ox, y = hzn_mid)) +
  geom_point()

cor.test(all_data$fe_ox, all_data$hzn_mid, method = "spearman") 

all_data %>% 
  drop_na(al_ox) %>% 
  ggplot(aes(x = al_ox, y = fe_ox)) +
  geom_point()

cor.test(all_data$al_ox, all_data$fe_ox, method = "spearman") 

### Prepare data for model w/o HLZ 
all_data_lm_depth <- all_data %>% 
  drop_na(al_ox, fe_ox) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005)) 

### Fit model
## Mox
# Random slope for depth
lmm_al_fe_depth_wo_HLZ <- lme(carbon ~ al_fe_ox + al_fe_ox*depth_cat, 
                              random = ~hzn_mid|ID,
                              data = all_data_lm_depth)

summary(lmm_al_fe_depth_wo_HLZ)
r.squaredGLMM(lmm_al_fe_depth_wo_HLZ)
anova(lmm_al_fe_depth_wo_HLZ)
sqrt(mean(residuals(lmm_al_fe_depth_wo_HLZ)^2))

LMM_df <- anova(lmm_al_fe_depth_wo_HLZ) %>% 
  base::as.data.frame()

LMM_df$predictor <- rownames(LMM_df)

LMM_df

write_csv(x = LMM_df, file = paste0("./Output/TableA9_", 
                                    Sys.Date(), ".csv"))

# Graphical exploration of residuals
plot(lmm_al_fe_depth_wo_HLZ)

F_Final <- fitted(lmm_al_fe_depth_wo_HLZ)
R_Final <- residuals(lmm_al_fe_depth_wo_HLZ, type = "pearson", scaled = TRUE) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ all_data_lm_depth$depth_cat)
plot(R_Final ~ all_data_lm_depth$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_al_fe_depth_wo_HLZ)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * all_data_lm_depth$al_fe_ox

{scatter.smooth(all_data_lm_depth$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*all_data_lm_depth$al_fe_ox ~ all_data_lm_depth$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous

# Check for spatial correlation
resid_norm <- resid(lmm_al_fe_depth_wo_HLZ, type = "normalized")
data_noHLZ <- data.frame(resid_norm, all_data_lm_depth$Longitude, 
                         all_data_lm_depth$Latitude)
coordinates(data_noHLZ) <- c("all_data_lm_depth.Longitude",
                            "all_data_lm_depth.Latitude")

bubble_p <- bubble(data_noHLZ, "resid_norm", col = c("black", "grey"),
                   main = "Normalised residuals",
                   xlab = "X-coordinates", ylab = "Y-coordinates")

jpeg(filename = paste0("./Output/FigureA2_", Sys.Date(), ".jpeg"),
     width = 880, height = 480)
plot(bubble_p)
dev.off()

Vario1 <- variogram(resid_norm ~ 1, data_noHLZ)
jpeg(filename = paste0("./Output/FigureA1_", Sys.Date(), ".jpeg"),
     width = 880, height = 480)
plot(Vario1)
dev.off()

### Post-Hoc analysis
## Mox
em_al_fe_depth <- emtrends(lmm_al_fe_depth_wo_HLZ, "depth_cat", var = "al_fe_ox")

pairs(em_al_fe_depth)
em_al_fe_depth

em_al_fe_depth %>% 
  as.data.frame() %>% 
  ggplot(aes(x = depth_cat, y = al_fe_ox.trend)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Slope: SOC ~ Mox", expand = c(0,0), limits = c(0,0.6)) +
  scale_x_discrete("Depth")

em_al_fe_depth_df <- em_al_fe_depth %>% 
  as.data.frame()

### Prepare data for model w/ HLZ
# Check for zones with too little data (<= 10)
rm_HLZ <- all_data %>% 
  group_by(DESC_filled, depth_cat) %>% 
  count() %>% 
  filter(n <= 10) %>% 
  ungroup() %>% 
  count(DESC_filled)

all_data_lm_HLZ <- all_data %>% 
  drop_na(al_ox, fe_ox) %>% 
  #remove HLZ with too little data
  filter(!DESC_filled %in% rm_HLZ$DESC_filled) %>% 
  # set 0s to smallest measured value/2 = 0.005
  dplyr::mutate(al_ox = log(al_ox + 0.005),
                carbon = log(carbon + 0.005),
                fe_ox = log(fe_ox + 0.005),
                al_fe_ox = log(al_fe_ox + 0.005)) 

skimr::skim_without_charts(all_data_lm_HLZ)

all_data_lm_HLZ %>% 
  count(HLZ_new, depth_cat) %>% view()

## Update Color Scheme for plotting
polar <- brewer.pal(4, "Greys")
boreal <- brewer.pal(4, "Blues")
cool_temp <- brewer.pal(6, "Purples")
warm_temp <- brewer.pal(5, "Greens")
subtrop <- brewer.pal(3, "Oranges")
trop <- brewer.pal(3, "RdPu")

hlz_color_red <- c(polar, boreal, cool_temp, warm_temp, subtrop, trop)

### Fit model
## Mox
# Random slope for depth
lmm_al_fe_HLZ_depth <- lme(carbon ~ al_fe_ox*HLZ_new + al_fe_ox*depth_cat + 
                             al_fe_ox*HLZ_new*depth_cat, 
                           random = ~hzn_mid|ID,
                           data = all_data_lm_HLZ)

summary(lmm_al_fe_HLZ_depth)
r.squaredGLMM(lmm_al_fe_HLZ_depth)
anova(lmm_al_fe_HLZ_depth)
sqrt(mean(residuals(lmm_al_fe_HLZ_depth)^2))

LMM_HLZ_df <- anova(lmm_al_fe_HLZ_depth) %>% 
  base::as.data.frame()

LMM_HLZ_df$predictor <- rownames(LMM_HLZ_df)

LMM_HLZ_df

write_csv(x = LMM_HLZ_df, file = paste0("./Output/TableA8_", 
                                        Sys.Date(), ".csv"))

# Graphical exploration of residuals
plot(lmm_al_fe_HLZ_depth)

F_Final <- fitted(lmm_al_fe_HLZ_depth)
R_Final <- residuals(lmm_al_fe_HLZ_depth) 
#no difference in graphs with pearson vs response

hist(R_Final)
boxplot(R_Final ~ all_data_lm_HLZ$HLZ_new)
boxplot(R_Final ~ all_data_lm_HLZ$depth_cat)
plot(R_Final ~ all_data_lm_HLZ$al_fe_ox)

# Partial residual plots
mox_c <- fixef(lmm_al_fe_HLZ_depth)[2] #predictor coefficient
mox_pr <- R_Final + mox_c * all_data_lm_HLZ$al_fe_ox

{scatter.smooth(all_data_lm_HLZ$al_fe_ox, mox_pr,
                lpars = list(col = "green", lwd = 3, lty = 3)) #residual loess
  abline(lm(mox_c*all_data_lm_HLZ$al_fe_ox ~ all_data_lm_HLZ$al_fe_ox), col = "red")}

#--> looks all good: similar variance across groups and generally homogeneous


# Check for spatial correlation
resid_norm <- resid(lmm_al_fe_HLZ_depth, type = "normalized")
data_HLZ <- data.frame(resid_norm, all_data_lm_HLZ$Longitude, 
                       all_data_lm_HLZ$Latitude)
coordinates(data_HLZ) <- c("all_data_lm_HLZ.Longitude",
                           "all_data_lm_HLZ.Latitude")
bubble(data_HLZ, "resid_norm", col = c("black", "grey"),
       main = "Normalised residuals",
       xlab = "X-coordinates", ylab = "Y-coordinates")

Vario1 <- variogram(resid_norm ~ 1, data_HLZ)
plot(Vario1)

### Post-Hoc analysis
## Mox
em_al_fe_HLZ <- emtrends(lmm_al_fe_HLZ_depth, pairwise ~ HLZ_new | depth_cat, 
                         var = "al_fe_ox")

em_al_fe_HLZ$emtrends
#trend - slope direction
em_al_fe_HLZ$contrasts
#significant differences between slopes

emip_al_fe_HLZ <- emmip(lmm_al_fe_HLZ_depth, HLZ_new ~ al_fe_ox | depth_cat, 
                        cov.reduce = range, 
                        plotit = FALSE, CIs = TRUE)

emip_al_fe_HLZ$HLZ_new <- as.ordered(emip_al_fe_HLZ$HLZ_new)

em_al_fe_HLZ_data <- emip_al_fe_HLZ %>% 
  as.data.frame() %>% 
  left_join(all_data_lm_HLZ %>% 
              droplevels() %>% 
              distinct(DESC_filled, ZONE_filled, HLZ_moist, HLZ_temp, HLZ_new), 
              multiple = "all")

em_al_fe_HLZ_data$HLZ_new <- factor(em_al_fe_HLZ_data$HLZ_new, 
                                      levels = unique(em_al_fe_HLZ_data$HLZ_new[
                                        order(em_al_fe_HLZ_data$ZONE_filled)]), 
                                      ordered = TRUE)

em_al_fe_HLZ_data %>%   
  ggplot(aes(x = al_fe_ox, y = yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = HLZ_new), alpha = 0.25) +
  geom_line(aes(color = HLZ_new), linewidth = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(vars(depth_cat)) +
  scale_color_manual("HLZ zones", values = hlz_color_red) +
  scale_fill_manual("HLZ zones", values = hlz_color_red) +
  scale_x_continuous(expression(paste("M"[ox])), expand = c(0,0)) +
  scale_y_continuous("Linear predicted SOC", expand = c(0,0))


em_al_fe_HLZ_trends <- em_al_fe_HLZ$emtrends %>% 
  data.frame()

em_al_fe_HLZ_trends$HLZ_new <- as.ordered(em_al_fe_HLZ_trends$HLZ_new)

em_al_fe_HLZ_trends %>% 
  left_join(all_data_lm_HLZ %>% 
              droplevels() %>% 
              distinct(DESC_filled, ZONE_filled, HLZ_moist, HLZ_temp, HLZ_new), 
            multiple = "all") %>% 
  ggplot(aes(x = HLZ_new, y = al_fe_ox.trend)) +
  geom_hline(data = em_al_fe_depth_df, mapping = aes(yintercept = al_fe_ox.trend)) +
  geom_hline(data = em_al_fe_depth_df, mapping = aes(yintercept = lower.CL),
             linetype = "dashed") +
  geom_hline(data = em_al_fe_depth_df, mapping = aes(yintercept = upper.CL),
             linetype = "dashed") +
  geom_vline(xintercept = c(4.5, 8.5, 14.5, 19.4, 22.6), linewidth = 0.25) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point(aes(fill = HLZ_new), size = 3, shape = 21) +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  facet_wrap(~depth_cat) +
  scale_y_continuous(expression(paste("Linear predicted slope: SOC ~ M"[ox])), expand = c(0,0),
                     limits = c(-0.4,1.2), breaks = seq(-0.4,1.2,0.4)) +
  scale_fill_manual("HLZ zones", values = hlz_color_red) +
  annotate(geom = "text", y = 1.1, x = 2.5, label = "(sub-)polar", size = 4,
           fontface = "italic") +
  annotate(geom = "text", y = 1.1, x = 6.5, label = "boreal", size = 4,
         fontface = "italic") +
  annotate(geom = "text", y = 1.1, x = 11.5, label = "cool temperate", size = 4,
           fontface = "italic") +
  annotate(geom = "text", y = 1.1, x = 17, label = "warm temperate", size = 4,
           fontface = "italic") +
  annotate(geom = "text", y = 1.1, x = 21, label = "subtropical", size = 4,
           fontface = "italic") +
  annotate(geom = "text", y = 1.1, x = 24, label = "tropical", size = 4,
           fontface = "italic")
ggsave(file = paste0("./Output/Figure5_", 
                     Sys.Date(), ".jpeg"), width = 14, height = 9)


