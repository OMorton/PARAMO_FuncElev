###############################
####  Model interpretation ####
###############################

#.libPaths("C:/Packages") ## Set up for working from home.

library(brms)
library(tidybayes)
library(bayestestR)
library(ggpubr)
library(tidyverse)

## functions script contains a bunch of convenience functions
source("Functions.R")

#### Data ####
SRGAM <- readRDS("Models/DB/SRGAMwfz_st.rds")
FRicGAM <- readRDS("Models/DB/FRicGAMwfz_st.rds")
SES.fricGAM <- readRDS("Models/DB/SESFRicGAMwfz_st.rds")
FOriGAM <- readRDS("Models/DB/FOriGAMwfz_st.rds")
FDisGAM <- readRDS("Models/DB/FDisGAMwfz_st.rds")
FSpeGAM <- readRDS("Models/DB/FSpeGAMwfz_st.rds")

FRic_all<- data.table::fread("Outputs/Summaries/DB/DB_FRic_fitting_data.csv")
FMulti_all <- data.table::fread("Outputs/Summaries/DB/DB_FMulti_fitting_data.csv")

#### First summary plot ####
## Simulate data range
new_fric_dat <- get_newdata(fitting_data = FRic_all)
new_multi_dat <- get_newdata(fitting_data = FMulti_all)

## Get fitted and summarize
SR_fit <- add_predicted_draws(SRGAM, newdata = new_fric_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15))

FRic1_fit <- add_predicted_draws(FRicGAM, newdata = new_fric_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15))

SESFRic1_fit <- add_predicted_draws(SES.fricGAM, newdata = new_fric_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15))

FOri_fit <- add_predicted_draws(FOriGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15)) 

FSpe_fit <- add_predicted_draws(FSpeGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15))

FDis_fit <- add_predicted_draws(FDisGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z, elev) %>% median_hdci(.prediction, .width = .9) %>%
  mutate(elev_plot = ifelse(habitat == "Forest", elev - 15, elev + 15))

## plot
SR_sum_plot <- ggplot(SR_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  geom_point(data = FRic_all, aes(ele_jaxa, sp_richn), alpha = .5, shape = 16) +
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  labs(y = "SR", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5))

fric_sum_plot <- ggplot(FRic1_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  geom_point(data = FRic_all, aes(ele_jaxa, fric), alpha = .5, shape = 16) +
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  labs(y = "FRic", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5))

sesfric_sum_plot <- ggplot(SESFRic1_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FRic_all, aes(ele_jaxa, SES.fric), alpha = .5, shape = 16) +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  labs(y = "SES.FRic", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5))

fori_sum_plot <- ggplot(FOri_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fori), alpha = .5, shape = 16) +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000, 3400),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000, 3400)) +
  labs(y = "FOri", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none",
        , axis.text.x = element_text(angle = 45, vjust = 0.5))

fspe_sum_plot <- ggplot(FSpe_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fspe), alpha = .5, shape = 16) +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000, 3400),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000, 3400)) +
  labs(y = "FSpe", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5))

fdis_sum_plot <- ggplot(FDis_fit, aes(elev_plot, .prediction, colour = habitat)) + 
  #geom_line(size = 1) +
  #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fdis), alpha = .5, shape = 16) +
  geom_point(size = 3, shape = 16) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = 1) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000, 3400),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000, 3400)) +
  labs(y = "FDis", x = "Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.5))

#### First derivatives ####

diff_sum_fric <- f_deriv(data = new_fric_dat, model = FRicGAM,
                    elev_raw = FRic_all$ele_jaxa, prediction = TRUE,
                    eps = 1e-4, summary = TRUE)

diff_sum_sesfric <- f_deriv(data = new_fric_dat, model = SES.fricGAM,
                    elev_raw = FRic_all$ele_jaxa, prediction = TRUE,
                    eps = 1e-4, summary = TRUE)

diff_sum_fori <- f_deriv(data = new_multi_dat, model = FOriGAM,
                            elev_raw = FMulti_all$ele_jaxa, prediction = TRUE,
                            eps = 1e-4, summary = TRUE)

diff_sum_fspe <- f_deriv(data = new_multi_dat, model = FSpeGAM,
                         elev_raw = FMulti_all$ele_jaxa, prediction = TRUE,
                         eps = 1e-4, summary = TRUE)

diff_sum_fdis <- f_deriv(data = new_multi_dat, model = FDisGAM,
                         elev_raw = FMulti_all$ele_jaxa, prediction = TRUE,
                         eps = 1e-4, summary = TRUE)

fric_smooth_plot <- ggplot(diff_sum_fric, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FRic delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

sesfric_smooth_plot <- ggplot(diff_sum_sesfric, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("SES.FRic delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

fori_smooth_plot <- ggplot(diff_sum_fori, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FOri delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

fspe_smooth_plot <- ggplot(diff_sum_fspe, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FSpe delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

fdis_smooth_plot <- ggplot(diff_sum_fdis, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FDis delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

#### Contrasts ####
SR_contr <- p_contrasts(pasture_data = filter(new_fric_dat, habitat == "Pasture"),
                          forest_data = filter(new_fric_dat, habitat == "Forest"), model = SRGAM)
fric_contr <- p_contrasts(pasture_data = filter(new_fric_dat, habitat == "Pasture"),
                          forest_data = filter(new_fric_dat, habitat == "Forest"), model = FRicGAM)
sesfric_contr <- p_contrasts(pasture_data = filter(new_fric_dat, habitat == "Pasture"),
                             forest_data = filter(new_fric_dat, habitat == "Forest"), model = SES.fricGAM)
fori_contr <- p_contrasts(pasture_data = filter(new_multi_dat, habitat == "Pasture"),
                          forest_data = filter(new_multi_dat, habitat == "Forest"), model = FOriGAM)
fspe_contr <- p_contrasts(pasture_data = filter(new_multi_dat, habitat == "Pasture"),
                          forest_data = filter(new_multi_dat, habitat == "Forest"), model = FSpeGAM)
fdis_contr <- p_contrasts(pasture_data = filter(new_multi_dat, habitat == "Pasture"),
                          forest_data = filter(new_multi_dat, habitat == "Forest"), model = FDisGAM)

write.csv(SR_contr, "Outputs/Summaries/Contrasts/DB_SR_diff.csv")
write.csv(fric_contr, "Outputs/Summaries/Contrasts/DB_FRic_diff.csv")
write.csv(sesfric_contr, "Outputs/Summaries/Contrasts/DB_SES.FRic_diff.csv")
write.csv(fori_contr, "Outputs/Summaries/Contrasts/DB_FOri_diff.csv")
write.csv(fspe_contr, "Outputs/Summaries/Contrasts/DB_FSpe_diff.csv")
write.csv(fdis_contr, "Outputs/Summaries/Contrasts/DB_FDis_diff.csv")

SR_contr_plot <- ggplot(SR_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200),
                     labels = c(1000, 1400, 1800, 2200)) +
  ylab("Difference in SR") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))

fric_contr_plot <- ggplot(fric_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200),
                     labels = c(1000, 1400, 1800, 2200)) +
  ylab("Difference in FRic") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
sesfric_contr_plot <- ggplot(sesfric_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200),
                     labels = c(1000, 1400, 1800, 2200)) +
  ylab("Difference in SES.FRic") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
fori_contr_plot <- ggplot(fori_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FOri") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
fspe_contr_plot <- ggplot(fspe_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FSpe") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
fdis_contr_plot <- ggplot(fdis_contr, aes(elev, Contr)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FDis") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
