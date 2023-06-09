###############################
####  Model interpretation ####
###############################

.libPaths("C:/Packages") ## Set up for working from home.

library(brms)
library(tidybayes)
library(bayestestR)
library(ggpubr)
library(tidyverse)

source("Functions.R")

#### Data ####
FRicGAM <- readRDS("Models/DB/FRicGAMwz.rds")
SES.fricGAM <- readRDS("Models/DB/SESFRicGAMwz.rds")
FOriGAM <- readRDS("Models/DB/FOriGAMwz.rds")
FDisGAM <- readRDS("Models/DB/FDisGAMwz.rds")
FSpeGAM <- readRDS("Models/DB/FSpeGAMwz.rds")

FRic_all<- data.table::fread("Outputs/Summaries/DB/DB_FRic_fitting_data.csv")
FMulti_all <- data.table::fread("Outputs/Summaries/DB/DB_FMulti_fitting_data.csv")

#### First summary plot ####
## Simulate data range
new_dat <- FRic_all %>% group_by(habitat) %>% 
  reframe(elev_z = seq(from = min(elev_z), to = max(elev_z),length.out = 25),
          Wood_Veg_z = mean(Wood_Veg_z),
          Cluster_dummy = 0, cluster = "cluster_MOF_1")

new_multi_dat <- FMulti_all %>% group_by(habitat) %>% 
  reframe(elev_z = seq(from = min(elev_z), to = max(elev_z),length.out = 25),
          Wood_Veg_z = mean(Wood_Veg_z),
          Cluster_dummy = 0, cluster = "cluster_MOF_1")

## Get fitted and summarize
FRic1_fit <- add_epred_draws(FRicGAM, newdata = new_dat) %>% 
  group_by(habitat, elev_z) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(elev = elev_z*sd(FRic_all$ele_jaxa) + mean(FRic_all$ele_jaxa))

SESFRic1_fit <- add_epred_draws(SES.fricGAM, newdata = new_dat) %>% 
  group_by(habitat, elev_z) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(elev = elev_z*sd(FRic_all$ele_jaxa) + mean(FRic_all$ele_jaxa))

FOri_fit <- add_epred_draws(FOriGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(elev = elev_z*sd(FMulti_all$ele_jaxa) + mean(FMulti_all$ele_jaxa))

FSpe_fit <- add_epred_draws(FSpeGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(elev = elev_z*sd(FMulti_all$ele_jaxa) + mean(FMulti_all$ele_jaxa))

FDis_fit <- add_epred_draws(FDisGAM, newdata = new_multi_dat) %>% 
  group_by(habitat, elev_z) %>% median_hdci(.epred, .width = .9) %>% 
  mutate(elev = elev_z*sd(FMulti_all$ele_jaxa) + mean(FMulti_all$ele_jaxa))

## plot
fric_sum_plot <- ggplot(FRic1_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FRic_all, aes(ele_jaxa, fric)) +
  labs(y = "FRic", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

sesfric_sum_plot <- ggplot(SESFRic1_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FRic_all, aes(ele_jaxa, SES.fric)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "SES.FRic", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fori_sum_plot <- ggplot(FOri_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fori)) +
  labs(y = "FOri", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fspe_sum_plot <- ggplot(FSpe_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fspe)) +
  labs(y = "FSpe", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fdis_sum_plot <- ggplot(FDis_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fdis)) +
  labs(y = "FDis", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

#### First derivatives ####

diff_sum_fric <- f_deriv(data = new_dat, model = FRicGAM,
                    elev_raw = FRic_all$ele_jaxa,
                    eps = 1e-4, summary = TRUE)

diff_sum_sesfric <- f_deriv(data = new_dat, model = SES.fricGAM,
                    elev_raw = FRic_all$ele_jaxa,
                    eps = 1e-4, summary = TRUE)

diff_sum_fori <- f_deriv(data = new_multi_dat, model = FOriGAM,
                            elev_raw = FMulti_all$ele_jaxa,
                            eps = 1e-4, summary = TRUE)

diff_sum_fspe <- f_deriv(data = new_multi_dat, model = FSpeGAM,
                         elev_raw = FMulti_all$ele_jaxa,
                         eps = 1e-4, summary = TRUE)

diff_sum_fdis <- f_deriv(data = new_multi_dat, model = FDisGAM,
                         elev_raw = FMulti_all$ele_jaxa,
                         eps = 1e-4, summary = TRUE)

fric_smooth_plot <- ggplot(diff_sum_fric, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FRic delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

sesfric_smooth_plot <- ggplot(diff_sum_sesfric, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("SES.FRic delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fori_smooth_plot <- ggplot(diff_sum_fori, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FOri delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fspe_smooth_plot <- ggplot(diff_sum_fspe, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FSpe delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

fdis_smooth_plot <- ggplot(diff_sum_fdis, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FDis delta (p/100m)") +
  xlab("Elevation (m)") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

#### Contrasts ####

fric_contr <- p_contrasts(fitting_data = FRic_all, model = FRicGAM)
sesfric_contr <- p_contrasts(fitting_data = FRic_all, model = SES.fricGAM)
fori_contr <- p_contrasts(fitting_data = FMulti_all, model = FOriGAM)
fspe_contr <- p_contrasts(fitting_data = FMulti_all, model = FSpeGAM)
fdis_contr <- p_contrasts(fitting_data = FMulti_all, model = FDisGAM)

fric_contr_plot <- ggplot(fric_contr, aes(elev, Contr)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FRic") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 10)

sesfric_contr_plot <- ggplot(sesfric_contr, aes(elev, Contr)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in SES.FRic") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 10)

fori_contr_plot <- ggplot(fori_contr, aes(elev, Contr)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FOri") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 10)

fspe_contr_plot <- ggplot(fspe_contr, aes(elev, Contr)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FSpe") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 10)

fdis_contr_plot <- ggplot(fdis_contr, aes(elev, Contr)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FDis") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 10)
