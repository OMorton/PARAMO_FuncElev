###############################
####  Model interpretation ####
###############################

.libPaths("C:/Packages") ## Set up for working from home.

library(brms)
library(tidybayes)
library(bayestestR)
library(ggpubr)

source("Functions.R")

FRicGAM <- readRDS("Models/DB/FRicGAMwz.rds")
SES.fricGAM <- readRDS("Models/DB/SESFRicGAMwz.rds")

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
  mutate(elev = elev_z*sd(FRic_all$ele_jaxa) + mean(FRic_all$ele_jaxa))

ggplot(FRic1_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FRic_all, aes(ele_jaxa, fric)) +
  labs(y = "FRic", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

ggplot(SESFRic1_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FRic_all, aes(ele_jaxa, SES.fric)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "SES.FRic", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

ggplot(FOri_fit, aes(elev, .epred, colour = habitat)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_point(data = FMulti_all, aes(ele_jaxa, fori)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "SES.FRic", x = "Elevation") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

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

ggplot(diff_sum, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FRic smooth slope") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

ggplot(diff_sum_sesfric, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FRic smooth slope") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

ggplot(diff_sum_fori, aes(elev, diff, colour = habitat)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA, linetype = "dashed") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  ylab("FRic smooth slope") +
  scale_colour_manual(values = c("olivedrab4", "tan3")) +
  theme_minimal(base_size = 14)

#### Contrasts ####

fric_contr <- p_contrasts(fitting_data = FRic_all, model = FRicGAM)
sesfric_contr <- p_contrasts(fitting_data = FRic_all, model = SES.fricGAM)
fori_contr <- p_contrasts(fitting_data = FMulti_all, model = FOriGAM)

ggplot(fric_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FRic") + xlab("Altitude") +
  theme_bw(base_size = 14)

ggplot(sesfric_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FRic") + xlab("Altitude") +
  theme_bw(base_size = 14)

ggplot(fori_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size =1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Difference in FRic") + xlab("Altitude") +
  theme_bw(base_size = 14)
