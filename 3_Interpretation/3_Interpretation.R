###############################
####  Model interpretation ####
###############################

.libPaths("C:/Packages") ## Set up for working from home.

library(brms)
library(tidybayes)
library(bayestestR)
library(ggpubr)
library(tidyverse)

## functions script contains a bunch of convenience functions
source("Functions.R")

#### Data ####
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
        axis.ticks = element_line(), legend.position = "none")

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
        axis.ticks = element_line(), legend.position = "none")

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
        axis.ticks = element_line(), legend.position = "none")

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
        axis.ticks = element_line(), legend.position = "none")

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
        axis.ticks = element_line(), legend.position = "none")

#### First derivatives ####

diff_sum_fric <- f_deriv(data = new_dat, model = FRicGAM,
                    elev_raw = FRic_all$ele_jaxa,
                    eps = 1e-4, summary = FALSE)

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

fric_contr_plot <- ggplot(fric_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200),
                     labels = c(1000, 1400, 1800, 2200)) +
  ylab("Difference in FRic") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line())
  
sesfric_contr_plot <- ggplot(sesfric_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200),
                     labels = c(1000, 1400, 1800, 2200)) +
  ylab("Difference in SES.FRic") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line())
  
fori_contr_plot <- ggplot(fori_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FOri") +
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line())
  
fspe_contr_plot <- ggplot(fspe_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FSpe") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line())
  
fdis_contr_plot <- ggplot(fdis_contr, aes(elev, Contr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
  scale_x_continuous(breaks = c(1000, 1400, 1800, 2200, 2600, 3000),
                     labels = c(1000, 1400, 1800, 2200, 2600, 3000)) +
  ylab("Difference in FDis") + 
  xlab("Elevation (m)") +
  theme_minimal(base_size = 12) +
  theme(axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.ticks = element_line())
  
#### Arrangement ####
library(ggpubr)
library(png)
library(grid)

## Get images
DB <-  rasterGrob(readPNG("Misc/DB_image.png"), interpolate = TRUE)

FR_arrange <- ggarrange(fric_sum_plot, fric_contr_plot, 
          sesfric_sum_plot, sesfric_contr_plot, ncol = 2, nrow = 2,
          labels = c("A.", "B.", "C.", "D."), widths = c(1, 0.5))

FMulti_arrange <- ggarrange(fori_sum_plot, fori_contr_plot, 
                            fspe_sum_plot, fspe_contr_plot, 
                            fdis_sum_plot, fdis_contr_plot, ncol = 2, nrow = 3,
                        labels = c("A.", "B.", "C.", "D.", "E.", "F."), 
                        widths = c(1, 0.5))


FR_arrange2 <- FR_arrange + annotation_custom(DB, xmin = 0.92, xmax = 0.99, ymin = 0.90, ymax = .97)
Fmulti_arrange2 <- FMulti_arrange + annotation_custom(DB, xmin = 0.92, xmax = 0.99, ymin = 0.90, ymax = .97)

ggsave(path = "Outputs/Figures/DB", FR_arrange2, filename = "FR_arrangement.png",  bg = "white",
       device = "png", width = 25, height = 15, units = "cm")

ggsave(path = "Outputs/Figures/DB", Fmulti_arrange2, filename = "FMulti_arrangement.png",  bg = "white",
       device = "png", width = 25, height = 20, units = "cm")
