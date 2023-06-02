#####################################
####  GAMs of functional metrics ####
#####################################

## Preliminary exploration and analysis of functional richness across both elevations and land-uses. 

## Set up packages for working from home
options(scipen=999)
.libPaths("C:/Packages") 

library(brms)
library(tidyverse)
library(tidybayes)

#### Data ####
## FRic is scaled by the convhull of the complete species pool so is bounded by
## 1 and 0.
FRic_all_raw <- read.csv("Outputs/Summaries/DB/DB_FRic_Point_all_NEW.csv") %>%
  select(point, fric, SES.fric)

FMulti_all_raw <- read.csv("Outputs/Summaries/DB/DB_Multi_Point_all_NEW.csv") %>%
  select(point, fdis, fori, fspe, SES.fdis, SES.fori, SES.fspe)

## Read in the distance to forest edge
Landscape_metrics <- read.csv("Data/Point_info/ECpts_landscapemetricsMASTER_v2.csv")

## bind and remove the two reference levels for paramo and forest at the high altitudes.
FRic_all <- left_join(FRic_all_raw, Landscape_metrics)
FMulti_all <- left_join(FMulti_all_raw, Landscape_metrics)

#### Preparation ####
## Standardize
FRic_all <- FRic_all %>% filter(habitat != "Paramo") %>% 
  mutate(elev_z = (ele_jaxa - mean(ele_jaxa)) / sd(ele_jaxa),
         For_Dist_z = (dist_forest_edge - mean(dist_forest_edge)) / sd(dist_forest_edge),
         Wood_Veg_z = (woodyveg_index - mean(woodyveg_index)) / sd(woodyveg_index),
         Cluster_dummy = 1)

FMulti_all <- FMulti_all %>% filter(habitat != "Paramo") %>% 
  mutate(elev_z = (ele_jaxa - mean(ele_jaxa)) / sd(ele_jaxa),
         For_Dist_z = (dist_forest_edge - mean(dist_forest_edge)) / sd(dist_forest_edge),
         Wood_Veg_z = (woodyveg_index - mean(woodyveg_index)) / sd(woodyveg_index),
         Cluster_dummy = 1)

ggplot(FRic_all, aes(elev_z, fric, colour = habitat)) + geom_point() +
  geom_smooth()
ggplot(FRic_all, aes(Wood_Veg_z, fric, colour = habitat)) + geom_point() +
  geom_smooth()
ggplot(FRic_all, aes(For_Dist_z, fric, colour = habitat)) + geom_point() +
  geom_smooth()

#### FRic Model fit ####

## basic model of focal variables and samplng structure
FRicGAM <- brm(fric ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                 s(cluster, bs="re", by= Cluster_dummy),
               family = Beta(),
               data = FRic_all,
               sample_prior = TRUE,  prior = c(
                 prior(normal(0,1), "b"),
                 prior(normal(0,1), "Intercept")),
               control = list(adapt_delta = .95),
               file = "Models/DB/FRicGAM.rds",
               chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

## more complete model incorporating that increasing wood in pasture sites 
## may covary with species composition
FRicGAM <- brm(fric ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                 Wood_Veg_z:habitat + Wood_Veg_z +
                 s(cluster, bs="re", by= Cluster_dummy),
               family = Beta(),
               data = FRic_all,
               sample_prior = TRUE,  prior = c(
                 prior(normal(0,1), "b"),
                 prior(normal(0,1), "Intercept")),
               control = list(adapt_delta = .95),
               file = "Models/DB/FRicGAMwz.rds",
               chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

model_check(data = FRic_all, model = FRicGAM)



#### SES FRic Model fit ####

## basic model of focal variables and samplng structure
SES.fricGAM <- brm(SES.fric ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                 s(cluster, bs="re", by= Cluster_dummy),
               family = gaussian(),
               data = FRic_all,
               sample_prior = TRUE,  prior = c(
                 prior(normal(0,1), "b"),
                 prior(normal(0,1), "Intercept")),
               control = list(adapt_delta = .95),
               file = "Models/DB/SESFRicGAM.rds",
               chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

## more complete model incorporating that increasing wood in pasture sites 
## may covary with species composition
SES.fricGAM <- brm(SES.fric ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                 Wood_Veg_z:habitat + Wood_Veg_z +
                 s(cluster, bs="re", by= Cluster_dummy),
               family = gaussian(),
               data = FRic_all,
               sample_prior = TRUE,  prior = c(
                 prior(normal(0,1), "b"),
                 prior(normal(0,1), "Intercept")),
               control = list(adapt_delta = .95),
               file = "Models/DB/SESFRicGAMwz.rds",
               chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

model_check(data = FRic_all, model = SES.fricGAM)

#### FOri Model fit ####
## basic model of focal variables and samplng structure
FOriGAM <- brm(fori ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                     s(cluster, bs="re", by= Cluster_dummy),
                   family = Beta(),
                   data = FMulti_all,
                   sample_prior = TRUE,  prior = c(
                     prior(normal(0,1), "b"),
                     prior(normal(0,1), "Intercept")),
                   control = list(adapt_delta = .95),
                   file = "Models/DB/FOriGAM.rds",
                   chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

## more complete model incorporating that increasing wood in pasture sites 
## may covary with species composition
FOriGAM <- brm(fori ~ habitat + s(elev_z, by = habitat, bs = "tp", k = 7) + 
                     Wood_Veg_z:habitat + Wood_Veg_z +
                     s(cluster, bs="re", by= Cluster_dummy),
                   family = Beta(),
                   data = FMulti_all,
                   sample_prior = TRUE,  prior = c(
                     prior(normal(0,1), "b"),
                     prior(normal(0,1), "Intercept")),
                   control = list(adapt_delta = .95),
                   file = "Models/DB/FOriGAMwz.rds",
                   chains = 4, iter = 600, thin = 1, cores = 4, warmup = 300)

model_check(data = FMulti_all, model = FOriGAM)
