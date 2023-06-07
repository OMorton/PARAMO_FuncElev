######################################
#### Targetted Functional Indices ####
######################################

## Summary functional values across land use and elevation.

## AIMS
## 1 - compute the quality of functional dendrogramm and multidimensional functional spaces.
## 2 - compute and illustrate multidimensional functional diversity indices for all bands.

## Set up packages for working from home
.libPaths("C:/Packages") 
library(tidyverse)
library(mFD)

#### Trait Preparation ####
## All species are coprophagus? Trait should be removed
DB_trait_raw <- read.csv("Data/output_allDBtraits_EC_MASTER.csv") %>% select(-X)

# FLA is 0.9 correlated with bodysize
cor(DB_trait_raw[,2:4][complete.cases(DB_trait_raw[,2:4]), ])

# Create dataframe of relevant traits
traits.all <- DB_trait_raw %>%
  mutate(med_bodysize_z = (med_bodysize - mean(med_bodysize))/sd(med_bodysize),
         med_legratio_z = (med_legratio - mean(med_legratio))/sd(med_legratio)) %>%
  # Select trait columns 
  select(scientificName, med_bodysize_z, med_legratio_z, nest_guild,
         necrophagus, rotting_fruit, rotting_fungi, activity) %>%
  mutate(activity = as.factor(activity), nest_guild = as.factor(nest_guild)) %>% 
  column_to_rownames(var = "scientificName") 

Trait_Descr <- data.frame(trait_name = c("med_bodysize_z", "med_legratio_z", "nest_guild",
                                         "necrophagus", "rotting_fruit", "rotting_fungi",
                                         "activity"), 
                          trait_type = c("Q", "Q", "N",
                                         "F", "F", "F",
                                         "N"), 
                          fuzzy_name = c(NA, NA, NA,
                                         "Diet", "Diet", "Diet", NA))

## Computing distances between species based on functional traits
sp_dist <- funct.dist(traits.all, Trait_Descr, metric = "gower")

## Compute multidimensional functional spaces and assess their quality
dist_qual <- quality.fspaces(sp_dist, maxdim_pcoa = 10,
                             deviation_weighting = "absolute", fdist_scaling = FALSE, fdendro = "average")

plt <- mFD::traits.faxes.cor(
  sp_tr  = traits.all, 
  sp_faxes_coord = dist_qual$details_fspaces$sp_pc_coord[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot = TRUE)
plt

## examine the quality of fspace
round(dist_qual$"quality_fspaces", 3) %>% as.data.frame()


#### Presence data ####
## Read in the point data.
Point_info_raw <- read.csv("Data/Point_info/output_ECRAWdata.ptsinfo.csv") %>% rename(habitat = habitat_all_points)
Point_info <- Point_info_raw %>% 
  select(point, habitat, elev_bands) %>%
  distinct()

## Read in the point abundance data.
## 108 species across 392 sites.
DB_Abund_raw <-  read.csv("Data/output_EC_fulldataDB_abunmatrix.csv") 

DB_Abund_long <- DB_Abund_raw %>% pivot_longer(-X) %>% 
  rename(species = 1, point = 2, abundance = 3) %>%
  mutate(presence = ifelse(abundance > 0, 1, 0)) %>%
  ## add env variables
  left_join(Point_info) %>%
  ## remove paramo
  filter(habitat != "Paramo") %>% 
  ## remove this species and site per Felicity (Chase why?)
  filter(species !="Canthon_sp._08H", point != "CCF1")

#### FRic ####

## Trade off between the number of axes and the avaialble data  

## Select axes, that minimize the MAD
coord_DB_3D <- dist_qual$details_fspaces$sp_pc_coord[, 1:3]

## 107 sp and 354 sites to start
length(unique(DB_Abund_long$species))
length(unique(DB_Abund_long$point))

## Remove sites that have no species and any species that never appear (0).
## Now 100 sp and 167 sites
DB_filter <- DB_Abund_long %>% 
  ## remove all sites with 2 or less species for the func metrics
  group_by(point) %>% filter(sum(presence) > 3) %>% 
  group_by(species) %>% filter(sum(presence) > 0)
length(unique(DB_filter$species))
length(unique(DB_filter$point))

## Make wide form again
DB_Ab_mat <- DB_filter %>% select(species, point, presence) %>% 
  pivot_wider(values_from = presence, names_from = point) %>% 
  column_to_rownames(var = "species")  %>% as.matrix()

## Transpose the matrix
DB_Ab_t <- t(DB_Ab_mat)
## Check the min number of species per site
rowSums(DB_Ab_t) %>% as.data.frame() %>% summarise(min(.))

## Functional alpha diversity indices in a multidimensional space
test <- alpha.fd.multidim(sp_faxes_coord = coord_DB_3D,
                          asb_sp_w = DB_Ab_t,
                          ind_vect = c("fric"),
                          scaling = TRUE, check_input = TRUE, details_returned = TRUE)


#### Standardised FRic ####
iter <- 400

## Make an empty object to store the data in
Storage <- data.frame(matrix(NA, nrow = 0, ncol = 0))

## For loop to get simulated expected communities.
for (i in 1:iter) {
  ## Prints an update to the console with each iteration
  cat(i, " out of ", iter, " simulated communities ", " (",i/iter*100,"%)\n", sep = "")
  
  ## Randomizes the species order of the occurence matrix so you get random
  ## communities from each country
  random_matrix <- picante::randomizeMatrix(DB_Ab_t, null.model = "richness")

  ## Calc fd
  expected_fd <- alpha.fd.multidim(
    sp_faxes_coord   = coord_DB_3D,
    asb_sp_w         = random_matrix,
    ind_vect         = c("fric"),
    scaling          = TRUE, verbose = FALSE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  ## tidying - add country column (exporter/importer/year etc) and add a collum
  ## that stores the simulation number
  metrics_iter <- expected_fd$functional_diversity_indices %>% 
    rownames_to_column(var = "point") %>% mutate(sim = i)
  
  ## save to storage
  Storage <- rbind(Storage, metrics_iter)
  
  ## issue that some combinations of species introduce coplanearity in the 
  ## geometry::convhull() function called internally by mfd helper functions.
  ## therefore iterate across communities until we have 250 sets of 
  ## communities with no coplanar incidences. 
  if (Storage %>% 
           group_by(sim) %>% filter(!is.na(fric)) %>% tally() %>%
           filter(n == 167) %>% nrow() == 250) break
}

## Summarise the simulated comms
expected_sum <- Storage %>% group_by(sim) %>% filter(all(!is.na(fric))) %>% 
  group_by(point) %>% 
  summarise(ave_expected = mean(fric), sd_expected = sd(fric))

FRic_DB_Points_NA <-  test$functional_diversity_indices %>% as.data.frame() %>% 
  rownames_to_column(var = "point") %>% left_join(expected_sum, by = "point") %>%
  mutate(SES.fric = (fric - ave_expected)/sd_expected)

#### Checking ####

## Add the elevation data
Point_elev_info <- Point_info_raw %>% 
  select(point, elev_ALOS_all_points, habitat) %>%
  distinct()
FRic_DB_Points_elev <- left_join(FRic_DB_Points_NA, Point_elev_info)

#### PLOTTING ####
ggplot(FRic_DB_Points_elev, aes(elev_ALOS_all_points, SES.fric, colour = habitat)) + 
  geom_point() + geom_smooth()

#### OUTPUTS ####

## write out band FRic
write.csv(FRic_DB_Points_elev, "Outputs/Summaries/DB/DB_FRic_Point_all_NEW.csv")


#### FMulti FSpe FOri FDis ####

## 107 sp and 354 sites to start
length(unique(DB_Abund_long$species))
length(unique(DB_Abund_long$point))

## Make long form again
DB_filter <- DB_Abund_long %>% 
  ## remove all sites with 0 or less species for the func metrics
  group_by(point) %>% filter(sum(presence) > 0) %>% 
  group_by(species) %>% filter(sum(presence) > 0)
length(unique(DB_filter$species)) ## 107
length(unique(DB_filter$point)) ## 280

## Make long form again
DB_Ab_mat_FOther <- DB_filter %>% select(species, point, presence) %>% 
  pivot_wider(values_from = presence, names_from = point) %>% 
  column_to_rownames(var = "species")  %>% as.matrix()

## Transpose the matrix
DB_Ab_t_FOther <- t(DB_Ab_mat_FOther)

## Check the min number of species per site
rowSums(DB_Ab_t_FOther) %>% as.data.frame() %>% summarise(min(.))

## Functional alpha diversity indices in a multidimensional space
test <- alpha.fd.multidim(sp_faxes_coord = coord_DB_3D,
                          asb_sp_w = DB_Ab_t_FOther,
                          ind_vect = c("fdis", "fori", "fspe"),
                          scaling = TRUE, check_input = TRUE, details_returned = TRUE)


#### Standardized FMulti FSpe FOri FDis ####
iter <- 250

## Make an empty object to store the data in
Storage <- data.frame(matrix(NA, nrow = 0, ncol = 0))

## For loop to get simulated expected communities.
for (i in 1:iter) {
  ## Prints an update to the console with each iteration
  cat(i, " out of ", iter, " simulated communities ", " (",i/iter*100,"%)\n", sep = "")
  
  ## Randomizes the species order of the occurence matrix so you get random
  ## communities from each country
  random_matrix <- picante::randomizeMatrix(DB_Ab_t_FOther, null.model = "richness")
  
  ## Calc fd
  expected_fd <- alpha.fd.multidim(
    sp_faxes_coord   = coord_DB_3D,
    asb_sp_w         = random_matrix,
    ind_vect         = c("fdis", "fori", "fspe"),
    scaling          = TRUE, verbose = FALSE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  ## tidying - add country column (exporter/importer/year etc) and add a collum
  ## that stores the simulation number
  metrics_iter <- expected_fd$functional_diversity_indices %>% 
    rownames_to_column(var = "point") %>% mutate(sim = i)
  
  ## save to storage
  Storage <- rbind(Storage, metrics_iter)
  
}

## Summarise the simulated comms
expected_sum_multi <- Storage %>% group_by(sim) %>% filter(all(!is.na(fdis))) %>% 
  group_by(point) %>% 
  summarise(ave_expected_fdis = mean(fdis), sd_expected_fdis = sd(fdis),
            ave_expected_fori = mean(fori), sd_expected_fori = sd(fori),
            ave_expected_fspe = mean(fspe), sd_expected_fspe = sd(fspe))

FMulti_DB_Points_NA <-  test$functional_diversity_indices %>% as.data.frame() %>% 
  rownames_to_column(var = "point") %>% left_join(expected_sum_multi, by = "point") %>%
  mutate(SES.fdis = (fdis - ave_expected_fdis)/sd_expected_fdis,
         SES.fori = (fori - ave_expected_fori)/sd_expected_fori,
         SES.fspe = (fspe - ave_expected_fspe)/sd_expected_fspe,)


#### Checking ####

## Add the elevation data
Point_elev_info <- Point_info_raw %>% 
  select(point, elev_ALOS_all_points, habitat) %>%
  distinct()
FMulti_DB_Points_elev <- left_join(FMulti_DB_Points_NA, Point_elev_info)

#### PLOTTING ####
ggplot(FMulti_DB_Points_elev, aes(elev_ALOS_all_points, fori, colour = habitat)) + 
  geom_point() + geom_smooth()

#### OUTPUTS ####

## write out band FRic
write.csv(FMulti_DB_Points_elev, "Outputs/Summaries/DB/DB_Multi_Point_all_NEW.csv")
