
#### Functions ####
library(ggpubr)
#### Model check ####

## Simple convenience function to do and present 4 posterior predictive checks.
## Will to updated to include more/flexibility
## model - model fitted
## data - the data used to fit the model
model_check <- function(model = model, data = data) {
  p1 <- pp_check(model) + theme_minimal() + theme(legend.position = "none")
  p2 <- pp_check(model, type = "ecdf_overlay") + theme_minimal() + theme(legend.position = "none")
  p3 <- pp_check(model, type = "stat_2d") + theme_minimal() + theme(legend.position = "none")

  p4 <- data %>%
    add_residual_draws(model) %>%
    ggplot(aes(x = .row, y = .residual)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
    stat_pointinterval(colour = "#011f4b", alpha = .25) + 
    theme_minimal()
  
  ch <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,
                  labels = c("A.", "B.", "C.", "D."))
  return(ch)
}

#### New data for models ####
## simple
get_newdata <- function(fitting_data = data) {
  
  constants <- fitting_data %>% group_by(habitat) %>% 
    reframe(For_Dist_z = median(For_Dist_z),
            Wood_Veg_z = median(Wood_Veg_z),
            Cluster_dummy = 0, cluster = "cluster_MOF_1")
  
  elev_range <- data.frame(habitat = c("Forest", "Pasture")) %>% group_by(habitat) %>%
    reframe(elev = seq(from = 800, to = 3600, by = 200)) %>% 
    mutate(elev_z = (elev - mean(fitting_data$ele_jaxa))/sd(fitting_data$ele_jaxa))
  
  min_forest_elev <- min(filter(fitting_data, habitat == "Forest")$ele_jaxa)
  max_forest_elev <- max(filter(fitting_data, habitat == "Forest")$ele_jaxa)
  min_pasture_elev <- min(filter(fitting_data, habitat == "Pasture")$ele_jaxa)
  max_pasture_elev <- max(filter(fitting_data, habitat == "Pasture")$ele_jaxa)
  
  elev_trim <- 
    elev_range %>% filter((habitat == "Forest" & elev >= min_forest_elev &
                           elev <= max_forest_elev ) | 
                        (habitat == "Pasture" & elev >= min_pasture_elev &
                           elev <= max_pasture_elev )
                          )
  
  new_dat <- left_join(constants, elev_trim)
  
  return(new_dat)
}
#### First derivatives finite differences ####

## data - Not raw data used to fit the model, a data frame of values to estimate
## the firest derivative for.
## model - model fitted 
## elev_raw - the raw vector of elevation data used to internally reverse the 
## standardisation of the elevation variable. Not strictly necessary but convenient.
## eps - the epsilon value. This is the adjustment difference to calculate the
## the difference over. Works best at smaller values.
## summary - logical. Either returns a summary (median and 90% hdci) or the 
## raw data.

f_deriv <- function(data = data , model = model, prediction = TRUE,
                    elev_raw = focal_raw,
                    eps = 1e-4, summary = TRUE) {
  
  ## make data
  first <- data
  second <- mutate(data, elev_z = elev_z + eps)
  
  if (prediction == TRUE) {
  ## get predictions
  first_preds <-
    add_predicted_draws(model, newdata = first) %>% ungroup() %>% rename(first = .prediction)
  second_preds <-
    add_predicted_draws(model, newdata = second) %>% ungroup() %>% select(.prediction) %>% rename(second = .prediction)
  } else {
  
  ## get expectation
  first_preds <-
    add_epred_draws(model, newdata = first) %>% ungroup() %>% rename(first = .epred)
  second_preds <-
    add_epred_draws(model, newdata = second) %>% ungroup() %>% select(.epred) %>% rename(second = .epred)
  }
  
  diff <- cbind(first_preds, second_preds)
  
  ## Calcualte the difference and divide by epsilon - this is analgous to the dx/dt 
  diff_sum <- diff %>% mutate(diff = (second - first)/eps) %>% 
    mutate(elev = elev_z*sd(elev_raw) + mean(elev_raw),
           ## put on the per 100 metre scale rather than the sd scale
           diff = (diff/sd(elev_raw))*100 ) %>% 
    group_by(elev, habitat) %>%
    mutate(pd = sum(sign(diff) == sign(median(diff)))/max(.draw))
  
  if( summary == TRUE) {
    diff_sum <- diff_sum %>% group_by(habitat, elev_z, elev, pd) %>% 
      median_hdci(diff, .width = .9)
  }
  return(diff_sum)
}


#### Posterior contrasts ####

## fitting_data - data used to fit the model.
## model - fitted model

p_contrasts <- function(forest_data = forest, pasture_data = pasture, model = FSpeGAM) {

  min_shared_elev <- max(c(min(forest_data$elev), min(pasture_data$elev)))
  max_shared_elev <- min(c(max(forest_data$elev), max(pasture_data$elev)))
  
  forest_trim <- forest_data %>% filter(elev >= min_shared_elev & elev <= max_shared_elev)
  pasture_trim <- pasture_data %>% filter(elev >= min_shared_elev & elev <= max_shared_elev)
  
  Forest_fit <- add_predicted_draws(model, newdata = forest_trim) %>% 
    rename("Forest" = ".prediction")
  
  Pasture_fit <- add_predicted_draws(model, newdata = pasture_trim) %>% 
    rename("Pasture" = ".prediction") %>% ungroup() %>% select(Pasture)
  
  ## Do the contrast
  Contr <- cbind(Forest_fit, Pasture_fit) %>% mutate(Contr = Forest - Pasture) %>%
    group_by(elev) %>% mutate(pd = (sum(sign(Contr) == sign(median(Contr)))/max(.draw))*100)
  
  Out <- Contr %>% group_by(elev, pd) %>% median_hdci(Contr, .width = .9) %>%
    mutate(pd_binary = ifelse(pd >= 97.5, "Yes", "No"))
  
  return(Out)
}


#### Bird habitat:elevation summary plot wrapper ####
plots_fx <- function(x) {
  
    ggplot(x, aes(elev_plot, value, color = habitat, group = habitat)) + 
    geom_point(size = 3, shape = 16) +
    geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, linewidth = 1) +
    scale_color_manual(values = c("olivedrab4", "tan3")) +
    scale_x_continuous(name = "Elevation (m)",
                       breaks = c(1000, 1400, 1800, 2200, 2600, 3000, 3400),
                       labels = c(1000, 1400, 1800, 2200, 2600, 3000, 3400)) +
    ylab(paste(unique(x$metric))) +
    theme_minimal(base_size = 12) +
    theme(axis.line.x.bottom = element_line(colour = "black"),
          axis.line.y.left = element_line(colour = "black"),
          axis.ticks = element_line(), legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5))
}

#### Bird habitat:elevation difference plot wrapper ####
plots_fx_diff <- function(y) {
  cordillera_diff_long %>% filter(., metric == y) %>% 
    ggplot(., aes(elev_ALOS, value)) +
    stat_summary(fun.data = median_hdci, 
                 fun.args = c(.width = .9),
                 geom = "pointrange") +
    geom_hline(yintercept = 0, linetype = "dashed", size = .75) +
    scale_x_discrete(name = "Elevation (m)",
                     breaks = c("1000", "1400", "1800", "2200", "2600", "3000", "3400")) +
    ylab(paste0("Difference in ", str_replace(y, "_diff", ""))) +
    theme_minimal(base_size = 12) +
    theme(axis.line.x.bottom = element_line(colour = "black"),
          axis.line.y.left = element_line(colour = "black"),
          axis.ticks = element_line(), axis.text.x = element_text(angle = 45, vjust = 0.5))
}
#### Bird metric summary ####
bird_sum <- function(x) {
  x1 <- x %>% group_by(elev_ALOS, point_type, habitat, metric) %>%
    median_hdci(value, .width = .9, na.rm = TRUE)  %>%
    mutate(elev_plot = ifelse(habitat == "Forest", 
                              as.numeric(as.character(elev_ALOS)) - 15, 
                              as.numeric(as.character(elev_ALOS)) + 15))
  return(x1)
}
