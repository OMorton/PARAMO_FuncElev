
#### Functions ####

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

f_deriv <- function(data = data , model = model, 
                    elev_raw = focal_raw,
                    eps = 1e-4, summary = TRUE) {
  
  ## make data
  first <- data
  second <- mutate(data, elev_z = elev_z + eps)
  
  ## get predictions
  first_preds <-
    add_epred_draws(model, newdata = first) %>% ungroup() %>% rename(first = .epred)
  second_preds <-
    add_epred_draws(model, newdata = second) %>% ungroup() %>% select(.epred) %>% rename(second = .epred)
  
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

p_contrasts <- function(fitting_data = FMulti_all, model = FSpeGAM) {

  ## Set up the raw data
  new_dat <- data.frame(elev = seq(from = 1000, to = 3000, by = 50)) %>% 
    mutate(elev_z = (elev - mean(fitting_data$ele_jaxa)) / sd(fitting_data$ele_jaxa),
           Cluster_dummy = 0, cluster = "cluster_MOF_1")

  Forest_fit <- add_epred_draws(model, newdata = 
                                  mutate(new_dat, habitat = "Forest",
                                         Wood_Veg_z = mean(filter(fitting_data, habitat == "Forest")$Wood_Veg_z))) %>% 
    rename("Forest" = ".epred")
  Pasture_fit <- add_epred_draws(model, newdata = 
                                   mutate(new_dat, habitat = "Pasture",
                                          Wood_Veg_z = mean(filter(fitting_data, habitat == "Pasture")$Wood_Veg_z))) %>% 
    rename("Pasture" = ".epred") %>% ungroup() %>% select(Pasture)
  
  ## Do the contrast
  Contr <- cbind(Forest_fit, Pasture_fit) %>% mutate(Contr = Forest - Pasture) %>%
    group_by(elev) %>% mutate(pd = sum(sign(Contr) == sign(median(Contr)))/max(.draw))
  
  Out <- Contr %>% group_by(elev, pd) %>% median_hdci(Contr, .width = .9) 
  
  return(Out)
}
