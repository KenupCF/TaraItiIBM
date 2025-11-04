# ------------------------------------------------------------------------------
# OVERVIEW ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# define functions to calculate:
# proportion of females breeding
# proportion females having second and third nesting attempt
# clutch size 
# define PVA

# Code  by Nicole Whitelock @ the University of Auckland and Finnbar Lee @ Cawthron,
# in 2024, finnbar.lee@cawthron.org.nz. 

# PVA originally coded by Thalassa Hamilton, with updates by Finnbar Lee here.

# ------------------------------------------------------------------------------
# load packages ----------------------------------------------------------------
# ------------------------------------------------------------------------------

library(tidyverse)
library(truncnorm)

# ------------------------------------------------------------------------------
# Define functions -------------------------------------------------------------
# ------------------------------------------------------------------------------

# calculate proportion of females breeding -------------------------------------
# this is a bit niggly and there is likely better ways to do this but it works:
# need to include sightings data to see how many birds there are in each year 

fun_prop_female_breed_og <- function(sightings_data,
                                     capture_morphometrics,
                                     egg_data,
                                     years) {
  
  # start with annual female sightings
  # remove NAs for joining 
  sightings_tidy <- sightings %>%
    filter(!is.na('C band number')) %>% 
    rename(band_number = `C band number`) # rename to make life easier 
  
  # determine the sex of each bird - this data is available in capture morphometrics (can join using band number)
  sexV <- capture_morphometrics %>%
    mutate(band_number = gsub("C", "", Band_Number)) %>% # change name for joining 
    dplyr::select(-Band_Number) %>% 
    dplyr::select(band_number, sex_verified) 
  
  # join sex verification with sightings
  sightings_sexV <- sightings_tidy %>% 
    left_join(sexV, by = "band_number", relationship = "many-to-many")
  
  # determine sightings each year 
  annual_sighting <- sightings_sexV %>%
    group_by(Year, sex_verified) %>%
    distinct(band_number) %>% # each bird counted once per year
    ungroup()
  
  # filter for female sightings only 
  female_sightings <- annual_sighting %>% 
    filter(sex_verified == "female") %>% 
    group_by(Year) %>%
    distinct(band_number) %>% 
    ungroup()
  
  # determine the number of female sightings in each year 
  female_summary <- female_sightings %>% 
    group_by(Year) %>% 
    summarise(total_female = n()) %>% 
    ungroup()
  
  #crop season - remove years with iffy/incomplete sightings data
  female_summary_crop <- female_summary %>%
    filter(Year %in% years)
  
  # now include annual females breeding
  female_breed <- egg_data %>%
    group_by(season) %>%
    summarise(female_breed = n_distinct(female_id)) %>% 
    ungroup()
  
  # change season to year and crop to match sightings data
  female_breed_crop <- female_breed %>% 
    mutate(Year = as.numeric(str_extract(season, "\\d{4}"))) %>% 
    filter(Year %in% years)
  
  
  # join female breeding and sighting dfs 
  breed_sightings <- left_join(female_breed_crop, female_summary_crop, by = "Year")
  
  # calculate how many females sighted each year are breeding
  summary_prop_female_breed <- breed_sightings %>%
    summarise(avg_breeding_proportion = mean(female_breed / total_female),
              min=min(female_breed / total_female),
              max=max(female_breed / total_female),
              sd_breeding_proportion = sd(female_breed / total_female))
  
  
  return(list(prop_breeding = round(summary_prop_female_breed$avg_breeding_proportion[1],3),
              sd_breeding = round(summary_prop_female_breed$sd_breeding_proportion[1],3)))
  
}

fun_prop_female_breed <- function(sightings_data,
                                  capture_morphometrics,
                                  egg_data,
                                  years,
                                  years_cov=NULL) {
  
  if(is.null(years_cov)){
    years_cov<-data.frame(Year=years,foo=T)
    
  }
  # start with annual female sightings
  # remove NAs for joining 
  sightings_tidy <- sightings %>%
    filter(!is.na('C band number')) %>% 
    rename(band_number = `C band number`) # rename to make life easier 
  
  # determine the sex of each bird - this data is available in capture morphometrics (can join using band number)
  sexV <- capture_morphometrics %>%
    mutate(band_number = gsub("C", "", Band_Number)) %>% # change name for joining 
    dplyr::select(-Band_Number) %>% 
    dplyr::select(band_number, sex_verified) 
  
  # join sex verification with sightings
  sightings_sexV <- sightings_tidy %>% 
    left_join(sexV, by = "band_number", relationship = "many-to-many")
  
  # determine sightings each year 
  annual_sighting <- sightings_sexV %>%
    group_by(Year, sex_verified) %>%
    distinct(band_number) %>% # each bird counted once per year
    ungroup()
  
  # filter for female sightings only 
  female_sightings <- annual_sighting %>% 
    filter(sex_verified == "female") %>% 
    group_by(Year) %>%
    distinct(band_number) %>% 
    ungroup()
  
  # determine the number of female sightings in each year 
  female_summary <- female_sightings %>% 
    group_by(Year) %>% 
    summarise(total_female = n()) %>% 
    ungroup()
  
  #crop season - remove years with iffy/incomplete sightings data
  female_summary_crop <- female_summary %>%
    filter(Year %in% years)
  
  # now include annual females breeding
  female_breed <- egg_data %>%
    group_by(season) %>%
    summarise(female_breed = n_distinct(female_id)) %>% 
    ungroup()
  
  # change season to year and crop to match sightings data
  female_breed_crop <- female_breed %>% 
    mutate(Year = as.numeric(str_extract(season, "\\d{4}"))) %>% 
    filter(Year %in% years)
  
  
  # join female breeding and sighting dfs 
  breed_sightings <- left_join(female_breed_crop, female_summary_crop, by = "Year")%>%
    dplyr::left_join(years_cov)%>%
    dplyr::mutate(prop= female_breed/total_female)
  
  
  mod<-glm(formula =   prop ~ NewMgmt,
           data = breed_sightings)
  
  newdata<- years_cov%>%filter(!duplicated(data.frame(NewMgmt)))%>%
    mutate(Year=NULL)
  
  pred<-predict(object = mod,
                newdata =newdata,
                se.fit=T)
  
  coeffs<-data.frame(estimate=mod$coefficients,se=sqrt(diag(vcov(mod))))
  resu<-newdata%>%
    mutate(av=pred$fit,se=pred$se.fit)
  
  return(list(pred=resu,coeffs=coeffs,mod=mod))
  
}

# Calculate proportion females having second and third nesting attempt --------

fun_prop_multi_nest <- function(nest_data,
                                years) {
  
  prop_multi_attempt <- nest_data %>%
    mutate(Year = as.numeric(str_extract(season, "\\d{4}"))) %>%
    filter(Year %in% years) %>%
    dplyr::select(Year, fate_nest, clutch_no, no_hatched, female_band) %>%
    group_by(Year, female_band) %>%
    # only get those cases where the next was not successful
    filter(sum(clutch_no == 1 & fate_nest %in% c("fledged", "hatched")) == 0) %>%
    group_by(Year, clutch_no) %>%
    reframe(n_atempts = table(clutch_no)) %>%
    filter(Year %in% years) %>%
    group_by(Year) %>%
    mutate(prop_attempt = n_atempts/ lag(n_atempts)) %>%
    filter(clutch_no %in% 2:3) %>%
    group_by(clutch_no) %>%
    summarise(mn_at = mean(prop_attempt),
              sd_at = sd(prop_attempt))
  
  return(list(prop_2nd_nest = round(prop_multi_attempt$mn_at[1],3),
              sd_2nd_nest = round(prop_multi_attempt$sd_at[1],3),
              prop_3rd_nest = round(prop_multi_attempt$mn_at[2],3),
              sd_3rd_nest = round(prop_multi_attempt$sd_at[2],3)))
}


# calculate mean clutch size ---------------------------------------------------
fun_clutch_size_og <- function(nest_data) {
  
  clutchsize <- nest_data %>% 
    summarise(mean_clutchSize = mean(total_eggs_clutch), 
              sd_clutchSize = sd(total_eggs_clutch)) %>% 
    ungroup()
  
  return(list(clutch_size = round(clutchsize$mean_clutchSize[1],3),
              sd_clutch_size = round(clutchsize$sd_clutchSize[1],3)))
}

fun_clutch_size <- function(nest_data,years,years_cov=NULL) {
  
  dat <- nest_data %>% 
    filter(year%in%years)%>%
    mutate(Year=year)%>%
    left_join(years_cov)
  
  mod<-glm(formula = total_eggs_clutch~NewMgmt,family = "poisson",data=dat)
  
  coeffs <- data.frame(
    estimate = coef(mod),
    se = sqrt(diag(vcov(mod)))
  )
  
  newdata<- years_cov%>%filter(!duplicated(data.frame(NewMgmt)))%>%
    mutate(Year=NULL)
  
  pred<-predict(object = mod,
                newdata =newdata,
                se.fit=T)
  

  resu<-newdata%>%
    mutate(av=pred$fit,se=pred$se.fit)
  
  
  clutchsize<-dat%>%
    summarise(mean_clutchSize = mean(total_eggs_clutch), 
              sd_clutchSize = sd(total_eggs_clutch)) %>% 
    ungroup()
  
  return(list(clutch_size = round(clutchsize$mean_clutchSize[1],3),
              sd_clutch_size = round(clutchsize$sd_clutchSize[1],3),
              coeffs=coeffs))
}

# ------------------------------------------------------------------------------
# PVA MODEL --------------------------------------------------------------------
# ------------------------------------------------------------------------------

# model as a function
# use 2022 abundances (juv, imm, adult)as default
fun_pva_model <- function(mod_params,
                          nsim = 1000,
                          initial_abund = c(2, 7, 11),
                          nyears = 50) {
  
  # helper function
  # translate mean and sd of survival into parameters of a beta distribution
  beta.params <- function(x_bar, sd.x){
    u <- x_bar * (1 - x_bar) / sd.x ^ 2 - 1
    alpha <- x_bar * u
    beta <- (1 - x_bar) * u
    return(list(alpha = alpha, beta = beta))
  }
  
  
  # SET UP PARAMETERS ------------------------------------------------------------
  
  # Number of years to run model
  nyears <- nyears
  
  # Pull parameter values
  
  # females breeding
  K <-   mod_params$mean_est[1]   # mean prop of females breeding
  K.e <- mod_params$sd_est[1]     # standard deviation
  
  # clutch size
  cl <-   mod_params$mean_est[2]   # mean clutch size 
  cl.e <- mod_params$sd_est[2]     # standard deviation
  
  # eggs hatching
  h <-   mod_params$mean_est[3]    # mean prob of eggs hatching
  h.e <- mod_params$se_est[3]      # std error
  h.t <- mod_params$temp_est[3]    # temporal variability
  
  # hatchlings fledging to juvenile
  fl <- mod_params$mean_est[4]    # mean prob of hatchlings fledging to juvenile 
  fl.e <- mod_params$se_est[4]    # std error
  fl.t <- mod_params$temp_est[4]  # temporal variability
  
  # Juvenile survival
  Sjuv <- mod_params$mean_est[5]   # mean juv survival
  Sjuv.e <- mod_params$se_est[5]   # std error
  Sjuv.t <- mod_params$temp_est[5] # temporal variability
  
  # immature survival
  Simm <- mod_params$mean_est[6]   # mean immature survival 
  Simm.e <- mod_params$se_est[6]   # std error
  Simm.t <- mod_params$temp_est[6] # temporal variability
  
  # Adult survival
  Sad <- mod_params$mean_est[7]    # mean adult survival
  Sad.e <- mod_params$se_est[7]   # std error
  Sad.t <- mod_params$temp_est[7] # temporal variability
  
  ## Carrying capacity
  Car <- mod_params$mean_est[8]
  
  # harvest - proportion eggs harvested for captive rearing, assume females
  # relay after harvest, and harvested eggs have the same survival  as wild individuals
  harvest <- mod_params$mean_est[9]
  #harvest <- 1
  
  # maximum harvest number
  max_harv <- mod_params$mean_est[10]
  max_harv <- 1000
  
  # SOMEWHERE TO STORE RESULTS ---------------------------------------------------
  
  # Define population matrix
  N <- array(0, dim = c(3, nyears + 1, nsim))
  
  # starting population vector
  Ni <- initial_abund
  N[,1,] <- Ni 
  
  # summary info
  alive <- matrix(0, nrow = nyears, ncol = nsim)
  r.sq <- matrix(0, nrow = nyears, ncol = nsim)
  mean.r <- lambda.sq <- numeric()
  sj.val <- si.val <- sa.val <- Fa.val <- K.val <- c.val <-
    h.val <- fl.val <- persist <- har.val <- max_h.val <- car.val <- numeric()
  
  
  # RUN REPLICATE MODELs -------------------------------------------------------
  
  # Loop over replicate simulations
  for (s in seq_len(nsim)){
    
    #s <- 1
    
    # Generate demographic parameters from beta and normal distributions 
    
    # females breeding
    if (K.e == 0) {
      K.sim <- K
    } else {
      K.sim <- rbeta(1, beta.params(K, K.e)$alpha, beta.params(K, K.e)$beta)
    }
    
    # clutch size
    if (cl.e != 0) {
      c.sim <- truncnorm::rtruncnorm(1, a = 0, b = 2, mean = cl, sd = cl.e) 
    } else {
      c.sim <- cl
    }
    
    # eggs hatching - convert to logit scale (SE already on logit scale)
    h.sim <-  rnorm(1, qlogis(h), h.e)
    
    # chicks fledging to juvenile - convert to logit scale (SE already on logit scale) 
    fl.sim <- rnorm(1, qlogis(fl), fl.e)
    
    # juvenile survival 
    if (Sjuv.e == 0) {
      sj.sim <- Sjuv
    } else {
      sj.sim <- rbeta(1, beta.params(Sjuv, Sjuv.e)$alpha, beta.params(Sjuv, Sjuv.e)$beta)
    }
    
    # immature survival
    if (Simm.e == 0) {
      si.sim <- Simm
    } else {
      si.sim <-  rbeta(1, beta.params(Simm, Simm.e)$alpha, beta.params(Simm, Simm.e)$beta)
    }
    
    # adult survival
    if (Sad.e == 0) {
      sa.sim <- Sad
    } else {
      sa.sim <- rbeta(1, beta.params(Sad, Sad.e)$alpha, beta.params(Sad, Sad.e)$beta)
    }
    
    # Generate annual demographic rates (subject to temporal variability)
    Ka <- rnorm(nyears, K.sim, 0) # no temporal var
    ca <- rnorm(nyears, c.sim, 0) # no temporal var
    ha <- plogis(rnorm(nyears, h.sim, h.t))
    fla <- plogis(rnorm(nyears, fl.sim, fl.t))
    sj <- plogis(rnorm(nyears, qlogis(sj.sim), Sjuv.t))
    si <- plogis(rnorm(nyears, qlogis(si.sim), Simm.t))
    sa <- plogis(rnorm(nyears, qlogis(sa.sim), Sad.t))
    
    
    # Project population
    for (t in seq_len(nyears)){ # Loop over years
      
      # t <- 1
      
      # harvests go to 0 if more than 50 adult females
      harv_t <- ifelse(N[3,t,s] < 51, harvest, 0) 
      
      # fecundity equation - number of female fledglings produced per adult female
      # Fa <- (Ka[t] * ca[t] * ha[t] * fla[t]) / 2
      
      Fa <- (Ka[t] * ca[t] * ((1 - harv_t) + (harv_t * Ka[t])) * ha[t] * fla[t]) / 2
      
      
      # wild juveniles corrected for carrying capacity
      new_j <- rpois(1, min(Car, N[3,t,s]) * sa[t] * Fa) 
      
      # prop females breeding * mean clutch size * number of adult females * prop harvest = number of eggs harvested
      # if this is more than max_harv, then only take max_harv
      # then number harvested * hatch rate * fledge rate / 2 for just females = total new juveniles from captivity
      cap_j <- rpois(1, (min(Ka[t] * ca[t] * min(Car, N[3,t,s]) * sa[t] * harv_t, max_harv) * ha[t] * fla[t]) / 2)
      
      
      
      # combine wild and captive reared new juveniles
      N[1,t+1,s] <- new_j + cap_j
      
      # update stages
      N[2,t+1,s] <- rbinom(1, (N[1,t,s]), sj[t])
      N[3,t+1,s] <- rbinom(1, (N[2,t,s]), si[t]) + rbinom(1, (N[3,t,s]), sa[t])
      
      # stop if extinct, or calculate intrinsic growth rate
      if (sum(N[,t+1,s]) == 0) {break} else {
        r.sq[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s]))
        alive[t,s] <- t
      } # else
      
      # persistence
      #persist[t] <- sum(N[3,t,]> 3)  / s 
      
    } # t
    
    # save param values 
    K.val[s] = K.sim
    c.val[s] = c.sim
    h.val[s] = plogis(h.sim)
    fl.val[s] = plogis(fl.sim)
    sj.val[s] = sj.sim
    si.val[s] = si.sim
    sa.val[s] = sa.sim  
    Fa.val[s] = Fa
    har.val[s] = harvest
    max_h.val[s] = max_harv
    car.val[s] = Car
    
    # get intrinsic growth rate and finite rate of increase
    mean.r[s] <- mean(r.sq[min(alive[,s], na.rm = TRUE):max(alive[,s], na.rm = TRUE),s])
    lambda.sq[s] <- exp(mean.r[s]) 
  } # s
  
  # results
  res <- data.frame(sim = 1:nsim,
                    r = mean.r,
                    lambda = lambda.sq,
                    n_50_af = N[3,nyears + 1,1:nsim],
                    pbreed = K.val,
                    clutch_size = c.val,
                    hatch = h.val,
                    fledge = fl.val,
                    sj = sj.val,
                    si = si.val,
                    sa = sa.val,
                    Fa = Fa.val,
                    harv = har.val,
                    max_h = max_h.val,
                    K = car.val)
  
  # prep data for plots
  N.mean <- apply(N,c(1,2),mean)
  N.median <- apply(N,c(1,2),median)
  N.10 <- apply(N,c(1,2),quantile, probs=c(.10))
  N.90 <- apply(N,c(1,2),quantile, probs=c(.90))
  N.2_5 <- apply(N,c(1,2),quantile, probs=c(.025))
  N.97_5 <- apply(N,c(1,2),quantile, probs=c(.975))
  
  
  n_time <- data.frame(years = 1:(nyears + 1),
                       N.mean = N.mean[3,],
                       N.median = N.median[3,],
                       N.10 = N.10[3,],
                       N.90 = N.90[3,],
                       N.2_5 = N.2_5[3,],
                       N.97_5 = N.97_5[3,])
  
  # get prop runs persist (N > 3 adult females) through time
  prop_persist <- apply(N, c(1, 2), function(x){
    sum(x > 3) / nsim
  })
  
  prop_persist <- data.frame(years = 1:(nyears + 1),
                             p.persist = prop_persist[3,])
  
  
  
  # calculate a few metrics ------------------------------------------------------
  # probability of adult female population size being > 2 throughout the 50 years
  prob_gr_2 <- 1 - sum(apply(N[3,,1:nsim], 2, min) < 3) / nsim # adults only
  
  # Quasi-extinction
  # probability of less than 3 birds after 50 years
  prob_less_3 <- sum(N[3, nyears + 1, ]<= 3) / nsim
  
  # mean pop size after 50 years  
  pop_size <- mean(N[3,nyears + 1,])
  
  # Quantiles - Thalassa had 0, 36
  pop_quant <- quantile(N[3,nyears,], c(0.1, 0.9))
  
  # Persistence probability (after nyears)
  prob_persist <- sum(alive[nyears,]!= 0) / nsim
  
  # probability of extinction after 50 years
  prob_extint <- sum(N[3, nyears + 1, ] == 0) / nsim
  
  res_summary <- data.frame(prob_gr_2 = prob_gr_2,
                            prob_less_3 = prob_less_3,
                            pop_size = pop_size,
                            pop_quant_10 = pop_quant[1],
                            pop_quant_90 = pop_quant[2],
                            prob_persist = prob_persist,
                            prob_extint = prob_extint)
  
  res_all <- list(results = res,
                  n_time = n_time,
                  res_summary = res_summary,
                  prop_persist = prop_persist)
  
  return(res_all)
  
}

# END --------------------------------------------------------------------------