run_model<-function(start_conditions,model_pars,idx=NA){

  # Bookkeeping and basic setup ----------------------------------------------
  start. <- now()                               # start time stamp
  n_years <- model_pars$sim$n_years             # number of simulated years
  
  pars <- list()                                # container for parameters
  
  # Biology and demography parameters -----------------------------------------
  pars$max_age        <- model_pars$bio$inherent$max_age
  pars$breeding_age   <- model_pars$bio$inherent$age_first_breed
  pars$sex_ratio      <- model_pars$bio$inherent$sex_ratio
  pars$max_brood_size <- model_pars$bio$inherent$max_prog_per_brood
  pars$av_clutch_size <- model_pars$bio$inherent$av_clutch_size
  
  # Genetics parameters --------------------------------------------------------
  pars$Fp               <- model_pars$bio$gen$starting_inbreeding
  pars$founder_kinship  <- model_pars$bio$gen$founder_kinship
  # pars$dip_leth_eq   <- model_pars$bio$gen$dip_leth_eq   # not used here
  
  # Environment and vital rate coefficients -----------------------------------
  pars$carr_capac_df <- model_pars$bio$carr_capac_df  # carrying capacity by site or year
  
  pars$phi_df <- model_pars$bio$surv_coeff            # survival coefficients
  pars$cs_df  <- model_pars$bio$cs_coeffs             # clutch size coefficients
  pars$hp_df  <- model_pars$bio$hatch_prob_coeff      # hatch probability coefficients
  pars$fp_df  <- model_pars$bio$fledge_prob_coeff     # fledge probability coefficients
  pars$fb_df  <- model_pars$bio$fb_coeffs             # fecundity or brood size coefficients
  
  # Assisted incubation or artificial incubation parameters -------------------
  pars$ai_hatch_prob  <- model_pars$bio$ai_hatch_prob
  pars$ai_fledge_prob <- model_pars$bio$ai_fledge_prob
  
  # Dispersal settings ---------------------------------------------------------
  pars$dispersalMat  <- model_pars$bio$dispersalMat
  pars$dispersalAges <- 2:model_pars$bio$inherent$max_age # ages that can disperse
  
  
  # Admixture release schedule and intensity ----------------------------------
  year_adm_rel <- model_pars$mgmt$year_admix_releases
  # Draw once whether an admixture program happens
  bin_adm_rel <- rbinom(n = 1, size = 1, prob = model_pars$mgmt$prob_admix_releases) *
  # this sets the draw to always be zero if admixture is not part of the strategy  
    model_pars$mgmt$admix_releases
  
  if (bin_adm_rel == 0) { year_adm_rel <- Inf }        # never occurs if draw is zero
  # Logical vector of years where admix releases occur
  
  years<-1:n_years
  candidate <- years >= year_adm_rel &
    ((years - year_adm_rel) %% model_pars$mgmt$admix_release_freq == 0)
  # Keep only the first 'max_number' TRUEs
  idx <- which(candidate)
  if (length(idx) > max_number) {
    idx <- idx[1:model_pars$mgmt$admix_total_releases]
  }
  release_vec <- rep(0, n_years)
  release_vec[idx] <- 1

  
  # pars$admix_release_years <- round(pmin(1, pmax(0, (1:n_years) - year_adm_rel + 1))) == 1
  pars$admix_release_years <- release_vec==1
  pars$admix_prop_released <- model_pars$mgmt$admix_prop_released
  pars$admix_age_released  <- model_pars$mgmt$admix_age_released
  pars$admix_no_released   <- model_pars$mgmt$admix_no_released
  
  # Initial population and ID bookkeeping -------------------------------------
  pop <- start_conditions$init_pop
  pars$all_ids <- pars$founder_ids <- unique(pop$id)
  
  # Management parameters ------------------------------------------------------
  pars$field <- model_pars$mgmt$field
  pars$first_clutch_harvest_prop <- model_pars$mgmt$egg_harvest_rate
  pars$prop_managed_clutches <- model_pars$mgmt$prop_managed %>%
    dplyr::filter(NewMgmt == "TRUE") %>%
    dplyr::pull(propManaged)
  
  # Set up initial inbreeding table at individual level -----------------------
  Fi_df <- pop %>%
    dplyr::arrange((is.na(Fi))) %>%         # ensures non NA come first
    dplyr::filter(!duplicated(id)) %>%      # one row per individual
    dplyr::select(id, Fi, nz_heritage)      # keep relevant genetic cols
  
  # Containers to store outputs across years ----------------------------------
  egg_fate   <- data.frame()                 # currently unused, left for future
  envir_stoch <- data.frame()                # to record stochastic draws
  
  j <- 1                                     # simulation year counter
  
  suppressMessages({
    while (j <= n_years) {
      
      # Keep ID registry updated with any newcomers
      pars$all_ids <- unique(pop$id)
      
      # Determine management releases for this year ---------------------------
      released <- releases(pars = pars, currentT = j)
      
      # Full population is previous plus current releases
      pars$full_pop <- plyr::rbind.fill(pop, released)
      
      # Select alive individuals at t = j - 1 and advance them to t = j -------
      currentPop0 <- pop %>%
        dplyr::filter(t == j - 1, alive) %>%
        dplyr::mutate(t = j)
      
      # Add releases into the current time step pool
      pars$all_ids <- unique(c(pars$all_ids, pars$full_pop$id))
      currentPop0 <- currentPop0 %>% plyr::rbind.fill(released)
      
      # If no individuals remain, stop early
      if (nrow(currentPop0) == 0) { break }
      
      # Survival process with environmental stochasticity ---------------------
      # mortality() returns a list with $alive and $stoch_q
      survival_outcome <- mortality(pop = currentPop0, currentT = j, pars = pars)
      
      # Record draw of survival enviromental stochastic component for diagnostics
      envir_stoch <- plyr::rbind.fill(
        envir_stoch,
        data.frame(par = "survival", q = survival_outcome$stoch_q, t = j)
      )
      
      currentPop1 <- survival_outcome$alive
      
      # Update pairings if one partner died -----------------------------------
      currentPop2 <- unpair_if_dead(pop = currentPop1, currentT = j)
      
      # Pairing process for breeding ------------------------------------------
      currentPop3 <- pairing(pop = currentPop2, currentT = j, pars = pars)
      
      # Recruitment and reproduction ------------------------------------------
      born_resu <- recruitment(pop = currentPop3, currentT = j, pars = pars)
      
      # Update ID registry with newly born individuals
      pars$all_ids <- unique(c(pars$all_ids, pars$full_pop$id, born_resu$born$id))
      
      # Pass harvested eggs info into artificial incubation step
      pars$harvested_eggs <- born_resu$harvested_eggs
      
      # Artificial incubation recruits (if any)
      ai_recruits <- ai_recruitment(pars = pars)
      
      # Clear harvested eggs to avoid reuse
      pars$harvested_eggs <- NULL
      
      # Combine adults and new recruits, then age everyone ---------------------
      currentPop4 <- plyr::rbind.fill(currentPop3, born_resu$born, ai_recruits) %>%
        aging(currentT = j, pars = pars)
      
      # Append this year to the main population time series and tidy ----------
      pop <- plyr::rbind.fill(pop, currentPop4) %>% tidy_pop_df()
      
      # Recompute inbreeding coefficients given updated pedigree --------------
      Fi_df <- calculate_inbreeding(pop_df = pop, pars = pars)
      
      # Merge fresh Fi back into population table ------------------------------
      pop <- pop %>%
        dplyr::mutate(Fi = NULL, nz_heritage = NULL) %>%
        dplyr::left_join(Fi_df, by = "id")
      
      # Advance time step
      j <- j + 1
    }
  })
  
  # Final summaries and returns -----------------------------------------------
  suppressMessages({
    # Snapshot of unique individuals with parentage and genetics
    full_ids <- pop %>%
      dplyr::filter(!duplicated(id)) %>%
      dplyr::select(id, father_id, mother_id, origin, sex) %>%
      dplyr::left_join(Fi_df, by = "id") %>%
      dplyr::mutate(i = idx)
    
    # kinship <- calculate_kinship(pop_df = full_ids, pars = pars, rm_non_breeders = FALSE)
  })
  
  # Assemble output list. Note: egg_fate and kinship are currently unused -----
  return(list(
    pop        = pop %>% dplyr::mutate(i = idx),
    full_ids   = full_ids,
    envir_stoch = envir_stoch,
    # kinship   = kinship,
    pars       = pars,
    # egg_fate = egg_fate %>% dplyr::mutate(i = idx),
    model_pars = model_pars,
    time_run   = difftime(time1 = now(), time2 = start., units = "secs")
  ))
}