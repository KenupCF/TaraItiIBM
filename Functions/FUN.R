# Function to initialize a population data frame
init_population <- function(pars, seed = 19) {
  require(dplyr)
  
  # Validate required parameters
  needed_pars <- c("StartN_df", "sex_ratio", "no_age_classes", "Fp", "age_structure")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Prepare containers
  pop <- data.frame()
  existing_ids <- character()
  
  # Remove zero rows from starting counts
  pars$StartN_df <- pars$StartN_df %>% filter(count > 0)
  
  # Build initial individuals per row of StartN_df
  for (i in seq_along(pars$StartN_df[, 1])) {
    # Create new rows for this stratum of age and sex
    new_indivs <- data.frame(
      id  = generate_unique_ids(n = pars$StartN_df$count[i], existing = existing_ids),
      age = rep(pars$StartN_df$age_class[i], pars$StartN_df$count[i]),
      sex = rep(pars$StartN_df$sex[i],       pars$StartN_df$count[i])
    )
    pop <- rbind(pop, new_indivs)
    existing_ids <- unique(pop$id)
  }
  
  # Add state and meta columns
  pop <- pop %>%
    dplyr::mutate(
      t = 0,
      alive = TRUE,
      pair = NA_character_,
      mother_id = NA_character_,
      father_id = NA_character_,
      origin = "Wild",
      nz_heritage = 1,
      sex = toupper(substr(sex, 1, 1)),
      Fi = pars$Fp,
      year_born = age - t,
      subpop = "A"
    )
  
  return(pop)
}

# Mortality used in your main loop (simpler, field and age specific effects)
mortality <- function(pop, currentT, pars, seed = 19) {
  
  # Validate required parameters
  needed_pars <- c("phi_df", "max_age", "full_pop")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Alive at current time
  current <- pop %>%
    dplyr::filter(t == currentT, alive) %>%
    mutate(field = pars$field) %>%
    ungroup()
  
  # Carrying capacity truncation terms
  N_df <- current %>%
    dplyr::group_by(subpop) %>%
    dplyr::summarise(N = sum(alive), .groups = "drop") %>%
    dplyr::left_join(pars$carr_capac_df) %>%
    dplyr::mutate(AboveC = N > C, surv_trunc = 1 / (N / C))
  
  # Shared stochastic term (zeroed below)
  q_sample <- runif(n = 1, min = 0, max = 1)
  
  # Intercept and optional variation
  surv_mean <- mean <- pars$phi_df$Intercept
  cyc_var <- 0
  sto_var <- sapply(seq_along(current$id), function(i) {
    qnorm(p = q_sample, mean = 0, sd = pars$cs_df$RE_time)
  })
  surv_int <- surv_mean + sto_var + cyc_var
  
  # Survival probability by age, field, heritage, and origin
  surv <- inv.logit(
    surv_int +
      (pars$phi_df$Beta_imm        * (current$age == 2)) +
      (pars$phi_df$Beta_adu        * (current$age >= 3)) +
      (pars$phi_df$Beta_newsq_juv  * (current$age == 1) * (current$field == "SQ")) +
      (pars$phi_df$Beta_newsq_imm  * (current$age == 2) * (current$field == "SQ")) +
      (pars$phi_df$Beta_newsq_adu  * (current$age >= 3) * (current$field == "SQ")) +
      (pars$phi_df$Beta_newsq_female * (current$sex == "F") * (current$field == "SQ")) +
      (pars$phi_df$Beta_newsq_male   * (current$sex == "M") * (current$field == "SQ")) +
      (pars$phi_df$Beta_admix_juv  * (current$age == 1) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) +
      (pars$phi_df$Beta_admix_imm  * (current$age == 2) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) +
      (pars$phi_df$Beta_admix_adu  * (current$age >= 3) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) +
      (pars$phi_df$Beta_ai_juv     * (current$age == 1) * (current$origin == "AI")) +
      (pars$phi_df$Beta_ai_imm     * (current$age == 2) * (current$origin == "AI")) +
      (pars$phi_df$Beta_ai_adu     * (current$age >= 3) * (current$origin == "AI")) +
      (pars$phi_df$Beta_field2_juv * (current$age == 1) * (current$field == "Field2")) +
      (pars$phi_df$Beta_field2_imm * (current$age == 2) * (current$field == "Field2")) +
      (pars$phi_df$Beta_field2_adu * (current$age >= 3) * (current$field == "Field2"))
  )
  
  current$surv <- surv
  
  # Apply carrying capacity truncation
  current <- current %>%
    dplyr::left_join(N_df) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(surv = dplyr::case_when(AboveC ~ min(surv_trunc, surv), TRUE ~ surv))
  
  # Validate probabilities
  if (any(is.na(surv)) | any(surv < 0) | any(surv > 1)) {
    stop("Invalid estimated survival probability")
  }
  
  # Realized survival and age limit
  alive_vec <- sapply(seq_along(current$id), function(i) {
    rbinom(n = 1, size = 1, prob = current$surv[i])
  }) == 1 & current$alive & current$age < pars$max_age
  
  # Update alive flag and clear pairs for deaths
  current$alive <- alive_vec
  current <- current %>%
    dplyr::mutate(
      pair = dplyr::case_when(!alive ~ NA, TRUE ~ pair)
      # Age increment occurs in aging()
    ) %>%
    tidy_pop_df()
  
  resu <- list(alive = current, stoch_q = q_sample)
  return(resu)
}


# Increment ages for individuals at current time
aging <- function(pop, currentT, pars, seed = 19) {
  
  # Validate parameters used downstream
  needed_pars <- c("phi_df", "max_age", "field", "full_pop")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Select current step and add 1 year for alive individuals
  current <- pop %>%
    dplyr::filter(t == currentT) %>%
    dplyr::mutate(age = dplyr::case_when(alive ~ age + 1, TRUE ~ age)) %>%
    tidy_pop_df()
  
  return(current)
}


# Pair breeding age adults by subpopulation
pairing <- function(pop, currentT, pars, seed = 19) {
  require(dplyr)
  
  # Need breeding_age
  needed_pars <- c("breeding_age")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Priority column used to keep latest updates when binding
  pop <- pop %>% dplyr::mutate(priority = 0)
  
  # Pool of available adults at current time without pairs
  adults <- pop %>%
    filter(alive, t == currentT, age >= pars$breeding_age) %>%
    dplyr::ungroup()
  
  pairs <- data.frame()
  subpops <- unique(adults$subpop)
  
  # Pair within each subpopulation
  for (sp in subpops) {
    males   <- adults %>% filter(sex == "M", subpop == sp, is.na(pair))
    females <- adults %>% filter(sex == "F", subpop == sp, is.na(pair))
    
    # Randomize order to avoid bias
    if (nrow(males) > 0)   males   <- males[sample(nrow(males)), ]
    if (nrow(females) > 0) females <- females[sample(nrow(females)), ]
    
    n_pairs <- min(nrow(males), nrow(females))
    if (n_pairs > 0) {
      pairs <- plyr::rbind.fill(
        pairs,
        (males[1:n_pairs, ])   %>% dplyr::mutate(pair = females$id[1:n_pairs], priority = 1),
        (females[1:n_pairs, ]) %>% dplyr::mutate(pair = males$id[1:n_pairs],   priority = 1)
      )
    }
  }
  
  # Merge paired individuals back and tidy
  adults_new <- plyr::rbind.fill(pop, pairs) %>% tidy_pop_df()
  return(adults_new)
}


# Remove pair references to dead partners at current time
unpair_if_dead <- function(pop, currentT) {
  require(dplyr)
  
  alive_ids <- pop %>%
    dplyr::filter(t == currentT, alive) %>%
    dplyr::pull(id)
  
  pop <- pop %>%
    dplyr::mutate(pair = dplyr::case_when(!(pair %in% alive_ids) ~ NA, TRUE ~ pair))
  
  return(pop)
}


# Reproduction pipeline: who breeds, eggs laid, hatch, fledge, offspring added
recruitment <- function(pop, currentT, pars, seed = 19, envir_stoch = TRUE, startAge = 0) {
  
  # Validate required parameters
  needed_pars <- c("breeding_age", "fb_df", "cs_df", "hp_df", "fp_df",
                   "sex_ratio", "max_brood_size",
                   "first_clutch_harvest_prop", "prop_managed_clutches",
                   "carr_capac_df", "all_ids")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Candidates of breeding pairs at current time
  temp <- pop %>%
    dplyr::filter(!is.na(pair), alive, t == currentT, age >= pars$breeding_age) %>%
    ungroup() %>%
    mutate(able_to_breed = TRUE)
  
  # Track breed capacity for both sexes to enforce pair feasibility
  male_breed_capacity   <- temp %>% filter(sex == "M") %>% pull_named(able_to_breed, id)
  female_breed_capacity <- temp %>% filter(sex == "F") %>% pull_named(able_to_breed, id)
  female_pair           <- temp %>% filter(sex == "F") %>% pull_named(pair, id)
  
  f_breed_df <- data.frame(
    f_able_to_breed = female_breed_capacity,
    id   = names(female_breed_capacity),
    pair = female_pair,
    m_able_to_breed = male_breed_capacity[female_pair]
  )
  
  # Eligible reproducing females with both partners able
  reproducing <- pop %>%
    dplyr::filter(sex == "F", !is.na(pair), alive, t == currentT, age >= pars$breeding_age) %>%
    dplyr::left_join(f_breed_df) %>%
    dplyr::mutate(able_to_breed = f_able_to_breed & m_able_to_breed) %>%
    dplyr::filter(able_to_breed)
  
  # Enforce carrying capacity by limiting number of breeding females
  fem_able_to_breed <- sample(reproducing$id,
                              size = min(length(reproducing$id), pars$carr_capac_df$C),
                              replace = FALSE)
  reproducing <- reproducing %>%
    dplyr::mutate(able_to_breed = id %in% fem_able_to_breed)
  
  # Probability of breeding given field and origin
  reproducing$breeding_prob <- pars$fb_df$Intercept +
    (pars$fb_df$Beta_newsq  * (pars$field == "SQ")) +
    (pars$fb_df$Beta_field2 * (pars$field == "Field2")) +
    (pars$fb_df$Beta_ai     * (reproducing$origin == "AI")) + 0
  
  # Bernoulli draw for breed or not
  reproducing$breed <- sapply(seq_along(reproducing$id), function(i) {
    rbinom(n = 1, size = 1, prob = reproducing$breeding_prob)
  }) == 1
  reproducing <- reproducing %>% dplyr::filter(breed)
  
  # Environmental stochastic quantile for random effects
  if (envir_stoch) { q_sample <- runif(n = 1, min = 0, max = 1) } else { q_sample <- .5 }
  
  if (nrow(reproducing) > 0) {
    females_breeding <- nrow(reproducing)
    
    # Mean clutch size on log scale
    cs_mean <- pars$cs_df$Intercept
    pars$cs_df$RE_time <- 0
    cs_st_var <- sapply(seq_along(reproducing$id), function(i) {
      qnorm(p = q_sample, mean = 0, sd = pars$cs_df$RE_time)
    })
    cs_int <- cs_mean + cs_st_var
    
    # Expected brood size per female, then draw Poisson and cap at max
    brood_size <- exp(
      cs_int +
        (pars$cs_df$Beta_newsq  * (pars$field == "SQ")) +
        (pars$cs_df$Beta_field2 * (pars$field == "Field2")) +
        (pars$cs_df$Beta_admix  * (reproducing$nz_heritage < 0.95 & reproducing$nz_heritage > 0.05)) +
        (pars$cs_df$Beta_ai     * (reproducing$origin == "AI")) + 0
    ) %>%
      pmax(0) %>%
      identity()
    
    reproducing$cs <- brood_size
    no_eggs <- pmin(sapply(reproducing$cs, function(x) rpois(n = 1, lambda = x)), pars$max_brood_size)
    
    # Expand mother and father ids for each egg
    motherIDs <- lapply(seq_along(reproducing$id), function(i) rep(reproducing$id[i],    no_eggs[i])) %>% unlist()
    fatherIDs <- lapply(seq_along(reproducing$id), function(i) rep(reproducing$pair[i], no_eggs[i])) %>% unlist()
    
    # Midparent heritage
    nz_heritage_vec <- (((temp %>% pull_named(nz_heritage, id))[fatherIDs]) +
                          ((temp %>% pull_named(nz_heritage, id))[motherIDs])) / 2
    
    # Offspring sex based on sex ratio
    newSex <- c("F", "M")[rbinom(n = sum(no_eggs), size = 1, prob = pars$sex_ratio) + 1]
    
    if (sum(no_eggs) > 0) {
      # Create egg rows
      eggs <- data.frame(
        id = generate_unique_ids(n = sum(no_eggs), existing = pars$all_ids),
        age = 0,
        sex = newSex,
        subpop = "A",
        t = currentT,
        alive = TRUE,
        pair = NA,
        mother_id = motherIDs,
        origin = "Wild",
        father_id = fatherIDs,
        nz_heritage = nz_heritage_vec,
        clutch_no = 1,
        year_born = currentT
      )
      
      # Add maternal covariates
      eggs$mother_origin <- (temp %>% pull_named(origin, id))[eggs$mother_id]
      eggs$mother_age    <- (temp %>% pull_named(age, id))[eggs$mother_id]
      
      # Harvest a proportion of first clutches for AI
      harvested_clutches_mom_id <- sample(
        unique(eggs$mother_id),
        size = round(length(unique(eggs$mother_id)) * pars$first_clutch_harvest_prop),
        replace = FALSE
      )
      harvested_eggs <- eggs %>% dplyr::filter(mother_id %in% harvested_clutches_mom_id)
      harvested_eggs$origin <- "AI"
      harvested_eggs$id <- NA_character_
      
      # Mothers whose first clutch was harvested can lay a second clutch
      eggs <- eggs %>%
        dplyr::mutate(
          clutch_no = dplyr::case_when(mother_id %in% harvested_clutches_mom_id ~ clutch_no + 1, TRUE ~ clutch_no)
        )
      
      # Assign managed clutches
      clutches_ids <- unique(eggs$mother_id)
      managed_clutches <- rbinom(n = length(clutches_ids), size = 1, prob = pars$prop_managed_clutches)
      managed_clutches_ids <- clutches_ids[managed_clutches == 1]
      eggs$managed <- eggs$mother_id %in% managed_clutches_ids
      
      # Hatch probability per egg
      hp_mean <- pars$hp_df$Intercept
      hp_st_var <- sapply(seq_along(eggs$id), function(i) qnorm(p = q_sample, mean = 0, sd = pars$hp_df$RE_time))
      hp_int <- hp_mean + hp_st_var
      
      hatch_prob <- inv.logit(
        hp_int +
          (pars$hp_df$Beta_f_age   * eggs$mother_age) +
          (pars$hp_df$Beta_newsq   * (pars$field == "SQ")) +
          (pars$hp_df$Beta_field2  * (pars$field == "Field2")) +
          (pars$hp_df$Beta_admix   * (eggs$nz_heritage < 0.95 & eggs$nz_heritage > 0.05)) +
          (pars$hp_df$Beta_ai      * (eggs$mother_origin == "AI")) +
          (pars$hp_df$Beta_managed * (eggs$managed)) +
          (pars$hp_df$Beta_clutch_no * eggs$clutch_no) + 0
      )
      eggs$hp <- hatch_prob
      
      # Realized hatching
      eggs$alive <- sapply(seq_along(eggs$id), function(i) rbinom(n = 1, size = 1, prob = eggs$hp[i])) == 1
      eggs_laid <- nrow(eggs)
      eggs_hatched <- sum(eggs$alive)
      
      if (eggs_hatched > 0) {
        hatchlings <- eggs %>% dplyr::filter(alive)
        
        # Fledge probability per hatchling
        fp_mean <- pars$fp_df$Intercept
        fp_st_var <- sapply(seq_along(hatchlings$id), function(i) qnorm(p = q_sample, mean = 0, sd = pars$fp_df$RE_time))
        fp_int <- fp_mean + fp_st_var
        
        fledge_prob <- inv.logit(
          fp_int +
            (pars$fp_df$Beta_f_age   * hatchlings$mother_age) +
            (pars$fp_df$Beta_newsq   * (pars$field == "SQ")) +
            (pars$fp_df$Beta_field2  * (pars$field == "Field2")) +
            (pars$fp_df$Beta_admix   * (hatchlings$nz_heritage < 0.95 & hatchlings$nz_heritage > 0.05)) +
            (pars$fp_df$Beta_ai      * (hatchlings$mother_origin == "AI")) +
            (pars$fp_df$Beta_managed * (hatchlings$managed)) +
            (pars$fp_df$Beta_clutch_no * hatchlings$clutch_no) + 0
        )
        hatchlings$fp <- fledge_prob
        
        # Realized fledging
        hatchlings$alive <- sapply(seq_along(hatchlings$id), function(i) rbinom(n = 1, size = 1, prob = hatchlings$fp[i])) == 1
        fledged <- sum(hatchlings$alive)
        
        # Newborn rows
        if (fledged > 0) {
          born <- hatchlings %>% filter(alive) %>% mutate(Fi = NA) %>% tidy_pop_df()
        } else {
          born <- data.frame()
        }
      } else {
        fledged <- 0
        born <- data.frame()
      }
    } else {
      born <- data.frame()
      harvested_eggs <- data.frame()
      eggs_laid <- 0
      eggs_hatched <- 0
      fledged <- 0
    }
    
  } else {
    # No reproduction possible this year
    born <- data.frame()
    harvested_eggs <- data.frame()
    females_breeding <- 0
    eggs_laid <- 0
    eggs_hatched <- 0
    fledged <- 0
  }
  
  # Return components for diagnostics and AI pipeline
  return(list(
    born = born,
    females_breeding = if (exists("females_breeding")) females_breeding else 0,
    harvested_eggs = if (exists("harvested_eggs")) harvested_eggs else data.frame(),
    eggs_laid = if (exists("eggs_laid")) eggs_laid else 0,
    eggs_hatched = if (exists("eggs_hatched")) eggs_hatched else 0,
    fledged = if (exists("fledged")) fledged else 0
  ))
}


# Convert harvested eggs via AI pathway into recruits
ai_recruitment <- function(pars, seed = 19, envir_stoch = TRUE, startAge = 0) {
  
  # Validate required fields on pars
  needed_pars <- c("ai_hatch_prob", "ai_fledge_prob", "harvested_eggs", "all_ids")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  eggs <- pars$harvested_eggs
  if (nrow(eggs) > 0) {
    # Assign fresh ids, then hatch, then fledge
    eggs$id <- generate_unique_ids(n = nrow(eggs), existing = pars$all_ids)
    
    eggs$alive <- sapply(seq_along(eggs$id), function(i) {
      rbinom(n = 1, size = 1, prob = pars$ai_hatch_prob)
    }) == 1
    hatchlings <- eggs %>% dplyr::filter(alive)
    
    if (nrow(hatchlings) > 0) {
      hatchlings$alive <- sapply(seq_along(hatchlings$id), function(i) {
        rbinom(n = 1, size = 1, prob = pars$ai_fledge_prob)
      }) == 1
      
      born <- hatchlings %>%
        mutate(Fi = NA) %>%
        tidy_pop_df() %>%
        dplyr::filter(alive)
    } else {
      born <- data.frame()
    }
  } else {
    born <- data.frame()
  }
  
  return(born)
}


# Move unpaired dispersers among subpops based on a probability matrix
dispersal <- function(pop, currentT, pars, seed = 19) {
  
  # Validate parameters
  needed_pars <- c("dispersalMat", "dispersalAges")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Priority column to keep updates
  pop <- pop %>% dplyr::mutate(priority = 0)
  
  # Label set of subpops (row names of movement matrix)
  subpops <- rownames(pars$dispersalMat)
  
  # Eligible dispersers are alive, current time, in dispersal ages, and unpaired
  dispersers <- pop %>%
    dplyr::filter(alive, t == currentT, age %in% pars$dispersalAges, is.na(pair)) %>%
    dplyr::mutate(priority = 1)
  
  # If none, return unchanged
  if (nrow(dispersers) == 0) return(pop)
  
  # Sample destination subpop by row-specific probabilities
  dispersal_destinations <- sapply(seq_along(dispersers$id), function(i) {
    subpops[which(rmultinom(n = 1, size = 1, prob = pars$dispersalMat[dispersers$subpop[i], ])[, 1] == 1)]
  })
  dispersers$subpop <- dispersal_destinations
  
  # Merge updated rows and tidy
  newpop <- plyr::rbind.fill(pop, dispersers) %>% tidy_pop_df()
  return(newpop)
}


# Zero truncated Poisson between bounds [min, max)
rtpois <- function(n, lambda, min = 0, max = Inf) {
  x <- rpois(n, lambda)
  while (any(x <= min | x >= max)) {
    x[x <= min | x >= max] <- rpois(sum(x <= min | x >= max), lambda)
  }
  return(x)
}


# Generate n unique IDs not in existing set
generate_unique_ids <- function(n, existing = character(0),
                                strlength = 8,
                                charset = c(0:9, LETTERS),
                                max_tries = 10 * n) {
  unique_ids <- character(0)
  tries <- 0
  
  while (length(unique_ids) < n && tries < max_tries) {
    needed <- n - length(unique_ids)
    candidates <- replicate(needed, paste0(sample(charset, strlength, replace = TRUE), collapse = ""))
    # Remove collisions with existing or already generated
    new_ids <- setdiff(candidates, c(existing, unique_ids))
    unique_ids <- unique(c(unique_ids, new_ids))
    tries <- tries + 1
  }
  
  if (length(unique_ids) < n) {
    stop("Failed to generate the required number of unique IDs after ", max_tries, " attempts.")
  }
  return(unique_ids)
}


# Keep latest rows by priority and time, standardize order and columns
tidy_pop_df <- function(df) {
  suppressWarnings({ if (is.null(df$priority)) { df$priority <- 0 } })
  
  resu <- df %>%
    dplyr::arrange(desc(priority), desc(t)) %>%          # prioritize updates
    dplyr::filter(!duplicated(data.frame(id, t))) %>%    # one record per id and time
    dplyr::select(
      id, subpop, sex, age, t, alive, pair,
      mother_id, father_id,
      Fi, nz_heritage,
      year_born,
      origin
    ) %>%
    dplyr::arrange(desc(alive), subpop, age, mother_id, origin)
  
  return(resu)
}


# Blend a probability with an odds ratio over a fraction of the period
adjust_probability <- function(prob, odds_ratio, prop = 1, agg_rule = "survival") {
  
  # Basic checks
  if (any(prob < 0) || any(prob > 1)) stop("Probability must be between 0 and 1 (inclusive).")
  if (any(odds_ratio <= 0))           stop("Odds ratio must be greater than 0.")
  
  # Missing prop means fully adjusted
  prop[is.na(prop)] <- 1
  
  # Convert to odds, multiply by OR, convert back
  odds <- prob / (1 - prob)
  adjusted_odds <- odds * odds_ratio
  adjusted_prob <- adjusted_odds / (1 + adjusted_odds)
  
  # Keep exact 0 or 1 unchanged
  adjusted_prob[prob %in% c(0, 1)] <- prob[prob %in% c(0, 1)]
  
  # Aggregate over the fraction of time exposed
  if (agg_rule == "survival") {
    # Geometric blend on survival scale
    resu <- (prob^(1 - prop)) * (adjusted_prob^prop)
  }
  # Alternative for success scale could be added if needed
  
  return(resu)
}


# Create releases in admixture years according to pars
releases <- function(pars, currentT) {
  
  # Validate required inputs
  needed_pars <- c("admix_release_years", "admix_prop_released", "admix_age_released",
                   "admix_no_released", "all_ids")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if (length(missing_pars) > 0) {
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  releaseYears <- which(pars$admix_release_years)
  
  if (currentT %in% releaseYears & pars$admix_no_released > 0) {
    # Build release cohort
    released <- data.frame(
      id   = generate_unique_ids(n = pars$admix_no_released, existing = pars$all_ids),
      age  = pars$admix_age_released,
      sex  = c("F", "M")[rbinom(n = pars$admix_no_released, size = 1, prob = 0.5) + 1],
      subpop = "A",
      t = currentT,
      alive = TRUE,
      pair = NA,
      mother_id = NA,
      origin = "Captivity",
      nz_heritage = sample(c(1, 0.5), replace = TRUE, size = pars$admix_no_released,
                           prob = c(1 - pars$admix_prop_released, pars$admix_prop_released)),
      father_id = NA,
      Fi = 0
    )
  } else {
    released <- data.frame()
  }
  
  return(released)
}


# Replace NA strings in character columns with empty strings
replace_na_characters <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col[is.na(col)] <- ""
    }
    return(col)
  })
  return(df)
}


# Greedy bin packing by a target column into n groups
split_evenly_by_col <- function(df, n_groups, target_col) {
  library(dplyr)
  
  # Sort by descending target to place large items first
  df <- df %>% arrange(desc(!!sym(target_col)))
  
  # Initialize group sums and assignment
  group_sums <- rep(0, n_groups)
  df$group <- NA_integer_
  
  for (i in 1:nrow(df)) {
    # Assign to the group with the smallest total so far
    min_group <- which.min(group_sums)
    df$group[i] <- min_group
    group_sums[min_group] <- group_sums[min_group] + df[[target_col]][i]
  }
  
  # Split into list by group
  split_list <- split(df, df$group)
  return(split_list)
}


# Upload result files to Google Drive and move to backup directory, one by one
output_clean_up <- function(rev = FALSE) {
  files <- list.files(path = "./Results", pattern = ".RData", full.names = TRUE)
  cat(paste0("\nFound ", length(files), " files\n"))
  
  if (rev) files <- rev(files)
  if (length(files) > 0) {
    for (f in files) {
      suppressMessages({
        # Upload each file
        gdrive_upload <- sapply(f, drive_upload, path = as_id("1FSqFfBrvyedxTOUuIFmI44u9tXxEMMU0"))
        uploaded_file <- attributes(gdrive_upload)$dimnames[[2]]
        # Move uploaded file into backup_dir
        file.rename(uploaded_file, file.path(backup_dir, basename(f)))
      })
    }
  }
}


# Zip all .RData results, upload once, then move both zip and sources to backup_dir
output_clean_up_zip <- function(rev = FALSE) {
  files <- list.files(path = "./Results", pattern = "\\.RData$", full.names = TRUE)
  cat(paste0("\nFound ", length(files), " files\n"))
  
  if (length(files) == 0) return(invisible(NULL))
  if (rev) files <- rev(files)
  
  # Build a unique zip name in tempdir
  zip_filename <- file.path(
    tempdir(),
    paste0("results_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
  )
  
  # Zip without directory structure
  zip(zipfile = zip_filename, files = files, flags = "-j")
  
  suppressMessages({
    # Upload zip to Drive
    gdrive_upload <- drive_upload(zip_filename, path = as_id("1FSqFfBrvyedxTOUuIFmI44u9tXxEMMU0"))
    # Move zip and sources to backup_dir
    file.rename(zip_filename, file.path(backup_dir, basename(zip_filename)))
    file.rename(files, file.path(backup_dir, basename(files)))
  })
}
