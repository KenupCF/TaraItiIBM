# ============================
# calculate_kinship
# ============================
# Builds a pedigree from pop_df and returns a kinship matrix.
# If rm_non_breeders = TRUE, removes individuals that never appear as a parent.
# Founder kinship for specified founders is overridden using pars.
calculate_kinship <- function(pop_df, pars, rm_non_breeders = TRUE) {
  
  require(kinship2)  # For pedigree() and kinship()
  
  # Parameters expected in 'pars'
  needed_pars <- c("founder_ids", "founder_kinship")
  # (No explicit check here; see notes below if you want to enforce)
  
  if (rm_non_breeders) {
    # Keep only individuals that are listed as a mother or father at least once
    pop_df <- pop_df %>%
      dplyr::filter(id %in% unique(c(pop_df$mother_id, pop_df$father_id)))
  }
  
  # Prepare unique rows for the pedigree
  ped_df <- pop_df %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::select(id, mother_id, father_id, sex)
  
  # Create pedigree object
  # kinship2 expects sex codes: 1 = male, 2 = female
  ped <- with(
    ped_df,
    pedigree(
      id    = id,
      momid = mother_id,
      dadid = father_id,
      sex   = ifelse(sex == "M", 1, 2)
    )
  )
  
  # Compute kinship matrix from the pedigree
  kin <- kinship(ped)
  
  # Override founder kinship on the diagonal block for listed founders that exist in ped_df
  f_ids <- pars$founder_ids[pars$founder_ids %in% ped_df$id]
  kin[f_ids, f_ids] <- pars$founder_kinship
  
  # Return full kinship matrix
  return(kin)
}

# ============================
# calculate_kinship_genlib
# ============================
# Builds a GENLIB genealogy and computes the kinship matrix.
# Returns a base R matrix. Founder kinship block is then overridden.
calculate_kinship_genlib <- function(pop_df, pars) {
  
  # Unique individuals with father, mother, and sex mapped to GENLIB field names
  ped_df <- pop_df %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::select(ind = id, father = father_id, mother = mother_id, sex)
  
  # Build a shared ID universe so we can recode to integer IDs as GENLIB expects
  id_levels <- c(ped_df$ind, ped_df$father, ped_df$mother) %>% unique()
  
  # Keep a mapping if you need to translate back later
  id_map <- data.frame(
    id_chr = id_levels,
    id_num = as.numeric(factor(id_levels, levels = id_levels)),
    stringsAsFactors = FALSE
  )
  
  # Recode IDs to integers; use 0 to represent unknown parents
  ped_df <- ped_df %>%
    dplyr::mutate(
      ind    = as.numeric(factor(ind,    levels = id_levels)),
      father = as.numeric(factor(father, levels = id_levels)),
      mother = as.numeric(factor(mother, levels = id_levels))
    ) %>%
    dplyr::mutate(
      father = dplyr::case_when(is.na(father) ~ 0, TRUE ~ father),
      mother = dplyr::case_when(is.na(mother) ~ 0, TRUE ~ mother),
      sex    = dplyr::case_when(sex == "M" ~ 1, sex == "F" ~ 2, TRUE ~ NA_real_)
    )
  
  # Create genealogy object
  gen <- GENLIB::gen.genealogy(ped = ped_df)
  
  # Compute kinship (phi) using GENLIB
  k <- GENLIB::gen.phi(gen = gen)
  kin <- as.matrix(k)
  
  # Override founder kinship for founders present in this matrix
  kin[pars$founder_ids, pars$founder_ids] <- pars$founder_kinship
  
  return(kin)
}

# ============================
# calculate_inbreeding
# ============================
# Fills missing individual inbreeding coefficients Fi.
# Rules:
# - Founders listed in pars$founder_ids get Fi = pars$Fp when Fi is NA
# - Individuals with both parents NA and not in founders get Fi = 0
# - Remaining NA Fi are estimated from parental Fi and parent kinship
calculate_inbreeding <- function(pop_df, pars) {
  
  # Parameters expected in 'pars'
  needed_pars <- c("founder_ids", "Fp")
  
  # Start with one row per individual, prioritize those with existing Fi at top
  temp <- pop_df %>%
    dplyr::arrange(is.na(Fi)) %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::mutate(
      # Impute founder Fi when missing
      Fi = dplyr::case_when(
        is.na(Fi) & id %in% pars$founder_ids ~ pars$Fp,
        # If both parents unknown and not a listed founder, set Fi = 0
        is.na(Fi) & is.na(mother_id) & is.na(father_id) ~ 0,
        TRUE ~ Fi
      ),
      priority = 1,  # Track original vs imputed
      t = NULL       # Drop any existing time variable
    )
  
  # Still-missing Fi need estimation via parental info
  missing_Fi <- temp %>%
    dplyr::filter(is.na(Fi)) %>%
    dplyr::mutate(priority = 2)
  
  if (nrow(missing_Fi) > 0) {
    
    # Optional: unique parent pairs that need kinship lookups
    id_filter <- missing_Fi %>%
      dplyr::filter(!duplicated(data.frame(mother_id, father_id))) %>%
      dplyr::select(mother_id, father_id)
    
    id_filter_vec <- c(id_filter$mother_id, id_filter$father_id) %>% unique()
    
    # Compute kinship matrix
    kin <- calculate_kinship(pop_df = pop_df, pars = pars)
    
    # Pull pairwise kinship between each mother and father
    parent_kin <- kin[missing_Fi$mother_id, missing_Fi$father_id]
    if (length(dim(parent_kin)) == 2) parent_kin <- diag(parent_kin)  # One value per child
    
    # Collect Fi for fathers in the same order as missing_Fi
    Fdad <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "M") %>%
      dplyr::pull(Fi)
    names(Fdad) <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "M") %>%
      dplyr::pull(id)
    Fdad <- Fdad[missing_Fi$father_id]
    
    # Collect Fi for mothers in the same order as missing_Fi
    Fmom <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "F") %>%
      dplyr::pull(Fi)
    names(Fmom) <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "F") %>%
      dplyr::pull(id)
    Fmom <- Fmom[missing_Fi$mother_id]
    
    # Standard inbreeding recursion:
    # Fi = (Fi_m + Fi_f)/2 + kin(m, f) * (1 - (Fi_m + Fi_f)/2)
    Fi <- ((Fdad + Fmom) / 2) + (parent_kin * (1 - ((Fdad + Fmom) / 2)))
    
    # Assign back
    missing_Fi$Fi <- Fi
  }
  
  # Merge original and estimated, mark as t = -1
  new <- plyr::rbind.fill(temp, missing_Fi) %>%
    dplyr::mutate(t = -1)
  
  # Tidy and keep only relevant fields
  resu <- tidy_pop_df(new) %>%
    dplyr::select(id, Fi, nz_heritage)
  
  if (any(is.na(resu$Fi))) stop("Some Fi's not calculated")
  
  return(resu)
}

# ============================
# get_kinship_pair
# ============================
# Returns the kinship coefficient between two individuals within a pedigree.
get_kinship_pair <- function(id1, id2, ped) {
  kin_mat <- kinship(ped, ids = c(id1, id2))
  return(kin_mat[id1, id2])
}

# ============================
# pull_named
# ============================
# dplyr::pull with optional names assignment taken from another column.
pull_named <- function(.data, col, names_col = NULL) {
  col <- rlang::enquo(col)
  names_col <- rlang::enquo(names_col)
  
  vec <- dplyr::pull(.data, !!col)
  
  if (!rlang::quo_is_null(names_col)) {
    vec_names <- dplyr::pull(.data, !!names_col)
    names(vec) <- vec_names
  }
  
  return(vec)
}


# =======================
# Processing function
# Loads one .RData result and returns a list of data frames ready to write
# =======================
process_result_file <- function(filename, folder_id){
  suppressMessages({
    load(filename)  # Loads an object called 'output' from the file
    
    # Simulation identifiers
    idx   <- unique(output$pop$i)
    label <- output$run_label
    
    # Keep a copy of full population to detect zero population time
    output$pop0 <- output$pop
    
    # Find time steps where total alive == 0
    zeroNs <- which(
      output$pop0 %>%
        dplyr::group_by(t) %>%
        dplyr::summarise(N = sum(alive)) %>%
        dplyr::pull(N) == 0
    )
    
    # First time step with zero population, if any
    suppressWarnings({
      zeroN <- min(zeroNs)
    })
    
    # Truncate time series to before extinction time
    output$pop <- output$pop %>% dplyr::filter(t < zeroN)
    # output$egg_fate <- output$egg_fate %>% dplyr::filter(t < zeroN)
    output$envir_stoch <- output$envir_stoch %>% dplyr::filter(t < zeroN)
    
    # Identify released individuals (not originating from wild)
    released_individuals <- output$pop %>%
      filter(origin != "Wild") %>%
      pull(id) %>%
      unique()
    
    # Block for computing release deaths is disabled by design (6==9)
    if (6 == 9) {
      if (length(released_individuals) > 0) {
        release_deaths_df <- output$pop %>%
          dplyr::filter(id %in% released_individuals, !alive) %>%
          dplyr::mutate(
            dead_before_first_june = case_when(
              age_release == 0 ~ age == 1,
              age_release  > 0 ~ tsr <= 1
            )
          ) %>%
          dplyr::left_join(output$model_pars$mgmt$release_year_cont) %>%
          dplyr::group_by(age_release) %>%
          dplyr::summarise(
            mortality = mean(dead_before_first_june),
            noDead = sum(dead_before_first_june),
            yr_duration = unique(yr_duration),
            expectedDeaths = n() * (mortality ^ unique(yr_duration))
          )
        release_deaths <- release_deaths_df %>% pull(expectedDeaths) %>% sum()
      } else {
        release_deaths <- 0
      }
    }
    
    # Compute final outcome metrics as a weighted average over the last 5 steps
    gen_final_outcome <- output$pop %>%
      dplyr::filter(t %in% (max(t):(max(t) - 5)), alive) %>%
      dplyr::group_by(t) %>%
      dplyr::summarise(
        N  = sum(alive),
        Kp = mean(nz_heritage),
        Fp = mean(Fi),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Nw = N / sum(N)) %>%
      dplyr::summarise(
        Kp = sum(Kp * Nw),
        Fp = sum(Fp * Nw)
      )
    
    # Full N time series by replicate i
    N_df <- output$pop %>%
      dplyr::group_by(t, i) %>%
      dplyr::summarise(
        N  = sum(alive),
        Fp = sum(Fi * alive) / sum(alive),
        Kp = sum(nz_heritage * alive) / sum(alive),
        .groups = "drop"
      )
    
    # Adult female population (age >= 3) by replicate i
    N_ad_fem_df <- output$pop %>%
      dplyr::filter(sex == "F", age >= 3) %>%
      dplyr::group_by(t, i) %>%
      dplyr::summarise(
        N  = sum(alive),
        Fp = sum(Fi * alive) / sum(alive),
        Kp = sum(nz_heritage * alive) / sum(alive),
        .groups = "drop"
      )
    
    # Calculate multiplicative trend lambda and its geometric mean
    Ns <- pull(N_ad_fem_df, N)
    Ns_2 <- Ns
    if (last(Ns_2) == 0) Ns_2 <- Ns_2[-length(Ns_2)]  # Avoid zero at the end for ratio calc
    lambda <- (Ns_2[-1]) / (Ns_2[-length(Ns_2)])
    trend <- prod(lambda) ^ (1 / length(lambda))
    
    # Extinction status and time
    finalN <- tail(Ns, 1)
    extinct <- finalN <= 2
    time_extinct <- ifelse(extinct, length(Ns), NA)
    
    # Optional egg fate metrics are disabled in this script
    # egg_fate <- output$egg_fate
    
    # Build summary row
    summ <- data.frame(
      extinct = extinct,
      finalN = finalN,
      trend = trend,
      time_extinct = time_extinct,
      # release_deaths = release_deaths,
      # total_eggs_released = total_eggs_released,
      # total_eggs_lost = total_eggs_lost,
      i = idx,
      label = label
    ) %>%
      merge(gen_final_outcome)
    
    # Management and run parameters
    run_pars <- output$run_pars
    # release_sch <- output$model_pars$mgmt$release_schedule
    
    # mgmt <- release_sch %>%
    #   cbind(SuppFeed = as.character(run_pars$SuppFeed)) %>%
    #   mutate(i = idx)
    
    mgmt <- output$model_pars$mgmt
  })
  
  # Return structured results
  resu <- list(
    summary  = summ,
    N_series = N_df,
    # egg_fate = egg_fate,
    mgmt     = mgmt,
    run_pars = as.data.frame(run_pars)
  )
  
  # Tag each output with provenance
  resu <- lapply(resu, function(x){
    x$folder   <- folder_id
    x$filename <- filename
    x$i        <- idx
    return(x)
  })
  
  return(resu)
}