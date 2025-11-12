# ============================
# calculate_kinship
# ============================
# Builds a pedigree from pop_df and returns a kinship matrix.
# If rm_non_breeders = TRUE, removes individuals that never appear as a parent.
# Founder kinship for specified founders is overridden using pars.
calculate_kinship <- function(pop_df, pars, rm_non_breeders = TRUE) {
  require(kinship2)  # pedigree(), kinship()
  
  # Expected in pars: founder_ids (character), founder_kinship (scalar)
  # You can enforce with a check if needed.
  
  if (rm_non_breeders) {
    # Keep only individuals that were listed at least once as a parent
    pop_df <- pop_df %>%
      dplyr::filter(id %in% unique(c(pop_df$mother_id, pop_df$father_id)))
  }
  
  # One row per individual with pedigree fields
  ped_df <- pop_df %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::select(id, mother_id, father_id, sex)
  
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
  
  # Pairwise kinship coefficients from the pedigree
  kin <- kinship(ped)
  
  # If supplied founders exist in this pedigree, set their mutual kinship
  f_ids <- pars$founder_ids[pars$founder_ids %in% ped_df$id]
  kin[f_ids, f_ids] <- pars$founder_kinship
  
  return(kin)
}


# ============================
# calculate_kinship_genlib
# ============================
# Builds a GENLIB genealogy and computes the kinship matrix.
# Returns a base R matrix. Founder kinship block is then overridden.
calculate_kinship_genlib <- function(pop_df, pars) {
  # Unique individuals mapped to GENLIB field names
  ped_df <- pop_df %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::select(ind = id, father = father_id, mother = mother_id, sex)
  
  # Create a shared ID universe (character -> integer codes)
  id_levels <- c(ped_df$ind, ped_df$father, ped_df$mother) %>% unique()
  
  # Optional: mapping table if you need to translate back after
  id_map <- data.frame(
    id_chr = id_levels,
    id_num = as.numeric(factor(id_levels, levels = id_levels)),
    stringsAsFactors = FALSE
  )
  
  # Recode IDs; GENLIB uses 0 for unknown parents
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
  
  # Build genealogy and compute kinship (phi)
  gen <- GENLIB::gen.genealogy(ped = ped_df)
  k <- GENLIB::gen.phi(gen = gen)
  kin <- as.matrix(k)
  
  # Override founders (note: indices must match integer recoding)
  kin[pars$founder_ids, pars$founder_ids] <- pars$founder_kinship
  
  return(kin)
}


# ============================
# calculate_inbreeding
# ============================
# Fills missing individual inbreeding coefficients Fi.
# Rules:
# - Founders in pars$founder_ids get Fi = pars$Fp when Fi is NA
# - If both parents are NA and not a founder, Fi := 0
# - Remaining NA Fi are computed via Fi = (Fi_m + Fi_f)/2 + kin(m,f)*(1 - (Fi_m + Fi_f)/2)
calculate_inbreeding <- function(pop_df, pars) {
  # Expected in pars: founder_ids, Fp
  
  # Start with unique individuals; prefer rows where Fi already present
  temp <- pop_df %>%
    dplyr::arrange(is.na(Fi)) %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::mutate(
      Fi = dplyr::case_when(
        is.na(Fi) & id %in% pars$founder_ids ~ pars$Fp,                       # founders: set Fi
        is.na(Fi) & is.na(mother_id) & is.na(father_id) ~ 0,                  # unknown parents: set 0
        TRUE ~ Fi
      ),
      priority = 1,  # mark original/imputed rows
      t = NULL       # drop time if present
    )
  
  # Identify those still missing Fi after founder/unknown-parent rules
  missing_Fi <- temp %>%
    dplyr::filter(is.na(Fi)) %>%
    dplyr::mutate(priority = 2)
  
  if (nrow(missing_Fi) > 0) {
    # Optional: identify unique parent pairs involved
    id_filter <- missing_Fi %>%
      dplyr::filter(!duplicated(data.frame(mother_id, father_id))) %>%
      dplyr::select(mother_id, father_id)
    id_filter_vec <- c(id_filter$mother_id, id_filter$father_id) %>% unique()
    
    # Compute kinship matrix for whole set
    kin <- calculate_kinship(pop_df = pop_df, pars = pars)
    
    # Extract pairwise kinship for each child’s parent pair
    parent_kin <- kin[missing_Fi$mother_id, missing_Fi$father_id]
    if (length(dim(parent_kin)) == 2) parent_kin <- diag(parent_kin)  # one value per child
    
    # Vector of fathers’ Fi aligned to missing_Fi rows
    Fdad <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "M") %>%
      dplyr::pull(Fi)
    names(Fdad) <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "M") %>%
      dplyr::pull(id)
    Fdad <- Fdad[missing_Fi$father_id]
    
    # Vector of mothers’ Fi aligned to missing_Fi rows
    Fmom <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "F") %>%
      dplyr::pull(Fi)
    names(Fmom) <- temp %>%
      dplyr::filter(!is.na(Fi), sex == "F") %>%
      dplyr::pull(id)
    Fmom <- Fmom[missing_Fi$mother_id]
    
    # Recursion formula for offspring inbreeding
    Fi <- ((Fdad + Fmom) / 2) + (parent_kin * (1 - ((Fdad + Fmom) / 2)))
    missing_Fi$Fi <- Fi
  }
  
  # Combine and tidy; mark t = -1 to indicate derived values
  new <- plyr::rbind.fill(temp, missing_Fi) %>%
    dplyr::mutate(t = -1)
  
  resu <- tidy_pop_df(new) %>%
    dplyr::select(id, Fi, nz_heritage)
  
  if (any(is.na(resu$Fi))) stop("Some Fi's not calculated")
  return(resu)
}


# ============================
# get_kinship_pair
# ============================
# Returns kinship coefficient between two IDs using a kinship2 pedigree.
get_kinship_pair <- function(id1, id2, ped) {
  kin_mat <- kinship(ped, ids = c(id1, id2))
  return(kin_mat[id1, id2])
}


# ============================
# pull_named
# ============================
# dplyr::pull with optional names from another column.
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
# process_result_file
# =======================
# Loads one .RData (expects object 'output') and returns tidy summaries and series.
# Adds provenance fields (folder, filename, i) to each returned data frame.
process_result_file <- function(filename, folder_id) {
  suppressMessages({
    load(filename)  # must load an object named 'output'
    
    # Identify replicate index and label
    idx   <- unique(output$pop$i)
    label <- output$run_label
    
    # Keep copy to find extinction time
    output$pop0 <- output$pop
    
    # Find time steps where total alive == 0
    zeroNs <- which(
      output$pop0 %>%
        dplyr::group_by(t) %>%
        dplyr::summarise(N = sum(alive), .groups = "drop") %>%
        dplyr::pull(N) == 0
    )
    
    # First time step with zero population (NA if never zero)
    suppressWarnings({
      zeroN <- min(zeroNs)
    })
    
    # Truncate series prior to extinction time (keeps NA -> no truncation)
    output$pop         <- output$pop         %>% dplyr::filter(t < zeroN)
    # output$egg_fate  <- output$egg_fate    %>% dplyr::filter(t < zeroN)
    output$envir_stoch <- output$envir_stoch %>% dplyr::filter(t < zeroN)
    
    # Identify released individuals (origin != "Wild")
    released_individuals <- output$pop %>%
      dplyr::filter(origin != "Wild") %>%
      dplyr::pull(id) %>%
      unique()
    
    # Block kept intentionally disabled (conditional always FALSE)
    if (6 == 9) {
      if (length(released_individuals) > 0) {
        release_deaths_df <- output$pop %>%
          dplyr::filter(id %in% released_individuals, !alive) %>%
          dplyr::mutate(
            dead_before_first_june = dplyr::case_when(
              age_release == 0 ~ age == 1,
              age_release  > 0 ~ tsr <= 1
            )
          ) %>%
          dplyr::left_join(output$model_pars$mgmt$release_year_cont) %>%
          dplyr::group_by(age_release) %>%
          dplyr::summarise(
            mortality     = mean(dead_before_first_june),
            noDead        = sum(dead_before_first_june),
            yr_duration   = unique(yr_duration),
            expectedDeaths = n() * (mortality ^ unique(yr_duration)),
            .groups = "drop"
          )
        release_deaths <- release_deaths_df %>% dplyr::pull(expectedDeaths) %>% sum()
      } else {
        release_deaths <- 0
      }
    }
    
    # Weighted means over last ~5 time steps for genetics outcomes
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
    
    # Total N series by replicate
    N_df <- output$pop %>%
      dplyr::group_by(t, i) %>%
      dplyr::summarise(
        N  = sum(alive),
        Fp = sum(Fi * alive) / sum(alive),
        Kp = sum(nz_heritage * alive) / sum(alive),
        .groups = "drop"
      )
    
    # Adult female series (age >= 3), used for trend & extinction metrics
    N_ad_fem_df <- output$pop %>%
      dplyr::filter(sex == "F", age >= 3) %>%
      dplyr::group_by(t, i) %>%
      dplyr::summarise(
        N  = sum(alive),
        Fp = sum(Fi * alive) / sum(alive),
        Kp = sum(nz_heritage * alive) / sum(alive),
        .groups = "drop"
      )
    
    # Multiplicative growth rate lambda and its geometric mean
    Ns <- dplyr::pull(N_ad_fem_df, N)
    Ns_2 <- Ns
    if (dplyr::last(Ns_2) == 0) Ns_2 <- Ns_2[-length(Ns_2)]  # avoid trailing zero
    lambda <- (Ns_2[-1]) / (Ns_2[-length(Ns_2)])
    trend  <- prod(lambda) ^ (1 / length(lambda))
    
    # Extinction flag and time (adult females <= 2)
    finalN <- tail(Ns, 1)
    extinct <- finalN <= 2
    time_extinct <- ifelse(extinct, length(Ns), NA)
    
    # Build compact summary row
    summ <- data.frame(
      extinct = extinct,
      finalN = finalN,
      trend = trend,
      time_extinct = time_extinct,
      i = idx,
      label = label
    ) %>%
      merge(gen_final_outcome)
    
    # Keep management and run parameters for reference
    run_pars <- output$run_pars
    mgmt     <- output$model_pars$mgmt
  })
  
  # Pack outputs
  resu <- list(
    summary  = summ,
    N_series = N_df,
    # egg_fate = egg_fate,
    mgmt     = mgmt,
    run_pars = as.data.frame(run_pars)
  )
  
  # Tag with provenance
  resu <- lapply(resu, function(x) {
    x$folder   <- folder_id
    x$filename <- filename
    x$i        <- idx
    x
  })
  
  return(resu)
}
