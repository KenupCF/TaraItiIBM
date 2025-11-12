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

