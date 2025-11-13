# ============================
# Management strategy setup
# ============================

# Build a full factorial grid of candidate management strategies.
# Each row defines one management alternative to be simulated.
model_pars$mgmt$strategies <- expand.grid(
  egg_harvest_rate = c(0, 0.1, 0.5, 1),  # proportion of eggs harvested for ex-situ or head-starting
  admix_releases   = c(FALSE, TRUE),      # whether admixed individuals (from external population) are released
  admix_release_freq = 3,
  admix_total_releases = 10,
  admix_prop_released = 0.2,              # proportion of released individuals that are admixed
  admix_age_released  = 2,                # age (years) at which admixed individuals are released
  admix_no_released   = 6,                # number of individuals released per event
  gen_mgmt = c(FALSE                      # toggle for active genetic management
               # , TRUE                    # (left commented: enable later if needed)
  ),
  field = c("SQ", "Field2")               # site/field conditions, “SQ” vs “Field2”
) %>%
  # ============================
# Logical filtering of infeasible or undesired combinations
# ============================
# Remove any combination where Field2 uses gen_mgmt
dplyr::filter(!(field == "Field2" & gen_mgmt)) %>%
  
  # Remove combinations where both admixture releases and genetic management occur simultaneously
  dplyr::filter(!(admix_releases & gen_mgmt)) %>%
  
  # Remove cases where genetic management is paired with high egg harvest (>= 0.5)
  dplyr::filter(!(gen_mgmt & egg_harvest_rate >= 0.5)) %>%
  
  # Remove cases where admixture releases coincide with high egg harvest (>= 0.5)
  dplyr::filter(!(admix_releases & egg_harvest_rate >= 0.5)) %>%
  
  # (optional filters left commented: could restrict more)
  # dplyr::filter(!(field == "Field2" & admix_releases))
  # dplyr::filter(!(field == "Field2" & egg_harvest_rate != min(egg_harvest_rate)))
  
  ungroup() %>%
  
  # Add an identifier column for each management alternative
  dplyr::mutate(alt = 1:n())

# Each row in model_pars$mgmt$strategies now represents one scenario combining:
# - field type (SQ / Field2)
# - egg harvest intensity
# - whether admixed individuals are released
# - release age, number, and proportion
# - whether genetic management is on/off
# - alt: unique integer ID for iteration labeling
