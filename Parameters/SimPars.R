# ============================
# Simulation controls
# ============================

model_pars$sim <- list()

model_pars$sim$use_genetics <- TRUE        # toggle genetics submodel on/off
model_pars$sim$sensitivity_analysis <- FALSE        # toggle genetics submodel on/off

model_pars$sim$n_years <- 50               # years per simulation run
model_pars$sim$n_iter  <- 1e3              # number of replicate runs

# Alternative quick test settings (kept for convenience)
# model_pars$sim$n_years <- 10
# model_pars$sim$n_iter  <- 2

model_pars$sim$parametric_uncertainty <- TRUE   # draw parameters from priors (TRUE) vs fixed medians (FALSE)
model_pars$sim$n_samples_quantile_function <- 1e4  # samples to build/approx quantile functions

model_pars$sim$idx_add <- 0                # optional offset to tag replicate indices

# Parallelization settings
model_pars$sim$parallel_across_runs <- FALSE
model_pars$sim$clusters_to_run      <- min(96, parallel::detectCores() - 2)  # safety margin: leave 2 cores idle
model_pars$sim$batching_clusters    <- 64                                    # batch size for parallel work scheduling

# Prior: start phase for cyclical effects (not used in this model)
model_pars$priors$start_cycle <- data.frame(min = -0.5, max = 4.5, dist = "unif")

# ============================
# Template for quantile draws
# ============================
# If parametric uncertainty is on, draw anywhere in [0,1];
# otherwise always take the median (0.5) of each distribution.
if (model_pars$sim$parametric_uncertainty) {
  qrunif_template <- data.frame(min = 0.0, max = 1.0, dist = "unif")
} else {
  qrunif_template <- data.frame(min = 0.5, max = 0.5, dist = "unif")
}

# Quantile function library (precomputed MC quantile functions)
qFUN <- agg_info$LP$mc_q_functions


sensitivity_analysis<-expand.grid(ad_fem_removed=c(0,2,4,8),
                                  ad_mal_removed=c(0,2,4,8))%>%
  dplyr::filter((ad_fem_removed==ad_mal_removed)| ad_fem_removed==0)%>%
  dplyr::mutate(q=1:n())

if(!model_pars$sim$sensitivity_analysis){
  
  sensitivity_analysis<-sensitivity_analysis[1,]
  
}

# ============================
# Structure of priors for Q-draws
# Each Q_* is a uniform on [min,max] over the probability (quantile) space,
# used to pick quantiles from the corresponding qFUN distribution.
# Suffixes:
#   *_gen    -> general baseline
#   *_ai_f1  -> AI pathway, first effect (e.g., F1)
#   *_ai_f2  -> AI pathway, second effect (e.g., F2)
# ============================

# Demographic/fitness components (field/wild context)
model_pars$priors$Q_juv_surv     <- qrunif_template
model_pars$priors$Q_imm_surv     <- qrunif_template
model_pars$priors$Q_adu_surv     <- qrunif_template
model_pars$priors$Q_hatch_prob   <- qrunif_template
model_pars$priors$Q_fledge_prob  <- qrunif_template
model_pars$priors$Q_prop_breed   <- qrunif_template
model_pars$priors$Q_clutch_size  <- qrunif_template

# General baselines (could be site-agnostic or reference-level effects)
model_pars$priors$Q_juv_surv_gen    <- qrunif_template
model_pars$priors$Q_imm_surv_gen    <- qrunif_template
model_pars$priors$Q_adu_surv_gen    <- qrunif_template
model_pars$priors$Q_hatch_prob_gen  <- qrunif_template
model_pars$priors$Q_clutch_size_gen <- qrunif_template
model_pars$priors$Q_fledge_prob_gen <- qrunif_template

# Assisted incubation (AI) pathway — first effect block
model_pars$priors$Q_hatch_prob_ai_f1  <- qrunif_template
model_pars$priors$Q_fledge_prob_ai_f1 <- qrunif_template
model_pars$priors$Q_juv_surv_ai_f1    <- qrunif_template
model_pars$priors$Q_imm_surv_ai_f1    <- qrunif_template
model_pars$priors$Q_adu_surv_ai_f1    <- qrunif_template
model_pars$priors$Q_clutch_size_ai_f1 <- qrunif_template
model_pars$priors$Q_prop_breed_ai_f1  <- qrunif_template

# Assisted incubation (AI) pathway — second effect block
model_pars$priors$Q_hatch_prob_ai_f2  <- qrunif_template
model_pars$priors$Q_fledge_prob_ai_f2 <- qrunif_template

# Release process priors (e.g., timing and probability of releases)
model_pars$priors$Q_release_year <- qrunif_template
model_pars$priors$Q_release_prob <- qrunif_template
