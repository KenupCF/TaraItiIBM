# Define the quantile sampling function from the linear-pooled or refitted distributions
# qFUN <- agg_info$LP$mc_q_functions 

# Number of samples to draw from each quantile function (used for Monte Carlo simulation)
n_samples_mc <- model_pars$sim$n_samples_quantile_function

# ============================
# Field 2 effects (vs SQ baseline)
# ============================
# NOTE: Effects are computed on the appropriate link scale:
#  - probabilities on logit scale
#  - clutch size on log scale
# This yields additive effects suitable for GLM-style coefficients.

prior_rng$field2_effect_juv_surv <-
  logit(qFUN$thali_juv_surv_field2(x = prior_rng$Q_juv_surv)) -
  logit(qFUN$thali_juv_surv_sq    (x = prior_rng$Q_juv_surv))

prior_rng$field2_effect_imm_surv <-
  logit(qFUN$thali_imm_surv_field2(x = prior_rng$Q_imm_surv)) -
  logit(qFUN$thali_imm_surv_sq    (x = prior_rng$Q_imm_surv))

prior_rng$field2_effect_adu_surv <-
  logit(qFUN$thali_adu_surv_field2(x = prior_rng$Q_adu_surv)) -
  logit(qFUN$thali_adu_surv_sq    (x = prior_rng$Q_adu_surv))

prior_rng$field2_effect_hatch_prob <-
  logit(qFUN$thali_hatch_prob_field2(x = prior_rng$Q_hatch_prob)) -
  logit(qFUN$thali_hatch_prob_sq    (x = prior_rng$Q_hatch_prob))

prior_rng$field2_effect_fledge_prob <-
  logit(qFUN$thali_fledge_prob_field2(x = prior_rng$Q_fledge_prob)) -
  logit(qFUN$thali_fledge_prob_sq    (x = prior_rng$Q_fledge_prob))

# Difference is on identity scale
prior_rng$field2_effect_prop_breed <-
  (qFUN$thali_prop_breed_field2(x = prior_rng$Q_prop_breed)) -
  (qFUN$thali_prop_breed_sq    (x = prior_rng$Q_prop_breed))

# Clutch size is positive and modeled on log scale
prior_rng$field2_effect_cs <-
  log(qFUN$thali_cs_field2(x = prior_rng$Q_clutch_size)) -
  log(qFUN$thali_cs_sq     (x = prior_rng$Q_clutch_size))

# ============================
# Effects of admixing (elicited distributions)
# ============================
# Mixed-offspring effects compared to elicitation anchors.
# Each qFUN call optionally draws n_samples_mc realizations to approximate a quantile.

prior_rng$admix_effect_juv_surv <-
  logit(qFUN$`juv_survival_Mixed offspring`(prior_rng$Q_juv_surv_gen, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_juv_surv_anchor        (prior_rng$Q_juv_surv_gen))

prior_rng$admix_effect_imm_surv <-
  logit(qFUN$`imm_survival_Mixed offspring`(prior_rng$Q_imm_surv_gen, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_imm_surv_anchor        (prior_rng$Q_imm_surv_gen))

# NOTE: This line uses *juvenile* mixed-offspring for adult survival effect.
# If unintended, consider qFUN$`adu_survival_Mixed offspring` here.
prior_rng$admix_effect_adu_surv <-
  logit(qFUN$`adu_survival_Mixed offspring`(prior_rng$Q_adu_surv_gen, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_adu_surv_anchor        (prior_rng$Q_adu_surv_gen))

prior_rng$admix_effect_cs <-
  log(qFUN$`clutch_size_Mixed offspring`(prior_rng$Q_clutch_size_gen, n_samples = n_samples_mc)) -
  log(qFUN$elicit_cs_anchor               (prior_rng$Q_clutch_size_gen))

prior_rng$admix_effect_hatch_prob <-
  logit(qFUN$`hatch_prob_Mixed offspring`(prior_rng$Q_hatch_prob_gen, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_hp_anchor             (prior_rng$Q_hatch_prob_gen))

prior_rng$admix_effect_fledge_prob <-
  logit(qFUN$`fledge_prob_Mixed offspring`(prior_rng$Q_fledge_prob_gen, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_fp_anchor             (prior_rng$Q_fledge_prob_gen))

# ============================
# Effects of artificial/captive rearing
# ============================

# Direct effect on adult female probability of breeding (identity scale)
prior_rng$ai_effect_prop_breed <-
  qFUN$`prop_breeding_AI females`(prior_rng$Q_prop_breed_ai_f1, n_samples = n_samples_mc) -
  qFUN$elicit_prob_breed_anchor (prior_rng$Q_prop_breed_ai_f1)

# Assisted incubation success probabilities (kept on probability scale; used downstream as-is)
prior_rng$ai_hatch_prob  <- qFUN$`hatch_prob_Artificial incubation`(prior_rng$Q_hatch_prob_ai_f1,  n_samples = n_samples_mc)
prior_rng$ai_fledge_prob <- qFUN$`fledge_prob_Artificial incubation`(prior_rng$Q_fledge_prob_ai_f1, n_samples = n_samples_mc)

# Survival effects for AI pathway on logit scale
prior_rng$ai_effect_juv_surv <-
  logit(qFUN$`juv_survival_Artificial incubation`(prior_rng$Q_juv_surv_ai_f1, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_juv_surv_anchor          (prior_rng$Q_juv_surv_ai_f1))

prior_rng$ai_effect_imm_surv <-
  logit(qFUN$`imm_survival_Artificial incubation`(prior_rng$Q_imm_surv_ai_f1, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_imm_surv_anchor            (prior_rng$Q_imm_surv_ai_f1))

prior_rng$ai_effect_adu_surv <-
  logit(qFUN$`adu_survival_Artificial incubation`(prior_rng$Q_adu_surv_ai_f1, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_adu_surv_anchor            (prior_rng$Q_adu_surv_ai_f1))

# Clutch size effect for AI females (log scale)
prior_rng$ai_effect_clutch_size <-
  log(qFUN$`clutch_size_AI females`(prior_rng$Q_clutch_size_ai_f1, n_samples = n_samples_mc)) -
  log(qFUN$elicit_cs_anchor         (prior_rng$Q_clutch_size_ai_f1))

# ============================
# Effects on offspring of AI females
# ============================

prior_rng$ai_effect_hatch_prob <-
  logit(qFUN$`hatch_prob_AI females`(prior_rng$Q_hatch_prob_ai_f2, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_hp_anchor       (prior_rng$Q_hatch_prob_ai_f2))

prior_rng$ai_effect_fledge_prob <-
  logit(qFUN$`fledge_prob_AI females`(prior_rng$Q_fledge_prob_ai_f2, n_samples = n_samples_mc)) -
  logit(qFUN$elicit_fp_anchor       (prior_rng$Q_fledge_prob_ai_f2))

# ============================
# Release timing and probability (quantile-based)
# ============================

prior_rng$year_admix_releases <- qFUN$Q2_1(prior_rng$Q_release_year, n_samples = n_samples_mc)  # e.g., year index
prior_rng$prob_admix_releases <- qFUN$Q1_1(prior_rng$Q_release_prob, n_samples = n_samples_mc)  # e.g., Bernoulli p

