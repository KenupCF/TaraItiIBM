# Initialize containers for biological parameter groups
model_pars$bio$inherent <- list()
model_pars$bio$gen <- list()

# Genetics starting points
model_pars$bio$gen$starting_inbreeding <- 0.68   # initial mean individual inbreeding (F)
model_pars$bio$gen$founder_kinship     <- 0.25   # average kinship among founders

# Starting population age-sex structure
model_pars$bio$StartN_df <- plyr::rbind.fill(list(
  data.frame(
    age_class = c(21, 17, 12, 11, 10, 8, 7, 6, 5, 3),
    sex  = "Male",
    count = c(1, 1, 1, 5, 1, 1, 2, 1, 1, 2)
  ),
  data.frame(
    age_class = c(13, 10, 9, 8, 7, 5, 4, 3),
    sex  = "Female",
    count = c(1, 1, 2, 1, 1, 1, 2, 2)
  ),
  # Extra age-classes for both sexes
  data.frame(
    age_class = c(1, 1, 2, 2),
    sex  = c("Female", "Male", "Female", "Male"),
    count = c(1, 0, 4, 5)
  )
))

### Results from Thali's thesis used to build quantile functions (PERT)
# Each function returns the quantile of a PERT distribution at p = x

qFUN$thali_juv_surv_sq      <- function(x) { qpert(p = x, min = .55, mode = .81, max = .93) }  # juvenile survival, site SQ
qFUN$thali_adu_surv_sq      <- function(x) { qpert(p = x, min = .86, mode = .92, max = .95) }  # adult survival, SQ
qFUN$thali_imm_surv_sq      <- function(x) { qpert(p = x, min = .68, mode = .93, max = .99) }  # immature survival, SQ
qFUN$thali_hatch_prob_sq    <- function(x) { qpert(p = x, min = .68, mode = .81, max = .89) }  # hatch prob, SQ
qFUN$thali_fledge_prob_sq   <- function(x) { qpert(p = x, min = .52, mode = .70, max = .84) }  # fledge prob, SQ
qFUN$thali_prop_breed_sq    <- function(x) { qpert(p = x, min = .64, mode = .72, max = .80) }  # proportion breeding, SQ
qFUN$thali_cs_sq            <- function(x) { qpert(p = x, min = 1.68, mode = 1.73, max = 1.78) } # clutch size, SQ

qFUN$thali_juv_surv_field2  <- function(x) { qpert(p = x, min = .53, mode = .82, max = .94) }  # juvenile survival, Field2
qFUN$thali_adu_surv_field2  <- function(x) { qpert(p = x, min = .85, mode = .92, max = .95) }  # adult survival, Field2
qFUN$thali_imm_surv_field2  <- function(x) { qpert(p = x, min = .64, mode = .93, max = .98) }  # immature survival, Field2
qFUN$thali_hatch_prob_field2<- function(x) { qpert(p = x, min = .71, mode = .89, max = .96) }  # hatch prob, Field2
qFUN$thali_fledge_prob_field2<-function(x) { qpert(p = x, min = .66, mode = .83, max = .92) }  # fledge prob, Field2
qFUN$thali_prop_breed_field2<- function(x) { qpert(p = x, min = .64, mode = .77, max = .83) }  # proportion breeding, Field2
qFUN$thali_cs_field2        <- function(x) { qpert(p = x, min = 1.70, mode = 1.79, max = 1.89) } # clutch size, Field2

### Anchors for expert elicitation baselines (again via PERT), with interval_size as ±k*SE/SD
interval_size <- 2

# Juvenile survival anchor
juv_surv_mean <- 0.687
juv_surv_se   <- 0.071
qFUN$elicit_juv_surv_anchor <- function(x) {
  qpert(p = x,
        min  = juv_surv_mean - (interval_size * juv_surv_se),
        mode = juv_surv_mean,
        max  = juv_surv_mean + (interval_size * juv_surv_se))
}

# Immature survival anchor (cap max at 0.99)
imm_surv_mean <- 0.922
imm_surv_se   <- 0.063
qFUN$elicit_imm_surv_anchor <- function(x) {
  qpert(p = x,
        min  = imm_surv_mean - (interval_size * imm_surv_se),
        mode = imm_surv_mean,
        max  = min(imm_surv_mean + (interval_size * imm_surv_se), .99))
}

# Adult survival anchor (cap max at 0.99)
adu_surv_mean <- 0.895
adu_surv_se   <- 0.026
qFUN$elicit_adu_surv_anchor <- function(x) {
  qpert(p = x,
        min  = adu_surv_mean - (interval_size * adu_surv_se),
        mode = adu_surv_mean,
        max  = min(adu_surv_mean + (interval_size * adu_surv_se), .99))
}

# Clutch size anchor
cs_mean <- 1.71
cs_sd   <- 0.45
qFUN$elicit_cs_anchor <- function(x) {
  qpert(p = x,
        min  = cs_mean - (interval_size * cs_sd),
        mode = cs_mean,
        max  = cs_mean + (interval_size * cs_sd))
}

# Hatch probability anchor on logit scale then mapped back
hp_mean_logit <- logit(.6)
hp_sd_logit   <- 0.47
qFUN$elicit_hp_anchor <- function(x) {
  qpert(p = x,
        min  = inv.logit(hp_mean_logit - (interval_size * hp_sd_logit)),
        mode = inv.logit(hp_mean_logit),
        max  = inv.logit(hp_mean_logit + (interval_size * hp_sd_logit)))
}

# Fledge probability anchor on logit scale
fp_mean_logit <- logit(.58)
fp_sd_logit   <- 0.22
qFUN$elicit_fp_anchor <- function(x) {
  qpert(p = x,
        min  = inv.logit(fp_mean_logit - (interval_size * fp_sd_logit)),
        mode = inv.logit(fp_mean_logit),
        max  = inv.logit(fp_mean_logit + (interval_size * fp_sd_logit)))
}

# Probability of breeding anchor
prob_breed_mean <- 0.568
prob_breed_sd   <- 0.135
qFUN$elicit_prob_breed_anchor <- function(x) {
  qpert(p = x,
        min  = prob_breed_mean - (interval_size * prob_breed_sd),
        mode = prob_breed_mean,
        max  = prob_breed_mean + (interval_size * prob_breed_sd))
}

# Structural settings
NoAgeClasses <- 3                 # number of age classes in model (used elsewhere)
NoSubPops    <- 1                 # number of subpopulations
subpops <- model_pars$bio$subpops <- LETTERS[1:NoSubPops]  # label subpops as A, B, ...

# Starting N per subpopulation (sum for check)
StartN <- c(65)
sum(StartN)                        # quick sanity check

# Carrying capacity per subpopulation
SubPopK <- c(100)
model_pars$bio$carr_capac_df <- data.frame(subpop = subpops, C = SubPopK)

# Dispersal matrix between subpops (here 1x1, so no movement beyond A)
dispersalMat <- matrix(c(100), nrow = NoSubPops, byrow = TRUE) / 100
rownames(dispersalMat) <- colnames(dispersalMat) <- subpops
model_pars$bio$dispersalMat <- dispersalMat

# Environmental correlation placeholders
cor_env_repro_surv <- 0  # correlation between reproduction and survival stochasticity
cor_env_among_pops <- 1  # correlation among subpopulations

#### Genetics
startK  <- 0.1  # initial kinship scaling placeholder
startF_p <- 0.1 # initial inbreeding proportion placeholder

# Prior for diploid lethal equivalents (redundant with priors above, but kept explicit)
diploidLethalEquivalents <- data.frame(min = 3, max = 15, dist = "unif")

# Life history invariants
model_pars$bio$inherent$age_first_breed   <- 3
model_pars$bio$inherent$age_last_breed    <- 24
model_pars$bio$inherent$max_age           <- 24
model_pars$bio$inherent$max_brood_p_year  <- 1
model_pars$bio$inherent$max_prog_per_brood<- 2
model_pars$bio$inherent$sex_ratio         <- 0.5


model_pars$bio$catastrophe_probability <- 1/100


# Survival coefficient random effect from analysis object
model_pars$bio$re_surv_temp <- surv_data_analysis$bio$survival_coeff_temp_re
