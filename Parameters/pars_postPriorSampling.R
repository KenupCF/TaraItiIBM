# all_iterations<-prior_rng

all_iterations<-all_iterations%>%
  dplyr::arrange(i)

model_pars$bio$fb_coeffs<-data.frame(
  "Intercept" = all_iterations$fb_intercept[i],
  "Beta_newsq" = all_iterations$fb_newsq[i],
  "Beta_ai" = all_iterations$ai_effect_prop_breed[i],
  "Beta_field2" = all_iterations$field2_effect_prop_breed[i]
)

model_pars$bio$cs_coeffs<-data.frame(
  "Intercept" = all_iterations$cs_intercept[i],
  "Beta_newsq" = all_iterations$cs_newsq[i],
  "Beta_admix" = all_iterations$admix_effect_cs[i],
  "Beta_ai" = all_iterations$ai_effect_clutch_size[i],
  "Beta_field2" = all_iterations$field2_effect_cs[i]
)

model_pars$bio$hatch_prob_coeff <- 
  
  # Start by creating a data frame for adult survival coefficients
  data.frame(
    "Intercept" = all_iterations$hatch_prob_coeff_intercept[i],             # effect of supplementary feeding
    "Beta_managed" = all_iterations$hatch_prob_coeff_managed[i],
    "Beta_f_age" = all_iterations$hatch_prob_coeff_f_age[i],
    "Beta_clutch_no" = all_iterations$hatch_prob_coeff_clutch_no[i],
    "Beta_newsq"  = all_iterations$hatch_prob_coeff_sq[i],
    "Beta_field2" = all_iterations$field2_effect_hatch_prob[i],
    "Beta_admix"  = all_iterations$admix_effect_hatch_prob[i],
    "Beta_ai"     = all_iterations$ai_effect_hatch_prob[i],
    'RE_time'     = model_pars$bio$hatch_prob_coeff_temp_re
  )


model_pars$bio$fledge_prob_coeff <- 
  
  # Start by creating a data frame for adult survival coefficients
  data.frame(
    "Intercept" = all_iterations$fledge_prob_coeff_intercept[i],             # effect of supplementary feeding
    "Beta_managed" = all_iterations$fledge_prob_coeff_managed[i],
    "Beta_f_age" = all_iterations$fledge_prob_coeff_f_age[i],
    "Beta_clutch_no" = all_iterations$fledge_prob_coeff_clutch_no[i],
    "Beta_newsq"  = all_iterations$fledge_prob_coeff_sq[i],
    "Beta_field2" = all_iterations$field2_effect_fledge_prob[i],
    "Beta_admix"  = all_iterations$admix_effect_fledge_prob[i],
    "Beta_ai"     = all_iterations$ai_effect_fledge_prob[i],
    'RE_time'     = model_pars$bio$fledge_prob_coeff_temp_re
  )


model_pars$bio$surv_coeff <- 
  
  # Start by creating a data frame for adult survival coefficients
  data.frame(
    "Intercept" = all_iterations$`Phi:(Intercept)`[i]
      + 0.6931472
    ,             # effect of supplementary feeding
    "Beta_imm" = all_iterations$`Phi:age_classimm`[i],               # effect of improved foraging
    "Beta_adu"  = all_iterations$`Phi:age_classad`[i],
    "Beta_male" = all_iterations$`Phi:sexmale`[i],
    "Beta_newsq_female" = all_iterations$`Phi:sexfemale:NewMgmtTRUE`[i],
    "Beta_newsq_male" = all_iterations$`Phi:sexmale:NewMgmtTRUE`[i],
    "Beta_newsq_juv" = all_iterations$`Phi:age_classjuv:NewMgmtTRUE`[i],
    "Beta_newsq_imm" = all_iterations$`Phi:age_classimm:NewMgmtTRUE`[i],
    "Beta_newsq_adu" = all_iterations$`Phi:age_classad:NewMgmtTRUE`[i],
    "Beta_admix_juv" = all_iterations$admix_effect_juv_surv[i],
    "Beta_admix_imm" = all_iterations$admix_effect_imm_surv[i],
    "Beta_admix_adu" = all_iterations$admix_effect_adu_surv[i],
    "Beta_ai_juv" = all_iterations$ai_effect_juv_surv[i],
    "Beta_ai_imm" = all_iterations$ai_effect_juv_surv[i],
    "Beta_ai_adu" = all_iterations$ai_effect_juv_surv[i],
    "Beta_field2_juv" = all_iterations$field2_effect_juv_surv[i],
    "Beta_field2_imm" = all_iterations$field2_effect_imm_surv[i],
    "Beta_field2_adu" = all_iterations$field2_effect_adu_surv[i],
    'RE_time' = model_pars$bio$re_surv_temp
  )

# cat()

model_pars$bio$ai_hatch_prob<-all_iterations$ai_hatch_prob[i]
model_pars$bio$ai_fledge_prob<-all_iterations$ai_fledge_prob[i]

model_pars$mgmt$year_admix_releases<-all_iterations$year_admix_releases[i]
model_pars$mgmt$prob_admix_releases<-all_iterations$prob_admix_releases[i]



model_pars$mgmt$field<-all_iterations$field[i]
model_pars$mgmt$egg_harvest_rate<-all_iterations$egg_harvest_rate[i]

model_pars$mgmt$admix_releases<-all_iterations$admix_releases[i]
model_pars$mgmt$admix_prop_released<-all_iterations$admix_prop_released[i]
model_pars$mgmt$admix_age_released<-all_iterations$admix_age_released[i]
model_pars$mgmt$admix_no_released<-all_iterations$admix_no_released[i]

model_pars$mgmt$admix_release_freq<-all_iterations$admix_release_freq[i]
model_pars$mgmt$admix_total_releases<-all_iterations$admix_total_releases[i]
model_pars$mgmt$gen_mgmt<-all_iterations$gen_mgmt[i]


init_pars$fems_to_rm<-all_iterations$ad_fem_removed[i]
init_pars$mals_to_rm<-all_iterations$ad_mal_removed[i]







