
model_pars$sim<-list()

model_pars$sim$use_genetics<-TRUE

model_pars$sim$n_years<-50
model_pars$sim$n_iter<-10

# model_pars$sim$n_years<-30
# model_pars$sim$n_iter<-5

model_pars$sim$parametric_uncertainty<-T
model_pars$sim$n_samples_quantile_function<-1e4

model_pars$sim$idx_add<-0

model_pars$sim$parallel_across_runs<-TRUE
model_pars$sim$clusters_to_run<-min(64,parallel::detectCores()-2)
model_pars$sim$batching_clusters<-64

model_pars$priors$start_cycle<-data.frame(min=-0.5,max=4.5,dist="unif")

# template for quantile extraction (0 to 1 means all quantiles are sampled, 0.5 means to always get the median)
if(model_pars$sim$parametric_uncertainty){
  qrunif_template<-data.frame(min=0.0,max=1.0,dist="unif")
}else{
  qrunif_template<-data.frame(min=0.5,max=0.5,dist="unif")
}



qFUN <- agg_info$LP$mc_q_functions





#### structure of variables to calculate coefficients

model_pars$priors$Q_juv_surv<-qrunif_template
model_pars$priors$Q_imm_surv<-qrunif_template
model_pars$priors$Q_adu_surv<-qrunif_template
model_pars$priors$Q_hatch_prob<-qrunif_template
model_pars$priors$Q_fledge_prob<-qrunif_template
model_pars$priors$Q_prop_breed<-qrunif_template
model_pars$priors$Q_clutch_size<-qrunif_template

model_pars$priors$Q_juv_surv_gen<-qrunif_template
model_pars$priors$Q_imm_surv_gen<-qrunif_template
model_pars$priors$Q_adu_surv_gen<-qrunif_template
model_pars$priors$Q_hatch_prob_gen<-qrunif_template
model_pars$priors$Q_clutch_size_gen<-qrunif_template
model_pars$priors$Q_fledge_prob_gen<-qrunif_template

model_pars$priors$Q_hatch_prob_ai_f1<-qrunif_template
model_pars$priors$Q_fledge_prob_ai_f1<-qrunif_template
model_pars$priors$Q_juv_surv_ai_f1<-qrunif_template
model_pars$priors$Q_imm_surv_ai_f1<-qrunif_template
model_pars$priors$Q_adu_surv_ai_f1<-qrunif_template
model_pars$priors$Q_clutch_size_ai_f1<-qrunif_template
model_pars$priors$Q_prop_breed_ai_f1<-qrunif_template


model_pars$priors$Q_hatch_prob_ai_f2<-qrunif_template
model_pars$priors$Q_fledge_prob_ai_f2<-qrunif_template


model_pars$priors$Q_release_year<-qrunif_template
model_pars$priors$Q_release_prob<-qrunif_template


