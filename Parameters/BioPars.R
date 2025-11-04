

model_pars$bio$inherent<-list()
model_pars$bio$gen<-list()

model_pars$bio$gen$starting_inbreeding<-0.68
model_pars$bio$gen$founder_kinship<-0.25


if(model_pars$sim$parametric_uncertainty){
  model_pars$priors$diploid_eq<-data.frame(min=3,max=15,dist="unif")
}else{
  model_pars$priors$diploid_eq<-data.frame(min=7,max=7,dist="unif")
}


model_pars$bio$StartN_df<-plyr::rbind.fill(list(
    data.frame(age_class=c(21,17,12,11,10,8,7,6,5,3),
                       sex="Male",
                           count=c(1,1,1,5,1,1,2,1,1,2)),
                      data.frame(age_class=c(13,10,9,8,7,5,4,3),
                                 sex="Female",
                                 count=c(1,1,2,1,1,1,2,2)),
    data.frame(
      age_class = c(1, 1, 2, 2),
      sex = c("Female", "Male", "Female", "Male"),
      count = c(1, 0, 4, 5)
    )))

### Results from Thali's thesis

qFUN$thali_juv_surv_sq<-function(x){qpert(p=x,min=.55,mode=.81,max=.93)}
qFUN$thali_adu_surv_sq<-function(x){qpert(p=x,min=.86,mode=.92,max=.95)}
qFUN$thali_imm_surv_sq<-function(x){qpert(p=x,min=.68,mode=.93,max=.99)}
qFUN$thali_hatch_prob_sq<-function(x){qpert(p=x,min=.68,mode=.81,max=.89)}
qFUN$thali_fledge_prob_sq<-function(x){qpert(p=x,min=.52,mode=.70,max=.84)}
qFUN$thali_prop_breed_sq<-function(x){qpert(p=x,min=.64,mode=.72,max=.80)}
qFUN$thali_cs_sq<-function(x){qpert(p=x,min=1.68,mode=1.73,max=1.78)}


qFUN$thali_juv_surv_field2<-function(x){qpert(p=x,min=.53,mode=.82,max=.94)}
qFUN$thali_adu_surv_field2<-function(x){qpert(p=x,min=.85,mode=.92,max=.95)}
qFUN$thali_imm_surv_field2<-function(x){qpert(p=x,min=.64,mode=.93,max=.98)}
qFUN$thali_hatch_prob_field2<-function(x){qpert(p=x,min=.71,mode=.89,max=.96)}
qFUN$thali_fledge_prob_field2<-function(x){qpert(p=x,min=.66,mode=.83,max=.92)}
qFUN$thali_prop_breed_field2<-function(x){qpert(p=x,min=.64,mode=.77,max=.83)}
qFUN$thali_cs_field2<-function(x){qpert(p=x,min=1.70,mode=1.79,max=1.89)}

### Estimates used as baselines for elicitation 

interval_size<-2

juv_surv_mean<-0.687
juv_surv_se<-0.071
qFUN$elicit_juv_surv_anchor<-function(x){qpert(p=x,min=juv_surv_mean-(interval_size*juv_surv_se),
                                               mode=juv_surv_mean,
                                               max=juv_surv_mean+(interval_size*juv_surv_se))}

imm_surv_mean<-0.922
imm_surv_se<-0.063
qFUN$elicit_imm_surv_anchor<-function(x){qpert(p=x,min=imm_surv_mean-(interval_size*imm_surv_se),
                                               mode=imm_surv_mean,
                                               max=min(imm_surv_mean+(interval_size*imm_surv_se),.99))}

adu_surv_mean<-0.895
adu_surv_se<-0.026
qFUN$elicit_adu_surv_anchor<-function(x){qpert(p=x,min=adu_surv_mean-(interval_size*adu_surv_se),
                                               mode=adu_surv_mean,
                                               max=min(adu_surv_mean+(interval_size*adu_surv_se),.99))}


cs_mean<-1.71
cs_sd<-0.45
qFUN$elicit_cs_anchor<-function(x){qpert(p=x,min=cs_mean-(interval_size*cs_sd),
                                               mode=cs_mean,
                                               max=cs_mean+(interval_size*cs_sd))}

hp_mean_logit<-logit(.6)
hp_sd_logit<-0.47
qFUN$elicit_hp_anchor<-function(x){qpert(p=x,min=inv.logit(hp_mean_logit-(interval_size*hp_sd_logit)),
                                         mode=inv.logit(hp_mean_logit),
                                         max=inv.logit(hp_mean_logit+(interval_size*hp_sd_logit)))}


fp_mean_logit<-logit(.58)
fp_sd_logit<-0.22
qFUN$elicit_fp_anchor<-function(x){qpert(p=x,min=inv.logit(fp_mean_logit-(interval_size*fp_sd_logit)),
                                         mode=inv.logit(fp_mean_logit),
                                         max=inv.logit(fp_mean_logit+(interval_size*fp_sd_logit)))}

prob_breed_mean<-0.568
prob_breed_sd<-0.135
qFUN$elicit_prob_breed_anchor<-function(x){qpert(p=x,
                               min=prob_breed_mean-(interval_size*prob_breed_sd),
                              mode=prob_breed_mean,
                              max=prob_breed_mean+(interval_size*prob_breed_sd))}



NoAgeClasses<-3

NoSubPops<-1
subpops<-model_pars$bio$subpops<-LETTERS[1:NoSubPops]

StartN<-c(65)
sum(StartN)


SubPopK<-c(100)

model_pars$bio$carr_capac_df<-data.frame(subpop=subpops,C=SubPopK)


dispersalMat <- matrix(c(
  100), nrow = NoSubPops, byrow = TRUE)/100

rownames(dispersalMat) <- colnames(dispersalMat) <- subpops

model_pars$bio$dispersalMat<-dispersalMat

cor_env_repro_surv <- 0
cor_env_among_pops <- 1

#### Genetics

startK<-0.1
startF_p<-0.1

diploidLethalEquivalents<-data.frame(min=3,max=15,dist="unif")


model_pars$bio$inherent$age_first_breed <- 3
model_pars$bio$inherent$age_last_breed <-	24
model_pars$bio$inherent$max_age <-	24
model_pars$bio$inherent$max_brood_p_year <-	1
model_pars$bio$inherent$max_prog_per_brood <- 2
model_pars$bio$inherent$sex_ratio <-	0.5


model_pars$bio$re_surv_temp<-0.001       
