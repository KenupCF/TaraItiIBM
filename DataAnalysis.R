# Set working directory based on available paths
wd<-"/models/TaraIti"

if(!dir.exists(wd)){
  wd<-"D:/03-Work/01-Science/00-Research Projects/Tara Iti/IBM"
  
}
if(!dir.exists(wd)){
  wd<-"C:/Users/Caio.Kenup/TaraItiIBM"
  
}
setwd(wd)

# Load required packages
source("packageLoader.R")
require(ggplot2)
require(duckdb)
require(scales)  # for alpha/lighten
require(RColorBrewer)
require(colorspace)

# Connect to DuckDB database
db_path <- "D:/03-Work/01-Science/00-Research Projects/Tara Iti/IBM/Results/bigGoMoreAlts/bigGoMoreAlts.duckdb"
con <- dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)

# Load tables from the DuckDB database
summary <- dbGetQuery(con, "SELECT * FROM summary")
run_pars <- dbGetQuery(con, "SELECT * FROM run_pars")


# hb_labels<-c("No habitat improvement","Habitat improved")
# sf_labels<-c("Continued","Temporary")
# 
# # Add habitat scenario labels
# run_pars <- run_pars %>%
#   dplyr::mutate(Habitat_Scenario = case_when(
#     prob_imp_for == 0 ~  hb_labels[1],
#     prob_imp_for == 1 ~  hb_labels[2]
#   ))

# Load management strategies
# mgmt <- dbGetQuery(con, "SELECT * FROM mgmt")

dbDisconnect(con, shutdown = TRUE)

mgmt<-run_pars%>%
  select(alt,field,admix_releases,egg_harvest_rate)%>%
  filter(!duplicated(alt))

resu<-left_join(summary,run_pars%>%select(i,alt,p))

resu_summary<-resu%>%
  group_by(alt)%>%
  summarise(probExt=mean(extinct),avTrend=mean(trend),
            avN=mean(finalN),
            avFp=mean(Fp),
            avKp=mean(Kp))%>%
  left_join(mgmt)

# Merge summary data with management and run indexes
dat <- summary %>%
  dplyr::left_join(mgmt_trim%>%select(i,alt))%>%
  dplyr::left_join(mgmt_summ)%>%
  dplyr::left_join(run_pars %>% dplyr::select(i, p, alt))


# Summarise results by strategy, habitat and feeding
resu <- summary %>%
  left_join(run_pars%>%select(alt,i))%>%
  dplyr::group_by(alt) %>%
  dplyr::summarise(
    noRuns=n(),
    avFinalN = mean(finalN * (!extinct)),
    lclFinalN = quantile(finalN * (!extinct), 0.025, na.rm = TRUE),
    uclFinalN = quantile(finalN * (!extinct), 0.975, na.rm = TRUE),
    probPersist = 1 - mean(extinct),
    probExtinct = 1 - probPersist,
    propPopulationsDeclining = mean(trend < 1, na.rm = TRUE),
    avTrend = mean(trend, na.rm = TRUE),
    lclTrend = quantile(trend, 0.025, na.rm = TRUE),
    uclTrend = quantile(trend, 0.975, na.rm = TRUE),
    avFp = mean(Fp),
    lclFp = quantile(Fp, 0.025, na.rm = TRUE),
    uclFp = quantile(Fp, 0.975, na.rm = TRUE),
    avKp = mean(Kp),
    lclKp = quantile(Kp, 0.025, na.rm = TRUE),
    uclKp = quantile(Kp, 0.975, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    reportTrend = sprintf("%.2f (%.2f-%.2f)", avTrend, lclTrend, uclTrend),
    reportFinalN = sprintf("%.2f (%.2f-%.2f)", avFinalN, lclFinalN, uclFinalN),
    reportFp = sprintf("%.2f (%.2f-%.2f)", avFp, lclFp, uclFp),
    reportKp = sprintf("%.2f (%.2f-%.2f)", avKp, lclKp, uclKp)
  )

# Join back to mgmt summary for contextual labels
resu <- mgmt %>%
  dplyr::left_join(resu)

write.csv(resu,file = "Model results summary.csv",row.names = F)

# Create table of vital rates for tara iti

pars2<-run_pars%>%filter(!duplicated(p))%>%
  mutate(
    
         Mal_adult_survival=`Phi:(Intercept)`+`Phi:age_classad`+`Phi:sexmale`+`Phi:age_classad:NewMgmtTRUE`+`Phi:sexmale:NewMgmtTRUE`,
         Fem_adult_survival=`Phi:(Intercept)`+`Phi:age_classad`+`Phi:age_classad:NewMgmtTRUE`+`Phi:sexfemale:NewMgmtTRUE`,
         Imm_survival=`Phi:(Intercept)`+`Phi:age_classimm`+`Phi:age_classimm:NewMgmtTRUE`,
         Juv_survival=`Phi:(Intercept)`+`Phi:age_classjuv:NewMgmtTRUE`,
         
         Mal_adult_survival_field2=`Phi:(Intercept)`+`Phi:age_classad`+`Phi:sexmale`+field2_effect_adu_surv,
         Fem_adult_survival_field2=`Phi:(Intercept)`+`Phi:age_classad`+field2_effect_adu_surv,
         Imm_survival_field2=`Phi:(Intercept)`+`Phi:age_classimm`+field2_effect_imm_surv,
         Juv_survival_field2=`Phi:(Intercept)`+field2_effect_juv_surv,
         
         Mal_adult_survival_admix=`Phi:(Intercept)`+`Phi:age_classad`+`Phi:sexmale`+ `Phi:age_classad:NewMgmtTRUE` + `Phi:sexmale:NewMgmtTRUE`+admix_effect_adu_surv,
         Fem_adult_survival_admix=`Phi:(Intercept)`+`Phi:age_classad`+ `Phi:age_classad:NewMgmtTRUE`+`Phi:sexfemale:NewMgmtTRUE`+admix_effect_adu_surv,
         Imm_survival_admix=`Phi:(Intercept)`+`Phi:age_classimm`+ `Phi:age_classimm:NewMgmtTRUE`+admix_effect_imm_surv,
         Juv_survival_admix=`Phi:(Intercept)`+`Phi:age_classjuv:NewMgmtTRUE`+admix_effect_juv_surv,
         
         Mal_adult_survival_ai=`Phi:(Intercept)`+`Phi:age_classad`+`Phi:sexmale`+ `Phi:age_classad:NewMgmtTRUE` + `Phi:sexmale:NewMgmtTRUE`+ai_effect_adu_surv,
         Fem_adult_survival_ai=`Phi:(Intercept)`+`Phi:age_classad`+ `Phi:age_classad:NewMgmtTRUE`+`Phi:sexfemale:NewMgmtTRUE`+ai_effect_adu_surv,
         Imm_survival_ai=`Phi:(Intercept)`+`Phi:age_classimm`+ `Phi:age_classimm:NewMgmtTRUE`+ai_effect_imm_surv,
         Juv_survival_ai=`Phi:(Intercept)`+`Phi:age_classjuv:NewMgmtTRUE`+ai_effect_juv_surv,
         
         Hatch_prob_managed_eggs=hatch_prob_coeff_intercept+hatch_prob_coeff_managed+(hatch_prob_coeff_clutch_no*1)+hatch_prob_coeff_sq,
         Hatch_prob_managed_eggs_field2=Hatch_prob_managed_eggs-hatch_prob_coeff_sq+field2_effect_hatch_prob,
         Hatch_prob_managed_eggs_admix=Hatch_prob_managed_eggs+admix_effect_hatch_prob,
         Hatch_prob_managed_eggs_ai_parents=Hatch_prob_managed_eggs+ai_effect_hatch_prob,
         
         Fledge_prob_managed_eggs=fledge_prob_coeff_intercept+fledge_prob_coeff_managed+fledge_prob_coeff_sq,
         Fledge_prob_managed_eggs_field2=Fledge_prob_managed_eggs-fledge_prob_coeff_sq+field2_effect_fledge_prob,
         Fledge_prob_managed_eggs_admix=Fledge_prob_managed_eggs+admix_effect_hatch_prob,
         Fledge_prob_managed_eggs_ai_parents=Fledge_prob_managed_eggs+ai_effect_fledge_prob,
         
         Prop_breeding = fb_intercept+fb_newsq,
         Prop_breeding_field2 = fb_intercept+field2_effect_prop_breed,
         Prob_breeding_ai_parents = fb_intercept+ai_effect_prop_breed,
         
         Clutch_size = cs_intercept+cs_newsq,
         Clutch_size_admix = Clutch_size + admix_effect_cs,
         Clutch_size_field2 = cs_intercept+field2_effect_cs,
         Clutch_size_ai_parents = Clutch_size + ai_effect_clutch_size,
         
         Hatch_prob_ai_eggs=ai_hatch_prob,
         Fledge_prob_ai_eggs=ai_fledge_prob,
         foo=1
         
         )

new_par_names<-colnames(pars2)[!colnames(pars2)%in%c(colnames(run_pars),"foo")]

logit_pars<-unique(c(str_subset(new_par_names,"survival"),
                     str_subset(new_par_names,"Fledge_prob"),
                     str_subset(new_par_names,"Hatch_prob")))

logit_pars<-logit_pars[!str_detect(logit_pars,"ai_eggs")]

log_pars<-unique(c(str_subset(new_par_names,"Clutch_size")))

pars2[,logit_pars]<-apply(pars2[,logit_pars],2,inv.logit)
pars2[,log_pars]<-apply(pars2[,log_pars],2,exp)

pars3<-pars2[,new_par_names]




summary_df <- pars3 %>%
  summarise(across(
    everything(),
    list(
      av  = mean,
      lcl = ~quantile(.x, 0.025),
      ucl = ~quantile(.x, 0.975)
    ),
    .names = "{.col}__{.fn}"   # use a separator that doesn't appear in names
  )) %>%
  pivot_longer(
    everything(),
    names_to   = c("variable", ".value"),
    names_sep  = "__"
  )
summary_df

write.csv(summary_df,file = "Parameter estimates.csv",row.names = F)
