# Set working directory based on available paths
wd<-"/models/TaraIti"

if(!dir.exists(wd)){
  wd<-"D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM"  
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
db_path <- "D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM/Results/bigRunV7/bigRunV7.duckdb"

# db_path <- "D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM/Results/Sensitivity_Analysis/Sensitivity_Analysis.duckdb"

con <- dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)

# Load tables from the DuckDB database
summary <- dbGetQuery(con, "SELECT * FROM summary")
run_pars <- dbGetQuery(con, "SELECT * FROM run_pars")
pop <- dbGetQuery(con, "SELECT * FROM N_series")
surv_rate <- dbGetQuery(con, "SELECT * FROM surv_rate")
reprod_pars <- dbGetQuery(con, "SELECT * FROM reprod_pars")


dbDisconnect(con, shutdown = TRUE)

mgmt <- run_pars %>%
  dplyr::select(alt, field, admix_releases, egg_harvest_rate, gen_mgmt) %>%
  dplyr::filter(!duplicated(alt)) %>%
  dplyr::mutate(
    egg_harvest_desc = paste0(
      dplyr::case_when(
        egg_harvest_rate == 0.0 ~ "no egg harvest",
        egg_harvest_rate == 0.1 ~ "10% egg harvest",
        egg_harvest_rate == 0.5 ~ "50% egg harvest",
        egg_harvest_rate == 1.0 ~ "100% egg harvest",
        TRUE ~ "Unknown"
      )
    ),
    egg_harvest_desc = factor(
      egg_harvest_desc,
      levels = c(
        "no egg harvest",
        "10% egg harvest",
        "50% egg harvest",
        "100% egg harvest",
        "Unknown"
      )
    )
  ) %>%
  dplyr::mutate(
    gen_mgmt_desc = dplyr::case_when(
      gen_mgmt ~ "Genetic management",
      !gen_mgmt ~ "No genetic management"
    ),
    admix_desc = dplyr::case_when(
      admix_releases ~ "Admixing",
      !admix_releases ~ "No admixing"
    )
  ) %>%
  dplyr::mutate(
    gen_desc = dplyr::case_when(
      admix_releases & gen_mgmt ~ paste0(gen_mgmt_desc, ", ", admix_desc),
      admix_releases & !gen_mgmt ~ "Admixing",
      !admix_releases & gen_mgmt ~ "Genetic management",
      !admix_releases & !gen_mgmt ~ "No genetic action"
    ),
    gen_desc = factor(
      gen_desc,
      levels = c(
        "No genetic action",
        "Genetic management",
        "Admixing",
        "Genetic management, Admixing"
      )
    )
  ) %>%
  dplyr::mutate(
    originalAlts =
      (field == "SQ" & egg_harvest_rate == .5 & !gen_mgmt & !admix_releases) |
      (field == "Field2" & egg_harvest_rate == .5 & !gen_mgmt & !admix_releases) |
      (field == "SQ" & egg_harvest_rate == .5 & !gen_mgmt & admix_releases) |
      (field == "SQ" & egg_harvest_rate == .5 & gen_mgmt & !admix_releases) |
      (field == "SQ" & egg_harvest_rate > .5 & !gen_mgmt & !admix_releases),
    sdmAlts =
      originalAlts |
      (field == "SQ" & egg_harvest_rate == 0 & !gen_mgmt & !admix_releases)
  )%>%
  dplyr::mutate(SDM_Strat_Name=case_when(
    !sdmAlts ~ NA_character_,
    (field == "SQ" & egg_harvest_rate == .5 & !gen_mgmt & !admix_releases) ~ "Status Quo",
    (field == "Field2" & egg_harvest_rate == .5 & !gen_mgmt & !admix_releases) ~ "Enhanced Field",
    (field == "SQ" & egg_harvest_rate == .5 & !gen_mgmt & admix_releases)  ~ "Genetic Rescue",
    (field == "SQ" & egg_harvest_rate == .5 & gen_mgmt & !admix_releases) ~ "Genetic Management",
    (field == "SQ" & egg_harvest_rate == 1 & !gen_mgmt & !admix_releases) ~ "Captive Upscaled",
    (field == "SQ" & egg_harvest_rate == 0 & !gen_mgmt & !admix_releases) ~ "No Captive")
   )%>%
  dplyr::mutate(SDM_Strat_Name=factor(SDM_Strat_Name,
                                      levels=c("Status Quo","No Captive","Captive Upscaled",
                                               "Enhanced Field","Genetic Management", "Genetic Rescue")))

sensit<-run_pars%>%
  dplyr::select(q,ad_fem_removed,ad_mal_removed)%>%
  dplyr::filter(!duplicated(q))

pop<-pop%>%
  dplyr::left_join(run_pars%>%select(i,alt,p,q))%>%
  dplyr::rename(TotalN=N)

surv_rate<-surv_rate%>%
  dplyr::left_join(run_pars%>%select(i,alt,p,q))

reprod_pars<-reprod_pars%>%
  dplyr::left_join(run_pars%>%select(i,alt,p,q))

pop_abund_names<-colnames(pop)
pop_abund_names<-pop_abund_names[!pop_abund_names%in%c("folder","filename","alt","p","t","i","q","Fp","Kp")]

fill_list<-lapply(pop_abund_names,function(x){0})
names(fill_list)<-pop_abund_names
fill_list$Kp<-NA
fill_list$Fp<-NA

filled_pop <- pop %>%
  dplyr::select(t, alt, p, q, dplyr::any_of(names(fill_list))) %>%
  tidyr::complete(t, alt, p, q, fill = fill_list)%>%
  # dplyr::left_join(run_pars%>%
  # select(alt,q,p,K))
  dplyr::mutate(Occupancy=BreedingPairs/100)

quasi_extinct_years<-filled_pop%>%
  dplyr::group_by(alt,q,p)%>%
  dplyr::summarise(value=mean(BreedingPairs<=2))%>%
  dplyr::group_by(alt,q)%>%
  dplyr::summarise(avQuasiExtinctionExposure=mean(value,na.rm=T),
                   # lcl50=quantile(value,1-.75,na.rm=T),ucl50=quantile(value,.75,na.rm=T),
                   # lcl90=quantile(value,1-.95,na.rm=T),ucl90=quantile(value,.95,na.rm=T),
                   avQuasiExtinctionExposurelcl95=quantile(value,1-.975,na.rm=T),
                   avQuasiExtinctionExposureucl95=quantile(value,.975,na.rm=T),
                   foo=0)%>%  
  dplyr::mutate(foo=NULL
                # ,Quantity="QuasiExtinctionExposure"
  )%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)

time_to_quasi_extinction<-filled_pop%>%
  dplyr::group_by(alt,q,p)%>%
  dplyr::arrange(t)%>%
  dplyr::summarise(value=min(which(BreedingPairs<=2)))%>%
  dplyr::mutate(value=case_when(value==Inf~NA,TRUE~value))%>%
  dplyr::group_by(alt,q)%>%
  dplyr::summarise(
    # median=quantile(value,.5,na.rm=T),
    avTimeToQuasiExtinction=mean(value,na.rm=T),
    # lcl50=quantile(value,1-.75,na.rm=T),ucl50=quantile(value,.75,na.rm=T),
    # lcl90=quantile(value,1-.95,na.rm=T),ucl90=quantile(value,.95,na.rm=T),
    TimeToQuasiExtinctionlcl95=quantile(value,1-.975,na.rm=T),
    TimeToQuasiExtinctionucl95=quantile(value,.975,na.rm=T),
    foo=0)%>%  
  dplyr::mutate(foo=NULL
                # ,Quantity="TimeToQuasiExtinction"
  )%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)

persistence_over_time<-filled_pop%>%
  dplyr::group_by(alt,q,t)%>%
  dplyr::summarise(Persistence=mean(BreedingPairs>2))%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)


persistence_df <- persistence_over_time %>%
  dplyr::rename(av = Persistence) %>%
  dplyr::mutate(
    median = av,
    lcl50 = NA_real_,
    ucl50 = NA_real_,
    lcl90 = NA_real_,
    ucl90 = NA_real_,
    lcl95 = NA_real_,
    ucl95 = NA_real_,
    Quantity = "Persistence"
  ) %>%
  dplyr::select(
    alt, t, q, Quantity, median, av,
    lcl50, ucl50, lcl90, ucl90, lcl95, ucl95,
    dplyr::everything()
  )


abundance_over_time<-filled_pop%>%
  tidyr::pivot_longer(
    cols = c(dplyr::all_of(names(fill_list)), Occupancy,Fp,Kp),
    names_to = "Quantity",
    values_to = "count"
  )%>%
  dplyr::filter(t>0)%>%
  dplyr::group_by(alt,t,q,Quantity)%>%
  dplyr::summarise(median=quantile(count,.5,na.rm=T),av=mean(count,na.rm=T),
                   lcl50=quantile(count,1-.75,na.rm=T),ucl50=quantile(count,.75,na.rm=T),
                   lcl90=quantile(count,1-.95,na.rm=T),ucl90=quantile(count,.95,na.rm=T),
                   lcl95=quantile(count,1-.975,na.rm=T),ucl95=quantile(count,.975,na.rm=T),
                   foo=0)%>%
  dplyr::mutate(foo=NULL)%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)

surv_rate_avg_time<-surv_rate%>%
  dplyr::group_by(alt,t,q,age_class)%>%
  dplyr::summarise(median=quantile(surv_rate,.5),av=mean(surv_rate),
                   lcl50=quantile(surv_rate,1-.75),ucl50=quantile(surv_rate,.75),
                   lcl90=quantile(surv_rate,1-.95),ucl90=quantile(surv_rate,.95),
                   lcl95=quantile(surv_rate,1-.975),ucl95=quantile(surv_rate,.975),
                   foo=0)%>%
  dplyr::mutate(foo=NULL)%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)%>%
  dplyr::mutate(Quantity=paste("Survival_",age_class,sep=""),age_class=NULL)

reprod_pars_avg_time<-reprod_pars%>%
  dplyr::mutate(egg_hatch_rate=eggs_hatched/eggs_laid,
                fledge_rate=fledged/eggs_hatched)%>%
  tidyr::pivot_longer(
    cols = c(eggs_laid,egg_hatch_rate,fledge_rate),
    names_to = "Quantity",
    values_to = "value"
  )%>%
  dplyr::group_by(alt,t,q,Quantity)%>%
  dplyr::summarise(median=quantile(value,.5,na.rm=T),av=mean(value,na.rm=T),
                   lcl50=quantile(value,1-.75,na.rm=T),ucl50=quantile(value,.75,na.rm=T),
                   lcl90=quantile(value,1-.95,na.rm=T),ucl90=quantile(value,.95,na.rm=T),
                   lcl95=quantile(value,1-.975,na.rm=T),ucl95=quantile(value,.975,na.rm=T),
                   foo=0)%>%
  dplyr::mutate(foo=NULL)%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)


big_tracker_df<-plyr::rbind.fill(reprod_pars_avg_time,
                                 surv_rate_avg_time,
                                 abundance_over_time
                                 ,persistence_df
                                 # ,time_to_quasi_extinction
                                 # ,quasi_extinct_years
)

big_tracker_df$Quantity <- factor(
  big_tracker_df$Quantity,
  levels = c(
    # Reproductive metrics
    "eggs_laid",
    "egg_hatch_rate",
    "fledge_rate",
    
    # Survival rates
    "Survival_juv",
    "Survival_imm",
    "Survival_ad_fem",
    "Survival_ad_mal",
    
    # Demographic groups
    "Juveniles",
    "Immatures",
    "UnpairedFemales",
    "UnpairedMales",
    "AdultFemales",
    "AdultMales",
    "BreedingPairs",
    
    # Population state metrics
    "Occupancy",
    "TotalN",
    
    # Add persistence
    "Persistence",
    
    # Genetics
    "Fp", "Kp"
  )
)


# Create lookup descriptions
variable_descriptions <- c(
  eggs_laid        = "Eggs laid (total)",
  egg_hatch_rate   = "Egg hatch success rate",
  fledge_rate      = "Fledging success rate",
  
  Survival_juv     = "Juvenile survival",
  Survival_imm     = "Immature survival",
  Survival_ad_fem  = "Adult female survival",
  Survival_ad_mal  = "Adult male survival",
  
  Juveniles        = "Number of juveniles",
  Immatures        = "Number of immatures",
  UnpairedFemales  = "Unpaired adult females",
  UnpairedMales    = "Unpaired adult males",
  AdultFemales     = "Adult females",
  AdultMales       = "Adult males",
  BreedingPairs    = "Breeding pairs",
  
  Occupancy        = "Breeding territories occupancy",
  TotalN           = "Total population size",
  
  Persistence      = "Persistence probability",
  
  Fp               = "Inbreeding coefficient (Fp)",
  Kp               = "Provenance coefficient (Kp)"
)

# Build map dataframe
variable_map <- data.frame(
  Quantity = factor(variable_levels,levels=levels(big_tracker_df$Quantity)),
  Quantity_label = variable_descriptions[variable_levels],
  row.names = NULL
)


big_tracker_df<-big_tracker_df%>%
  dplyr::arrange(Quantity,q,alt,t)%>%
  dplyr::left_join(variable_map)

last_points_qt <- big_tracker_df %>%
  dplyr::ungroup()%>%
  dplyr::filter(t == max(t))

classes<-levels(big_tracker_df$Quantity)


plots_all_strats<-list()
plots_sdm_strats<-list()

for(cl in classes){
  
  tempDat<-big_tracker_df%>%
    dplyr::filter(q==1,Quantity==cl)
  
  tempFinalPoints<-last_points_qt%>%
    dplyr::filter(q==1,Quantity==cl)
  
  ylabel<- unique(tempDat$Quantity_label)
  
  plots_all_strats[[cl]]<-ggplot(tempDat, aes(x = t, y = median,
                                                                group = gen_desc,
                                                                color = gen_desc)) +
    geom_ribbon(
      aes(ymin = lcl50, ymax = ucl50, fill = gen_desc),
      alpha = 0.2,       # transparency of the ribbon
      color = NA          # no border line on ribbon
    ) +
    geom_line(alpha = 1,lwd=1) +
    geom_text_repel(
      data = tempFinalPoints,
      aes(label = gen_desc),
      color = "black",
      size = 3,
      direction = "y",      # keeps labels spread vertically
      hjust = 0,            # align to the right of points
      nudge_x = 1,          # move labels a bit to the right
      segment.color = "grey50",  # line color
      segment.size = 0.3,        # line thickness
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    facet_grid(egg_harvest_desc~field,scales = "fixed")+
    labs(
      # title = "Population Size over Time by Alternative",
      x = "Time (t)",
      y = ylabel,
      # color = "Genetic action"
    ) +
    theme_minimal() +
    scale_color_viridis_d(name = "Genetic action", option = "D") +
    scale_fill_viridis_d(name = "Genetic action", option = "D")+
    guides(
      fill = guide_legend("Genetic action"),
      color = guide_legend("Genetic action")
    )+
    theme(
      # plot.margin = margin(5.5, 100, 5.5, 5.5), # increase right margin
      legend.position = "bottom"
    ) +
    coord_cartesian(clip = "off")  # ensures labels beyond plot area are visible
  
  plots_sdm_strats[[cl]]<-ggplot(tempDat, aes(x = t, y = median,
                                              group = SDM_Strat_Name,
                                              color = SDM_Strat_Name)) +
    geom_ribbon(
      aes(ymin = lcl50, ymax = ucl50, fill = SDM_Strat_Name),
      alpha = 0.2,       # transparency of the ribbon
      color = NA          # no border line on ribbon
    ) +
    geom_line(alpha = 1,lwd=1) +
    geom_text_repel(
      data = tempFinalPoints,
      aes(label = SDM_Strat_Name),
      color = "black",
      size = 3,
      direction = "y",      # keeps labels spread vertically
      hjust = 0,            # align to the right of points
      nudge_x = 1,          # move labels a bit to the right
      segment.color = "grey50",  # line color
      segment.size = 0.3,        # line thickness
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    # facet_grid(egg_harvest_desc~field,scales = "fixed")+
    labs(
      # title = "Population Size over Time by Alternative",
      x = "Time (t)",
      y = ylabel,
      # color = "Genetic action"
    ) +
    theme_minimal() +
    scale_color_viridis_d(name = "Strategy", option = "B", na.translate = FALSE,direction=-1) +
    scale_fill_viridis_d(name = "Strategy", option = "B", na.translate = FALSE,direction=-1)+
    guides(
      fill = guide_legend("Strategy"),
      color = guide_legend("Strategy")
    )+
    theme(
      # plot.margin = margin(5.5, 100, 5.5, 5.5), # increase right margin
      legend.position = "bottom"
    ) +
    coord_cartesian(clip = "off")  # ensures labels beyond plot area are visible
  
  
  
}

lapply(names(plots_sdm_strats),function(nm){
  w<-720*3
  h<-480*3
  ggsave(plot=plots_sdm_strats[[nm]],filename = paste0("SDM_",nm,".png"),width = w,height = h,units = "px");
  ggsave(plot=plots_all_strats[[nm]],filename = paste0("All_",nm,".png"),width = w,height = h,units = "px");
  
  })
zip("SDM_outputs_figures.zip", files = list.files(pattern = "^(SDM_|All_).+\\.png$"))



resu<-left_join(summary,run_pars%>%select(i,alt,p,q))

resu_summary<-resu%>%
  dplyr::group_by(alt,q)%>%
  dplyr::summarise(probExt=mean(extinct),avTrend=mean(trend),
                   avN=mean(finalN),
                   sdN=sd(finalN),
                   avFp=mean(Fp),
                   avKp=mean(Kp))%>%
  dplyr::left_join(mgmt)%>%
  dplyr::left_join(sensit)

# Summarise results by strategy, habitat and feeding
resu <- summary %>%
  left_join(run_pars%>%
              dplyr::select(alt,i,q))%>%
  dplyr::group_by(alt,q) %>%
  dplyr::summarise(
    noRuns=n(),
    avFinalN = mean(finalN * (!extinct)),
    lclFinalN = pmax(0,avFinalN - 1.96*sd(finalN* (!extinct))),
    uclFinalN = pmax(0,avFinalN + 1.96*sd(finalN* (!extinct))),
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

# Join quasi-ext variables into resu
resu <- resu %>%
  dplyr::left_join(time_to_quasi_extinction) %>%
  dplyr::left_join(quasi_extinct_years) %>%
  dplyr::mutate(
    reportQuasiExtExposure = sprintf(
      "%.3f (%.3f–%.3f)",
      avQuasiExtinctionExposure,
      avQuasiExtinctionExposurelcl95,
      avQuasiExtinctionExposureucl95
    ),
    reportTimeToQuasiExt = sprintf(
      "%.1f (%.1f–%.1f)",
      avTimeToQuasiExtinction,
      TimeToQuasiExtinctionlcl95,
      TimeToQuasiExtinctionucl95
    )
  )

# Join back to mgmt summary for contextual labels
resu <- mgmt %>%
  dplyr::left_join(resu)

write.csv(resu,file = "Model results summary.csv",row.names = F)

# Consequence table 

# resu is the  biological side of it






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
         Hatch_prob_managed_eggs_field2=hatch_prob_coeff_intercept+hatch_prob_coeff_managed+(hatch_prob_coeff_clutch_no*1)+field2_effect_hatch_prob,
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
