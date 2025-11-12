# ------------------------------------------------------------------------------
# OVERVIEW ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# generate the sightings records needed by RMark/MARK.

# Code  by Nicole Whitelock @ the University of Auckland and Finnbar Lee @ Cawthron,
# in 2024, finnbar.lee@cawthron.org.nz. 


surv_data_analysis<-list(priors=list(),bio=list())


# ------------------------------------------------------------------------------
# Define start of dataset ------------------------------------------------------
# ------------------------------------------------------------------------------

# min_date <- "1990-01-01"
min_date <- "2010-01-01"

# ------------------------------------------------------------------------------
# load data --------------------------------------------------------------------
# ------------------------------------------------------------------------------
sightings <- read_xlsx(path = "./Data/Sightings_data_1991_2024.xlsx")
capture_morphometrics <- read.csv("./Data/captures_morphometrics_data_1990_2025.csv")
egg_data <- read.csv("./Data/egg_level_data_2000_2024.csv")


# find birds that fledged to juvenile

# from the egg_level_data.csv I have:

# 137 where fledged == 1
# of these 122 have a fledgling_id

# remove birds fledged in most recent season (2022/2023) - no chance to resight in the following Q4
# now 113 birds 

# captive reared birds - using old captive rearing method - these all died quickly
# now 112 birds
# cap_b <- read_rds("captive_vec.rds")
# cap_b <- as.integer(gsub("C", "", cap_b))

cap_b <- egg_data %>% filter(chick_captive_reared == 1)
cap_b <- gsub(x=cap_b$fledgling_id,"C","")
cap_b <- cap_b[!is.na(cap_b)]

# odd date issues with a couple birds, removed these
# 69018, 75815, 75897, 75804
# now 108 birds

# get fledging year/quarter
fledg_birds <- egg_data %>%
  filter(fledged == 1) %>%
  # filter(season != "2022/23") %>%
  mutate(fledgling_id = stringi::stri_trans_general(fledgling_id, "latin-ascii"),
         fledgling_id = gsub("C-69066", "C69066'", fledgling_id),
         fledgling_id = gsub("C69090)", "C69090", fledgling_id),
         fledgling_id = gsub("C69054�", "C69054", fledgling_id)) %>%
  mutate(fledgling_id = gsub("C", "", fledgling_id),
         fledgling_id = as.integer(fledgling_id),
         date  = as.Date(estimated_lay , tryFormats = c("%d/%m/%Y"))) %>%
  dplyr::select(fledgling_id,  date) %>%
  filter(!is.na(fledgling_id)) %>%
  filter(!fledgling_id %in% cap_b) %>%
  filter(!fledgling_id %in% c(69018, 75815, 75897, 75804)) %>%
  dplyr::rename(band_number = fledgling_id)


# ------------------------------------------------------------------------------
# Sightings of parents at nests  -----------------------------------------------
# ------------------------------------------------------------------------------
# inferred based on egg lay date, date doesn't need to be perfect, as we are
# binning observations into yearly quarters

parent_sight <- egg_data %>%
  mutate(female_id = gsub("C", "", female_id),
         male_id = gsub("C", "", male_id),
         female_id = as.integer(female_id),
         male_id = as.integer(male_id),
         date  = as.Date(estimated_lay , tryFormats = c("%d/%m/%Y"))) %>%
  dplyr::select(date, female_id, male_id) %>%
  distinct() %>%
  pivot_longer(cols = female_id:male_id, values_to = "band_number") %>%
  dplyr::select(-name) %>%
  filter(band_number %in% fledg_birds$band_number)

# ------------------------------------------------------------------------------
# Sightings data ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# direct observations of birds

# captures
# ignore warning - character to int
sight2 <- sightings %>%
  clean_names() %>%
  dplyr::select(c_band_number, sight_date) %>%
  dplyr::rename(band_number = c_band_number,
                date = sight_date) %>%
  mutate(band_number = as.integer(band_number)) %>%
  drop_na %>%
  filter(band_number %in% fledg_birds$band_number) %>%
  rbind(fledg_birds) %>%
  filter(date >= min_date)

# ------------------------------------------------------------------------------
# Combine and organise data ----------------------------------------------------
# ------------------------------------------------------------------------------

# combine
sight2 <- rbind(fledg_birds, parent_sight, sight2) %>%
  na.omit()

# create a quarterly sequence from earliest record to latest
sight3 <- data.frame(date = seq.Date(min(sight2$date),
                                     max(sight2$date),
                                     by = "quarter"),
                     band_number = 1) %>%
  rbind(sight2) %>%
  mutate(date_qtr = zoo::as.yearqtr(date),
         date_qtr2 = lubridate::quarter(date)) %>%
  filter(date_qtr2 == 4) %>%
  dplyr::select(-date) %>%
  distinct()


# convert to wide format
sight5 <- sight3 %>%
  arrange(date_qtr) %>%
  dplyr::select(-date_qtr2) %>%
  distinct() %>%
  mutate(filler = 1) %>%
  pivot_wider(names_from = date_qtr , values_from = filler) %>%
  filter(band_number != 1)

# replace NA with 0
sight5[is.na(sight5)] <- 0

# final sightings data set
sight_record <- data.frame(id = sight5$band_number,
                           ch = do.call(paste, c(sight5[2:ncol(sight5)], sep="")))


# get sex info
sex_df <- capture_morphometrics %>%
  clean_names() %>%
  mutate(band_number  = gsub("C", "", band_number ),
         band_number  = as.integer(band_number )) %>%
  dplyr::select(band_number, sex_verified) %>%
  filter(band_number %in% sight_record$id)

# combine
sight_record <- sight_record %>%
  left_join(sex_df, by = c("id" = "band_number"))

# tidy
colnames(sight_record)[3] <- "sex"
sight_record$sex[is.na(sight_record$sex)] <- "unknown"

# save
# write.csv(sight_record,
          # file = "./data_RDA/sightings_history.csv",
          # row.names = F)


last_year<-2025

bird_ch<-sight_record
# # post 2010 dataset
# pva_params <- read.csv("./data/PVA/created/param_estimates_post_2010.csv")
# bird_ch <- read.csv("./data/PVA/raw/sightings_history_2010.csv",
#                     colClasses = c("integer", "character", "factor"))


# remove unknown sex
bird_ch <- bird_ch %>%
  filter(sex != "unknown")

# ------------------------------------------------------------------------------
# Check for GOF ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# perform goodness-of-fit test for for the Cormack-Jolly-Seber model. 

# Get frequencies
bird_ch_tst <- bird_ch %>%
  group_by(ch, sex) %>%
  summarise(Freq = n())

# convert to matrix
bird.hist <- matrix(as.numeric(unlist(strsplit(as.character(bird_ch_tst$ch),""))),
                    nrow = length(bird_ch_tst$ch),
                    byrow = T)
bird.freq <- bird_ch_tst$Freq
bird.group <- bird_ch_tst$sex

# prep for tests
mask <- (bird.group == "female")
bird.fem.hist <- bird.hist[mask,]
bird.fem.freq <- bird.freq[mask]

mask <- (bird.group == "male")
bird.mal.hist <- bird.hist[mask,]
bird.mal.freq <- bird.freq[mask]

# design data
hade.process <- process.data(bird_ch, model = "CJS", groups = "sex")
hade.ddl <- make.design.data(hade.process)

# run tests - not significant - all good
R2ucare::overall_CJS(bird.fem.hist, bird.fem.freq)

# ------------------------------------------------------------------------------
# Set up models ----------------------------------------------------------------
# ------------------------------------------------------------------------------

# add age classes to P and Phi
hade.ddl2 <- add.design.data(hade.process,
                             hade.ddl,
                             parameter = "Phi",
                             type = "age",
                             bins = c(0,1,2,40),
                             right = F, 
                             name = "age_class",
                             replace = T)

levels(hade.ddl2$Phi$age_class) <- c("juv", "imm", "ad")

hade.ddl2 <- add.design.data(hade.process,
                             hade.ddl2,
                             parameter = "p",
                             type = "age",
                             bins = c(0,1,2,40),
                             right = F,
                             name = "age_class",
                             replace = T)

levels(hade.ddl2$p$age_class) <- c("juv", "imm", "ad")

years_cov<-data.frame(year=rev(seq(from=last_year,by=-1,
                                   length.out=length(unique(hade.ddl2$Phi$time)))),
                      time=unique(hade.ddl2$Phi$time))%>%
  mutate(NewMgmt=year%in%2021:2025)

hade.ddl2$Phi<-left_join(hade.ddl2$Phi,years_cov)
# ------------------------------------------------------------------------------
# Run models -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# define models to run - was having issues so have created the simplest set of models
run_models <- function() {
  
  # ----------------------------------------------------------------------------
  # 1990 dataset ---------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # note initially run with all model combinations - then rerun
  # with models that contribute 95% of Akike weights
  
  # survival
  #Phi.dot<-list(formula = ~1)
  #Phi.sex<-list(formula = ~sex)
  Phi.age_class<-list(formula = ~age_class)
  Phi.age_class_sex<-list(formula = ~age_class+sex)
  Phi.age_class_sex_sq<-list(formula = ~age_class+sex+NewMgmt:sex+NewMgmt:age_class)
  Phi.age_class_sq<-list(formula = ~age_class+NewMgmt:age_class)
  
  #Phi.time <- list(formula = ~time)
  #Phi.time.sex <- list(formula = ~time+sex)
  #Phi.time.age_class <- list(formula = ~time+age_class)
  #Phi.age_class.time_sex <- list(formula = ~age_class+time+sex)
  
  # recapture
  #p.dot<-list(formula = ~1)
  #p.age_class<-list(formula = ~age_class)
  
  #p.time <- list(formula = ~time)
  
  p.time.age_class <- list(formula = ~time+age_class)
  
  cml <- create.model.list("CJS")
  
  model.list <- mark.wrapper(cml,
                             data = hade.process,
                             ddl = hade.ddl2,
                             output = F,
                             silent = T)
  
  return(model.list)
}

# run models
mod_res <- run_models()

# check over dispersion
mod_test <- release.gof(hade.process)

# overdispersion factor (chat)
mod_overdisp <- mod_test[3,1] / mod_test[3,2]

mod_overdisp

# # adjust model weighting if overdispersion present
mod_res <- adjust.chat(mod_overdisp, mod_res)


## model averaging the coefficients

m<-mod_res$Phi.age_class.p.time.age_class
phi_coeffs<-lapply(mod_res,function(m){
  rn<-rownames(m$results$beta)
  rn<-rn[str_detect(rn,"Phi")]
  return(rn)
})%>%
  unlist%>%
  unique

coeffs_full_set<-lapply(names(mod_res),function(nm){
  m<-mod_res[[nm]]
  if("mark"%in%class(m)){
    coeffs<-m$results$beta
    coeffs<-coeffs[str_detect(rownames(coeffs),"Phi"),]
    zero_coeffs<-setdiff(y = rownames(coeffs),x = phi_coeffs)
    # Create a data frame with the same number of rows as zero_coeffs
    newcoeffs <- data.frame(estimate = rep(0, length(zero_coeffs)),
                            se = rep(0, length(zero_coeffs)),
                            lcl = rep(0, length(zero_coeffs)),
                            ucl = rep(0, length(zero_coeffs)))
    
    # Assign rownames only if there are any
    if (length(zero_coeffs) > 0) {
      rownames(newcoeffs) <- zero_coeffs
      coeffs<-rbind(coeffs,newcoeffs)    
      
    }
    coeffs$model<-m$model.name
    
    coeffs$weight<-mod_res$model.table%>%
      filter(model==m$model.name)%>%
      pull(weight)
    resu<-coeffs
    resu$par<-rownames(resu)
  } else{
    resu<-NULL
  }
  return(resu)
})%>%
  plyr::rbind.fill()%>%
  mutate(var=se^2)

estimate_avg<-coeffs_full_set%>%
  group_by(par)%>%
  summarise(av_estimate=sum(estimate*weight))

estimate_se<-coeffs_full_set%>%
  left_join(estimate_avg)%>%
  group_by(par)%>%
  summarise(av_var=sum(weight*(var+((estimate-estimate)^2))))%>%
  mutate(av_se=sqrt(av_var))

surv_simulated_re<-estimate_se$av_se[1]/2

estimate_se$av_se[1]<-estimate_se$av_se[1]-surv_simulated_re
estimate_se$av_var[1]<-estimate_se$av_se[1]^2

mod_av_coeffs<-estimate_avg%>%
  left_join(estimate_se)

# average models 
mod_av <- model.average(mod_res, "Phi",vcv=T)


# tidy results
mod_params <- mod_av$estimate %>%
  filter(age%in%0:2)%>%
  left_join(years_cov)%>%
  filter(!duplicated(data.frame(sex,age,NewMgmt,estimate)))%>%
  arrange(group,age,NewMgmt)%>%
  left_join(data.frame(stage=factor(c("juvenile", "immature", "adult"),
                                    levels = c("juvenile", "immature", "adult")),
                       age=factor(0:2)))

surv_data_analysis$bio$survival_coeff_temp_re<-surv_simulated_re

surv_data_analysis$bio$survival<-mod_params
surv_data_analysis$bio$survival_coeffs<-mod_av_coeffs

for(i in 1:nrow(surv_data_analysis$bio$survival_coeffs)) {
  
  parName<-surv_data_analysis$bio$survival_coeffs$par[i]
  
  dist_mean <- surv_data_analysis$bio$survival_coeffs$av_estimate[i]
  dist_sd <- surv_data_analysis$bio$survival_coeffs$av_se[i]
  dist_name <- "norm"
  
  surv_data_analysis$priors[[parName]]<-data.frame(mean=dist_mean,sd=dist_sd,dist="norm")

}




# plot
gg<-ggplot(mod_params,
       aes(sex, estimate, 
           ymin = lcl,
           ymax = ucl,
           colour = stage)) +
  geom_pointrange(size = 1, position = position_dodge(width = 0.5)) +
  # scale_colour_manual(values = CawthronColours::get_pal("caw_cat_1")) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Sex",
       y = "Survival probability",
       colour = "Stage") +
  theme_bw() 

# save figure
# ggsave("./figures/PVA/rmark_surival_non_captive.png",
#        width = 14,
#        height = 8,
#        dpi = 300,
#        units = "cm")

cleanup(ask = FALSE)

save(surv_data_analysis,file = "./Data/SurvDataAnalysis.RData")

# Remove all objects creates on this script
rm(mod_test,parent_sight,gg,sightings,sex_df,sight_record,sight2,sight3,bird_ch,bird_ch_tst,
   bird.fem.freq,bird.fem.hist,bird.freq,bird.group,bird.hist,bird.mal.freq,bird.mal.hist,capture_morphometrics,
   coeffs_full_set,egg_data,estimate_avg,estimate_se,fledg_birds,hade.ddl,hade.ddl2,hade.process,m,
   mod_av,mod_av_coeffs,mod_params,mod_res,sight5,sightings,years_cov,cap_b,last_year,mask,min_date,mod_overdisp,
   phi_coeffs)

# END --------------------------------------------------------------------------
