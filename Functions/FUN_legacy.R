# Function to initialize a population data frame
init_population <- function(pars, seed=19) {
  require(dplyr)
  
  # Check that required parameters are present
  needed_pars <- c("StartN_df", "sex_ratio", "no_age_classes","Fp","age_structure")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse=", ")))
  }
  
  # Define subpopulation labels
  pop <- data.frame()
  existing_ids<-character()
  
  pars$StartN_df<-pars$StartN_df%>%
    filter(count>0)
  # Loop through each subpopulation
  for (i in seq_along(pars$StartN_df[,1])) {
    # set.seed(i + seed)  # Ensure reproducibility with seed offset
    
    # Generate new individuals with IDs, random ages and sexes
    new_indivs <- data.frame(
      id = generate_unique_ids(n = pars$StartN_df$count[i],existing = existing_ids),
      age = rep(pars$StartN_df$age_class[i],pars$StartN_df$count[i]),
      sex = rep(pars$StartN_df$sex[i],pars$StartN_df$count[i])
    )
    
    # Append new individuals to the population
    pop <- rbind(pop, new_indivs)
    
    existing_ids<-unique(pop$id)
  }
  
  # Add time (t), alive, pair, and parental info
  pop <- pop %>%
    dplyr::mutate(t = 0, alive = TRUE, 
                  pair = NA_character_, mother_id = NA_character_, father_id = NA_character_,
                  origin="Wild",nz_heritage=1,sex=toupper(substr(sex,1,1)),
                  
                  Fi=pars$Fp,year_born=age-t,subpop="A") #tsr = time since release (years)
  
  
  return(pop)
}

# Function to simulate mortality for individuals at a given time step
mortlity_agng <- function(pop, currentT, pars,seed=19) {
  
  # Check that required mortality parameter is provided
  needed_pars <- c("phi_df","max_age","improved_foraging",
                   "start_cycle",
                   "dip_leth_eq",
                   "carr_capac_df",
                   "acc_period_df","release_year_cont","supp_feeding_current","full_pop")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse=", ")))
  }
  
  # Filter population for individuals at the current time step
  
  current <- pop %>%
    dplyr::filter(t == currentT,alive)
  
  
  current<-current%>%
    
    dplyr::left_join(pars$phi_df)%>%
    dplyr::left_join(pars$acc_period_df)%>%
    dplyr::left_join(pars$release_year_cont)%>%
    dplyr::left_join(pars$supp_feeding_current)%>%
    
    dplyr::mutate(yr_duration=case_when(
      tsr>0~1,
      TRUE~yr_duration
    ))%>%
    # dplyr::group_by(id)%>%
    dplyr::mutate(prop_year_acc=pmin(pmax(acc_period-tsr,0),yr_duration)/yr_duration)
  
  N_df<-current%>%
    dplyr::group_by(subpop)%>%
    dplyr::summarise(N=sum(alive))%>%
    dplyr::left_join(pars$carr_capac_df)%>%
    dplyr::mutate(AboveC=N>C,surv_trunc=1/(N/C))
  
  # ib_df<-calculate_inbreeding(pop = pars$full_pop,pars=pars)
  
  # Assign constant survival probability
  
  # environmental stochasticity (random increases throughout all parameters)
  # cyclcical variation
  # set.seed(seed+1)
  q_sample<-runif(n = 1,min = 0,max = 1)
  
  surv_mean<-mean<-current$Beta0
  cyc_var <- current$Beta0_sd_cyc*sin(pi*(currentT+pars$start_cycle)/2.5)
  sto_var<- sapply(seq_along(current$id),function(i){
    qnorm(p = q_sample,mean = 0,sd=current$Beta0_sd_sto[i])})
  
  surv_int<- surv_mean + sto_var + cyc_var
  
  
  surv <-
    ### Calculate survival based on habitat conditions
    inv.logit(surv_int + 
                current$Beta_sf * current$sf +
                current$Beta_if * pars$improved_foraging[currentT] +
                current$Beta_b  * current$sf * pars$improved_foraging[currentT]
              + 0)
  
  
  # add inbreeding depressio on survival if it is age class 1
  surv <- exp(log(surv) - current$Fi*(pars$dip_leth_eq/2)*(current$age==1) )
  
  #### Adjust by how much of the year individual spent in the wild
  surv<-surv ^ current$yr_duration
  
  ### adjust probability of surviving using odds-ratio, given the proportion of the period individual was under acclimation    
  surv<-adjust_probability(surv,current$OR_release,prop=current$prop_year_acc)
  
  current$surv<-surv
  
  current<-current%>%
    dplyr::left_join(N_df)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(surv=case_when(AboveC~min(surv_trunc,surv),
                                 TRUE~surv))
  
  if(any(is.na(surv)) | any(surv < 0) | any(surv > 1)){stop("Invalid estimated survival probability")}
  
  # set.seed(seed)
  # Simulate survival for each individual using Bernoulli trial
  alive_vec <- sapply(seq_along(current$id), FUN = function(i) {
    rbinom(n = 1, size = 1, prob = current$surv[i])
  }) == 1 & 
    current$alive & 
    current$age < pars$max_age
  
  # Update alive status based on survival outcome
  current$alive <- alive_vec
  
  current<-current%>%
    dplyr::mutate(
      ## time since release
      tsr=case_when(tsr==0~tsr+yr_duration, # if 0, add the year duration from release to "new year"
                    TRUE~tsr+1), # if bigger than 0, add a full year
      
      pair=case_when(
        !alive ~ NA,
        TRUE ~ pair
      ),
      
      age=case_when(
        alive ~ age +1,
        TRUE ~ age
      ))%>%
    tidy_pop_df()
  
  resu<-list(alive=current,stoch_q=q_sample)  
  return(resu)
}

mortality <- function(pop, currentT, pars,seed=19) {
  
  # Check that required mortality parameter is provided
  needed_pars <- c("phi_df","max_age",
                   # "improved_foraging",
                   # "start_cycle",
                   # "dip_leth_eq",
                   # "carr_capac_df",
                   # "acc_period_df",
                   # "release_year_cont",
                   # "supp_feeding_current",
                   "full_pop")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse=", ")))
  }
  
  # Filter population for individuals at the current time step
  
  current <- pop %>%
    dplyr::filter(t == currentT,alive)
  
  # current<-current%>%
  #   dplyr::filter(!is.na(age_release))
  current<-current%>%
    
    # dplyr::left_join(pars$phi_df)%>%
    # dplyr::left_join(pars$acc_period_df)%>%
    # dplyr::left_join(pars$release_year_cont)%>%
    # dplyr::left_join(pars$supp_feeding_current)%>%
    # dplyr::mutate(yr_duration=case_when(
    # tsr>0~1,
    # TRUE~yr_duration
    # ))%>%
    # dplyr::group_by(id)%>%
    # dplyr::mutate(prop_year_acc=pmin(pmax(acc_period-tsr,0),yr_duration)/yr_duration)currentPop0
    mutate(field=pars$field)%>%
    ungroup()
  
  N_df<-current%>%
    dplyr::group_by(subpop)%>%
    dplyr::summarise(N=sum(alive))%>%
    dplyr::left_join(pars$carr_capac_df)%>%
    dplyr::mutate(AboveC=N>C,surv_trunc=1/(N/C))
  
  # ib_df<-calculate_inbreeding(pop = pars$full_pop,pars=pars)
  
  # Assign constant survival probability
  
  # environmental stochasticity (random increases throughout all parameters)
  # cyclcical variation
  # set.seed(seed+1)
  q_sample<-runif(n = 1,min = 0,max = 1)
  
  surv_mean<-mean<-pars$phi_df$Intercept
  
  # cyc_var <- 0 * sin(pi*(currentT+pars$start_cycle)/2.5)
  cyc_var <- 0
  
  sto_var<- sapply(seq_along(current$id),function(i){
    qnorm(p = q_sample,mean = 0,sd=pars$phi_df$RE_time)})
  
  sto_var <- 0
  
  surv_int<- surv_mean + sto_var + cyc_var
  
  
  surv <-
    ### Calculate survival based on habitat conditions
    inv.logit(
      surv_int +
        (pars$phi_df$Beta_imm * (current$age == 2)) + 
        (pars$phi_df$Beta_adu * (current$age >= 3)) + 
        (pars$phi_df$Beta_newsq_juv * (current$age == 1) * (current$field == "SQ")) + 
        (pars$phi_df$Beta_newsq_imm * (current$age == 2) * (current$field == "SQ")) +
        (pars$phi_df$Beta_newsq_adu * (current$age >= 3) * (current$field == "SQ")) +
        (pars$phi_df$Beta_newsq_female * (current$sex == "F") * (current$field == "SQ")) + 
        (pars$phi_df$Beta_newsq_male * (current$sex == "M") * (current$field == "SQ")) + 
        (pars$phi_df$Beta_admix_juv * (current$age == 1) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) + 
        (pars$phi_df$Beta_admix_imm * (current$age == 2) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) + 
        (pars$phi_df$Beta_admix_adu * (current$age >= 3) * (current$nz_heritage < .95) * (current$nz_heritage > .05)) + 
        (pars$phi_df$Beta_ai_juv * (current$age == 1) * (current$origin == "AI")) + 
        (pars$phi_df$Beta_ai_imm * (current$age == 2) * (current$origin == "AI")) + 
        (pars$phi_df$Beta_ai_adu * (current$age >= 3) * (current$origin == "AI")) +
        (pars$phi_df$Beta_field2_juv * (current$age == 1) * (current$field == "Field2")) + 
        (pars$phi_df$Beta_field2_imm * (current$age == 2) * (current$field == "Field2")) + 
        (pars$phi_df$Beta_field2_adu * (current$age >= 3) * (current$field == "Field2"))
    )
  
  
  # add inbreeding depressio on survival if it is age class 1
  # surv <- exp(log(surv) - current$Fi*(pars$dip_leth_eq/2)*(current$age==1) )
  
  #### Adjust by how much of the year individual spent in the wild
  # surv<-surv ^ current$yr_duration
  
  ### adjust probability of surviving using odds-ratio, given the proportion of the period individual was under acclimation    
  # surv<-adjust_probability(surv,current$OR_release,prop=current$prop_year_acc)
  
  current$surv<-surv
  
  current<-current%>%
    dplyr::left_join(N_df)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(surv=case_when(AboveC~min(surv_trunc,surv),
                                 TRUE~surv))
  
  if(any(is.na(surv)) | any(surv < 0) | any(surv > 1)){stop("Invalid estimated survival probability")}
  
  # set.seed(seed)
  # Simulate survival for each individual using Bernoulli trial
  alive_vec <- sapply(seq_along(current$id), FUN = function(i) {
    rbinom(n = 1, size = 1, prob = current$surv[i])
  }) == 1 & 
    current$alive & 
    current$age < pars$max_age
  
  # Update alive status based on survival outcome
  current$alive <- alive_vec
  
  current<-current%>%
    dplyr::mutate(
      ## time since release
      # tsr=case_when(tsr==0~tsr+yr_duration, # if 0, add the year duration from release to "new year"
      # TRUE~tsr+1), # if bigger than 0, add a full year
      
      pair=case_when(
        !alive ~ NA,
        TRUE ~ pair
      ),
      
      # age=case_when(
      # alive ~ age +1,
      # TRUE ~ age
      # )
    )%>%
    tidy_pop_df()
  
  resu<-list(alive=current,stoch_q=q_sample)  
  return(resu)
}

aging <- function(pop, currentT, pars,seed=19) {
  
  # Check that required mortality parameter is provided
  needed_pars <- c("phi_df","max_age",
                   "field",
                   # "improved_foraging",
                   # "start_cycle",
                   # "dip_leth_eq",
                   # "carr_capac_df",
                   # "acc_period_df","release_year_cont",
                   # "supp_feeding_current",
                   "full_pop")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse=", ")))
  }
  
  # Filter population for individuals at the current time step
  
  current <- pop %>%
    dplyr::filter(t == currentT)
  
  
  
  # Assign constant survival probability
  
  current<-current%>%
    dplyr::mutate(
      age=case_when(
        alive ~ age +1,
        TRUE ~ age
      )
    )%>%
    tidy_pop_df()
  
  resu<-current
  return(resu)
}

#--------------------------#
# Function: pairing        #
#--------------------------#
# Pairs breeding-age adults within each subpopulation.
# Args:
#   pop: Population dataframe.
#   currentT: Current time step.
#   pars: A list with "breeding_age".
#   seed: Random seed for reproducibility.
# Returns:
#   Updated adult population dataframe with pair IDs assigned.
pairing <- function(pop, currentT, pars, seed=19) {
  require(dplyr)
  
  needed_pars <- c("breeding_age")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse=", ")))
  }
  
  # Add priority to manage replacement of existing pairs
  pop <- pop %>%
    dplyr::mutate(priority = 0)
  
  adults <- pop %>%
    filter(alive, t == currentT, age >= pars$breeding_age)%>%
    dplyr::ungroup()
  
  pairs <- data.frame()
  subpops <- unique(adults$subpop)
  
  # Pairing process for each subpopulation
  for (sp in subpops) {
    males <- adults %>% filter(sex == "M", subpop == sp, is.na(pair))
    females <- adults %>% filter(sex == "F", subpop == sp, is.na(pair))
    
    # set.seed(seed + currentT)
    males <- males[sample(nrow(males)), ]
    females <- females[sample(nrow(females)), ]
    
    n_pairs <- min(nrow(males), nrow(females))
    if(n_pairs>0){
      pairs <- plyr::rbind.fill(pairs,
                                (males[1:n_pairs, ]) %>%
                                  dplyr::mutate(pair = females$id[1:n_pairs], priority = 1),
                                (females[1:n_pairs, ]) %>%
                                  dplyr::mutate(pair = males$id[1:n_pairs], priority = 1)
      )
    }
  }
  
  # Prioritize newly formed pairs
  adults_new <- plyr::rbind.fill(pop, pairs) %>%
    tidy_pop_df()
  
  return(adults_new)
}

#-----------------------------------#
# Function: unpair_if_dead          #
#-----------------------------------#
# Removes pair references to individuals who died in currentT.
# Args:
#   pop: Population dataframe (includes both dead and alive).
#   currentT: Time step.
# Returns:
#   Updated dataframe with invalid pairings removed.
unpair_if_dead <- function(pop, currentT) {
  require(dplyr)
  
  alive_ids <- pop %>%
    dplyr::filter(t == currentT, alive) %>%
    dplyr::pull(id)
  
  pop <- pop %>%
    dplyr::mutate(pair = case_when(
      (!pair %in% alive_ids) ~ NA,
      TRUE ~ pair
    ))
  
  return(pop)
}


#----------------------------#
# Function: recruitment      #
#----------------------------#
# Adds new offspring to the population at a given time step.
# Only females that are alive, paired, of breeding age, and present at time `currentT` can reproduce.
# Offspring number per female is drawn from a truncated Poisson distribution if nesting is successful.
#
# Args:
#   pop: Current population dataframe.
#   currentT: Current time step.
#   pars: A list containing:
#     - breeding_age: Minimum age required for reproduction.
#     - nesting_success: Probability that a female successfully nests.
#     - av_brood_size: Mean number of offspring per nesting attempt.
#     - sex_ratio: Probability of female among offspring.
#     - max_brood_size: Maximum allowed number of offspring per female.
#     - all_ids: Character vector of all previously used IDs.
#   seed: Random seed for reproducibility.
# Returns:
#   Updated population dataframe with new individuals appended.
recruitment <- function(pop, currentT, pars, seed = 19,envir_stoch=TRUE,startAge=0) {
  
  # Ensure required parameters are present
  needed_pars <- c("breeding_age",
                   "fb_df",
                   "cs_df",
                   "hp_df",
                   "fp_df",
                   "sex_ratio",
                   "max_brood_size",
                   "first_clutch_harvest_prop",
                   "prop_managed_clutches",
                   "carr_capac_df",
                   "all_ids")
  
  par_names <- names(pars)
  
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  
  temp<-pop%>%
    dplyr::filter(
      #                 sex == "F",
      !is.na(pair),
      alive,
      t == currentT,
      age >= pars$breeding_age)%>%
    ungroup()%>%
    mutate(able_to_breed=T)
  
  male_breed_capacity <- temp%>%filter(sex=="M")%>%
    pull_named(able_to_breed,id)
  
  female_breed_capacity <- temp%>%filter(sex=="F")%>%
    pull_named(able_to_breed,id)
  
  female_pair <- temp%>%filter(sex=="F")%>%
    pull_named(pair,id)
  
  f_breed_df<-data.frame(f_able_to_breed=female_breed_capacity,
                         id=names(female_breed_capacity),
                         pair=female_pair,
                         m_able_to_breed=male_breed_capacity[female_pair])
  
  # Filter eligible reproducing females
  reproducing <- pop %>%
    dplyr::filter(sex == "F",
                  !is.na(pair),
                  alive,
                  t == currentT,
                  age >= pars$breeding_age)%>%
    dplyr::left_join(f_breed_df)%>%
    dplyr::mutate(able_to_breed = f_able_to_breed & m_able_to_breed)%>%
    dplyr::filter(able_to_breed)
  
  ### This line adds carrying capacity
  fem_able_to_breed<-sample(reproducing$id,size = min(length(reproducing$id),pars$carr_capac_df$C),replace = F)
  
  reproducing<-reproducing%>%
    dplyr::mutate(able_to_breed=id%in%fem_able_to_breed)
  
  reproducing$breeding_prob<- pars$fb_df$Intercept + 
    (pars$fb_df$Beta_newsq  * (pars$field=="SQ")) + 
    (pars$fb_df$Beta_field2 * (pars$field=="Field2")) + 
    (pars$fb_df$Beta_ai     * (reproducing$origin=="AI")) +
    0
  
  reproducing$breed <- sapply(seq_along(reproducing$id),function(i){
    rbinom(n = 1,size = 1,prob = reproducing$breeding_prob)
  })==1
  
  reproducing<-reproducing%>%
    dplyr::filter(breed)
  
  
  # # reproducing$male_breed<-male_breed_capacity[reproducing$pair]
  # reproducing$female_breed<-reproducing$able_to_breed
  # reproducing$able_to_breed<-reproducing$female_breed & reproducing$male_breed
  
  # set.seed(seed+1)
  if(envir_stoch){q_sample<-runif(n = 1,min = 0,max = 1)}else{q_sample<-.5}
  
  if(nrow(reproducing)>0){
    
    # Assign per-female nesting success and expected brood size
    females_breeding<-nrow(reproducing)
    
    cs_mean<-pars$cs_df$Intercept
    
    pars$cs_df$RE_time<-0
    cs_st_var<-sapply(seq_along(reproducing$id),function(i){
      qnorm(p = q_sample,mean = 0,sd=pars$cs_df$RE_time)
    })
    
    cs_int <- cs_mean + cs_st_var
    # bs_int <- pmax(bs_int,0)
    # cs_int <- exp(cs_int)
    
    brood_size <- exp(cs_int+
                        (pars$cs_df$Beta_newsq * (pars$field == "SQ"))+
                        (pars$cs_df$Beta_field2 * (pars$field == "Field2"))+
                        (pars$cs_df$Beta_admix * (reproducing$nz_heritage<0.95 & reproducing$nz_heritage>0.05))+
                        (pars$cs_df$Beta_ai * (reproducing$origin == "AI")) +
                        0)%>%
      pmax(0)%>%
      # pmin(pars$max_brood_size)%>%
      identity()
    
    reproducing$cs<-brood_size
    no_eggs<-pmin(sapply(reproducing$cs,function(x){rpois(n = 1,lambda = x)}),pars$max_brood_size)
    
    
    # Determine subpopulation for each offspring
    motherIDs <- lapply(seq_along(reproducing$id), function(i) {
      tempId<-reproducing$id[i]
      rep(tempId, no_eggs[i])
    }) %>%
      unlist() 
    
    # Determine subpopulation for each offspring
    fatherIDs <- lapply(seq_along(reproducing$id), function(i) {
      tempId<-reproducing$pair[i]
      rep(tempId, no_eggs[i])
    }) %>%
      unlist()  
    
    
    nz_heritage_vec<-(((temp%>%pull_named(nz_heritage,id))[fatherIDs])+
                        ((temp%>%pull_named(nz_heritage,id))[motherIDs]))/2
    
    # Assign sexes to offspring
    newSex <- c("F", "M")[rbinom(n = sum(no_eggs), size = 1, prob = pars$sex_ratio) + 1]
    
    if(sum(no_eggs)>0){
      # Generate unique IDs for all new offspring
      eggs <- data.frame(
        id = generate_unique_ids(n = sum(no_eggs), existing = pars$all_ids),
        age = 0,
        sex = newSex,
        subpop = "A",
        t = currentT,
        # t = currentT+1,
        alive = TRUE,
        pair = NA,
        mother_id=motherIDs,
        origin="Wild",
        father_id=fatherIDs,
        nz_heritage = nz_heritage_vec,
        clutch_no = 1,
        year_born=currentT
      )
      
      eggs$mother_origin<-(temp%>%pull_named(origin,id))[eggs$mother_id]
      eggs$mother_age<-(temp%>%pull_named(age,id))[eggs$mother_id]
      
      harvested_clutches_mom_id<-sample(unique(eggs$mother_id),
                                        size = round(length(unique(eggs$mother_id))*pars$first_clutch_harvest_prop),
                                        replace=F)
      
      harvested_eggs<-eggs%>%
        dplyr::filter(mother_id%in%harvested_clutches_mom_id)
      harvested_eggs$origin<-rep("AI",nrow(harvested_eggs))
      harvested_eggs$id<-rep(NA_character_,nrow(harvested_eggs))
      
      eggs<-eggs%>%
        dplyr::mutate(clutch_no=case_when(
          mother_id%in%harvested_clutches_mom_id~clutch_no+1,
          TRUE~clutch_no
        ))
      
      clutches_ids<-unique(eggs$mother_id)
      managed_clutches<-rbinom(n = length(clutches_ids),size=1,prob=pars$prop_managed_clutches)
      
      managed_clutches_ids<-clutches_ids[managed_clutches==1]
      
      eggs$managed<-eggs$mother_id%in%managed_clutches_ids
      
      hp_mean<-pars$hp_df$Intercept
      
      hp_st_var<-sapply(seq_along(eggs$id),function(i){
        qnorm(p = q_sample,mean = 0,sd=pars$hp_df$RE_time)
      })
      
      hp_int <- hp_mean + hp_st_var
      
      hatch_prob<-inv.logit(hp_int+
                              (pars$hp_df$Beta_f_age * eggs$mother_age) +
                              (pars$hp_df$Beta_newsq * (pars$field=="SQ"))+
                              (pars$hp_df$Beta_field2 * (pars$field=="Field2"))+
                              (pars$hp_df$Beta_admix * (eggs$nz_heritage<0.95 & eggs$nz_heritage>0.05))+
                              (pars$hp_df$Beta_ai * (eggs$mother_origin=="AI"))+
                              (pars$hp_df$Beta_managed * (eggs$managed))+
                              (pars$hp_df$Beta_clutch_no * eggs$clutch_no) + 
                              0)
      
      eggs$hp<-hatch_prob
      
      # cat("Sampling wether eggs survive\n")
      eggs$alive<-sapply(seq_along(eggs$id),function(i){rbinom(n = 1,size = 1,prob = eggs$hp[i])})==1
      
      eggs_laid<-nrow(eggs)
      eggs_hatched<-sum(eggs$alive)
      
      if(eggs_hatched>0){
        
        hatchlings<-eggs%>%
          dplyr::filter(alive)
        
        fp_mean<-pars$fp_df$Intercept
        
        fp_st_var<-sapply(seq_along(hatchlings$id),function(i){
          qnorm(p = q_sample,mean = 0,sd=pars$fp_df$RE_time)
        })
        
        fp_int <- fp_mean + fp_st_var
        
        fledge_prob<-inv.logit(fp_int+
                                 (pars$fp_df$Beta_f_age * hatchlings$mother_age) +
                                 (pars$fp_df$Beta_newsq * (pars$field=="SQ"))+
                                 (pars$fp_df$Beta_field2 * (pars$field=="Field2"))+
                                 (pars$fp_df$Beta_admix * (hatchlings$nz_heritage<0.95 & hatchlings$nz_heritage>0.05))+
                                 (pars$fp_df$Beta_ai * (hatchlings$mother_origin=="AI"))+
                                 (pars$fp_df$Beta_managed * (hatchlings$managed))+
                                 (pars$fp_df$Beta_clutch_no * hatchlings$clutch_no) + 
                                 0)
        
        hatchlings$fp<-fledge_prob
        
        # cat("Sampling wether hatchlings survive\n")
        
        hatchlings$alive<-sapply(seq_along(hatchlings$id),function(i){rbinom(n = 1,size = 1,prob = hatchlings$fp[i])})==1
        
        fledged<-sum(hatchlings$alive)
        
        if(fledged>0){
          born<-hatchlings%>%
            filter(alive)%>%
            mutate(Fi=NA)%>%
            tidy_pop_df()}else{born<-data.frame()}
        
      }else{
        
        fledged<-0
        born<-data.frame()}
      
      
    }else{
      
      
      
      born<-data.frame()
      harvested_eggs<-data.frame()
      
      eggs_laid<-0
      eggs_hatched<-0
      fledged<-0
    }
    
    # if no reproduction
  }else{
    
    born<-data.frame()
    harvested_eggs<-data.frame()
    females_breeding<-0
    eggs_laid<-0
    eggs_hatched<-0
    fledged<-0
    
  }
  
  # Append new offspring to existing population
  newpop <- plyr::rbind.fill(pop, born)
  
  return(list(born=born,
              females_breeding=females_breeding,
              harvested_eggs=harvested_eggs,
              eggs_laid=eggs_laid,eggs_hatched=eggs_hatched,fledged=fledged))
}


ai_recruitment <- function(pars,seed=19,envir_stoch=TRUE,startAge=0){
  
  # Ensure required parameters are present
  needed_pars <- c(
    "ai_hatch_prob",
    "ai_fledge_prob",
    "harvested_eggs",
    "all_ids")
  
  par_names <- names(pars)
  
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  eggs<-pars$harvested_eggs
  if(nrow(eggs)>0){
    
    eggs$id<-generate_unique_ids(n = nrow(eggs),existing = pars$all_ids)
    
    eggs$alive<-sapply(seq_along(eggs$id),function(i){
      rbinom(n = 1,size = 1,prob = pars$ai_hatch_prob)})==1
    
    hatchlings<-eggs%>%
      dplyr::filter(alive)
    
    if(nrow(hatchlings)>0){
      
      hatchlings$alive<-sapply(seq_along(hatchlings$id),function(i){
        rbinom(n = 1,size = 1,prob = pars$ai_fledge_prob)})==1
      
      born<-hatchlings%>%
        mutate(Fi=NA)%>%
        tidy_pop_df()%>%
        dplyr::filter(alive)
    }else{
      born<-data.frame()} 
  }else{
    born<-data.frame()
  }
  
  return(born)
}

#----------------------------#
# Function: dispersal        #
#----------------------------#
# Simulates the dispersal of individuals across subpopulations.
# Only alive individuals at time `currentT` and in allowed dispersal ages can disperse.
# Dispersal is based on a matrix of probabilities for moving between subpopulations.
#
# Args:
#   pop: Population dataframe.
#   currentT: Current time step.
#   pars: A list containing:
#     - dispersalMat: A square matrix of movement probabilities between subpops.
#       Rows = origin subpop, Columns = destination subpop.
#     - dispersalAges: Vector of ages at which dispersal is allowed.
#   seed: Random seed for reproducibility.
# Returns:
#   Updated population dataframe with new subpopulation assignments for dispersers.
dispersal <- function(pop, currentT, pars, seed = 19) {
  
  # Check required parameters
  needed_pars <- c("dispersalMat", "dispersalAges")
  par_names <- names(pars)
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  # Set default priority to 0 for everyone
  pop <- pop %>%
    dplyr::mutate(priority = 0)
  
  # Identify all subpop labels
  subpops <- rownames(pars$dispersalMat)
  
  # Select eligible dispersers: alive, at current time, and within dispersal age class
  dispersers <- pop %>%
    dplyr::filter(alive, t == currentT, age %in% pars$dispersalAges,is.na(pair)) %>%
    dplyr::mutate(priority = 1)
  
  # set.seed(seed)
  
  # If there are individuals eligible to disperse
  if(nrow(dispersers) > 0) {
    
    # Sample destination subpopulation based on dispersal matrix
    dispersal_destinations <- sapply(seq_along(dispersers$id), function(i) {
      subpops[which(
        rmultinom(n = 1, size = 1, 
                  prob = pars$dispersalMat[dispersers$subpop[i], ])[, 1] == 1
      )]
    })
    
    # Update subpopulation assignment for dispersers
    dispersers$subpop <- dispersal_destinations
    
    # Combine updated dispersers and original population
    newpop <- plyr::rbind.fill(pop, dispersers) %>%
      tidy_pop_df()
    
  }else{
    newpop<-pop
  }
  
  # Return updated population (with dispersal applied if needed)
  return(newpop)
}

# Function to simulate zero-truncated Poisson
rtpois <- function(n, lambda,min=0,max=Inf) {
  x <- rpois(n, lambda)
  while(any(x <= min | x >= max)) {
    x[x <= min | x >= max] <- rpois(sum(x <= min | x >= max), lambda)
  }
  return(x)
}

#----------------------------#
# Function: generate_unique_ids
#----------------------------#
# Generates multiple unique 8-character strings not found in the given vector.
# Args:
#   n: Number of unique IDs to generate.
#   existing: A character vector of already-used strings.
#   charset: Characters to sample from (default: letters and digits).
#   max_tries: Maximum attempts to find unique IDs (default: 10 * n).
# Returns:
#   A character vector of `n` unique strings.
generate_unique_ids <- function(n, existing = character(0),
                                strlength = 8,
                                charset = c(0:9, LETTERS), 
                                max_tries = 10 * n) {
  unique_ids <- character(0)
  tries <- 0
  
  while (length(unique_ids) < n && tries < max_tries) {
    needed <- n - length(unique_ids)
    candidates <- replicate(needed, paste0(sample(charset, strlength, replace = TRUE), collapse = ""))
    # Filter out existing and duplicate values
    new_ids <- setdiff(candidates, c(existing, unique_ids))
    unique_ids <- unique(c(unique_ids, new_ids))
    tries <- tries + 1
  }
  
  if (length(unique_ids) < n) {
    stop("Failed to generate the required number of unique IDs after ", max_tries, " attempts.")
  }
  
  return(unique_ids)
}



tidy_pop_df <- function(df){
  
  suppressWarnings({if(is.null(df$priority)){df$priority<-0}})
  
  resu<-df%>%
    dplyr::arrange(desc(priority),desc(t)) %>%           # Prioritize updates
    dplyr::filter(!duplicated(data.frame(id,t))) %>%           # Keep only one entry per individual
    dplyr::select(id, subpop, sex, age, t, alive, pair,
                  mother_id,father_id,
                  Fi,nz_heritage,
                  year_born,
                  origin)%>%
    dplyr::arrange(desc(alive),subpop,age,mother_id,origin)
  
  return(resu)
}



adjust_probability <- function(prob, odds_ratio,prop=1,agg_rule="survival") {
  
  if (any(prob < 0) || any(prob > 1)) {
    stop("Probability must be between 0 and 1 (inclusive).")
  }
  
  if (any(odds_ratio <= 0)) {
    stop("Odds ratio must be greater than 0.")
  }
  
  prop[is.na(prop)] <- 1
  
  odds <- prob / (1 - prob)
  
  adjusted_odds <- odds * odds_ratio
  adjusted_prob <- adjusted_odds / (1 + adjusted_odds)
  
  adjusted_prob[prob%in%c(0,1)]<-prob[prob%in%c(0,1)]
  
  if(agg_rule=="survival"){
    resu<- (prob^(1-prop)) * (adjusted_prob^prop)
  }
  
  # if(agg_rule=="success"){
  #   fail <- ((1-prob)^(1-prop)) * ((1 - adjusted_prob)^prop)
  #   resu <- 1-fail
  # }
  # 
  
  return(resu)
}


releases<-function(pars,currentT){
  
  # Ensure required parameters are present
  needed_pars <- c("admix_release_years",
                   "admix_prop_released",
                   "admix_age_released",
                   "admix_no_released",
                   "all_ids")
  
  par_names <- names(pars)
  
  missing_pars <- needed_pars[which(!needed_pars %in% par_names)]
  
  if(length(missing_pars) > 0){
    stop(paste0("Error: Missing parameters are ", paste0(missing_pars, collapse = ", ")))
  }
  
  
  releaseYears<-which(pars$admix_release_years)
  
  if(currentT%in%releaseYears & pars$admix_no_released>0 ){
    
    # Generate unique IDs for all new offspring
    released <- data.frame(
      id = generate_unique_ids(n = pars$admix_no_released, existing = pars$all_ids),
      age = pars$admix_age_released,
      sex = c("F", "M")[rbinom(n = pars$admix_no_released, size = 1, prob = 0.5) + 1],
      subpop = "A",
      t = currentT,
      # t = currentT+1,
      alive = TRUE,
      pair = NA,
      mother_id=NA,
      origin="Captivity",
      nz_heritage=sample(c(1,0.5),replace = T,size=pars$admix_no_released,
                         prob = c(1-pars$admix_prop_released,pars$admix_prop_released)),
      father_id=NA,
      Fi=0
    )      
  }else{
    
    released<-data.frame()
  }
  
  return(released)
  
}


replace_na_characters <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col[is.na(col)] <- ""
    }
    return(col)
  })
  return(df)
}

split_evenly_by_col <- function(df, n_groups, target_col) {
  library(dplyr)
  
  df <- df %>%
    arrange(desc(!!sym(target_col)))  # Sort by descending target column value
  
  # Initialize groups
  group_sums <- rep(0, n_groups)
  df$group <- NA  # Column to store group assignment
  
  for (i in 1:nrow(df)) {
    # Find the group with the smallest current sum
    min_group <- which.min(group_sums)
    
    # Assign the current row to that group
    df$group[i] <- min_group
    
    # Update the sum of that group
    group_sums[min_group] <- group_sums[min_group] + df[[target_col]][i]
  }
  
  # Split into a list of sub-dataframes
  split_list <- split(df, df$group)
  
  return(split_list)
}



output_clean_up<-function(rev=FALSE){
  ### Get files in folder
  files<-list.files(path = "./Results",pattern = ".RData",full.names = T)
  
  cat(paste0("\nFound ",length(files)," files\n"))
  
  if(rev){files<-rev(files)}
  if(length(files)>0){
    #  Upload to Gdrive
    for(f in files){
      
      suppressMessages({
        gdrive_upload<-sapply(f,drive_upload, path = as_id("1FSqFfBrvyedxTOUuIFmI44u9tXxEMMU0"))
        uploaded_file<-attributes(gdrive_upload)$dimnames[[2]]
        # Move files to the backup directory
        file.rename(uploaded_file, file.path(backup_dir, basename(f)))
      })
    }
  }
}



output_clean_up_zip <- function(rev = FALSE) {
  # Get files in folder
  files <- list.files(path = "./Results", pattern = "\\.RData$", full.names = TRUE)
  
  cat(paste0("\nFound ", length(files), " files\n"))
  
  if (length(files) == 0) return(invisible(NULL))
  if (rev) files <- rev(files)
  
  # Define zip filename
  zip_filename <- file.path(tempdir(), paste0("results_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"))
  
  # Create zip file containing all files
  zip(zipfile = zip_filename, files = files, flags = "-j")  # -j removes directory structure
  
  suppressMessages({
    # Upload the zip file to Google Drive
    gdrive_upload <- drive_upload(zip_filename, path = as_id("1FSqFfBrvyedxTOUuIFmI44u9tXxEMMU0"))
    
    # Move the zip file to backup directory
    file.rename(zip_filename, file.path(backup_dir, basename(zip_filename)))
    
    # Move original .RData files to backup directory
    file.rename(files, file.path(backup_dir, basename(files)))
  })
}












