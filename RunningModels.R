runLabel<-"bigRunV2"
get_runs_from_gsheet<-FALSE
replace_runs_gsheet<-FALSE
prior_rng_seed<-1e3
runAnalysis<-F
wd<-"~/TaraItiIBM"

if(!dir.exists(wd)){
  wd<-"D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM"
  
}
if(!dir.exists(wd)){
  wd<-"C:/Users/Caio.Kenup/TaraItiIBM"
  
}
setwd(wd)


source("packageLoader.R")
source("functionLoader.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")


sheet_url <- "https://docs.google.com/spreadsheets/d/1dCIkkofz0h2s9MWOtZfNqAN4DKY59Z2IqNIWlM3isMY/edit?gid=903185379#gid=903185379"
gs4_auth(path = "./.tokens/fresh-replica-344321-0e0618a3b5de.json")

### Acessing Google Drive
token <- readRDS(".tokens/token.rds")
drive_auth(token = token)

### Get device name
source(".tokens/setDeviceName.R")


### Parameter setting
model_pars<-list(priors=list(),sim=list(),bio=list(),mgmt=list())


### Import expert elicited info
load("./Data/Tara_Iti_Aggregated.RData")
source("./Parameters/SimPars.R")

if(runAnalysis){
source("./DataAnalysis/DataAnalysisNesting.R")
source("./DataAnalysis/DataAnalysisEggs.R")
source("./DataAnalysis/DataAnalysisSurvival.R")
}else{
}


datafiles <- list.files(path = "./Data/", pattern = "\\.RData$", full.names = TRUE)
for (f in datafiles) load(f, envir = .GlobalEnv)

model_pars$priors<-c(model_pars$priors,
                     egg_data_analysis$priors,
                     nest_data_analysis$priors,
                     surv_data_analysis$priors)

model_pars$bio<-c(model_pars$bio,
                    egg_data_analysis$bio,
                    nest_data_analysis$bio,
                    surv_data_analysis$bio)

model_pars$mgmt<-c(model_pars$mgmt,
                  egg_data_analysis$mgmt,
                  nest_data_analysis$mgmt,
                  surv_data_analysis$mgmt)

source("./Parameters/BioPars.R")
source("./Parameters/MgmtPars.R")


init_pars<-list(StartN_df=model_pars$bio$StartN_df,
           sex_ratio=0.5,
           breeding_age=3,
           no_age_classes=model_pars$bio$inherent$max_age,
           max_age=model_pars$bio$inherent$max_age,
            age_structure=model_pars$bio$inherent$age_structure,
           Fp=model_pars$bio$gen$starting_inbreeding)


### Sample uncertain parameters, given the description of their distribution
prior_rng<-priorSampling(model_pars$priors,
                         seed = 26051991,
                         # size=1e3)%>%
                         size=model_pars$sim$n_iter)%>%
  dplyr::mutate(p = 1:n())

if(!model_pars$sim$use_genetics){
prior_rng$diploid_eq<-0  # manually turn off genetics
}

source("./Parameters/priorHandling.R")


### ADJUSTMENTS

### weighing population towards adults

# model_pars$bio$inherent$age_structure[1:2]<-model_pars$bio$inherent$age_structure[1:2]*.5
# model_pars$bio$inherent$age_structure<-model_pars$bio$inherent$age_structure/sum(model_pars$bio$inherent$age_structure)

all_iterations<-merge(prior_rng,model_pars$mgmt$strategies)%>%
  dplyr::arrange(p,alt)%>%
  dplyr::mutate(i=(1:n())+model_pars$sim$idx_add,Label=runLabel)

if(get_runs_from_gsheet | replace_runs_gsheet){
 
  runs0<- read_sheet(sheet_url,sheet="Runs")
  runs0$Label<-as.character(runs0$Label)
  runs<-runs0%>%
    replace_na_characters()%>%
    dplyr::filter(Label==runLabel)
  
  runs2<-left_join(
                   all_iterations%>%
                     dplyr::mutate(Iteration=i)%>%
                     dplyr::select(Label,alt,Iteration,p),runs)%>%
    dplyr::arrange(Label,Iteration)%>%
    dplyr::mutate(ID=paste(Label,zero_pad(Iteration,5),sep="_"))
  
  if(replace_runs_gsheet){
    write_sheet(runs2, sheet_url, sheet = "Runs")
  }
  
}

#### Get iterations to run on this script run
if(get_runs_from_gsheet){
  
  n_p_used<-1e3
  
  p_summ<-runs2%>%
    filter(Label==runLabel,Scheduled==TRUE)%>%
    dplyr::group_by(p)%>%
    dplyr::summarise(n=n())%>%
    dplyr::ungroup()%>%
    dplyr::arrange(n)
  
  iterations_to_run<-runs2%>%
    filter(Label==runLabel,Scheduled==TRUE,DeviceToRun==device_name)%>%
    dplyr::arrange(desc(p))%>%
    dplyr::group_by(p)%>%
    dplyr::arrange(desc(alt))%>%
    dplyr::ungroup()%>%
    dplyr::rename(i=Iteration)%>%
    dplyr::filter(p%in%p_summ$p[1:n_p_used])%>%
    dplyr::select(alt,p,i)%>%
    # left_join(alternatives%>%dplyr::select(alt,ExpectedTimeMins))%>%
    dplyr::ungroup()
  # pull(Iteration)
  
}else{
  
  n_p_used<-1e3   
  
  p_summ<-all_iterations%>%
    # filter(Label==runLabel,Scheduled==TRUE)%>%
    dplyr::group_by(p)%>%
    dplyr::summarise(n=n())%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(n))
  
  iterations_to_run<-all_iterations%>%
    dplyr::arrange(desc(p))%>%
    dplyr::group_by(p)%>%
    dplyr::arrange(desc(alt))%>%
    dplyr::filter(p%in%p_summ$p[1:n_p_used])%>%
    dplyr::select(alt,p,i)%>%
    # left_join(alternatives%>%dplyr::select(alt,ExpectedTimeMins))%>%
    dplyr::ungroup()
  # %>%
  # dplyr::pull(i)
}

if(nrow(iterations_to_run)<model_pars$sim$clusters_to_run & model_pars$sim$parallel_across_runs){
  stop("More clusters than iterations")
  # iterations_to_run<-iterations_to_run[nrow(length(iterations_to_run),model_pars$sim$clusters_to_run),]
}


# iterations_to_run<-iterations_to_run[1:14,]
cat(paste0("\n running a model for ",model_pars$sim$n_years,
           " time steps",
           # length(spatial_info_list[[1]]$small$valid_cells),
           ", for ",nrow(iterations_to_run)," runs,",
           ifelse(model_pars$sim$parallel_across_runs,paste0(" in parallel using ",model_pars$sim$clusters_to_run," nodes"),
                  "sequentially"),
           "\n"))

if(model_pars$sim$parallel_across_runs){
# if(FALSE){
    
  require(parallel)
  
  model_pars$sim$print_crit<-Inf # remove printing outcomes - wont show up anyway
  
  # Define batch size (number of runs per node before restarting)
  batch_size <- model_pars$sim$batching_clusters * model_pars$sim$clusters_to_run
  
  iterations_to_run<-iterations_to_run%>%
    dplyr::group_by(p)%>%
    dplyr::mutate(n=n(),dummy=1)%>%
    dplyr::ungroup()%>%
    dplyr::arrange(n,p,alt)
    # dplyr::mutate(batch=rep_len(1:batch_size,length.out=n()))%>%
    # dplyr::group_by(batch)%>%
    
  
  # Split iterations into batches of batch_size
  
  iterations_to_run_bkp<-iterations_to_run
  
  all_batches<-split_evenly_by_col(df = iterations_to_run,
                                   n_groups = ceiling(nrow(iterations_to_run_bkp)/model_pars$sim$batching_clusters),
                                   target_col = "dummy")
  
  batch_summ<-plyr::rbind.fill(all_batches)%>%
    # dplyr::group_by(group)%>%
    # dplyr::summarise(time=sum(ExpectedTimeMins))%>%
    # dplyr::ungroup()%>%
    arrange((group),p,alt)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(batch=NA)
  batchLoop<-1
  for(i in 1:nrow(batch_summ)){
    
    batch_summ$batch[i]<-batchLoop
    
    if(sum(batch_summ$batch==batchLoop,na.rm = T)==batch_size){batchLoop<-batchLoop+1}
    
  }
  
  iteration_batches<-split(batch_summ,batch_summ$batch)
  iteration_batches<-lapply(iteration_batches,function(x){split(x,ceiling(seq_along(x[,1])/(batch_size/model_pars$sim$clusters_to_run)))})

  for (batch_index in seq_along(iteration_batches)) {
  # for (batch_index in 1:1) {
      
    # batch <- iteration_batches[[batch_index]]
    # 
    
    sub_batches<-iteration_batches[[batch_index]]
    cluster_size<-min(model_pars$sim$clusters_to_run, length(sub_batches))
  
    iteration_chunks <- lapply(sub_batches,function(x){sample(x$i)})
    
    # Split current batch among nodes
    cat(paste0("\nStarting batch ", batch_index, " at ", lubridate::now(), "\n"))
    
    # Start a new cluster for each batch
    cl <- makeCluster(cluster_size)
    
    # Initialize cluster with required functions and data
    invisible(clusterEvalQ(cl, {
      
      source("packageLoader.R")
      source("functionLoader.R")
      
      devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")
      devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")
      
      token <- readRDS(".tokens/token.rds")
      drive_auth(token = token)
      
      con <- file(nullfile(), open = "w")
      sink(con, type = "message")    
    }))
    
    # Assign objects to cluster
    cat("\n Passing objects to nodes\n")
    clusterExport(cl, varlist = c("init_pars",
                                  "model_pars","prior_rng_seed","iterations_to_run",
                                  "all_iterations"))
    
    
    cat("\nStarting parallel execution for batch", batch_index, "\n")
    
    # Run in parallel
    parLapply(cl = cl, iteration_chunks, fun = function(chunk_items) {
      
      j <- 1
      for (cc in seq_along(chunk_items)) {
        
        i <- chunk_items[cc]
        
        idx<-i  # clone iteration index 
        p <- iterations_to_run%>%filter(i==idx)%>%pull(p) # get parameter sampling index
        
        all_iterations<-all_iterations%>%dplyr::arrange(i)
        suppressMessages({
          source("./Parameters/pars_postPriorSampling.R",local = T)
        })
        
        set.seed(prior_rng_seed+p)
        init_pop<-init_population(pars=init_pars)
        set.seed(prior_rng_seed+p)
        init_pop<-pairing(pop=init_pop,currentT = 0,pars=init_pars)
        
        start_conditions<-list(init_pop=init_pop)
        
        # model_pars$mgmt$supp_feeding_df<-model_pars$mgmt$supp_feeding_master[[as.character(all_iterations$SuppFeed[i])]]
        # model_pars$mgmt$release_schedule<-model_pars$mgmt$release_schedule_master%>%
        # dplyr::filter(r==all_iterations$ReleaseStrat[i])
        
        # model_pars$mgmt$acc_period_df$acc_period<-1
        set.seed(p+50)
        output<-run_model(start_conditions=start_conditions,model_pars=model_pars,idx = i)
        output$run_pars <- all_iterations[i,]
        output$run_label <- all_iterations$Label[i]
        
        # Save output
        filename_output <- paste0("./Results/", output$run_label, "_Resu_", zero_pad(i, 5), ".RData")
        save(output, file = filename_output)    
        
        j <- j + 1
      }
      return(p)
    })
    
    # Stop the cluster after each batch
    cat("\nStopping cluster for batch", batch_index, "\n")
    stopCluster(cl)
    
    # Garbage collection to free memory before next batch
    gc()
  }
  
  cat("\nAll batches completed at ", lubridate::now(), "\n")
  
  
}else{
  
  
  # big_release<-model_pars$mgmt$release_schedule_master%>%dplyr::filter(release_size==24,age_release==2)%>%pull(r)
  # set.seed(19)
  # i<-all_iterations%>%filter(ReleaseStrat%in%big_release)%>%pull(i)%>%sample(1)
  i<-1
  # iterations_to_run<-iterations_to_run%>%
  #   dplyr::filter(i%in%(iteration_chunks%>%unlist%>%head(30)) )
  for(i in iterations_to_run$i){
    
    idx<-i  # clone iteration index 
    p <- iterations_to_run%>%filter(i==idx)%>%pull(p) # get parameter sampling index
    
    all_iterations<-all_iterations%>%dplyr::arrange(i)
    suppressMessages({
    source("./Parameters/pars_postPriorSampling.R",local = T)
    })
    
    set.seed(prior_rng_seed+p)
    init_pop<-init_population(pars=init_pars)
    set.seed(prior_rng_seed+p)
    init_pop<-pairing(pop=init_pop,currentT = 0,pars=init_pars)
    
    start_conditions<-list(init_pop=init_pop)
    
    # model_pars$mgmt$supp_feeding_df<-model_pars$mgmt$supp_feeding_master[[as.character(all_iterations$SuppFeed[i])]]
    # model_pars$mgmt$release_schedule<-model_pars$mgmt$release_schedule_master%>%
      # dplyr::filter(r==all_iterations$ReleaseStrat[i])
    
    # model_pars$mgmt$acc_period_df$acc_period<-1
    set.seed(p+50)
    output<-run_model(start_conditions=start_conditions,model_pars=model_pars,idx = i)
    output$run_pars <- all_iterations[i,]
    output$run_label <- all_iterations$Label[i]
    
    # Save output
    filename_output <- paste0("./Results/", output$run_label, "_Resu_", zero_pad(i, 5), ".RData")
    save(output, file = filename_output)    
  
}

}
