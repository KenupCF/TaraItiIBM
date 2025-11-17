# ============================
# Master run script
# ============================

runLabel <- "bigRunV6"                    # label for this run (used in outputs and sheet)
# runLabel <- "foo"                    # label for this run (used in outputs and sheet)
get_runs_from_gsheet <- FALSE             # pull scheduled runs from Google Sheet
replace_runs_gsheet   <- FALSE            # overwrite Google Sheet "Runs" tab with computed runs
prior_rng_seed <- 1e3                     # base seed for prior sampling per p
runAnalysis <- FALSE                      # whether to re-run the analysis scripts
wd <- "/kenup/TaraItiIBM"                      # preferred working directory

# Fallback working directories for different machines
if (!dir.exists(wd)) { wd <- "D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM" }
if (!dir.exists(wd)) { wd <- "C:/Users/Caio.Kenup/TaraItiIBM" }
setwd(wd)

# Load packages and functions
source("packageLoader.R")
source("functionLoader.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")

# Google services auth
sheet_url <- "https://docs.google.com/spreadsheets/d/1dCIkkofz0h2s9MWOtZfNqAN4DKY59Z2IqNIWlM3isMY/edit?gid=903185379#gid=903185379"
gs4_auth(path = "./.tokens/fresh-replica-344321-0e0618a3b5de.json")

# Google Drive auth
token <- readRDS(".tokens/token.rds")
drive_auth(token = token)

# Identify device (used to assign scheduled runs to machines)
source(".tokens/setDeviceName.R")

# Root parameter container
model_pars <- list(priors = list(), sim = list(), bio = list(), mgmt = list())

# Bring in aggregated expert elicitation and sim parameter scaffolding
load("./Data/Tara_Iti_Aggregated.RData")
source("./Parameters/SimPars.R")

# Optionally re-run data analysis scripts to refresh priors and bio blocks
if (runAnalysis) {
  source("./DataAnalysis/DataAnalysisNesting.R")
  source("./DataAnalysis/DataAnalysisEggs.R")
  source("./DataAnalysis/DataAnalysisSurvival.R")
} else {
  # no-op
}

# Load any additional .RData files under ./Data into the global env
datafiles <- list.files(path = "./Data/", pattern = "\\.RData$", full.names = TRUE)
for (f in datafiles) load(f, envir = .GlobalEnv)

# Merge priors, bio, and mgmt blocks from each analysis domain into model_pars
model_pars$priors <- c(model_pars$priors,
                       egg_data_analysis$priors,
                       nest_data_analysis$priors,
                       surv_data_analysis$priors)

model_pars$bio <- c(model_pars$bio,
                    egg_data_analysis$bio,
                    nest_data_analysis$bio,
                    surv_data_analysis$bio)

model_pars$mgmt <- c(model_pars$mgmt,
                     egg_data_analysis$mgmt,
                     nest_data_analysis$mgmt,
                     surv_data_analysis$mgmt)

# Load fixed biological and management parameter definitions
source("./Parameters/BioPars.R")
source("./Parameters/MgmtPars.R")

# Initial population construction parameters for init_population()
init_pars <- list(
  StartN_df     = model_pars$bio$StartN_df,
  sex_ratio     = 0.5,
  breeding_age  = 3,
  no_age_classes= model_pars$bio$inherent$max_age,
  max_age       = model_pars$bio$inherent$max_age,
  age_structure = model_pars$bio$inherent$age_structure,
  Fp            = model_pars$bio$gen$starting_inbreeding
)

# Sample uncertain parameters according to priors
prior_rng <- priorSampling(
  model_pars$priors,
  seed = 26051991,
  size = model_pars$sim$n_iter
) %>%
  dplyr::mutate(p = 1:n())                 # index per parameter draw

# Apply prior transforms and compute effects (writes into prior_rng and model_pars)
source("./Parameters/priorHandling.R")

# Optional adjustments to age structure (kept commented)
# model_pars$bio$inherent$age_structure[1:2] <- model_pars$bio$inherent$age_structure[1:2] * .5
# model_pars$bio$inherent$age_structure <- model_pars$bio$inherent$age_structure / sum(model_pars$bio$inherent$age_structure)

# Cross product of parameter draws p with management strategies alt
all_iterations <- merge(prior_rng, model_pars$mgmt$strategies) %>%
  merge(sensitivity_analysis)%>%
  dplyr::arrange(p, alt, q) %>%
  dplyr::mutate(
    i = (1:n()) + model_pars$sim$idx_add,  # unique iteration id
    Label = runLabel
  )

# Optional integration with Google Sheet to schedule or replace runs
if (get_runs_from_gsheet | replace_runs_gsheet) {
  runs0 <- read_sheet(sheet_url, sheet = "Runs")
  runs0$Label <- as.character(runs0$Label)
  runs <- runs0 %>%
    replace_na_characters() %>%
    dplyr::filter(Label == runLabel)
  
  # Left join adds p and alt from current all_iterations to existing sheet schedule
  runs2 <- left_join(
    all_iterations %>%
      dplyr::mutate(Iteration = i) %>%
      dplyr::select(Label, alt, Iteration, p),
    runs
  ) %>%
    dplyr::arrange(Label, Iteration) %>%
    dplyr::mutate(ID = paste(Label, zero_pad(Iteration, 5), sep = "_"))
  
  # If requested, overwrite the Runs tab with updated plan
  if (replace_runs_gsheet) {
    write_sheet(runs2, sheet_url, sheet = "Runs")
  }
}

# Choose which iterations to run in this session
if (get_runs_from_gsheet) {
  n_p_used <- 1e3
  
  # Count how many alts per parameter draw p are scheduled
  p_summ <- runs2 %>%
    dplyr::filter(Label == runLabel, Scheduled == TRUE) %>%
    dplyr::group_by(p) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::arrange(n)
  
  # Take scheduled runs for this device, keeping only a subset of p values
  iterations_to_run <- runs2 %>%
    dplyr::filter(Label == runLabel, Scheduled == TRUE, DeviceToRun == device_name) %>%
    dplyr::arrange(dplyr::desc(p)) %>%
    dplyr::group_by(p) %>%
    dplyr::arrange(dplyr::desc(alt)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(i = Iteration) %>%
    dplyr::filter(p %in% p_summ$p[1:n_p_used]) %>%
    dplyr::select(alt, p, i) %>%
    dplyr::ungroup()
  
} else {
  n_p_used <- 1e3
  
  # Count alts per p in full cross product
  p_summ <- all_iterations %>%
    dplyr::group_by(p) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n))
  
  # Prioritize higher p, then higher alt, limit to first n_p_used p values
  iterations_to_run <- all_iterations %>%
    dplyr::arrange(dplyr::desc(p)) %>%
    dplyr::group_by(p) %>%
    dplyr::arrange(dplyr::desc(alt)) %>%
    dplyr::filter(p %in% p_summ$p[1:n_p_used]) %>%
    dplyr::select(alt, p, i) %>%
    dplyr::ungroup()
}

# Safety check: number of parallel workers must not exceed iterations
if (nrow(iterations_to_run) < model_pars$sim$clusters_to_run && model_pars$sim$parallel_across_runs) {
  stop("More clusters than iterations")
}

# Run summary banner
cat(paste0(
  "\n running a model for ", model_pars$sim$n_years, " time steps",
  ", for ", nrow(iterations_to_run), " runs,",
  ifelse(model_pars$sim$parallel_across_runs,
         paste0(" in parallel using ", model_pars$sim$clusters_to_run, " nodes"),
         " sequentially"),
  "\n"
))

# ============================
# Parallel execution path
# ============================
if (model_pars$sim$parallel_across_runs) {
  require(parallel)
  
  model_pars$sim$print_crit <- Inf   # disable in-loop printing
  
  # Batch size = how many runs to send through before recycling the cluster
  batch_size <- model_pars$sim$batching_clusters * model_pars$sim$clusters_to_run
  
  # Compute helper columns for balanced batching, then sort
  iterations_to_run <- iterations_to_run %>%
    dplyr::group_by(p) %>%
    dplyr::mutate(n = n(), dummy = 1) %>%  # dummy used for split_evenly_by_col
    dplyr::ungroup() %>%
    dplyr::arrange(n, p, alt)
  
  # Keep a backup copy for later checks
  iterations_to_run_bkp <- iterations_to_run
  
  # Split into roughly even groups by count using a greedy bin packer
  all_batches <- split_evenly_by_col(
    df        = iterations_to_run,
    n_groups  = ceiling(nrow(iterations_to_run_bkp) / model_pars$sim$batching_clusters),
    target_col= "dummy"
  )
  
  # Flatten batches into a single df and assign batch ids of length batch_size
  batch_summ <- plyr::rbind.fill(all_batches) %>%
    dplyr::arrange((group), p, alt) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(batch = NA_integer_)
  
  batchLoop <- 1
  for (i in 1:nrow(batch_summ)) {
    batch_summ$batch[i] <- batchLoop
    if (sum(batch_summ$batch == batchLoop, na.rm = TRUE) == batch_size) batchLoop <- batchLoop + 1
  }
  
  # Two-level batching: batch -> sub-batches per cluster worker
  iteration_batches <- split(batch_summ, batch_summ$batch)
  iteration_batches <- lapply(
    iteration_batches,
    function(x) split(x, ceiling(seq_along(x[, 1]) / (batch_size / model_pars$sim$clusters_to_run)))
  )
  
  # Iterate batches sequentially, each on a fresh cluster
  for (batch_index in seq_along(iteration_batches)) {
    sub_batches <- iteration_batches[[batch_index]]
    cluster_size <- min(model_pars$sim$clusters_to_run, length(sub_batches))
    
    # Extract iteration id vectors for each worker and randomize within
    iteration_chunks <- lapply(sub_batches, function(x) { sample(x$i) })
    
    cat(paste0("\nStarting batch ", batch_index, " at ", lubridate::now(), "\n"))
    
    # Start cluster for this batch
    cl <- makeCluster(cluster_size)
    
    # Initialize each worker with code and secrets, silence messages
    invisible(clusterEvalQ(cl, {
      source("packageLoader.R")
      source("functionLoader.R")
      devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")
      devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")
      
      token <- readRDS(".tokens/token.rds")
      drive_auth(token = token)
      
      con <- file(nullfile(), open = "w")
      sink(con, type = "message")    # suppress messages on workers
    }))
    
    # Send shared objects to workers
    cat("\n Passing objects to nodes\n")
    clusterExport(cl, varlist = c("init_pars", "model_pars", "prior_rng_seed", "iterations_to_run", "all_iterations"))
    
    cat("\nStarting parallel execution for batch ", batch_index, "\n")
    
    # Main parallel loop: each worker handles a chunk of iteration ids
    parLapply(cl = cl, iteration_chunks, fun = function(chunk_items) {
      j <- 1
      for (cc in seq_along(chunk_items)) {
        i <- chunk_items[cc]                          # iteration id
        idx <- i                                      # alias
        
        # Find the parameter draw index p associated with this iteration
        p <- iterations_to_run %>% dplyr::filter(i == idx) %>% dplyr::pull(p)
        
        # Prepare per-iteration parameters based on p and alt
        all_iterations <- all_iterations %>% dplyr::arrange(i)
        suppressMessages({
          source("./Parameters/pars_postPriorSampling.R", local = TRUE)
        })
        
        # Build initial population for this iteration (deterministic under seed)
        set.seed(prior_rng_seed + p)
        init_pop <- init_population(pars = init_pars)
        set.seed(prior_rng_seed + p)
        init_pop <- pairing(pop = init_pop, currentT = 0, pars = init_pars)
        start_conditions <- list(init_pop = init_pop)
        
        # Run the simulation
        set.seed(p + 50)
        output <- run_model(start_conditions = start_conditions, model_pars = model_pars, idx = i)
        output$run_pars  <- all_iterations[i, ]
        output$run_label <- all_iterations$Label[i]
        
        # Persist result to disk
        filename_output <- paste0("./Results/", output$run_label, "_Resu_", zero_pad(i, 6), ".RData")
        save(output, file = filename_output)
        
        j <- j + 1
      }
      return(p)
    })
    
    # Cleanup for this batch
    cat("\nStopping cluster for batch ", batch_index, "\n")
    stopCluster(cl)
    gc()
  }
  
  cat("\nAll batches completed at ", lubridate::now(), "\n")
  
  # ============================
  # Sequential execution path
  # ============================
} else {
  # Example: choose one specific iteration or loop all
  # i <- 1
  for (i in iterations_to_run$i) {
    idx <- i
    p <- iterations_to_run %>% dplyr::filter(i == idx) %>% dplyr::pull(p)
    
    all_iterations <- all_iterations %>% dplyr::arrange(i)
    suppressMessages({
      source("./Parameters/pars_postPriorSampling.R", local = TRUE)
    })
    
    
    set.seed(prior_rng_seed + p)
    init_pop <- init_population(pars = init_pars)
    set.seed(prior_rng_seed + p)
    init_pop <- pairing(pop = init_pop, currentT = 0, pars = init_pars)
    start_conditions <- list(init_pop = init_pop)
    
    set.seed(prior_rng_seed + p)
    output <- run_model(start_conditions = start_conditions, model_pars = model_pars, idx = i)
    output$run_pars  <- all_iterations[i, ]
    output$run_label <- all_iterations$Label[i]
    
    filename_output <- paste0("./Results/", output$run_label, "_Resu_", zero_pad(i, 5), ".RData")
    save(output, file = filename_output)
  }
}

