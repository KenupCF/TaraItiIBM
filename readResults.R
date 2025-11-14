# Load required packages
library(dplyr)
library(progress)
library(duckdb)
library(stringr)
library(lubridate)

# Turn warnings into errors so failures are caught early
options(warn = 2)


source(".\\Functions\\FUN.R")
# =======================
# SETUP
# =======================
folder_extr <- "D:\\03-Work\\01-Science\\00-Research Projects\\Tara Iti\\TaraItiIBM\\Results\\bigRunV4"
folderID    <- gsub(x = folder_extr, "^.*/", "")   # Extract last path segment as a folder label
loopSize    <- 12e3                                  # Max number of files to attempt per run
time_limit_secs <- 60*60*(20/60)                  # time limit in seconds
buffer_size <- 25                                   # Number of records to batch-write into DuckDB

# =======================
# Discover result files and dedupe by replicate index i
# =======================
r_files <- list.files(path = folder_extr, pattern = ".RData", full.names = TRUE)

# Parse replicate index from filenames
idx <- r_files %>%
  gsub(".RData", "", .) %>%
  gsub(" \\(.*", "", .) %>%
  substr(nchar(.) - 5, nchar(.)) %>%
  gsub("_", "", .) %>%
  as.numeric()

# Data frame of files and parsed indices
files_df <- data.frame(file = r_files, i = idx) %>%
  mutate(nchar = nchar(file)) %>%
  arrange(desc(nchar))

# Remove duplicate indices, keeping the longest filename entry
duplicated_files <- files_df %>%
  filter(duplicated(i, fromLast = TRUE)) %>%
  pull(file)
sapply(duplicated_files, file.remove)   # Delete duplicate files on disk

files_df <- files_df %>% filter(!duplicated(i))

# =======================
# Connect to DuckDB
# =======================

db_path <- paste0(folder_extr, "\\bigRunV4.duckdb")
con <- dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)

# Find which i have already been imported for this folder
if (!"summary" %in% dbListTables(con)) {
  already_imported <- numeric(0)
} else {
  already_imported <- dbGetQuery(con, "SELECT * FROM summary") %>%
    dplyr::filter(folder == folderID) %>%
    pull(i)
}

# Keep only new files
files_df <- files_df %>% filter(!i %in% already_imported)

cat(paste0(nrow(files_df), " files still left to read."))

# Limit to loopSize
files_df <- files_df[1:min(loopSize, nrow(files_df)), ]

# =======================
# Progress bar setup
# =======================
pb <- progress_bar$new(
  total = nrow(files_df),
  format = "  Processing [:bar] :percent ETA: :eta",
  clear = FALSE,
  width = 60
)
counter <- 0

# =======================
# Prepare write buffers
# Each element will become a list of data frames to bind and append
# =======================
buffer <- list(
  summary  = vector("list", length = buffer_size),
  N_series = vector("list", length = buffer_size),
  egg_fate = vector("list", length = buffer_size),
  mgmt     = vector("list", length = buffer_size),
  run_pars = vector("list", length = buffer_size)
)

start_time <- Sys.time()

# =======================
# Main loop over files
# =======================
if(nrow(files_df)>0){
for (f in seq_len(nrow(files_df))) {
  
  # Process one file with error handling
  result <- tryCatch({
    process_result_file(files_df$file[f], folder_id = folderID)
  }, error = function(e) {
    message(sprintf("Error processing file %s: %s", files_df$file[f], e$message))
    return(NULL)
  })
  
  # Skip if failed
  if (is.null(result)) next
  
  # Collect into buffers
  # buffer$mgmt[[f]]      <- result$mgmt
  buffer$summary[[f]]   <- result$summary
  buffer$N_series[[f]]  <- result$N_series
  buffer$egg_fate[[f]]  <- result$egg_fate
  buffer$run_pars[[f]]  <- result$run_pars
  
  # Periodic writes by buffer size or at the end
  if (f %% buffer_size == 0 || f == nrow(files_df)) {
    # dbWriteTable(con, "mgmt",      bind_rows(buffer$mgmt), append = TRUE)
    # dbWriteTable(con, "egg_fate",  bind_rows(buffer$egg_fate), append = TRUE)
    dbWriteTable(con, "summary",   bind_rows(buffer$summary), append = TRUE)
    dbWriteTable(con, "N_series",  bind_rows(buffer$N_series), append = TRUE)
    dbWriteTable(con, "run_pars",  bind_rows(buffer$run_pars), append = TRUE)
    
    # Reset buffers after writing
    buffer <- list(summary = list(), N_series = list(), egg_fate = list(), mgmt = list(), run_pars = list())
  }
  
  # Update timing and check time limit
  counter <- counter + 1
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (elapsed_time > time_limit_secs) {
    message("Time limit reached. Stopping loop.")
    break
  }
  
  # Tick progress
  pb$tick()
}
}
# Report throughput
time_diff <- difftime(now(), start_time)
cat(paste0("Imported ", counter, " entries in ", round(time_diff, 2), " ", attr(time_diff, "units")))

# =======================
# Cleanup
# =======================
dbDisconnect(con, shutdown = TRUE)
