wd <- "/kenup/TaraItiIBM"                      # preferred working directory

# Fallback working directories for different machines
if (!dir.exists(wd)) { wd <- "D:/03-Work/01-Science/00-Research Projects/Tara Iti/TaraItiIBM" }
if (!dir.exists(wd)) { wd <- "C:/Users/Caio.Kenup/TaraItiIBM" }
setwd(wd)

source("RunningModels.R")