# ------------------------------------------------------------------------------
# OVERVIEW ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# calculate the proportion of females breeding and the average clutch size 

# Code authored 2024, by Finnbar Lee, finnbar.lee@cawthron.org.nz. and
# Nicole Whitelock @ the University of Auckland
# Modified by Caio Kenup @ Zoological Society of London


nest_data_analysis<-list(priors=list(),bio=list())
# ------------------------------------------------------------------------------
# load data --------------------------------------------------------------------
# ------------------------------------------------------------------------------

sightings <- read_xlsx(path = "./Data/Sightings_data_1991_2024.xlsx")
capture_morphometrics <- read.csv("./Data/captures_morphometrics_data_1990_2025.csv")
egg_data <- read.csv("./Data/egg_level_data_2000_2024.csv")
nest_data <- read.csv("./Data/nest_level_data_2000_2024_og.csv")


# ------------------------------------------------------------------------------
# Define start of dataset ------------------------------------------------------
# ------------------------------------------------------------------------------

# default - unreliable records before 2008 - Nicole
# min_year <- 2008
min_year <- 2010

nest_data <- nest_data %>%
  dplyr::mutate(year = as.numeric(str_extract(season, "\\d{4}"))) %>% 
  dplyr::filter(year >= min_year)

# ------------------------------------------------------------------------------
# prop females breeding --------------------------------------------------------
# ------------------------------------------------------------------------------

years_cov<-data.frame(Year=min_year:2025)%>%
  dplyr::mutate(NewMgmt=Year%in%2021:2025)

# prop females breeding
fb <- fun_prop_female_breed(sightings_data = sightings,
                            egg_data = egg_data,
                            capture_morphometrics = capture_morphometrics,
                            years = min_year:2025,
                            years_cov=years_cov)

# clutch size - nest averaged
cs <- fun_clutch_size(nest_data = nest_data,years = min_year:2025,years_cov = years_cov)

# nest_data_analysis$priors$prop_fem_breeding<-data.frame(mean=fb$av[2],sd=fb$se[2],dist="norm")

nest_data_analysis$priors$cs_intercept = data.frame(mean=cs$coeffs$estimate[1],sd=cs$coeffs$se[1],dist="norm")
nest_data_analysis$priors$cs_newsq = data.frame(mean=cs$coeffs$estimate[2],sd=cs$coeffs$se[2],dist="norm")

nest_data_analysis$priors$fb_intercept = data.frame(mean=fb$coeffs$estimate[1],sd=fb$coeffs$se[1],dist="norm")
nest_data_analysis$priors$fb_newsq = data.frame(mean=fb$coeffs$estimate[2],sd=fb$coeffs$se[2],dist="norm")

# nest_data_analysis$priors$clutch_size<-data.frame(mean=cs$clutch_size,sd=cs$sd_clutch_size,dist="norm")

save(nest_data_analysis,file = "./Data/NestDataAnalysis.RData")

rm(cs,fb,years_cov,sightings,capture_morphometrics,egg_data,nest_data,min_year)
