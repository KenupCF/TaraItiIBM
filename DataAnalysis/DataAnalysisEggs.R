# ------------------------------------------------------------------------------
# OVERVIEW ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Calculate egg hatching and hatchling survival rates. Here I have ignored fertile
# versus infertile eggs and just have probability of any eggs hatching. "We still
# don’t really have a good handle on this (some are early embryo death)" - Ilse Corkery

# Code originally from Thalassa Hamilton 2021. Subsequently updated by
# Finnbar Lee in 2024, finnbar.lee@cawthron.org.nz. 

# ------------------------------------------------------------------------------
# load data --------------------------------------------------------------------
# ------------------------------------------------------------------------------


egg_data_analysis<-list(priors=list(),bio=list())



# load egg monitoring data
NZFT <- read.csv("./Data/egg_level_data_2000_2024.csv")


# ------------------------------------------------------------------------------
# Define start of dataset ------------------------------------------------------
# ------------------------------------------------------------------------------

# default
# min_year <- 1990

# interest in seeing effect of only including 2010 and onwards data 
min_year <- 2010

years_cov<-data.frame(season=min_year:2025)%>%
  mutate(NewMgmt=factor(season%in%2021:2025),year=season)


# ------------------------------------------------------------------------------
# Data cleaning ----------------------------------------------------------------
# ------------------------------------------------------------------------------

NZFT <- NZFT %>%
  mutate(year = as.numeric(str_extract(season, "\\d{4}"))) %>% 
  filter(year >= min_year)

NZFT <-  NZFT %>%
  dplyr::select(season, nest_id, egg_id, clutch_no, total_eggs_clutch, fertile,
         egg_managed, f_age, female_id, hatched, fledged, release_site,
         natalsite, foster_site, chick_captive_reared)

# unsure natal sites - Thalia doesn't mention what she did with these
NZFT$natalsite[NZFT$natalsite == "Waipu?"] = "Waipu"
NZFT$natalsite[NZFT$natalsite == "Waip?"] = "Waipu"

# ignore warnings - converting text to NA
# remove nests 287 and 288, where an infertile egg was accidentally returned from the zoo - Thalia
# remove egg management == "add" or "shift" for model convergence - Thalia - doesn't cause an issue now
# remove clutch_no  == 4 for model convergence - Thalia - doesn't cause an issue now
# scale age - Thalia

# rescale age to continuous
f_age = scale(as.integer(NZFT$f_age))

NZFT_2 <- NZFT %>%
  mutate(fertile = ifelse(fertile == "unknown", NA, fertile),
         fertile = as.integer(fertile)) %>%
  filter(!nest_id %in% c("288", "287")) %>%
  mutate(clutch_no = as.numeric(clutch_no),
         nest_id = as.factor(as.integer(nest_id)),
         f_age = f_age,
         man = as.factor(ifelse(egg_managed == "none", "n", "y")),
         season = as.integer(str_extract(season, "[^/]+")),
         female_id = as.factor(ifelse(female_id == "unbanded female", NA, female_id))) 

# stuff needed for model
NZFT_3 <- NZFT_2 %>%
  dplyr::select(female_id, season, natalsite, f_age, clutch_no, man, hatched) # %>%
#filter(!female_id %in% cap_birds)


NZFT_4<-NZFT_3%>%
  left_join(years_cov)

prop_managed<-NZFT_4%>%
  group_by(NewMgmt)%>%
  summarise(propManaged=mean(man=="y"))

# ------------------------------------------------------------------------------
# Probability of any egg hatching  ---------------------------------------------
# ------------------------------------------------------------------------------

# start with female_id, season and breeding site as random effects
mod_res <- glmer(hatched  ~ 1 + clutch_no + f_age + man +
                   NewMgmt +
                   (1|female_id) + (1|season) + (1|natalsite),
                 data = NZFT_4,
                 family = 'binomial',
                 control = glmerControl(optimizer = "bobyqa"))

# no model convergence issues
summary(mod_res)


# # remove natal site for model convergence issue
# mod_res <- glmer(hatched  ~ 1 + clutch_no + f_age + man +
#                    (1|female_id) + (1|season),
#                  data = NZFT_3,
#                  family = 'binomial',
#                  control = glmerControl(optimizer = "bobyqa"))
# 
# # no model convergence issues
# summary(mod_res)


# diagnostic checks
simulationOutput <- simulateResiduals(fittedModel = mod_res,
                                      n = 500,
                                      plot = T)
# 
# plot(simulationOutput) 

# extract random effect for season
mod_re_stv <- broom.mixed::tidy(mod_res, effects = "ran_pars")
# mod_re_stv$estimate[2]

## test for overdispersion - ok
testDispersion(simulationOutput)

# plot model fixed effects
# plot_model(mod_res, type="est",
#            transform = "plogis",
#            show.intercept = T)

# pull estimates from plotting
probs1 <- ggpredict(mod_res, terms = c("f_age", "man"))
probs1 <- as.data.frame(probs1)

# ggplot() +
#   geom_line(data = probs1, aes(x, predicted, colour = group))


# extract mean values from the model - logit scale. Note SE can't be directly transformed
# have to use "delta method". PVA takes logit as input so all g. 
suv_est <- summary(emmeans::emmeans(mod_res, c("f_age")))
# suv_est <- summary(emmeans::emmeans(mod_res, c("f_age","NewMgmt")))
# suv_est <- summary(emmeans::emmeans(mod_res, c("f_age","man","NewMgmt")))


hatch_prob_mod<-mod_res

# Combine them in a tidy data frame
hatch_prob_coeffs<-data.frame(est = fixef(hatch_prob_mod), se = sqrt(diag(vcov(hatch_prob_mod))))
hatch_prob_coeffs$par<-rownames(hatch_prob_coeffs)


egg_data_analysis$priors$hatch_prob_coeff_intercept<-data.frame(
  mean=hatch_prob_coeffs%>%
  filter(par=="(Intercept)")%>%pull(est),
  sd=hatch_prob_coeffs%>%
    filter(par=="(Intercept)")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$hatch_prob_coeff_managed<-data.frame(
  mean=hatch_prob_coeffs%>%
    filter(par=="many")%>%pull(est),
  sd=hatch_prob_coeffs%>%
    filter(par=="many")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$hatch_prob_coeff_sq<-data.frame(
  mean=hatch_prob_coeffs%>%
    filter(par=="NewMgmtTRUE")%>%pull(est),
  sd=hatch_prob_coeffs%>%
    filter(par=="NewMgmtTRUE")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$hatch_prob_coeff_clutch_no<-data.frame(
  mean=hatch_prob_coeffs%>%
    filter(par=="clutch_no")%>%pull(est),
  sd=hatch_prob_coeffs%>%
    filter(par=="clutch_no")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$hatch_prob_coeff_f_age<-data.frame(
  mean=hatch_prob_coeffs%>%
    filter(par=="f_age")%>%pull(est),
  sd=hatch_prob_coeffs%>%
    filter(par=="f_age")%>%pull(se),
  dist="norm")
egg_data_analysis$bio$hatch_prob_coeff_temp_re<-mod_re_stv$estimate[2]


# mean and SE
# pva_params$mean_est[pva_params$parameter == "prob egg hatch"] <- round(plogis(suv_est$emmean), 3)
# pva_params$se_est[pva_params$parameter == "prob egg hatch"] <- round(suv_est$SE, 3)
# pva_params$temp_est[pva_params$parameter == "prob egg hatch"] <- round(mod_re_stv$estimate[2],3)
# pva_params$comments[pva_params$parameter == "prob egg hatch"] <- "SE & temp_est on logit scale"

# ------------------------------------------------------------------------------
# probability hatchling fledging -----------------------------------------------
# ------------------------------------------------------------------------------

# Note some of the GLMM's have issues, so need to remove some variables

# stuff needed
NZFT_fl <- NZFT_2 %>%
  filter(hatched == 1) %>%
  dplyr::select(female_id, season, natalsite, f_age, clutch_no, man, fledged) %>%
  na.omit()%>%
  left_join(years_cov)

# start with female_id, season and breeding site as random effects
# singularity issue

# fledge_mod <- glmer(fledged ~ 1 + clutch_no + f_age + man + 
                      # NewMgmt +
                      # (1|female_id) + (1|season) + (1|natalsite),
                    # data = NZFT_fl,
                    # family = 'binomial',
                    # control=glmerControl(optimizer = "bobyqa"))

# remove female_id and natal site 
fledge_mod <- glmer(fledged  ~ 1 + clutch_no + f_age + man + 
                      NewMgmt+
                      (1|season) ,
                    family = 'binomial',
                    control=glmerControl(optimizer = "bobyqa"),
                    data = NZFT_fl)

# summary(fledge_mod)

# extract random effect for season
mod_re_stv <- broom.mixed::tidy(fledge_mod, effects = "ran_pars")
# mod_re_stv$estimate[1]

# model diagnostics - ok
simulationOutput <- simulateResiduals(fittedModel = fledge_mod,
                                      n = 500,
                                      seed = 123)
# plot(simulationOutput)

## test for overdispersion
testDispersion(simulationOutput)

# extract probability for plotting
probs1 <- ggpredict(fledge_mod, terms = c("f_age"))
probs1 <- as.data.frame(probs1)

gg<-ggplot(probs1) +
  geom_ribbon(aes(x = x,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = group),
              alpha = 0.3) +
  geom_line(aes(x,
                predicted,
                colour = group)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Female age",
       y = "Prob. fledge",
       colour = "Managed",
       fill = "Managed") +
  theme_bw()

# extract mean values from the model - logit scale. Note SE can't be directly transformed
# have to use "delta method". PVA takes logit as input so all g. 
# suv_est <- summary(emmeans::emmeans(fledge_mod, "clutch_no"))
# suv_est <- summary(emmeans::emmeans(fledge_mod, c("clutch_no","man")))
suv_est <- summary(emmeans::emmeans(fledge_mod, c("clutch_no","man","NewMgmt")))

fledge_mod

# Combine them in a tidy data frame
fledge_prob_coeffs<-data.frame(est = fixef(fledge_mod), se = sqrt(diag(vcov(hatch_prob_mod))))
fledge_prob_coeffs$par<-rownames(fledge_prob_coeffs)


egg_data_analysis$priors$fledge_prob_coeff_intercept<-data.frame(
  mean=fledge_prob_coeffs%>%
    filter(par=="(Intercept)")%>%pull(est),
  sd=fledge_prob_coeffs%>%
    filter(par=="(Intercept)")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$fledge_prob_coeff_managed<-data.frame(
  mean=fledge_prob_coeffs%>%
    filter(par=="many")%>%pull(est),
  sd=fledge_prob_coeffs%>%
    filter(par=="many")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$fledge_prob_coeff_sq<-data.frame(
  mean=fledge_prob_coeffs%>%
    filter(par=="NewMgmtTRUE")%>%pull(est),
  sd=fledge_prob_coeffs%>%
    filter(par=="NewMgmtTRUE")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$fledge_prob_coeff_clutch_no<-data.frame(
  mean=fledge_prob_coeffs%>%
    filter(par=="clutch_no")%>%pull(est),
  sd=fledge_prob_coeffs%>%
    filter(par=="clutch_no")%>%pull(se),
  dist="norm")
egg_data_analysis$priors$fledge_prob_coeff_f_age<-data.frame(
  mean=fledge_prob_coeffs%>%
    filter(par=="f_age")%>%pull(est),
  sd=fledge_prob_coeffs%>%
    filter(par=="f_age")%>%pull(se),
  dist="norm")

egg_data_analysis$bio$fledge_prob_coeff_temp_re<-mod_re_stv$estimate[1]
egg_data_analysis$mgmt$prop_managed<-prop_managed

save(egg_data_analysis,file = "./Data/EggDataAnalysis.RData")

rm(NZFT_3,NZFT_4,NZFT_fl,probs1,simulationOutput,suv_est,years_cov,min_year,
   fledge_mod,fledge_prob_coeffs,gg,hatch_prob_mod,hatch_prob_coeffs,mod_re_stv,mod_res,
   NZFT,NZFT_2,f_age)

# END --------------------------------------------------------------------------

