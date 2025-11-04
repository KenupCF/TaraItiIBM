
# model_pars$mgmt$field_opts<-
# 
# model_pars$mgmt$egg_harvest_rate<-

# model_pars$mgmt$admix_releases<-
# model_pars$mgmt$gen_mgmt<-c(FALSE,TRUE)

model_pars$mgmt$strategies<-expand.grid(
    egg_harvest_rate=c(0,.1
                       ,.5,1
                       ),
    admix_releases=c(FALSE
                     ,TRUE
                     ),
    admix_prop_released=0.2,
    admix_age_released=2,
    admix_no_released=6,
    gen_mgmt=c(FALSE
               # ,TRUE
               ),
    field=c("SQ","Field2"))%>%
  # dplyr::filter(!(field=="Field2" & admix_releases))%>%
  dplyr::filter(!(field=="Field2" & gen_mgmt))%>%
  dplyr::filter(!(admix_releases & gen_mgmt))%>%
  # dplyr::filter(!(field=="Field2" & egg_harvest_rate!=min(egg_harvest_rate)))%>%
  dplyr::filter(!(gen_mgmt & egg_harvest_rate>=0.5))%>%
  dplyr::filter(!(admix_releases & egg_harvest_rate>=0.5))%>%
  ungroup()%>%
  mutate(alt=1:n())
