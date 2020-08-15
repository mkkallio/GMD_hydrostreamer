# script for the hydrostreamer v1.0 model description paper submitted to 
# Geoscientific Model Development.
# The script allows for reproduction of Figures 3, 4 and 5 and getting data 
# shown in Table 4, and Appendix A tables A1-A4.

# total run time of the entire script is approximately 1 hour for a machine
# equipped with Intel i7 processor


# requires the following packages installed:
# raster, tidyverse, sf, lubridate, hydrostreamer, patchwork, ggplot2, readr,
# glue, inlmisc, stringr, MCMCpack, ggrepel, units

# tested and works with the following session:
# > sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
# [3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] patchwork_1.0.1     hydrostreamer_1.0.0 lubridate_1.7.9     sf_0.9-5            forcats_0.5.0      
# [6] stringr_1.4.0       dplyr_1.0.1         purrr_0.3.4         readr_1.3.1         tidyr_1.1.1        
# [11] tibble_3.0.3        ggplot2_3.3.2       tidyverse_1.3.0     raster_3.3-13       sp_1.4-2           
# 
# loaded via a namespace (and not attached):
#     [1] tidyselect_1.1.0   haven_2.3.1        lattice_0.20-41    colorspace_1.4-1   vctrs_0.3.2       
# [6] generics_0.0.2     blob_1.2.1         rlang_0.4.7        e1071_1.7-3        pillar_1.4.6      
# [11] glue_1.4.1         withr_2.2.0        DBI_1.1.0          dbplyr_1.4.4       modelr_0.1.8      
# [16] readxl_1.3.1       lifecycle_0.2.0    munsell_0.5.0      gtable_0.3.0       cellranger_1.1.0  
# [21] rvest_0.3.6        codetools_0.2-16   class_7.3-17       fansi_0.4.1        broom_0.7.0       
# [26] Rcpp_1.0.5         KernSmooth_2.23-17 scales_1.1.1       backports_1.1.8    classInt_0.4-3    
# [31] jsonlite_1.7.0     fs_1.5.0           hms_0.5.3          stringi_1.4.6      grid_4.0.2        
# [36] quadprog_1.5-8     cli_2.0.2          tools_4.0.2        magrittr_1.5       crayon_1.3.4      
# [41] pkgconfig_2.0.3    Matrix_1.2-18      ellipsis_0.3.1     xml2_1.3.2         reprex_0.3.0      
# [46] assertthat_0.2.1   httr_1.4.2         rstudioapi_0.11    R6_2.4.1           units_0.6-7       
# [51] compiler_4.0.2    


library(raster)
library(tidyverse)
library(sf)
library(lubridate)
library(hydrostreamer)
library(patchwork)
library(ggplot2)

# ------------------------------------------------------------------------------
#  LOAD DATA AND DO SOME PREPROCESSING

system.time({
    rasters <- list.files("data/runoff/", full.names = TRUE)
    river <- "Data/Hydrosheds 3S rivers.gpkg"
    basin <- "Data/3S basin.gpkg"
    voronoi <- "Data/HS voronoi.gpkg"
    delbas <- "Data/HS_basins_dune.gpkg"
    
    river <- read_sf(river) %>%
        rename(riverID = ARCID)
    basin <- read_sf(basin)
    voronoi <- read_sf(voronoi) %>% st_transform(4326) %>%
        rename(riverID = ARCID)
    delbas <- read_sf(delbas)  %>%
        rename(riverID = X3S_drdir) %>%
        mutate(invDUNE = 1/DUNE)
    raster <- brick(rasters[1])
    
    # remove WATCH-forced model runs and .aux.xml files if there are any
    remove <- grepl("_watch_", rasters)
    rasters <- rasters[!remove]
    remove <- grepl(".aux.", rasters)
    rasters <- rasters[!remove]
    
    # get model names, and the starting date for each of them
    fcast_names <- vector("character", length(rasters))
    fcast_startdate <- vector("character", length(rasters))
    for (i in seq_along(rasters)) {
        temp <- basename(rasters[[i]])
        fcast_names[[i]] <- word(temp,c(1,2),sep = "_") %>% 
            glue::glue_collapse(sep="_")
        fcast_startdate[[i]] <- paste0(word(temp,9,sep = "_"),"-01-01" )
    }
    
    
    # load station locations and observation timeseries
    HYMOS <- read_sf("Data/HYMOS_stations_at_GHM_locations.gpkg") 
    
    obs <- readr::read_csv("data/observations/3S_station_discharge.csv") %>%
        rename(Date = DATE) %>%
        mutate(Date = lubridate::dmy(Date))
    obs[obs == 0] <- NA
    obs <- obs %>%
        mutate(year = year(Date),
               month = month(Date)) %>%
        group_by(year,month) %>%
        summarise_all(mean, na.rm=TRUE) %>%
        mutate(Date = ymd(paste(year,month,"01", sep="/"))) %>%
        ungroup() %>%
        select(-year, -month)
    
    temp <- colnames(obs)[-1] %>% match(HYMOS$Code)
    temp2 <- HYMOS$ARCID[temp]
    colnames(obs) <- c("Date",as.character(HYMOS$Study_name)[temp])
    
    
    # colour schemes used
    colors1 <- inlmisc::GetColors(11)
    colors2 <- c('#f5cac3', #Sekong d
                 '#e49689', #Sekong m
                 '#ce6153', #Sekong u
                 '#e4dcf6', #Sesan d
                 '#c8baec', #Sesan m
                 '#8c79d8', #Sesan u-n
                 '#ab99e2', #Sesan u-e
                 '#e2dec2', #Srepok d
                 '#c3be87', #Srepok m
                 '#a29e4d', #Srepok u
                 "#FFFFFF", #ungauged
                 "#FAEBD7") #ungauged)
    
    # station order for ggplot
    statorder <- c("Sekong Downstream",
                   "Sekong Midstream",
                   "Sekong Upstream",
                   "",
                   "Sesan Downstream",
                   "Sesan Midstream",
                   "Sesan Upstream-N",
                   "Sesan Upstream-E",
                   "Srepok Downstream",
                   "Srepok Midstream",
                   "Srepok Upstream",
                   " ")

    
    
    # GHM discharge files and model names
    files <- list.files("Data/Discharge/", 
                        full.names = TRUE)
    remove <- grepl("_watch_", files)
    files <- files[!remove]
    remove <- remove <- grepl(".aux.", files)
    files <- files[!remove]
    
    dis_names <- vector("character", length(files))
    startdate <- vector("character", length(files))
    for (i in seq_along(files)) {
        temp <- basename(files[[i]])
        dis_names[[i]] <- word(temp,c(1,2),sep = "_") %>% 
            glue::glue_collapse(sep="_")
        startdate[[i]] <- paste0(word(temp,-2,sep = "_"),"-01-01" )
    }
        
}) # 3 seconds



# ------------------------------------------------------------------------------
#  RUN HYDROSTREAMER

system.time({
    ##############
    # HSgrid
    
    # process only those runoff rasters which also have discharge counterparts
    keep <- which(fcast_names %in% dis_names)
    
    # process rasters to source zones
    system.time({
        HSgrid <- raster_to_HS(rasters[keep], 
                               unit = "mm/s",
                               date = lubridate::ymd(fcast_startdate[keep]),
                               timestep="month",
                               aoi = basin,
                               names = fcast_names[keep],
                               verbose = TRUE)
    }) # 89 seconds if global rasters, 19 secods when already cropped
    
    # create routing information for the river network
    river <- river_network(river, verbose = TRUE)

    
    # prepare observations
    temp <- colnames(obs)[-1] %>% match(HYMOS$Study_name)
    riverIDs <- HYMOS$ARCID[temp]
    names(riverIDs) <- as.character(HYMOS$Study_name)[temp]
    
    ############
    # Run hydrostreamer for the four downscaling methods
    
    # using river segments alone, without catchment estimation
    system.time({
        HS_length <- interpolate_runoff(HSgrid,
                                        river,
                                        verbose = TRUE) %>%
            add_observations(obs, "m3/s", riverIDs) %>%
            accumulate_runoff(verbose=TRUE)
    }) # 22 seconds
    
    # Estimate catchments using segment Voronoi diagram (Thiessen polygons)
    system.time({
        HS_voronoi <- interpolate_runoff(HSgrid,
                                         river, 
                                         basins = voronoi,
                                         verbose = TRUE) %>%
            add_observations(obs,"m3/s", riverIDs) %>%
            accumulate_runoff(verbose=TRUE)
    }) # 47 seconds
    
    # catchments delineated from the flow direction raster
    system.time({
        HS_basin <- interpolate_runoff(HSgrid,
                                       river, 
                                       basins = delbas,
                                       verbose = TRUE) %>%
            add_observations(obs,"m3/s", riverIDs) %>%
            accumulate_runoff(verbose=TRUE)
    }) # 50 seconds
    
    # dasymetric mapping using DUNE as ancillary variable
    system.time({
        HS_bdasy <- interpolate_runoff(HSgrid,
                                       river,
                                       basins = delbas,
                                       dasymetric = "invDUNE",
                                       verbose = TRUE) %>%
            add_observations(obs, "m3/s", riverIDs) %>%
            accumulate_runoff(verbose=TRUE)
    }) # 46 seconds
    
}) # overall this hydrostreamer part here takes 167 seconds



# ------------------------------------------------------------------------------
# Process Discharge predictions (hacky way to create a HS object for this)

system.time({
    
    # coordinates for the monitoring stations for GHM at 0.5 degree resolution
    points <- st_coordinates(HYMOS)
    
    # keep only those discharge predictions which are included in runoff
    keep <- which(dis_names %in% fcast_names)
    
    
    # extract discharge raster timeseries for station locations
    discharge <- list()
    for (i in keep) {
        temp <- raster::brick(files[[i]])
        temp2 <- cellFromXY(temp,points)
        temp3 <- raster::extract(temp,temp2)
        dates <- seq(lubridate::ymd(startdate[[i]]),
                     length.out = ncol(temp3), 
                     by = 'month')
        temp3 <- t(temp3)
        temp3 <- data.frame(dates, temp3)
        colnames(temp3) <- c("Date", as.character(HYMOS$Study_name))
        
        discharge[[ dis_names[i] ]] <- temp3
    }
    
    
    # Ensure that each timeseries has equal dates
    alldates <- lapply(discharge, function(x) x$Date) %>%
        unlist %>% 
        unique %>%
        as_date %>%
        enframe(name = NULL) %>%
        rename(Date = value)
    
    for(i in seq_along(discharge)) {
        discharge[[i]] <- left_join(alldates, discharge[[i]], by="Date")
    }
    temp <- hydrostreamer:::spread_listc(discharge)
    
    
    
    # create a HS object for discharge
    discharge <- HS_length
    sel <- discharge$observation_station %>% unlist %in% names(temp)
    discharge <- discharge[sel,]
    for (i in seq_along(temp)) {
        row <- which(discharge$observation_station %>% 
                         unlist == names(temp)[i])
        
        check <- temp[[i]]
        dates <- discharge$observation_ts[[row]]
        dates <- dates[complete.cases(dates),1] %>% unlist %>% 
            as_date %>% unname
        check <- check[check$Date %in% dates,]
        
        testi <- lapply(check, is.na)
        testi <- lapply(testi, all) %>%
            unlist
        testi <- which(!testi)
        discharge$discharge_ts[[row]] <- check[,testi]
    }
    discharge <- discharge %>%
        select(-runoff_ts)
    
}) # 28 secs 

# Add benchmark data to the discharge object
system.time({
    glofas <- readr::read_csv("data/benchmark/glofas_monthly_discharge_1979-2019.csv") %>% 
        select(-X1)
    
    grades <- readr::read_csv("data/benchmark/grades_monthly_discharge_1979-2013.csv") %>% 
        select(-X1)
    
    for(i in 1:nrow(discharge)) {
        stat <- discharge$observation_station[i]
        gd <- glofas[,c("Date", stat)]
        colnames(gd)[2] <- "GLOFAS"
        discharge$discharge_ts[[i]] <- left_join(discharge$discharge_ts[[i]],
                                                 gd,
                                                 by="Date")
        
        gd <- grades[,c("Date", stat)]
        colnames(gd)[2] <- "GLOFAS"
        discharge$discharge_ts[[i]] <- left_join(discharge$discharge_ts[[i]],
                                                 gd,
                                                 by="Date")
    }
}) # < 1 sec

# Create a HS object with only those river segments with a monitoring station
# This is to speed up some of the computations later on and to reduce 
# memory consumption
keep <- which(!is.na(HS_length$observation_station))
HS_stations <- HS_length[keep,]


# remove unused variables
rm("rasters","remove", "temp", "temp2", "temp3", "dates", "alldates",
   "files", "points", "raster", "sel", "row", 
   "fcast_startdate", "startdate", "i", 
   "testi", "check")










# ------------------------------------------------------------------------------
# Evaluate goodness-of-fit for all predictions

system.time({
    # difference between three downscaling methods
    length_gof <- flow_gof(HS_length)
    voronoi_gof <- flow_gof(HS_voronoi)
    basin_gof <- flow_gof(HS_basin)
    dasy_gof <- flow_gof(HS_bdasy)
    lvb_gof <- bind_rows(add_column(length_gof, Type = "Length", .before=1),
                         add_column(voronoi_gof, Type = "Voronoi", .before=1),
                         add_column(basin_gof, Type = "Basin", .before=1),
                         add_column(dasy_gof, Type = "Dasymetric", .before=1)) %>%
        mutate(Impact_Model = word(Prediction, 1, sep="_"),
               Climate_Forcing = word(Prediction,2,sep="_")) %>%
        as_tibble()
    lvb_gof$Station <- factor(lvb_gof$Station, levels=statorder)
    
    discharge_gof <- flow_gof(discharge)
}) # 5 seconds
    

## performance rank of downscaling methods by ensemble members
lvb_gof %>%
    dplyr::select(Type, Prediction, MAE, RMSE, `PBIAS %`, NSE, KGE, R2) %>%
    group_by(Type, Prediction) %>%
    summarise_all(mean) %>%
    arrange(Prediction, desc(KGE)) %>%
    print(n=Inf) %>% {
        # print rank by ensemble member
        temp <- rbind(Basin = which(.$Type == "Basin") -
                          seq(0,by=4,length.out=15),
                      Voronoi = which(.$Type == "Voronoi") -
                          seq(0,by=4,length.out=15),
                      Length =  which(.$Type == "Length") -
                          seq(0,by=4,length.out=15),
                      Dasymetric =  which(.$Type == "Dasymetric") -
                          seq(0,by=4,length.out=15))
        colnames(temp) <- unique(.$Prediction)
        print(as_tibble(temp))
        
        
        # print average rank
        temp <- c(Basin = (which(.$Type == "Basin") -
                               seq(0,by=4,length.out=15)) %>% mean,
                  Voronoi = (which(.$Type == "Voronoi") -
                                 seq(0,by=4,length.out=15)) %>% mean,
                  Length =  (which(.$Type == "Length") -
                                 seq(0,by=4,length.out=15)) %>% mean,
                  Dasymetric =  (which(.$Type == "Dasymetric") -
                                     seq(0,by=4,length.out=15)) %>% mean)
        print(temp)
    }

# performance rank of downscaling methods by monitoring station
lvb_gof %>%
    dplyr::select(Type, Station, MAE, RMSE, `PBIAS %`, NSE, KGE, R2) %>%
    group_by(Type, Station) %>%
    summarise_all(mean) %>%
    arrange(Station, desc(KGE)) %>%
    print(n=Inf) %>% {
        # print rank by monitoring station
        temp <- rbind(Basin = which(.$Type == "Basin") -
                          seq(0,by=4,length.out=10),
                      Voronoi = which(.$Type == "Voronoi") -
                          seq(0,by=4,length.out=10),
                      Length =  which(.$Type == "Length") -
                          seq(0,by=4,length.out=10),
                      Dasymetric =  which(.$Type == "Dasymetric") -
                          seq(0,by=4,length.out=10))
        colnames(temp) <- unique(.$Station)
        print(temp)
        
        # print average rank
        temp <- c(Basin = (which(.$Type == "Basin") -
                               seq(0,by=4,length.out=10)) %>% mean,
                  Voronoi = (which(.$Type == "Voronoi") -
                                 seq(0,by=4,length.out=10)) %>% mean,
                  Length =  (which(.$Type == "Length") -
                                 seq(0,by=4,length.out=10)) %>% mean,
                  Dasymetric =  (which(.$Type == "Dasymetric") -
                                     seq(0,by=4,length.out=10)) %>% mean)
        print(temp)
    }
    
    
    
    lvb_gof %>% 
        select(Type, Prediction, Station, RMSE, `PBIAS %`, NSE, R2, KGE) %>%
        group_by(Station, Prediction) %>%
        arrange(RMSE)

    lvb_gof %>%
        mutate(model = word(Prediction, 1, sep="_"),
               cf = word(Prediction, 2, sep="_")) %>%
        select(Type, model, cf, Station, RMSE, `PBIAS %`, NSE, R2, KGE) %>%
        arrange(Station, model, cf, RMSE)

   
    ####### Practically no difference at the observation stations - the difference 
    # is in the small streams where we have no observations available -- continue 
    # with only length since it is the most simple one.
    
    

    ### EXPLORE AREA-TO-LINE PREDICTIONS AND PRODUCE APPENDIX TABLE
    
    
    # Appendix A Table A1
    lvb_gof %>%
        mutate(Impact_Model = word(Prediction, 1, sep="_"),
               Climate_Forcing = word(Prediction,2,sep="_")) %>%
        # group_by(Type, Impact_Model, Climate_Forcing) %>%
        group_by(Type, Station) %>%
        summarise(KGE = mean(KGE),
                  NSE = mean(NSE),
                  R2 = mean(R2),
                  RMSE = mean(RMSE),
                  `PBIAS %` = mean(`PBIAS %`)) %>%
        select(Type, Station, 
               RMSE, `PBIAS %`, NSE, R2, KGE) %>% 
        arrange(Station,desc(KGE))
    
    # Appendix A Table A2
    length_gof %>%
        group_by(Prediction) %>%
        summarise(KGE = mean(KGE),
                  NSE = mean(NSE),
                  R2 = mean(R2),
                  RMSE = mean(RMSE),
                  `PBIAS %` = mean(`PBIAS %`)) %>%
        select(Prediction, RMSE, `PBIAS %`, NSE, R2, KGE) %>% 
        arrange(desc(KGE)) %>%
        print(n=Inf) 
    
    
    # Appendix A Table A3
    length_gof %>%
        mutate(Impact_Model = word(Prediction, 1, sep="_"),
               Climate_Forcing = word(Prediction,2,sep="_")) %>%
        group_by(Climate_Forcing) %>%
        summarise(KGE = mean(KGE),
                  NSE = mean(NSE),
                  R2 = mean(R2),
                  RMSE = mean(RMSE),
                  `PBIAS %` = mean(`PBIAS %`)) %>%
        select(Climate_Forcing, RMSE, `PBIAS %`, NSE, R2, KGE) %>% 
        arrange(desc(KGE))
    
    
    # Appendix A Table A4
    length_gof %>%
        mutate(Impact_Model = word(Prediction, 1, sep="_"),
               Climate_Forcing = word(Prediction,2,sep="_")) %>%
        group_by(Impact_Model) %>%
        summarise(KGE = mean(KGE),
                  NSE = mean(NSE),
                  R2 = mean(R2),
                  RMSE = mean(RMSE),
                  `PBIAS %` = mean(`PBIAS %`)) %>%
        select(Impact_Model, RMSE, `PBIAS %`, NSE, R2, KGE) %>% 
        arrange(desc(KGE)) 
    



    # ------------------------------------------------------------------------------
    # Compute ensembles

    # use hydrostreamer's ensemble summary function to compute min, mean, 
    # median, and max for the ensemble
    # for the GHM data, we need to get rid of GRADES and GLOFAS first
    dis_ensemble <- discharge
    dis_ensemble$discharge_ts <- lapply(dis_ensemble$discharge_ts,
                                        function(x) {
                                            x <- dplyr::select(x,
                                                               -GLOFAS,
                                                               -GRADES)
                                        })
    dis_ensemble <- ensemble_summary(dis_ensemble, 
                                     summarise_over_timeseries = FALSE,
                                     drop = TRUE)
    pred_ensemble <- ensemble_summary(HS_stations, 
                                      summarise_over_timeseries = FALSE,
                                      drop = TRUE)
    
    pred_ensemble_gof <- flow_gof(pred_ensemble)
    dis_ensemble_gof <- flow_gof(dis_ensemble)
    

    
    
    
    
    
    # ------------------------------------------------------------------------------
    # Random weights and goodness of fit



system.time({
    set.seed(32648) # for reproducibility
    # staged sampling of
    # 1. number of models for the combination. 
    # 2. which models to pick
    # 3. sample weights from dirichlet distribution
    nm <- 15
    nmodels <- sample(2:nm,1e4,replace=T)
    random_weights <- lapply(nmodels,function(n){
        weights<-rep(0,nm)
        names(weights) <- names(HS_stations$discharge_ts[[1]][-1])
        nonzero_weights<-sample(nm,n)
        weights[nonzero_weights]<-MCMCpack::rdirichlet(1,rep(1,n))
        weights
    }) %>% do.call(rbind,.)
    
    # remove duplicates if any
    random_weights <- random_weights[!duplicated(random_weights),]
    
    # process them so that they can be used in combine_runoff()
    weights <- list()
    for(j in 1:nrow(random_weights)) {
        weights[[ as.character(j) ]] <- random_weights[j,]
    }
    
    # create timeseries from the random weights. Remove runoff_ts column 
    # since we only want to combine the predicted discharge (speed up)
    HS_random <- combine_runoff(HS_stations %>% select(-runoff_ts),
                                weights=weights,
                                #intercept = as.list(rep(0,length(weights))),
                                #bias =  as.list(rep(0,length(weights))),
                                drop=TRUE)
    
}) # ~ 1000 sec

# save the outputs for future so computation does not need to be repeated
# save(HS_random, file = "HSrandom_20200813.RData")
# save(random_weights, file = "random_weights_20200813.RData")

 
# ### this takes a long time
system.time({
    random_gof <- flow_gof(HS_random, verbose=TRUE)
}) # 505 sec

# mean performance of random combinations across all monitoring stations
random_gof %>%
    group_by(Prediction) %>%
    summarise_all(mean) %>%
    select(Prediction, ME, `PBIAS %`, NSE, KGE, R2) %>%
    arrange(desc(KGE))

# comapre to discharge
discharge_gof %>%
    group_by(Prediction) %>%
    summarise_all(mean) %>%
    select(Prediction, ME, `PBIAS %`, NSE, KGE, R2) %>%
    arrange(desc(KGE)) 

# compare to area-to-line predictions
length_gof %>%
    group_by(Prediction) %>%
    summarise_all(mean) %>%
    select(Prediction, ME, `PBIAS %`, NSE, KGE, R2) %>%
    arrange(desc(KGE)) 


## Best random combinations fare much better than any individual GHM




# check average goodness-of-fit for different numbers of models in 
# random combinations
iter <- unique(nmodels) %>% sort
nmod_gof <- vector()
for(i in iter) {
    inds <- which(nmodels == i)
    nmod_gof <- rbind(nmod_gof,
                          random_gof %>%
                              filter(Prediction %in% inds) %>%
                              select(RMSE, `PBIAS %`, NSE, R2, KGE) %>%
                              mutate(n = nrow(.)) %>%
                              summarise_all(mean) %>%
                              add_column(nmodels = i, .before=1))
}
nmod_gof


# check average goodness-of-fit by each ensemble member in
# random combinations
mod_gof <- vector()
for(i in fcast_names) {
    inds <- which(colnames(random_weights) == i)
    if(length(inds) == 0) next
    inds <- which(random_weights[,inds] > 0)
    mod_gof <- rbind(mod_gof, random_gof %>%
                             filter(Prediction %in% inds) %>%
                             select(RMSE, `PBIAS %`, NSE, R2, KGE) %>%
                             mutate(n = nrow(.)) %>%
                             summarise_all(mean) %>%
                             add_column(model = i, .before=1))
}
arrange(mod_gof, RMSE)


# check goodness of fit of all ensemble members, in all ensemble sizes (2-15)
# at all monitoring stations. This produces a large table of 2100 rows.
modn_gof <- vector()
for(j in iter) {
    inds <- nmodels == j
    for(i in fcast_names) {
        inds2 <- which(colnames(random_weights) == i)
        if(length(inds2) == 0) next
        inds2 <- random_weights[,inds2] > 0
        inds3 <- which(inds & inds2)
        modn_gof <- rbind(modn_gof, random_gof %>%
                                  filter(Prediction %in% inds3) %>%
                                  select(Station, RMSE, `PBIAS %`, 
                                         NSE, R2, KGE) %>%
                                  mutate(n = nrow(.)) %>%
                                  group_by(Station) %>%
                                  summarise_all(mean) %>%
                                  add_column(model = i, 
                                             nmod = j,
                                             .before=1))
    }
}
arrange(modn_gof, RMSE)



# add individual ensemble member goodness-of-fit to modn_gof
temp <- length_gof %>%
    select(Prediction, Station, RMSE, `PBIAS %`, NSE, R2, KGE) %>%
    group_by(Prediction, Station) %>%
    summarise_all(mean) %>%
    mutate(nmod = 1, n = 10) %>%
    rename(model = Prediction)
modn_gof <- bind_rows(modn_gof, temp)
modn_gof$Station <- factor(modn_gof$Station, levels=statorder)



## make Figure 5 
labels <- bind_rows(modn_gof %>% 
                        filter(nmod == 1) %>%
                        rename(Prediction = model) %>%
                        select(Prediction, nmod, Station, KGE) %>%
                        mutate(Prediction = rep(LETTERS[1:15], each=10),
                               shape = "Ensemble mean",
                               size=3), # to get the round shape
                    discharge_gof %>% 
                        filter(Prediction %in% c("GRADES", "GLOFAS")) %>%
                        mutate(nmod = 15,
                               shape = Prediction,
                               # size so that pdf export does not make dots small
                               size = 1.5) %>%
                        select(Prediction, nmod, Station, KGE, shape, size),
                    pred_ensemble_gof %>%
                        filter(Prediction == "mean") %>%
                        mutate(nmod = 15,
                               Prediction = "Ensemble mean",
                               shape = Prediction,
                               size = 3) %>%
                        select(Prediction, nmod, Station, KGE, shape, size))
labels$Station <- factor(labels$Station, levels=statorder)

breaks_y <- seq(-1,1,by=0.2)
breaks_x <- 1:15
modn_gof %>%
    filter(!model %in% c("GRUN", "LORA")) %>%
    mutate(Prediction = model,
           model = word(Prediction,1,sep="_"),
           forcing = word(Prediction, 2, sep="_")) %>%
    group_by(nmod, Station) %>%
    mutate(KGEmin = min(KGE),
           KGEmax = max(KGE)) %>%
    ggplot() + 
    geom_ribbon(aes(x=nmod, ymin=KGEmin, ymax=KGEmax), fill="grey90") +
    geom_line(aes(x=nmod, y=KGE, color=model, linetype=forcing), 
              size=0.5) +
    # geom_smooth(aes(x=nmod, y=KGE, color=model, linetype=forcing),
    #             se=FALSE) +
    #geom_jitter(aes(x=nmod, y=KGE), color = "grey40", width = 0.1) +
    geom_point(data = labels, aes(x=nmod, y=KGE, shape = shape), 
               size=1) +
    geom_hline(yintercept = 0.5) +
    ggrepel::geom_text_repel(data = labels,
                             aes(x=nmod, y=KGE, label=Prediction),
                             size = 3) +
    facet_wrap(~Station, scales="free_y", drop = FALSE) +
    theme_bw() +
    ylim(c(1,1)) +
    theme(text=element_text(size=12)) +
    scale_x_continuous(breaks = breaks_x) +
    scale_y_continuous(breaks = breaks_y) +
    scale_color_manual(values = rep('grey50', 6)) +
    labs(y = "Mean KGE",
         x = "Number of models in random combination")








# ------------------------------------------------------------------------------
# Optimal combinations and their GOF


# define function to do combinations and stores optimal weights, gof stats,
# and timesteps which were used for the training

do_comb <- function(HS, n, optim_method, combination, 
                    sampling, train, name,
                    bias_correction) {
    pb <- txtProgressBar(min=0,max=n, style = 3)
    
    wgt <- list()
    gof <- list()
    trained <- list()
    
    for (i in 1:n) {
        tempcomb <- optimise_point(HS,
                                   optim_method = optim_method,
                                   combination = combination,
                                   sampling = sampling,
                                   train = train)
        
        if(combination %in% c("ts", "ann", "timeseries", "annual")) {
            
            weights <- lapply(seq_along(tempcomb),
                              function(x) { 
                                  temp <- tempcomb[[x]]$Weights 
                                  temp <- data.frame(Station = 
                                                         names(tempcomb)[[x]],
                                                     Comb = paste0(name, "_",i),
                                                     t(temp))
                                  temp <- as_tibble(temp)
                                  colnames(temp) <- gsub("\\.", "-", 
                                                         colnames(temp))
                                  return(temp)})
            
            
            
            weights <- suppressWarnings(do.call(bind_rows, weights))
            
            wgt[[ paste0(name, "_", i) ]] <- weights
            
            gofs <- lapply(seq_along(tempcomb),
                           function(x)  tempcomb[[x]]$Goodness_of_fit %>% 
                               t() %>%
                               as.data.frame() %>%
                               tibble::rownames_to_column() %>% 
                               as_tibble() %>%
                               rename(Period = rowname) %>%
                               tibble::add_column(Station = 
                                                      names(tempcomb)[[x]],
                                                  Comb = paste0(name, "_",i),
                                                  .before=1))
            
            gofs <- do.call(bind_rows, gofs)
            
            gof[[ paste0(name, "_", i) ]] <- gofs
            
            trains <- lapply(seq_along(tempcomb),
                             function(x) tempcomb[[x]]$Optimized_ts$Date[
                                 tempcomb[[x]]$Optimized_ts$Train_or_test == "Train"] )
            
            
            names(trains) <- names(tempcomb)
            
            trained[[ paste0(name, "_", i) ]] <- trains
            
            setTxtProgressBar(pb, i)
        } else if (combination %in% c("mon", "monthly")) {
            
            weights <- lapply(seq_along(tempcomb),
                              function(x) { 
                                  temp <- tempcomb[[x]]$Weights 
                                  temp <- data.frame(Station = 
                                                         names(tempcomb)[[x]],
                                                     Comb = paste0(name, "_",i),
                                                     temp)
                                  temp <- as_tibble(temp)
                                  colnames(temp) <- gsub("\\.", "-", 
                                                         colnames(temp))
                                  return(temp)})
            
            
            
            weights <- suppressWarnings(do.call(bind_rows, weights))
            
            wgt[[ paste0(name, "_", i) ]] <- weights
            
            gofs <- lapply(seq_along(tempcomb),
                           function(x)tempcomb[[x]]$Goodness_of_fit$`Entire timeseries` %>% 
                               t() %>%
                               as.data.frame() %>%
                               tibble::rownames_to_column() %>% 
                               as_tibble() %>%
                               rename(Period = rowname) %>%
                               tibble::add_column(Station = names(tempcomb)[[x]],
                                                  Comb = paste0(name, "_",i),
                                                  .before=1))
            
            gofs <- do.call(bind_rows, gofs)
            
            gof[[ paste0(name, "_", i) ]] <- gofs
            
            
            trains <- lapply(seq_along(tempcomb),
                             function(x) tempcomb[[x]]$Optimized_ts$Date[
                                 tempcomb[[x]]$Optimized_ts$Train_or_test == "Train"] )
            
            
            names(trains) <- names(tempcomb)
            
            trained[[ paste0(name, "_", i) ]] <- trains
            
            setTxtProgressBar(pb, i)
        } else stop("wrong combination")
    }
    close(pb)
    
    output <- list(weights = suppressWarnings(do.call(bind_rows, wgt)),
                   gofs = suppressWarnings(do.call(bind_rows, gof)),
                   train = trained)
    return(output)
}


# set seed for reproducibility

####################################
# OPTIMISE USING TIMESERIES
####################################

system.time({
    set.seed(2019)
    
    down_ts_comb <- do_comb(HS_stations, 
                            n = 100,
                            optim_method = "CLS",
                            combination = "ts",
                            sampling = "random",
                            train = 0.5,
                            name = "down_ts",
                            bias_correction = FALSE)
    
}) # 54 sec

####################################
# OPTIMISE USING MONTHLY COMBINATIONS
####################################
system.time({
    
    down_mon_comb <- do_comb(HS_stations, 
                             n = 50,
                             optim_method = "CLS",
                             combination = "mon",
                             sampling = "random",
                             train = 0.5,
                             name = "down_mon",
                             bias_correction = FALSE)     
}) # 303 sec

####################################
# OPTIMISE USING FULL CALENDAR YEARS
####################################
system.time({
    
    down_ann_comb <- do_comb(HS_stations, 
                             n = 50,
                             optim_method = "CLS",
                             combination = "ann",
                             sampling = "random",
                             train = 0.5,
                             name = "down_ann",
                             bias_correction = FALSE)     
}) # 21 sec





#-------------------------------------------------------------------------------
# Gather optimal weights and combine flows using them

system.time({
    temp <- lapply(seq_along(1:nrow(down_ts_comb$weights)), 
                   function(x) {temp <- down_ts_comb$weights[x,3:17] %>% 
                       unlist()
                   })
    names(temp) <- paste(down_ts_comb$weights$Station,
                         down_ts_comb$weights$Comb,sep='_')
    HSf_ts <- combine_runoff(HS_stations, weights = temp, 
                             intercept = as.list(rep(0,length(temp))),
                             bias = as.list(rep(0,length(temp))),
                             drop=TRUE)
    fitted_gofs <- flow_gof(HSf_ts, verbose=TRUE)
    
    ##
    temp <- lapply(seq_along(1:nrow(down_ann_comb$weights)), 
                   function(x) {temp <- down_ann_comb$weights[x,3:17] %>% 
                       unlist()
                   })
    names(temp) <- paste(down_ann_comb$weights$Station,
                         down_ann_comb$weights$Comb,sep='_')
    HSf_ann <- combine_runoff(HS_stations, weights = temp, 
                              intercept = as.list(rep(0,length(temp))),
                              bias = as.list(rep(0,length(temp))),
                              drop=TRUE)
    fitted_gofs <- bind_rows(fitted_gofs, flow_gof(HSf_ann, verbose=TRUE))
    
    ##
    seq <- seq(1, nrow(down_mon_comb$weights)-11, by=12)
    temp <- lapply(seq_along(seq), function(x) {
        rows <- x:(x+11)
        temp <- down_mon_comb$weights[rows,4:18] 
        temp[is.na(temp)] <- 0
        temp <- as.matrix(temp)
        return(temp)
    })
    names(temp) <- paste(down_mon_comb$weights$Station[seq],
                         down_mon_comb$weights$Comb[seq],
                         sep='_')
    HSf_mon <- combine_runoff(HS_stations, weights = temp, 
                              intercept = rep(list(rep(0,12)), length(temp)),
                              bias = rep(list(rep(0,12)), length(temp)),
                              drop=TRUE, monthly=TRUE)
    fitted_gofs <- bind_rows(fitted_gofs,flow_gof(HSf_mon,verbose=TRUE))
}) # 429 sec





# ------------------------------------------------------------------------------
# get data for Table 4

# change assess to any of
# * discharge_gof to get GHM data
# * to get area_to_line downscaling data
# * random_gof to get random combination data
assess <- length_gof
temp <- assess %>% 
    select(Prediction, Station, R2, NSE, `PBIAS %`, KGE) %>%
    filter(!Prediction %in% c("GRADES", "GLOFAS")) # remove benchmarks 
testi <- apply(temp[,-c(1:2)], 1, function(x) {
    t1 <- x[1] >0.6
    t2 <- x[2] >0.5
    t3 <- abs(x[3]) < 15
    t4 <- x[4] >0.5
    t5 <- sum(t1,t2,t3,t4)
    return(data.frame(R2_ = t1, NSE_ = t2, PBIAS_ = t3, KGE = t4, sum = t5))
})
testi <- do.call("rbind", testi) %>% as_tibble() 
testi %>% summary # TRUE means the criteria of Moriarsi is met


# change assess to any of 
# * down_ts_comb$gofs for test period gof for timeseries combinations
# * down_ann_comb$gofs for test period gof for annual combinations
# * down_mon_comb$gofs for test period gof for monthly combinations
assess <- down_ts_comb$gofs
temp <- assess %>% 
    select(Comb, Station, Period, R2, NSE, `PBIAS %`, KGE) %>%
    filter(Period == "Test")
testi <- apply(temp[,-c(1:3)], 1, function(x) {
    t1 <- x[1] >0.6
    t2 <- x[2] >0.5
    t3 <- abs(x[3]) < 15
    t4 <- x[4] >0.5
    t5 <- sum(t1,t2,t3,t4)
    return(data.frame(R2_ = t1, NSE_ = t2, PBIAS_ = t3, KGE = t4, sum = t5))
})
testi <- do.call("rbind", testi) %>% as_tibble() %>%
    add_column(Prediction = temp$Comb,
               Station = temp$Station,
               .before=1)
testi %>% summary # TRUE means the criteria of Moriarsi is met



# change assess to any of
# * 'ts' for regionalized timeseries combination performance
# * 'ann' for regionalized annual combination performance
# * 'mon' for regionalized monthly combination performance
assess <- 'ts'
temp <- fitted_gofs %>% 
    select(Prediction, Station, R2, NSE, `PBIAS %`, KGE) %>%
    mutate(fit_stat = word(Prediction,1,sep="_"),
           fit_river = word(fit_stat, 1),
           stat_river = word(Station, 1),
           Type = word(Prediction, 3, sep="_")) %>%
    filter(Station != fit_stat,
           Type == assess,
           fit_river == stat_river)
testi <- apply(temp[,c(3:6)], 1, function(x) {
    t1 <- x[1] >0.6
    t2 <- x[2] >0.5
    t3 <- abs(x[3]) < 15
    t4 <- x[4] >0.5
    t5 <- sum(t1,t2,t3,t4)
    return(data.frame(R2_ = t1, NSE_ = t2, PBIAS_ = t3, KGE = t4, sum = t5))
})
testi <- do.call("rbind", testi) %>% as_tibble() %>%
    add_column(Prediction = temp$Prediction,
               Station = temp$Station,
               .before=1)
testi %>% summary # TRUE means the criteria of Moriarsi is met





# get data for table 4 for the optimized ensembles

# get the ensembles for the stations they are trained at
ts_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_ts[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    ind <- which(grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    ts_ensemble_gof <- rbind(ts_ensemble_gof,
                                 flow_gof(hs))
    
}
ts_ensemble_gof <- ts_ensemble_gof[complete.cases(ts_ensemble_gof),]


ann_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_ann[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    ind <- which(grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    ann_ensemble_gof <- rbind(ann_ensemble_gof,
                                  flow_gof(hs))
    
}
ann_ensemble_gof <- ann_ensemble_gof[complete.cases(ann_ensemble_gof),]

mon_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_mon[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    ind <- which(grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    mon_ensemble_gof <- rbind(mon_ensemble_gof,
                                  flow_gof(hs))
    
}
mon_ensemble_gof <- mon_ensemble_gof[complete.cases(mon_ensemble_gof),]

random_ensemble_gof <- ensemble_summary(HS_random, 
                                        summarise_over_timeseries = FALSE,
                                        funs = c("mean", "median"), 
                                        drop = TRUE)

# change assess to one of
# * random_ensemble_gof for mean of random combinations
# * dis_ensemble_gof for GHM prediction means
# * pred_ensemble_gof for area-to-line downscaled mean
# * ts_ensemble_gof
# * ann_ensemble_gof
# * mon_ensemble_gof
assess <- random_ensemble_gof
temp <- assess %>% 
    select(Prediction, Station, R2, NSE, `PBIAS %`, KGE) %>%
    filter(Prediction == "mean")
testi <- apply(temp[,-c(1:2)], 1, function(x) {
    t1 <- x[1] >0.6
    t2 <- x[2] >0.5
    t3 <- abs(x[3]) < 15
    t4 <- x[4] >0.5
    t5 <- sum(t1,t2,t3,t4)
    return(data.frame(R2_ = t1, NSE_ = t2, PBIAS_ = t3, KGE = t4, sum = t5))
})
testi <- do.call("rbind", testi) %>% as_tibble() %>%
    add_column(Prediction = temp$Prediction,
               Station = temp$Station,
               .before=1)
testi %>% summary # TRUE means the criteria of Moriarsi is met





# get data for table 4 for the regionalized ensembles

# get the ensembles for the stations they are NOT trained at
regts_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_ts[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    inds <- grepl(word(stat,1),names(hs$discharge_ts[[1]]))
    inds2 <- !grepl(stat, names(hs$discharge_ts[[1]]))
    ind <- which(inds & inds2)
    #ind <- which(!grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    regts_ensemble_gof <- rbind(regts_ensemble_gof,
                                    flow_gof(hs))
    
}
regts_ensemble_gof <- regts_ensemble_gof[complete.cases(regts_ensemble_gof),]


regann_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_ann[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    inds <- grepl(word(stat,1),names(hs$discharge_ts[[1]]))
    inds2 <- !grepl(stat, names(hs$discharge_ts[[1]]))
    ind <- which(inds & inds2)
    #ind <- which(!grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    regann_ensemble_gof <- rbind(regann_ensemble_gof,
                                     flow_gof(hs))
    
}
regann_ensemble_gof <- regann_ensemble_gof[complete.cases(regann_ensemble_gof),]

regmon_ensemble_gof <- vector()
for(i in 1:10) {
    hs <- HSf_mon[i,] %>% select(-runoff_ts)
    stat <- hs$observation_station
    inds <- grepl(word(stat,1),names(hs$discharge_ts[[1]]))
    inds2 <- !grepl(stat, names(hs$discharge_ts[[1]]))
    ind <- which(inds & inds2)
    #ind <- which(!grepl(stat, names(hs$discharge_ts[[1]])))
    hs$discharge_ts[[1]] <- hs$discharge_ts[[1]][,c(1,ind)]
    
    hs <- ensemble_summary(hs, summarise_over_timeseries = FALSE,
                           funs = c("mean", "median"), drop = TRUE)
    
    regmon_ensemble_gof <- rbind(regmon_ensemble_gof,
                                     flow_gof(hs))
    
}
regmon_ensemble_gof <- regmon_ensemble_gof[complete.cases(regmon_ensemble_gof),]


# change assess to one of
# * regts_ensemble_gof
# * regann_ensemble_gof
# * regmon_ensemble_gof
assess <- regts_ensemble_gof
temp <- assess %>% 
    select(Prediction, Station, R2, NSE, `PBIAS %`, KGE) %>%
    filter(Prediction == "mean")
testi <- apply(temp[,-c(1:2)], 1, function(x) {
    t1 <- x[1] >0.6
    t2 <- x[2] >0.5
    t3 <- abs(x[3]) < 15
    t4 <- x[4] >0.5
    t5 <- sum(t1,t2,t3,t4)
    return(data.frame(R2_ = t1, NSE_ = t2, PBIAS_ = t3, KGE = t4, sum = t5))
})
testi <- do.call("rbind", testi) %>% as_tibble() %>%
    add_column(Prediction = temp$Prediction,
               Station = temp$Station,
               .before=1)
testi %>% summary







# ------------------------------------------------------------------------------
# create figure 4
plotgofs <- bind_rows(discharge_gof %>% 
                          add_column(Type = "A Global model") %>% 
                          filter(!Prediction %in% c("GLOFAS", "GRADES")),
                      length_gof %>% 
                          add_column(Type = "B Downscaled model"))
temp <- discharge_gof %>%
    filter(!Prediction %in% c("GRADES", "GLOFAS")) %>% 
    group_by(Prediction) %>%
    summarise_all(mean) %>%
    #select(Prediction, ME, `PBIAS %`, NSE, KGE, R2) %>%
    arrange(Prediction)

temp2 <- length_gof %>%
    group_by(Prediction) %>%
    summarise_all(mean) %>%
    #select(Prediction, ME, `PBIAS %`, NSE, KGE, R2) %>%
    arrange(Prediction)
table <- data.frame(Prediction = temp$Prediction,
           RMSE_GM = temp$RMSE,
           RMSE_D = temp2$RMSE,
           PBIAS_GM = temp$`PBIAS %`,
           PBIAS_D = temp2$`PBIAS %`,
           NSE_GM = temp$NSE,
           NSE_D = temp2$NSE,
           KGE_GM = temp$KGE,
           KGE_D = temp2$KGE,
           R2_GM = temp$R2,
           R2_D = temp2$R2) %>%
    gather(var,val,-Prediction) %>%
    mutate(variable = word(var,1,sep="_"),
           model = word(var,2,sep="_")) %>%
    select(-var)%>%
    spread(model,val)

KGE <- ggplot(table %>% filter(variable == "KGE")) +
    geom_segment(aes(x= Prediction, xend=Prediction,
                     y= GM, yend = D),
                 arrow = arrow(length=unit(0.02, "npc"), type="closed")) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = GM, yend = GM)) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = D, yend = D),
                 color = 'red') +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(title = "KGE",
         y = "Average Kling-Gupta Efficiency",
         x = "Prediction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))

NSE <- ggplot(table %>% filter(variable == "NSE")) +
    geom_segment(aes(x= Prediction, xend=Prediction,
                     y= GM, yend = D),
                 arrow = arrow(length=unit(0.02, "npc"), type="closed")) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = GM, yend = GM)) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = D, yend = D),
                 color = 'red') +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(title = "NSE",
         y = "Average Nash-Sutcliffe Efficiency",
         x = "Prediction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))

PBIAS <- ggplot(table %>% filter(variable == "PBIAS")) +
    geom_segment(aes(x= Prediction, xend=Prediction,
                     y= GM, yend = D),
                 arrow = arrow(length=unit(0.02, "npc"), type="closed")) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = GM, yend = GM)) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = D, yend = D),
                 color = 'red') +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "PBIAS",
         y = "Average Bias %",
         x = "Prediction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))

R2 <- ggplot(table %>% filter(variable == "R2")) +
    geom_segment(aes(x= Prediction, xend=Prediction,
                     y= GM, yend = D),
                 arrow = arrow(length=unit(0.02, "npc"), type="closed")) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = GM, yend = GM)) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = D, yend = D),
                 color = 'red') +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(title = "R2",
         y = "Average Coefficient of Determination",
         x = "Prediction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))

RMSE <- ggplot(table %>% filter(variable == "RMSE")) +
    geom_segment(aes(x= Prediction, xend=Prediction,
                     y= GM, yend = D),
                 arrow = arrow(length=unit(0.02, "npc"), type="closed")) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = GM, yend = GM)) +
    geom_segment(aes(x = as.numeric(Prediction)-0.2,
                     xend = as.numeric(Prediction)+0.2,
                     y = D, yend = D),
                 color = 'red') +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "RMSE",
         y = "Average Root-Mean-Square-of-Error",
         x = "Prediction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))

(RMSE | PBIAS | NSE) / (plot_spacer() | R2 | KGE)
#ggsave(filename = "comparison_global_downscaled.pdf")







# ------------------------------------------------------------------------------
# create figure 3

gof_combined <- bind_rows(
    discharge_gof %>% add_column(Type = "A Global model") %>%
        filter(!Prediction %in% c("GRADES", "GLOFAS")),
    length_gof %>% add_column(Type = "B Downscaled model"),
    random_gof %>% add_column(Type = "C Random combination"),
    #dis_ts_comb$gofs %>% add_column(Type = "B Discharge_ts"),
    down_ts_comb$gofs %>% add_column(Type = "D Timeseries combination"),
    #dis_ann_comb$gofs %>% add_column(Type = "C DEischarge_ann"),
    down_ann_comb$gofs %>% add_column(Type = "E Annual combinations"),
    #dis_mon_comb$gofs %>% add_column(Type = "D Discharge_mon"),
    down_mon_comb$gofs %>% add_column(Type = "F Monthly combinations"),
    fitted_gofs %>%
        mutate(fit_stat = word(Prediction,1,sep="_"),
               Type = word(Prediction, 3, sep="_"),
               Type = case_when(
                   Type == "ts" ~ "G Regionalized timeseries",
                   Type == "ann" ~ "H Regionalized annual",
                   Type == "mon" ~ "I Regionalized monthly"
               )) %>%
        filter(Station != fit_stat)  %>%
        mutate(eval_riv = stringr::word(Station, 1),
               w_riv = stringr::word(Prediction, 1)) %>%
        filter(eval_riv == w_riv) %>%
        select(-eval_riv, -w_riv))
gof_combined$Station <- factor(gof_combined$Station,
                               levels = statorder)

point_labels <- discharge_gof %>% 
    filter(Prediction %in% c("GRADES", "GLOFAS")) %>%
    mutate(Type = "A Global model",
           size = 2) %>%
    select(Type, Prediction, Station, KGE, size) %>%
    bind_rows(pred_ensemble_gof %>% 
                  filter(Prediction %in% c("mean")) %>%
                  mutate(Type = "B Downscaled model",
                         Prediction = "Ensemble mean",
                         size = 3.5) %>%
                  select(Type, Prediction, Station, KGE, size))
point_labels$Station <- factor(point_labels$Station,
                               levels = statorder)

gof_combined %>%
    filter(Period == "Test" || is.na(Period)) %>%
    select(Station, Type, Comb, ME, `PBIAS %`, NSE, KGE, R2) %>%
    ggplot() +
    geom_violin(aes(Type,KGE, fill=Type),
                position='dodge',
                trim=TRUE,
                draw_quantiles = TRUE,
                scale="width",
                adjust=0.75) +
    # geom_jitter(aes(Type,KGE),
    #             size=0.05,
    #             width=0.5) +
    facet_wrap(~Station, drop=FALSE) +
    coord_cartesian(ylim=c(0, 1)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    geom_hline(yintercept = 0.75, linetype = 2, size = 0.5) +
    geom_hline(yintercept = 0.5, linetype = 2, size = 0.5) +
    geom_point(data = point_labels,
               aes(Type, KGE, shape=Prediction), 
               size=point_labels$size) +
    ggrepel::geom_text_repel(data = point_labels %>% 
                                 filter(KGE < -0.05) %>%
                                 mutate(label = paste(Prediction, KGE)),
                             aes(Type, KGE, label=label)) +
    scale_fill_manual(values=colors1[c(1:2,11,8:10,4:6)]) +
    theme(axis.text.x = element_blank(),
          text = element_text(size=12))

#ggsave("fig_optimized.pdf", height = 200, width=300, unit="mm")
