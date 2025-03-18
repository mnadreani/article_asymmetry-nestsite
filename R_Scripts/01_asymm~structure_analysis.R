# Script for the analysis of the manuscript:
# Choice or constraint? Interaction between the nest and external structures predicts nest's architecture

# Author: Nicolas M. Adreani
# Goal: Does having a structure in contact with the nest have an effect on it's asymmetric architecture?


##### 01 Require packages #####
sapply(c('dplyr','maps', 'sf', 'mapdata',
         'ggplot2', 'arm', 'gridExtra','DHARMa'),
       function(x) suppressPackageStartupMessages(
         require(x ,character.only = TRUE, quietly = TRUE)))


#### 02 Import data and organize variables

  db <- read.csv(file = "Data/database_asymm_structures.csv")
  
  db$ASYMM_SIDE <- as.factor(db$ASYMM_SIDE) 
# convert to factor
  db$ASYMM_SIDE <- as.numeric(db$ASYMM_SIDE) 
# convert to numeric: 1 corresponds to 'left' nest asymmetry and 2 correspond to 'right' nest asymmetry.
  db$ASYMM_SIDE[which(db$ASYMM_SIDE==2)] <- 0 
# replace 2 by 0, to have a binomial variable (0-1). for interpretation: 0=righ nest asymmetry \ 1=left nest asymmetry
  db$STRUCTURE_SIDE <- as.factor(db$STRUCTURE_SIDE) 
# convert to factor
  db$match_mismatch <- factor(db$match_mismatch, 
                                      levels = c("match", "mismatch")) 
  # convert to factor
  
  db$match_mismatch <- as.numeric(db$match_mismatch)
  db$match_mismatch[which(db$match_mismatch==2)] <- 0 #replace 2 for 0 so it is a binomial variable
  db$ASYMM_SIDE[which(db$ASYMM_SIDE==0)] <- "right" 
  db$ASYMM_SIDE[which(db$ASYMM_SIDE==1)] <- "left"
  db$ASYMM_SIDE <- as.factor( db$ASYMM_SIDE)
  db <- subset(db, db$RELIABILITY=="2") # only keep cases of high reliavility
  db$REP <- as.factor(db$REP) 
  db$REP <- droplevels( db$REP)
  
  
  # Sample sizes
  table(db$ASYMM_SIDE,db$REP)
  
  
  # number of simulations
  nsim=10000
  
  
  ############
  # 1) Model & Summary
  # 2) posterior distributions & infor for statistical tables
  # 3) Preeliminary plot
  # 4) final plot of figure
  ############
  
  
  #### 1. Model  ####
  mod_04 <- glm(match_mismatch ~ ASYMM_SIDE * REP, 
                family = binomial(link = "logit"),
                data = db)
  summary((mod_04))
  
  # Residual Check - ( package DHARMa)
  simulationOutput0 <- simulateResiduals(fittedModel = mod_04, plot = F)
  plot(simulationOutput0)
  
  ### Simulations and posterior distributions MODEL 4 ####
  # ASYMM_SIDE * REP
  
  bsim_04 <- sim(mod_04, n.sim =nsim)
  
  newdat_04 <- expand.grid(ASYMM_SIDE = levels(db$ASYMM_SIDE),
                           REP = levels(db$REP))
  
  # NOTA NICO: AGREGUE INTERACCION
  Xmat_04 <- model.matrix(~ASYMM_SIDE * REP, newdat_04)
  fitmat_04 <- matrix(nrow = nrow(newdat_04), ncol=nsim)
  for(i in 1:nsim) fitmat_04[,i] <- plogis(Xmat_04 %*% bsim_04@coef[i,])
  newdat_04$lwr <-apply(fitmat_04, 1, quantile, prob=0.025) # upper interval
  newdat_04$upr <-apply(fitmat_04, 1, quantile, prob=0.975) # lower interval
  newdat_04$fit <-plogis(Xmat_04 %*% coef(mod_04))   #fitter value
  
  # 2. Estimates and 95% CrI for Table 1
  ##To get the mean estimate of the fixed effects parameters apply(bsim@fixef, 2, mean)
  apply((bsim_04@coef), 2, mean)
  ##To get the 95% credible Interval estimate of the fixed effects parameters apply(bsim_mod_testo@coef, 2, quantile, prob=c(0.025, 0.975))
  apply((bsim_04@coef), 2, quantile, prob=c(0.025, 0.975))
  
  
  ### 3. Raw data plot.  ####  
  
  raw_plot_04 <- ggplot(data = db, aes(y =match_mismatch, x = ASYMM_SIDE,
                                               group=REP, shape=REP))+
    geom_jitter(data = db, size = 0.7, alpha = 0.7,color= "#0072B2",
                 position = position_jitterdodge(jitter.width =0.4, jitter.height = 0.09, dodge.width = 0.9)) +
     scale_x_discrete(breaks = c("left", "right"), labels = c("Left", "Right")) +
    # scale_y_continuous(breaks = c(0,0.25,0.5,0.75, 1))+
    scale_y_continuous(breaks = c(0, 0.5,0.75,1))+
    labs(x = "Nest Asymmetry",y = "P(Match)")+
    scale_shape_manual(values = c(16, 1))+
    theme_classic()+
    facet_grid(.~REP)+
    geom_hline(yintercept = c(0.75), linetype = "dotted", color="grey69")+  # Ensure this is placed correctly
    geom_hline(yintercept = c(0.5), color="grey69")  # Ensure this is placed correctly
  
  raw_plot_04

  # NICO: I change the variable name 'fit' into 'ASYMM_SIDE' so that the aes() defined in raw_plot_04 
  # can be found in the final plot and when merging the posterior probabilities with the raw data (newdat and db_vik_sub) there are no 
  # errors.
  
  names(newdat_04)[5] <- "match_mismatch"
  
  #### 4. Final plot for Figure 2B
  
  final_plot_04 <- raw_plot_04+ 
    geom_point(data = newdat_04, 
               aes(x = ASYMM_SIDE, y = match_mismatch), 
               color = "black", size = 3,
               position = position_dodge(width = 0.8)) +
    geom_errorbar(data = newdat_04,
                  aes(x=ASYMM_SIDE, ymin= lwr, ymax= upr), 
                  color="black", width= 0,linewidth = 1,
                  position = position_dodge(width = 0.8))
  
  final_plot_04
  
  
  
  
  ### Map Figure 2A ####
  # to visualize the distribution of matching-mismatching nests we generated a map
  # using the packages sf, maps, and mapdata
  
  # Database of the app to extract GPS coordinates
  db_app <- read.csv(file = "Data/db_app.csv", header = T)
  db_app <- db_app[,-c(2:4,7,9:11,13)]
  db_app <- subset(db_app, !is.na(db_app$LONGITUDE))
  db_app <- subset(db_app, !is.na(db_app$LATITUDE))
  
  #merge asymmetry database with hornero app pictures database to find the nests with GPS coordinates
  db_gps <- merge(db, db_app,by = "ID")
  
  db_gps_map <- st_as_sf(db_gps, coords = c("LONGITUDE", "LATITUDE"), 
                             crs = 4326)
  world_map <- map_data("world")
  db_gps_map$match_mismatch <- factor(db_gps_map$match_mismatch, 
                                          levels = c("match", "mismatch"))
  
  # Map for figure 2A
  map <- ggplot(data = db_gps_map, aes(colour = match_mismatch)) +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "white", color = "#0C0E1D", linewidth=0.2) +
    geom_sf(data = db_gps_map,  size = 1.2, alpha=0.5) +
    coord_sf(xlim = c(-80, -45), ylim = c(-55, -10), expand = FALSE) +
    theme_minimal() +
    scale_color_manual(labels= c("match","mismatch"),
                       values = c("mismatch" = "#009E73","match" = "#CC79A7"))+
    # guides(color='none')+
    labs(
      x = "",
      y = ""
    ) +
    theme_classic()
  
    map
  

    
  