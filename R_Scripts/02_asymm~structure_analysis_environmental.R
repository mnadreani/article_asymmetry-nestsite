

##### 01 Require packages #####
sapply(c('dplyr','maps', 'sf', 'mapdata',
         'ggplot2', 'arm', 'gridExtra','DHARMa'),
       function(x) suppressPackageStartupMessages(
         require(x ,character.only = TRUE, quietly = TRUE)))





#### Figure MAP








# GLMER w/User Random
############
# 1) Model & Summary
# 2) posterior distributions & infor for statistical tables
# 3) Preeliminary plot
# 4) final plot of figure
############

final_db <- read.csv(file = "Data/final_db_env.csv",header = T)

final_db$match_mismatch <- as.factor(final_db$match_mismatch)
final_db$match_mismatch <- as.numeric(final_db$match_mismatch)
final_db$match_mismatch[which(final_db$match_mismatch==2)] <- 0 #replace 2 for 0 so it is a binomial variable
# match = 1
# mismatch = 0
names(final_db)[27] <- "temp"
names(final_db)[28] <- "rain"

final_db = final_db[!is.na(final_db$LATITUDE), ]
final_db = final_db[!is.na(final_db$rain), ]
final_db = final_db[!is.na(final_db$altitude), ]
final_db = final_db[!is.na(final_db$temp), ]

final_db$s_temp <- scale(final_db$temp)
final_db$s_lat <- scale(final_db$LATITUDE)
final_db$s_long <- scale(final_db$LONGITUDE)
final_db$s_rain <- scale(final_db$rain)
final_db$s_alt <- scale(final_db$altitude)

final_db$REP <- as.factor(final_db$REP)

#### 2. Final Model
mod_07 <- glmer(match_mismatch ~ REP + 
                  s_long   + s_lat +
                  s_alt + s_temp + s_rain + (1|user), 
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa", 
                                       optCtrl = list(maxfun = 1e5)),
                data = final_db)
summary(mod_07)

# Residual Check - ( package DHARMa)
simulationOutput0 <- simulateResiduals(fittedModel = mod_07, plot = F)
plot(simulationOutput0)

### Simulations and posterior distributions MODEL 4 ####
# ASYMM_SIDE * REP
nsim<- 10000
bsim_07 <- sim(mod_07, n.sim =nsim)
newdat_07 <- expand.grid(REP = levels(final_db$REP),
                         s_lat = mean(final_db$s_lat),
                         s_long = mean(final_db$s_long),
                         s_temp = mean(final_db$s_temp),
                         s_rain = mean(final_db$s_rain),
                         s_alt = mean(final_db$s_alt))

Xmat_07 <- model.matrix(~ REP + 
                          s_long   + s_lat +
                          s_alt + s_temp + s_rain, newdat_07)

fitmat_07 <- matrix(nrow = nrow(newdat_07), ncol=nsim)
for(i in 1:nsim) fitmat_07[,i] <- plogis(Xmat_07 %*% bsim_07@fixef[i,])
newdat_07$lwr <-apply(fitmat_07, 1, quantile, prob=0.025) # upper interval
newdat_07$upr <-apply(fitmat_07, 1, quantile, prob=0.975) # lower interval
newdat_07$fit <-plogis(Xmat_07 %*% fixef(mod_07))   #fitter value

# 2. Estimates and 95% CrI for Table 1
##To get the mean estimate of the fixed effects parameters apply(bsim@fixef, 2, mean)
apply((bsim_07@fixef), 2, mean)
##To get the 95% credible Interval estimate of the fixed effects parameters apply(bsim_mod_testo@coef, 2, quantile, prob=c(0.025, 0.975))
apply((bsim_07@fixef), 2, quantile, prob=c(0.025, 0.975))


### 3. Raw data plot.  ####  

raw_plot_07 <- ggplot(data = final_db, aes(y =match_mismatch, x = REP,
                                           shape=REP))+
  geom_jitter(data = final_db, size = 0.7, alpha = 0.5,
              position = position_jitterdodge(jitter.width =0.8, jitter.height = 0.08, dodge.width = 0.9)) +
  scale_x_discrete(breaks = c("left", "right"), labels = c("Left", "Right")) +
  # scale_y_continuous(breaks = c(0,0.25,0.5,0.75, 1))+
  scale_y_continuous(breaks = c(0, 0.5,0.75,1))+
  scale_shape_manual(values = c(16, 1))+
  theme_classic()+
  geom_hline(yintercept = c(0.75), linetype = "dotted", color="grey69")+  # Ensure this is placed correctly
  geom_hline(yintercept = c(0.5), color="grey69")  # Ensure this is placed correctly

raw_plot_07



# NICO: I change the variable name 'fit' into 'match_mismatch' so that the aes() defined in raw_plot_05 
# can be found in the final plot and when merging the posterior probabilities with the raw data (newdat and db_vik_sub) there are no 
# errors.

names(newdat_07)[9] <- "match_mismatch"

#### 4. Final plot for Figure 2B

final_plot_07 <- raw_plot_07+ 
  geom_point(data = newdat_07, 
             aes(x = REP, y = match_mismatch), 
             color = "black", size = 3,
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = newdat_07,
                aes(x=REP, ymin= lwr, ymax= upr), 
                color="black", width= 0,linewidth = 1,
                position = position_dodge(width = 0.8))

final_plot_07


#### 2. Supplementary model w asymmetry

mod_08 <- glmer(match_mismatch ~ REP *ASYMM_SIDE +
                  s_long   + s_lat +
                  s_alt + s_temp + s_rain + (1|user), 
                family = binomial(link = "logit"),
                control = glmerControl(optimizer = "bobyqa", 
                                       optCtrl = list(maxfun = 1e5)),
                data = final_db)
summary(mod_08)

# Residual Check - ( package DHARMa)
simulationOutput0 <- simulateResiduals(fittedModel = mod_08, plot = F)
plot(simulationOutput0)

### Simulations and posterior distributions MODEL 4 ####
# ASYMM_SIDE * REP
nsim<- 10000
bsim_08 <- sim(mod_08, n.sim =nsim)
newdat_08 <- expand.grid(REP = levels(final_db$REP),
                         ASYMM_SIDE = levels(final_db$ASYMM_SIDE),
                         s_lat = mean(final_db$s_lat),
                         s_long = mean(final_db$s_long),
                         s_temp = mean(final_db$s_temp),
                         s_rain = mean(final_db$s_rain),
                         s_alt = mean(final_db$s_alt))

Xmat_08 <- model.matrix(~ REP*ASYMM_SIDE +
                          s_long   + s_lat +
                          s_alt + s_temp + s_rain, newdat_08)

fitmat_08 <- matrix(nrow = nrow(newdat_08), ncol=nsim)
for(i in 1:nsim) fitmat_08[,i] <- plogis(Xmat_08 %*% bsim_08@fixef[i,])
newdat_08$lwr <-apply(fitmat_08, 1, quantile, prob=0.025) # upper interval
newdat_08$upr <-apply(fitmat_08, 1, quantile, prob=0.975) # lower interval
newdat_08$fit <-plogis(Xmat_08 %*% fixef(mod_08))   #fitter value

# 2. Estimates and 95% CrI for Table 1
##To get the mean estimate of the fixed effects parameters apply(bsim@fixef, 2, mean)
apply((bsim_08@fixef), 2, mean)
##To get the 95% credible Interval estimate of the fixed effects parameters apply(bsim_mod_testo@coef, 2, quantile, prob=c(0.025, 0.975))
apply((bsim_08@fixef), 2, quantile, prob=c(0.025, 0.975))


### 3. Raw data plot.  ####  

raw_plot_08 <- ggplot(data = final_db, aes(y =match_mismatch, x = REP,
                                           shape=REP))+
  geom_jitter(data = final_db, size = 0.7, alpha = 0.5,
              position = position_jitterdodge(jitter.width =0.8, jitter.height = 0.08, dodge.width = 0.9)) +
  scale_x_discrete(breaks = c("left", "right"), labels = c("Left", "Right")) +
  # scale_y_continuous(breaks = c(0,0.25,0.5,0.75, 1))+
  scale_y_continuous(breaks = c(0, 0.5,0.75,1))+
  scale_shape_manual(values = c(16, 1))+
  theme_classic()+
  geom_hline(yintercept = c(0.75), linetype = "dotted", color="grey69")+  # Ensure this is placed correctly
  geom_hline(yintercept = c(0.5), color="grey69") +  # Ensure this is placed correctly
  facet_grid(.~ASYMM_SIDE)
raw_plot_08



# NICO: I change the variable name 'fit' into 'match_mismatch' so that the aes() defined in raw_plot_05 
# can be found in the final plot and when merging the posterior probabilities with the raw data (newdat and db_vik_sub) there are no 
# errors.

names(newdat_08)[10] <- "match_mismatch"

#### 4. Final plot for Figure 2B

final_plot_08 <- raw_plot_08+ 
  geom_point(data = newdat_08, 
             aes(x = REP, y = match_mismatch), 
             color = "black", size = 3,
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = newdat_08,
                aes(x=REP, ymin= lwr, ymax= upr), 
                color="black", width= 0,linewidth = 1,
                position = position_dodge(width = 0.8))

final_plot_08


