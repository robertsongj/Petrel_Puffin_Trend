
# Sarah having a play at modifying these simulations to actually have only 2-6 surveys per "colony" 
# and have 65% of surveys with SE
# and produce output for incoporation into ms

my_packs = c('tidyverse','readxl','RColorBrewer','viridis','jagsUI','mgcv','ggrepel','scales','ggthemes')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

library(viridis)
library(ggplot2)
# Plot theme ----

# Custom theme for plotting


CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(panel.grid.major = element_line(colour = 'transparent'),
                            panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        size = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "grey90",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0),size=14),
                            axis.title.x = element_text(margin = margin(10,0,0,0),size=14),
                            axis.text.y = element_text(size=12),
                            axis.text.x = element_text(size=12),
                            panel.background = element_rect(fill = "white"),
                            legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
                            legend.box="vertical",
                            legend.box.background = element_rect(color = "black",fill = "white"),
                            legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
                            legend.background = element_rect(color = "white"),
                            legend.text = element_text(size=16),
                            legend.title = element_text(size=14))



# 1: Simulate an entire 50-year time series for each of 15 colonies ----

start_time <- Sys.time()

trend_results <- data.frame()

for (run in 1:250){
  
  if (run %in% trend_results$run) next
  
  set.seed(run)
  
  ncolony <- 15
  nyears <- 50
  
  N_matrix <- matrix(NA,nrow=ncolony,ncol = nyears)
  
  # Initial abundance
  N_matrix[,1] <- rlnorm(ncolony, meanlog = log(10000), sdlog = 0.5) %>% sort()
  
  # Set a colony-specific variance in annual growth rates
  process_sd = runif(ncolony,0,0.3)
  
  # Simulate trajectories at each colony
  for (i in 1:ncolony){
    
    
    # Code to generate log-linear trajectories
    
    # 
    # trend <- rnorm(1,0,0.03)
    # 
    # # Generate a random trajectory at the colony (random walk)
    # for (t in 2:nyears) N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + trend)
    # 
    
    # Code to generate complex random walks, where some colonies 'switch' their average tendency way through the simulation
    
    
    # switch dynamics for some colonies (adds complexity to population dynamics)
    switch = sample(0:1,1)
    
    r_mean = runif(1,-0.02,0.01)
    r_switch <- runif(1,-0.1,0.05)
    
    # Generate a random trajectory at the colony (random walk)
    for (t in 2:nyears){
      N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + rnorm(1,mean=r_mean, sd = process_sd[i]))
      if (switch == 1 & (t > nyears/2))N_matrix[i,t] <- exp(log(N_matrix[i,t-1]) + rnorm(1,mean=r_mean+r_switch, sd = process_sd[i]))
    }
    
  }

# Convert to dataframe (to plot with ggplot)
N_df <- reshape2::melt(N_matrix) %>% 
  rename(Colony = Var1,Year=Var2,N=value)

# Plot dynamics (on log10 scale)
ggplot()+
  geom_line(aes(x = 1:nyears, y = colSums(N_matrix), col = "Regional Sum"),linewidth = 1)+
  geom_line(data = N_df, aes(x = Year, y = N, col = factor(Colony)))+
  scale_y_continuous(labels = comma, trans = "log10")+
  ggtitle("Simulated trajectories at each of 9 colonies")+
  ylab("Abundance")+
  xlab("Year")


# 2: Simulate intermittent surveys (2-6 surveys at each colony) ----


N_df$Count <- N_df$lambda_obs <- N_df$log_observation_SE <- NA

for (i in 1:ncolony){
  
  # Rows in N_df corresponding to observations at this colony
  colony_rows <- which(N_df$Colony == i)
  
  # Range of survey error at this colony among years (some years are worse than others)
  colony_obs_SE_lower <- 0.05 # observation SE is at least 0.05 (log scale)
  colony_obs_SE_upper <- runif(1,0.3,0.5) # in some years, obs SE is anywhere from 0.3 to 0.5
  
  # Magnitude of observation error varies among surveys
  N_df$log_observation_SE[colony_rows] <- runif(length(colony_rows), colony_obs_SE_lower,colony_obs_SE_upper)
  
  N_df$lambda_obs[colony_rows] <- rlnorm(length(colony_rows), meanlog = log(N_df$N)[colony_rows], sdlog = N_df$log_observation_SE[colony_rows])
  

  # Which surveys were actually conducted?
  
  survey_years <- c()
  
  # Simulate one count in first 10 years of surveys
  survey_years <- c(survey_years, sample(1:10,1))
  
  # Simulate one count in final 10 years of surveys
  survey_years <- c(survey_years, sample(nyears:(nyears-10),1))
  
  # Simulate 0-4 additional surveys
  survey_years <- c(survey_years, sample(11:(nyears-11),sample(0:4,1)))
  
  # Poisson observations
  N_df$Count[(N_df$Colony == i) & 
                     (N_df$Year %in% survey_years)] <- rpois(n = length(survey_years), 
                                                             lambda = N_df$lambda_obs[(N_df$Colony == i) & (N_df$Year %in% survey_years)])
  
}

# Confidence intervals on counts
N_df$Count_lci <- exp(log(N_df$Count) - 1.96*N_df$log_observation_SE)
N_df$Count_uci <- exp(log(N_df$Count) + 1.96*N_df$log_observation_SE)

# Standard errors on counts
N_df$SE <- ((N_df$Count_uci-N_df$Count)+(N_df$Count-N_df$Count_lci)/2)/(2*1.96)

# and keep standard errors for only 65% of counts to imitate real datasets
# index to where the SEs are not NA
rows.w.SE <- which(!is.na(N_df$SE))

# now randomly pick which 35% of the rows to drop
SEs.to.drop <- sample(rows.w.SE, size = floor(0.35 * length(rows.w.SE)))
# now change those to NA in the dataframe
N_df$SE[SEs.to.drop] <- NA

# how many counts did we get?
length(which(!is.na(N_df$Count)))
# and how many now have SEs?
length(which(!is.na(N_df$SE)))

# Plot survey counts at each of the colonies
ggplot()+
  geom_line(data = N_df, aes(x = Year, y = N), col = "gray80")+
  geom_errorbar(data = N_df, aes(x = Year, ymin = Count_lci, ymax = Count_uci), col = "black", width = 0)+
  geom_point(data = N_df, aes(x = Year, y = Count), col = "black")+
  facet_wrap(Colony~., scales = "free_y")+
  scale_y_continuous(labels = comma)+
  ggtitle("Simulated surveys at colonies")


# 3: Fit model to simulated survey data ----


spdat <- N_df %>%
  filter(!is.na(Count)) %>% 
  dplyr::select(Colony,Year,Count,SE)
summary(as.factor(spdat$Colony))
str(spdat)

year_table = data.frame(Year = min(spdat$Year):max(spdat$Year))
year_table$year_numeric = 1:nrow(year_table)

# Data for import into jags
nknots = 6
year <- spdat$Year - min(year_table$Year) + 1
ymax <- max(year_table$year_numeric)
nyears = length(1:ymax)
colony = spdat$Colony
ncolony <- max(colony)
count <- round(spdat$Count) # Must be integer
ncounts = length(count)
C0prior <- 0 # approximate mean of exp(N), gets model closer to the intercept to start

# Use jagam to prepare basis functions
nyearspred = length(1:ymax)
preddat = data.frame(yrs = 1:ymax,count = 1)
form = as.formula(paste("count ~ s(yrs,k =",nknots,")"))
gamprep = jagam(formula = form,
                data = preddat,
                file = "code/tempgam.txt",
                centred = T)

# Package data into a list for JAGS
jags_data = list(X = gamprep$jags.data$X,
                 S1 = gamprep$jags.data$S1,
                 zero = gamprep$jags.data$zero,
                 colony = colony,
                 ncounts = ncounts,
                 ncolony = ncolony,
                 count = count,
                 nknots = nknots,
                 nyearspred = nyearspred,
                 year = year,
                 survey_count = spdat$Count,
                 survey_SE = spdat$SE,
                 C0prior = C0prior)

# Fit model using JAGS
parameters.to.save = c("sdnoise","sdbeta","C","beta.X","population_index")

  # NOTE: CONSIDER DIFFERENT NUMBERS OF KNOTS (K VALUES)
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter =  550000,
              n.burnin = 50000,
              n.thin = 2500,
              model.file = "code/Seabird_Model-no int.jags",
              n.chains = 3,
              parallel = TRUE)
  
out$mcmc.info$elapsed.mins


# 4: Summarize predictions and compare to true (i.e., simulated) trajectories ----


# Extract predictions in dataframe format
fit_samples = reshape2::melt(out$sims.list$population_index) %>%
  rename(samp = Var1, Year = Var2, Colony = Var3, N_pred = value)

N_summary_colony = fit_samples %>% 
  group_by(Colony, Year) %>%
  summarize(q025 = quantile(N_pred,0.025),
            q50 = quantile(N_pred,0.500),
            mean = mean(N_pred),
            q975 = quantile(N_pred,0.975)) %>%
  
  # Join with true values
  full_join(N_df)

# Plot estimates
ggplot()+
  geom_ribbon(data = N_summary_colony, aes(x = Year, ymin = q025, ymax = q975), alpha = 0.2, fill = "dodgerblue")+
  geom_line(data = N_summary_colony, aes(x = Year, y = q50, col = "Estimate"))+
  
  geom_line(data = N_df, aes(x = Year, y = N, col = "True Trajectory"))+
  geom_errorbar(data = N_df, aes(x = Year, ymin = Count_lci, ymax = Count_uci), col = "black", width = 0)+
  geom_point(data = N_df, aes(x = Year, y = Count, col = "Observed Count"))+
  facet_wrap(Colony~.)+
  scale_y_continuous(trans="log10", labels = comma)+
  ggtitle("Simulated trajectories at each of 9 colonies")+
  scale_color_manual(values = c("dodgerblue","black","red"), name = "")+
  ylab("Index of abundance")


# 5: Calculate regional total ----


# True regional total
regional_df <- N_df %>%
  group_by(Year) %>%
  summarize(N = sum(N))

# Estimated regional total
regional_samples = fit_samples %>% 
  group_by(Year,samp) %>%
  summarize(N_pred = sum(N_pred)) 

# Summary (mean and 95% CI)
N_summary_regional <- regional_samples %>%
  group_by(Year) %>%
  summarize(q025 = quantile(N_pred,0.025),
            q50 = quantile(N_pred,0.500),
            mean = mean(N_pred),
            q975 = quantile(N_pred,0.975)) %>%
  
  # Join with true values
  full_join(regional_df)

# Trend estimate
baseline_year <- 1
trend_true <- 100*((regional_df$N[regional_df$Year == nyears]/regional_df$N[regional_df$Year == baseline_year])^(1/(nyears-baseline_year))-1)

trend_est <- 100*((regional_samples$N_pred[regional_samples$Year == nyears]/regional_samples$N_pred[regional_samples$Year == baseline_year])^(1/(nyears-baseline_year))-1)
trend_est <- quantile(trend_est,c(0.025,0.5,0.975))

ggplot()+
  geom_ribbon(data = N_summary_regional, aes(x = Year, ymin = q025, ymax = q975), alpha = 0.2, fill = "dodgerblue")+
  geom_line(data = N_summary_regional, aes(x = Year, y = q50, col = "Estimate"))+
  
  geom_line(data = regional_df, aes(x = Year, y = N, col = "True Trajectory"))+
  scale_y_continuous(labels = comma)+
  ggtitle("Regional trajectory")+
  scale_color_manual(values = c("dodgerblue","black","red"), name = "")+
  ylab("Index of abundance")+
  geom_text(aes(x = 0, 
                y = max(c(N_summary_regional$q975,N_summary_regional$N))), 
            label = paste0("True trend = ",round(trend_true,2),"% per year\nEst trend = ",round(trend_est[2],2),"% (",round(trend_est[1],2)," to ",round(trend_est[3],2),")"), hjust=0)


# 6: Append results for this simulation run to dataframe ----


trend_results <- rbind(trend_results,data.frame(run = run,
                                                trend_true = trend_true,
                                                trend_est_q025 = trend_est[1],
                                                trend_est_q500 = trend_est[2],
                                                trend_est_q975 = trend_est[3],
                                                cov = trend_true > trend_est[1] & trend_true < trend_est[3],
                                                max_Rhat = max(out$Rhat$population_index),
                                                n_counts = length(which(!is.na(spdat$Count))),
                                                n_counts_w_SE = length(which(!is.na(spdat$SE)))))
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time
cat("Elapsed time:", elapsed_time, "\n")

# save the trend_results dataframe
#write.csv(trend_results, "output/LESP_250simulations_allSE_trendresults.csv", row.names = FALSE)
#write.csv(trend_results, "output/LESP_250simulations_0.65SE_trendresults.csv", row.names = FALSE)
write.csv(trend_results, "output/LESP_250simulations_0.65SE_15colonies_trendresults.csv", row.names = FALSE)
# Elapsed time: 21.46644 

# 7:Summarize and plot results across repeated simulations ----
# load saved .csv from here if just plotting results
trend_results <- read.csv("output/LESP_250simulations_0.65SE_15colonies_trendresults.csv")

mean(trend_results$cov) # credible interval coverage (proportion of estimated trend credible intervals containing the "true" trend)
mean(trend_results$trend_est_q500 - trend_results$trend_true) # mean bias or accuracy (mean estimated trend is within x of true trend)
sd(trend_results$trend_est_q500 - trend_results$trend_true)
summary(trend_results$trend_est_q500 - trend_results$trend_true)
mean(trend_results$trend_est_q975 - trend_results$trend_est_q025) # precision (mean size of 95% credible interval)
# what proportion of trend estimates are the wrong sign based on the credible interval being entirely positive or negative?
(sum((trend_results$trend_true > 0 & trend_results$trend_est_q975 < 0) |
      (trend_results$trend_true < 0 & trend_results$trend_est_q025 > 0)))/nrow(trend_results)*100


# make coverage a factor
trend_results$cov_f <- factor(trend_results$cov, labels = c("No", "Yes"))
# fix so that no plots on top
trend_results$cov_f <- factor(trend_results$cov_f, levels = c("Yes", "No"))


col_pal <- viridis::plasma(9)
col_pal <- col_pal[c(2,6)]

trend_plot <- ggplot(data = trend_results, aes(x = trend_true, y = trend_est_q500, 
                                               ymin = trend_est_q025, ymax = trend_est_q975,col=cov_f))+
  geom_abline(intercept=0,slope=1,col="black", linetype="dashed")+
  geom_vline(xintercept=0,col="black")+
  geom_hline(yintercept=0,col="black")+
  geom_errorbar(width=0)+
  geom_point()+
  geom_point(data = subset(trend_results, cov_f=="Yes"), aes(x = trend_true, y = trend_est_q500, 
                                       ymin = trend_est_q025, ymax = trend_est_q975), color=col_pal[1])+
  geom_errorbar(data = subset(trend_results, cov_f=="Yes"), aes(x = trend_true, y = trend_est_q500, 
                                                             ymin = trend_est_q025, ymax = trend_est_q975), color=col_pal[1],width=0)+
  geom_point(data = subset(trend_results, cov_f=="Yes"), aes(x = trend_true, y = trend_est_q500, 
                                                            ymin = trend_est_q025, ymax = trend_est_q975),color=col_pal[1])+
  geom_point(data = subset(trend_results, cov_f=="No"), aes(x = trend_true, y = trend_est_q500, 
                                                            ymin = trend_est_q025, ymax = trend_est_q975),color=col_pal[2])+
  geom_errorbar(data = subset(trend_results, cov_f=="No"), aes(x = trend_true, y = trend_est_q500, 
                                                                ymin = trend_est_q025, ymax = trend_est_q975), color=col_pal[2],width=0)+
  coord_cartesian(ylim=c(-4,7),xlim=c(-4,7))+
  theme_bw()+
  xlab("True (Simulated) Regional Trend")+
  ylab("Estimated Regional Trend")+
  scale_color_manual(values=col_pal, name = "Credible\nInterval\nCoverage") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.y = element_text(margin = margin(0,10,0,0),size=14),
        axis.title.x = element_text(margin = margin(10,0,0,0),size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
        legend.box="vertical",
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        legend.position=c(0.15,0.8))
trend_plot

ggsave(filename="output/figures/goodness_of_fit/LESP_250simulations_0.65SE_15colonies_trendresults.png", plot=trend_plot, 
       device="png", dpi=300, units="cm", width=20, height=20)

# is there a relationship between true trend and estimated trend accuracy?
trend_results$accuracy <- trend_results$trend_est_q500 - trend_results$trend_true
ggplot(data=subset(trend_results, trend_true<8), aes(x=abs(accuracy), y=abs(trend_true))) + 
  geom_point() + 
  geom_smooth(method="lm") +
  coord_cartesian(ylim = c(0,7), xlim = c(0,3))
summary(lm(abs(trend_true) ~ abs(accuracy), data = subset(trend_results, trend_true<8)))
# positive relationship with slope of 0.97 and adjusted r-squared of 0.27, probably worth reporting?

# and what about the simulated datasets themselves?
summary(trend_results$n_counts) #median of 60 counts, min of 42 and max of 75
