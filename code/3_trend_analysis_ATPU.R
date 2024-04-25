# ------------------------------------------------
# Load libraries
# ------------------------------------------------

library(tidyverse) # for data formatting and plotting
library(readxl)    # for importing xlsx
library(jagsUI)    # for analysis
library(mgcv)      # for creating jagam object (bayesian gams)
library(scales)    # for plotting

#setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Seabirds/Petrel_Puffin_Trend/")

rm(list=ls())

# ------------------------------------------------
# Custom theme for plotting
# ------------------------------------------------

CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA), 
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'transparent'),
                            panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        size = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))

# ------------------------------------------------
# Load dataset
# ------------------------------------------------

spdat = read_xlsx("data/ATPU_trend_data_SEGupdate.xlsx", sheet = 1) %>%
  dplyr::rename(Count = `Mature individuals`)

spdat <- subset(spdat, Year >= 1965) # making this 1965 for ATPU vs 1970 for LESP because quite a few good counts from 1965

colonies_to_include = spdat %>%
  group_by(Colony) %>%
  summarize(mean_count = mean(Count),
            first_survey = min(Year),
            last_survey = max(Year),
            n_surveys = length(unique(Year)))

# Omit the smallest LESP colonies that are strongly influencing uncertainty
colonies_to_include = subset(colonies_to_include, mean_count > 1000) # keeping this from LESP but should compare with and without
spdat = subset(spdat, Colony %in% colonies_to_include$Colony)

# Tables that link colony numbers to colony names
spdat$colony_numeric <- as.integer(factor(spdat$Colony))
colony_name_table = unique(spdat[,c("Colony","colony_numeric")])

# Tables that link year index to actual year
spdat$year_numeric <- spdat$Year - min(spdat$Year) + 1
year_table = data.frame(Year = min(spdat$Year):max(spdat$Year))
year_table$year_numeric = 1:nrow(year_table)

# ------------------------------------------------
# Plot raw data
# ------------------------------------------------

ggplot(data = spdat,aes(x = Year, y = Count))+
  geom_point()+
  geom_line(linetype = 2)+
  facet_wrap(Colony~., scales = "free_y")+
  scale_y_continuous(labels = comma)

# ------------------------------------------------
# Empirical relationship between log CV and log count
# Note: this should probably go in an appendix
# ------------------------------------------------

ggplot(data = spdat, aes(x = log(Count), y = log(SE)))+
  geom_point()+
  ggtitle("Empirical relationship between log(Count) and log(SE)")

# ------------------------------------------------
# Fit model in JAGS and save
# ------------------------------------------------

# Data for import into jags
nknots = 6
year <- spdat$Year - min(year_table$Year) + 1
ymax <- max(year_table$year_numeric)
nyears = length(1:ymax)
colony = spdat$colony_numeric
ncolony <- max(colony)
count <- round(spdat$Count) # Must be integer
ncounts = length(count)

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
                 survey_SE = spdat$SE)

# Fit model using JAGS
parameters.to.save = c(
  
  # Estimate of annual process variance
  "ProcVar_sd",  
  
  # Parameters for estimating standard error for surveys that don't have SE estimates available
  "intercept_SE",
  "slope_SE",
  "sd_SE",
  
  # Magnitude of observation error for each survey
  "survey_SE",
  
  # Variance in GAM parameters
  "sdbeta",
  
  # GAM components
  "C",
  "beta.X",
  
  # Annual population indices
  "population_index",
  
  # Discrepancy measures for posterior predictive checks (goodness-of-fit testing)
  "RMSE_actual",
  "RMSE_simulated"
)

if (!file.exists("output/ATPU_fitted.rds")){
  
  # NOTE: CONSIDER DIFFERENT NUMBERS OF KNOTS (K VALUES)
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter =  5500000,
              n.burnin = 500000,
              n.thin = 2500,
              model.file = "code/Seabird_Model.jags",
              n.chains = 3,
              parallel = TRUE)
  
  # Save fitted model
  saveRDS(out, file = "output/ATPU_fitted.rds")
}

# Load fitted model
out <- readRDS(file = "output/ATPU_fitted.rds")

# How long it took to fit the model
out$mcmc.info$elapsed.mins 

# ************************************************
# Model convergence
# ************************************************

# Which parameters have not converged?

# Which parameters have effective sample sizes less than 1000?
n.eff <- unlist(out$n.eff)
n.eff <- n.eff[-which(n.eff == 1)]
n.eff[which(n.eff < 1000)]

# ************************************************
# Goodness-of-fit assessment
# ************************************************

# ~Posterior predictive check~
# Calculate Bayesian p-value by comparing discrepancy measures from simulated datasets
# under a 'perfectly specified model' to discrepancy measures from the empirical data
# (Note: Bayesian p-values close to 0.5 imply that the model fits the data reasonably well
#  i.e., it produces simulated datasets that "look like" the empirical dataset)

Bayesian_pval <- mean(out$sims.list$RMSE_actual > out$sims.list$RMSE_simulated) %>% round(2)

# Plot results of Posterior predictive check: We want to see values clustered along the 1:1 line
lim <- range(c(out$sims.list$RMSE_actual,out$sims.list$RMSE_simulated))
RMSE_df = data.frame(actual = out$sims.list$RMSE_actual, simulated = out$sims.list$RMSE_simulated)
plot1 <- ggplot(data = RMSE_df,
                aes(x = actual, y = simulated )) +
  geom_hex(binwidth = diff(lim)/50) +
  scale_fill_gradientn(colors = c("gray95","darkblue"), name = "Number of\nsimulated\ndatasets") +
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim,ylim=lim)+
  ggtitle(paste0("Posterior predictive check: \n\nBayesian p-value = ",pval_adult))+
  xlab("RMSE (actual dataset)")+
  ylab("RMSE (simulated dataset)")+
  theme_bw()
plot1

# ************************************************
# SUMMARIZE AND PLOT RESULTS
# ************************************************

# ----------------------------------
# Summarize and plot colony-level trajectories and trends
# ----------------------------------

# Extract predictions of annual indices in dataframe format
fit_samples_colony = reshape2::melt(out$sims.list$population_index) %>%
  rename(samp = Var1, year_numeric = Var2, colony_numeric = Var3, N_pred = value) %>%
  full_join(colony_name_table) %>% full_join(year_table)

# Estimates of SE for each survey (fills in missing ones)
spdat$SE_est = out$mean$survey_SE

annual_summary_colony = fit_samples_colony %>% 
  group_by(Colony, Year) %>%
  summarize(q025 = quantile(N_pred,0.025),
            q50 = quantile(N_pred,0.500),
            mean = mean(N_pred),
            q975 = quantile(N_pred,0.975))

# Plot results for each colony
colony_plot_freeaxis = ggplot() +
  
  # Full time series
  geom_ribbon(data = annual_summary_colony, aes(x = Year, ymin = q025, ymax = q975), fill = "dodgerblue",col = "transparent", alpha = 0.3)+
  geom_line(data = annual_summary_colony, aes(x = Year, y = q50), col = "dodgerblue", linewidth = 1, alpha = 0.7)+
  
  geom_errorbar(data = spdat, aes(x = Year, ymin = Count - 1.96*SE_est, ymax = Count + 1.96*SE_est), width = 0)+
  geom_point(data = spdat,aes(x = Year, y = Count))+
  xlab("Year")+
  ylab("Population index")+
  
  facet_wrap(Colony~., scales = "free_y")+
  scale_y_continuous(labels = comma)

colony_plot_freeaxis

# ----------------------------------
# Plot regional trajectory
# ----------------------------------

fit_samples_regional = fit_samples_colony %>% 
  group_by(samp,Year) %>%
  summarize(N_pred = sum(N_pred))

# Summarize the median and 95% CI
annual_summary_regional <- fit_samples_regional %>%
  group_by(Year) %>%
  summarize(N_med = median(N_pred),
            N_lcl = quantile(N_pred,0.025),
            N_ucl = quantile(N_pred,0.975))

# Choose 1000 samples from the posterior to visualize on plot
samples_to_plot <- sample(unique(fit_samples_regional$samp),1000)

# Choose years that are considered "the most reliable" for summarizing trends
# (i.e., the year range where surveys are available at the largest colonies)
t_start <- 1984
t_end <- 2023

# Choose y limits for plot
ylim = annual_summary_regional %>%
  subset(Year >= t_start & Year <= t_end) %>%
  summarize(min = 0,
            max = max(N_ucl))

ylim = fit_samples_regional %>%
  ungroup() %>%
  subset(samp %in% samples_to_plot & 
           Year >= t_start & 
           Year <= t_end) %>%
  
  summarize(min = 0,
            max = max(N_pred))

# Create plot of regional trajectory
regional_trajectory_plot <- ggplot()+
  
  # Entire trajectory (shows the estimates outside the "reliable" window)
  geom_ribbon(data = annual_summary_regional, 
              aes(x = Year, ymin=N_lcl, ymax = N_ucl),
              alpha = 0.2, fill = "gray80", col = "gray60", 
              linetype = 2, linewidth = 0.5)+
  
  geom_line(data = annual_summary_regional, 
            aes(x = Year, y=N_med), linewidth = 1, col = "gray60")+
  
  # "Reliable" survey window
  geom_ribbon(data = subset(annual_summary_regional, 
                            Year >= t_start & Year <= t_start), 
              aes(x = Year, ymin=N_lcl, ymax = N_ucl),
              alpha = 0.1, fill = "dodgerblue", col = "black", linetype = 2, linewidth = 0.5)+
  
  # Plot 1000 trajectories from Bayesian posterior
  geom_line(data = subset(fit_samples_regional, 
                          samp %in% samples_to_plot & 
                            Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
  
  # Thick black line for posterior median
  geom_line(data = subset(annual_summary_regional, 
                          Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y= N_med), linewidth = 1, col = "black")+
  
  # Thick black dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_regional,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_lcl, ymax = N_ucl),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 1)+
  
  coord_cartesian(ylim = c(ylim$min,ylim$max))+
  
  scale_y_continuous(labels = comma)+
  
  scale_color_manual(values=rep("blue",length(unique(fit_samples_regional$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance")
regional_trajectory_plot

png(paste0("output/figures/trajectory_plots/LESP_trajectory_regional.png"), width = 6, height = 4, units = "in", res = 600)
print(regional_trajectory_plot)
dev.off()

# ----------------------------------
# Calculate regional trend
# ----------------------------------

regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples)
quantile(regional_trend_samples,c(0.025,0.975))

# Probability trend is negative
mean(regional_trend_samples <= 0)

# Probability trend is less than -1% per year
mean(regional_trend_samples <= -1)

# Probability trend is less than -2% per year
mean(regional_trend_samples <= -2)

# Probability trend is less than -3% per year
mean(regional_trend_samples <= -3)

# Histogram illustrating posterior trend estimate
hist(regional_trend_samples, breaks = 500, 
     col = "dodgerblue", border = "dodgerblue",
     main = "Posterior estimate of regional population trend",
     xlab = "Trend (% change per year)",
     ylab = "Probability Density",
     freq = FALSE)
abline(v = 0)
abline(v = median(regional_trend_samples), col = "blue", lwd=3)
abline(v = quantile(regional_trend_samples,c(0.025,0.975)), col = "blue", lwd=3, lty =2)

# Violin plot to visualize posterior trend estimate

trend_plot <- ggplot()+
  geom_hline(yintercept = 0, col = "gray80", linewidth = 2)+
  geom_violin(aes(y = regional_trend_samples, x = "LESP"),
              fill = "dodgerblue",col = "blue",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend (% change per year)")+
  ggtitle("Posterior estimate of regional population trend")+
  coord_cartesian(ylim=c(-4,4))
trend_plot
