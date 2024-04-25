# Trying out the 2a_LESP on ATPU

# ------------------------------------------------
# Load libraries
# ------------------------------------------------

library(tidyverse) # for data formatting and plotting
library(readxl)    # for importing xlsx
library(jagsUI)    # for analysis
library(mgcv)      # for creating jagam object (bayesian gams)
library(scales)    # for plotting
library(here)      # for sanity

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

load(here("input", "Atlantic ATPU colonies clean.RData"))

# align variables with original code

names(spdat) <- tools::toTitleCase(names(spdat))

spdat <- subset(spdat, Year >= 1965) # was 1970 for LESP

colonies_to_include = spdat %>%
  group_by(Colony) %>%
  summarize(mean_count = mean(Count),
            first_survey = min(Year),
            last_survey = max(Year),
            n_surveys = length(unique(Year)))

# Omit the smallest ATPU colonies that are strongly influencing uncertainty
# colonies_to_include = subset(colonies_to_include, mean_count > 1000)
# spdat = subset(spdat, Colony %in% colonies_to_include$Colony)

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
  "B.X",
  
  # Annual population indices
  "population_index",
  
  # Discrepancy measures for posterior predictive checks (goodness-of-fit testing)
  "RMSE_actual",
  "RMSE_simulated",
  
  # Penalty measures
  'lambda'
)

if (!file.exists("output/ATPU_fitted.rds")){
  
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
unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)]

# Which parameters have effective sample sizes less than 1000?
n.eff <- unlist(out$n.eff)
n.eff <- n.eff[-which(n.eff == 1)]
n.eff[which(n.eff < 1000)]

# ************************************************
# Goodness-of-fit assessments
# ************************************************

# ~Posterior predictive check~
# Calculate Bayesian p-value by comparing discrepancy measures from simulated datasets
# under a 'perfectly specified model' to discrepancy measures from the empirical data
# (Note: Bayesian p-values close between 0.3 and 0.7 imply that the model fits the data reasonably well
#  i.e., it produces simulated datasets that "look like" the empirical dataset)

Bayesian_pval <- mean(out$sims.list$RMSE_actual > out$sims.list$RMSE_simulated) %>% round(2)

# Plot results of Posterior predictive check: We want to see values clustered along the 1:1 line
lim <- range(c(out$sims.list$RMSE_actual,out$sims.list$RMSE_simulated))
RMSE_df = data.frame(actual = out$sims.list$RMSE_actual, simulated = out$sims.list$RMSE_simulated)
PostPredCheckPlot <- ggplot(data = RMSE_df,
                            aes(x = actual, y = simulated )) +
  geom_hex(binwidth = diff(lim)/50) +
  scale_fill_gradientn(colors = c("gray95","darkblue"), name = "Number of\nsimulated\ndatasets") +
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim,ylim=lim)+
  ggtitle(paste0("Posterior predictive check: \n\nBayesian p-value = ",Bayesian_pval))+
  xlab("RMSE (actual datasets)")+
  ylab("RMSE (simulated datasets)")+
  theme_bw()

PostPredCheckPlot

png(paste0("output/figures/goodness_of_fit/ATPU_PPC.png"), width = 6, height = 4, units = "in", res = 600)
print(PostPredCheckPlot)
dev.off()

# ----------------------------------
# Plot observed counts versus estimated annual indices at each colony
# Do they track the 1:1 line?
# ----------------------------------

# Estimates of SE for each survey (fills in missing ones)
spdat$SE_est = out$mean$survey_SE

# Extract predictions of annual indices in dataframe format
fit_samples_colony = reshape2::melt(out$sims.list$population_index) %>%
  rename(samp = Var1, year_numeric = Var2, colony_numeric = Var3, N_pred = value) %>%
  full_join(colony_name_table) %>% full_join(year_table)

# Summarize
annual_summary_colony = fit_samples_colony %>% 
  group_by(Colony, Year) %>%
  summarize(q025 = quantile(N_pred,0.025),
            q50 = quantile(N_pred,0.500),
            mean = mean(N_pred),
            q975 = quantile(N_pred,0.975))

# Join with observed data
annual_summary_colony <- annual_summary_colony %>% full_join(spdat)

ObsPredPlot <- ggplot(annual_summary_colony, 
                      aes(x = Count, xmin = Count - 1.96*SE_est, xmax = Count + 1.96*SE_est, 
                          y = q50, ymin = q025, ymax = q975, 
                          col = Colony))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbar(width=0)+
  geom_errorbarh(height=0)+
  xlab("Observed Count")+
  ylab("Estimated Population Index")+
  scale_y_continuous(trans="log10", labels = comma)+
  scale_x_continuous(trans="log10", labels = comma)+
  ggtitle("Predicted annual indices vs observed counts\n\n(Note that axes are on a log scale)")

ObsPredPlot



png(paste0("output/figures/goodness_of_fit/ATPU_Obs_vs_Pred.png"), width = 8, height = 4, units = "in", res = 600)
print(ObsPredPlot)
dev.off()

# ************************************************
# SUMMARIZE AND PLOT RESULTS
# ************************************************

# ----------------------------------
# Summarize and plot colony-level trajectories and trends
# ----------------------------------


spdat$SE_est = out$mean$survey_SE

# Extract predictions of annual indices in dataframe format
fit_samples_colony = reshape2::melt(out$sims.list$population_index) %>%
  rename(samp = Var1, year_numeric = Var2, colony_numeric = Var3, N_pred = value) %>%
  full_join(colony_name_table) %>% full_join(year_table)

annual_summary_colony = fit_samples_colony %>% 
  group_by(Colony, Year) %>%
  summarize(N_med = quantile(N_pred,0.500),
            N_q025 = quantile(N_pred,0.025),
            N_q975 = quantile(N_pred,0.975))

# Choose 1000 samples from the posterior to visualize on plot
samples_to_plot <- sample(unique(fit_samples_colony$samp),600)

# Choose years that are considered "the most reliable" for summarizing trends
# (i.e., the year range where surveys are available at the largest colonies)
t_start <- 1984
t_end <- 2023

# Create plot of colony trajectories, overlaid with raw data
colony_trajectory_plot <- ggplot()+
  
  # Entire trajectory (shows the estimates outside the "reliable" window)
  geom_ribbon(data = annual_summary_colony, 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              alpha = 0.2, fill = "gray80", col = "gray60", 
              linetype = 2, linewidth = 0.5)+
  
  geom_line(data = annual_summary_colony, 
            aes(x = Year, y=N_med), linewidth = 1, col = "gray60")+
  
  # Plot 1000 trajectories from Bayesian posterior
  geom_line(data = subset(fit_samples_colony, 
                          samp %in% samples_to_plot & 
                            Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
  
  # Thick line for posterior median
  geom_line(data = subset(annual_summary_colony, 
                          Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y= N_med), linewidth = 1, col = "darkblue")+
  
  # Thick dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_colony,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              fill = "transparent", col = "darkblue", 
              linetype = 2, linewidth = 0.5)+
  
  # Observed counts
  geom_point(data = spdat, aes(x = Year, y = Count), size = 2)+
  geom_errorbar(data = spdat, aes(x = Year, ymin = Count-1.96*SE_est, ymax = Count+1.96*SE_est), width = 0, linewidth = 1)+
  
  scale_y_continuous(labels = comma)+
  
  scale_color_manual(values=rep("dodgerblue",length(unique(fit_samples_colony$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance") +
  facet_wrap(Colony~., scales = "free_y")

colony_trajectory_plot

png(paste0("output/figures/trajectory_and_trend_plots/ATPU_trajectory_colony.png"), width = 12, height = 6, units = "in", res = 600)
print(colony_trajectory_plot)
dev.off()

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
            N_q025 = quantile(N_pred,0.025),
            N_q975 = quantile(N_pred,0.975))

# Choose 1000 samples from the posterior to visualize on plot
samples_to_plot <- sample(unique(fit_samples_regional$samp),500)

# Choose years that are considered "the most reliable" for summarizing trends
# (i.e., the year range where surveys are available at the largest colonies)
t_start <- 1984
t_end <- 2023

# Choose y limits for plot
ylim = annual_summary_regional %>%
  subset(Year >= t_start & Year <= t_end) %>%
  summarize(min = 0,
            max = max(N_q975))

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
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              alpha = 0.2, fill = "gray80", col = "gray60", 
              linetype = 2, linewidth = 0.5)+
  
  geom_line(data = annual_summary_regional, 
            aes(x = Year, y=N_med), linewidth = 1, col = "gray60")+
  
  # Plot 1000 trajectories from Bayesian posterior
  geom_line(data = subset(fit_samples_regional, 
                          samp %in% samples_to_plot & 
                            Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
  
  # Thick darkblue line for posterior median
  geom_line(data = subset(annual_summary_regional, 
                          Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y= N_med), linewidth = 1, col = "darkblue")+
  
  # Thick darkblue dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_regional,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              fill = "transparent", col = "darkblue", 
              linetype = 2, linewidth = 0.5)+
  
  coord_cartesian(ylim = c(ylim$min,ylim$max))+
  
  scale_y_continuous(labels = comma)+
  
  scale_color_manual(values=rep("blue",length(unique(fit_samples_regional$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance")

regional_trajectory_plot

png(paste0("output/figures/trajectory_and_trend_plots/ATPU_trajectory_regional.png"), width = 6, height = 4, units = "in", res = 600)
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
png(paste0("output/figures/trajectory_and_trend_plots/ATPU_trend_histogram.png"), width = 6, height = 4, units = "in", res = 600)

hist(regional_trend_samples, breaks = 100, 
     col = "dodgerblue", border = "dodgerblue",
     main = "Posterior estimate of regional population trend",
     xlab = "Trend (% change per year)",
     ylab = "Probability Density",
     freq = FALSE)
abline(v = 0)
abline(v = median(regional_trend_samples), col = "blue", lwd=3)
abline(v = quantile(regional_trend_samples,c(0.025,0.975)), col = "blue", lwd=3, lty =2)
dev.off()

# Violin plot to visualize posterior trend estimate

trend_violin_plot <- ggplot()+
  geom_hline(yintercept = 0, col = "gray80", linewidth = 2)+
  geom_violin(aes(y = regional_trend_samples, x = "ATPU"),
              fill = "dodgerblue",col = "blue",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend (% change per year)")+
  ggtitle("Posterior estimate of\nregional population trend")+
  coord_cartesian(ylim=c(-4,4))
trend_violin_plot

png(paste0("output/figures/trajectory_and_trend_plots/ATPU_trend_violin.png"), width = 4, height = 4, units = "in", res = 600)
print(trend_violin_plot)
dev.off()


# Wed Feb 14 10:55:02 2024 ------------------------------

# some other ways to look at the data

out$mean$C
out$mean$beta.X
out$mean$B.X
out$mean$sdbeta
out$mean$lambda
out$mean$ProcVar_sd
hist(out$sims.list$ProcVar_sd)

out$sims.list$population_index
mean(rgamma(1000,0.05, 0.005))
