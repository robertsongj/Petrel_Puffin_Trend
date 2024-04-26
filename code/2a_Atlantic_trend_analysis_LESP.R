# Tue Feb 13 15:14:09 2024 ------------------------------

# GJR stole (i.e., forked and branched) this repo from Dave Iles to make a mess
# of things and see if he can get it to work

# then SEG messed around with reporting/visualizing results (and cleaning up code structure to have defined chunks)


# Load libraries ----


library(tidyverse) # for data formatting and plotting
library(readxl)    # for importing xlsx
library(jagsUI)    # for analysis
library(mgcv)      # for creating jagam object (bayesian gams)
library(scales)    # for plotting
library(here)      # for sanity
library(patchwork) # for stacking plots real nice
library(mapview) # for quick looks at spatial data

rm(list=ls())

# Plot theme ----

# Custom theme for plotting


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
                                                            fill = "grey90",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0),size=14),
                            axis.title.x = element_text(margin = margin(10,0,0,0),size=14),
                            axis.text.y = element_text(size=12),
                            axis.text.x = element_text(size=12),
                            panel.background = element_rect(fill = "white"))


# Load and manipulate dataset ----


load(here("input", "Atlantic LHSP colonies clean.RData"))

# align variables with original code

spdat <- mydf

names(spdat) <- tools::toTitleCase(names(spdat))

names(spdat)[names(spdat)=="Sd"] <-  "SE"

mycountry <- "Canada"

spdat <- spdat[spdat$Country == mycountry,]

spdat <- subset(spdat, Year >= 1970)

colonies_to_include = spdat %>%
  group_by(Colony) %>%
  summarize(mean_count = mean(Count),
            recent_count = Count[which.max(Year)],
            first_survey = min(Year),
            last_survey = max(Year),
            n_surveys = length(unique(Year)))

# Omit the smallest LESP colonies that are strongly influencing uncertainty
colonies_to_include = subset(colonies_to_include, mean_count > 1000)
spdat = subset(spdat, Colony %in% colonies_to_include$Colony)

# Tables that link colony numbers to colony names
spdat$colony_numeric <- as.integer(factor(spdat$Colony))
colony_name_table = unique(spdat[,c("Colony","colony_numeric")])

# Tables that link year index to actual year
spdat$year_numeric <- spdat$Year - min(spdat$Year) + 1
year_table = data.frame(Year = min(spdat$Year):max(spdat$Year))
year_table$year_numeric = 1:nrow(year_table)


# Plot raw data ----


ggplot(data = spdat,aes(x = Year, y = Count))+
  geom_point()+
  geom_line(linetype = 2)+
  facet_wrap(Colony~., scales = "free_y")+
  scale_y_continuous(labels = comma)


# Empirical relationship between log CV and log count
# Note: this should probably go in an appendix


ggplot(data = spdat, aes(x = log(Count), y = log(SE)))+
  geom_point()+
  ggtitle("Empirical relationship between log(Count) and log(SE)") + 
  CustomTheme


# Fit model in JAGS and save ----


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

if (!file.exists("output/LESP_fitted.rds")){
  
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
  saveRDS(out, file = "output/LESP_fitted.rds")
}

# Load fitted model ----
out <- readRDS(file = "output/LESP_fitted.rds")

# How long it took to fit the model
out$mcmc.info$elapsed.mins

# ************************************************
# Model convergence ----
# ************************************************

# Which parameters have not converged?
unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)]
# and how bad are we off?
mean(unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)])#1.16
hist(unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)])
# just one parameter >1.3, and two >1.2, the rest are between 1.1-1.2

# Which parameters have effective sample sizes less than 1000?
n.eff <- unlist(out$n.eff)
n.eff <- n.eff[-which(n.eff == 1)]
n.eff[which(n.eff < 1000)]
summary(n.eff[which(n.eff < 1000)])
hist(n.eff[which(n.eff < 1000)])
# a few very small n.eff, but mostly >550


# ************************************************
# Goodness-of-fit assessments ----
# ************************************************

# ~Posterior predictive check~
# Calculate Bayesian p-value by comparing discrepancy measures from simulated datasets
# under a 'perfectly specified model' to discrepancy measures from the empirical data
# (Note: Bayesian p-values close between 0.3 and 0.7 imply that the model fits the data reasonably well
#  i.e., it produces simulated datasets that "look like" the empirical dataset)

Bayesian_pval <- mean(out$sims.list$RMSE_actual > out$sims.list$RMSE_simulated) %>% round(2)
# not terrible?

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

png(paste0("output/figures/goodness_of_fit/LESP_PPC.png"), width = 6, height = 4, units = "in", res = 600)
print(PostPredCheckPlot)
dev.off()


# Plot observed counts versus estimated annual indices at each colony
# Do they track the 1:1 line?


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

png(paste0("output/figures/goodness_of_fit/LESP_Obs_vs_Pred.png"), width = 8, height = 4, units = "in", res = 600)
print(ObsPredPlot)
dev.off()

# ************************************************
# results and plotting ----
# ************************************************

# trajectories ----

# Summarize and plot colony-level trajectories and trends

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
summary(fit_samples_colony$samp)
samples_to_plot <- sample(unique(fit_samples_colony$samp),600)
summary(samples_to_plot)
temp<- subset(fit_samples_colony, samp %in% samples_to_plot)
# I don't think this does what we thought... this subset is just the entire fit_samples_colony dataframe again?

# we're going to try showing all years
# Choose years that are considered "the most reliable" for summarizing trends
# (i.e., the year range where surveys are available at the largest colonies)
t_start <- min(annual_summary_colony$Year)
t_end <- 2023

# we should re-order the colonies by latitude rather than alphabetical so that
# colonies geographically closer together are shown together
# we should also add province to help the reader orient relative to the map

# bring in the known colony coordinates
coords = read_xlsx("data/LESP_colony_coordinates.xlsx", sheet = 1) %>%
  subset(!is.na(Latitude)&!is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),crs=4326, remove = FALSE) %>%
  subset(Country %in% c("Canada")) %>%
  st_transform(AEA_proj) %>%
  dplyr::rename(Count = `Estimated no. of mature individuals`)

# add the coordinates to the colonies_to_include dataframe
colonies_to_include$lat <- NA
colonies_to_include$lat[1] <- coords$Latitude[coords$Colony=="Baccalieu Island"]
colonies_to_include$lat[2] <- coords$Latitude[coords$Colony=="Bon Portage Island"]
colonies_to_include$lat[3] <- coords$Latitude[coords$Colony=="Coleman Island, Wadham Islands"]
colonies_to_include$lat[4] <- coords$Latitude[coords$Colony=="Great Island (Witless Bay)"]
colonies_to_include$lat[5] <- coords$Latitude[coords$Colony=="Green Island (Fortune Bay)"]
colonies_to_include$lat[6] <- coords$Latitude[coords$Colony=="Gull Island (Witless Bay)"]
colonies_to_include$lat[7] <- coords$Latitude[coords$Colony=="Kent Island (Grand Manan Archipelago)"]
colonies_to_include$lat[8] <- coords$Latitude[coords$Colony=="Middle Lawn Island"]
colonies_to_include$lat[9] <- coords$Latitude[coords$Colony=="Small Island, Wadham Islands"]
colonies_to_include$lat[10] <- coords$Latitude[coords$Colony=="Penguin Island, South"]
colonies_to_include$lon <- NA
colonies_to_include$lon[1] <- coords$Longitude[coords$Colony=="Baccalieu Island"]
colonies_to_include$lon[2] <- coords$Longitude[coords$Colony=="Bon Portage Island"]
colonies_to_include$lon[3] <- coords$Longitude[coords$Colony=="Coleman Island, Wadham Islands"]
colonies_to_include$lon[4] <- coords$Longitude[coords$Colony=="Great Island (Witless Bay)"]
colonies_to_include$lon[5] <- coords$Longitude[coords$Colony=="Green Island (Fortune Bay)"]
colonies_to_include$lon[6] <- coords$Longitude[coords$Colony=="Gull Island (Witless Bay)"]
colonies_to_include$lon[7] <- coords$Longitude[coords$Colony=="Kent Island (Grand Manan Archipelago)"]
colonies_to_include$lon[8] <- coords$Longitude[coords$Colony=="Middle Lawn Island"]
colonies_to_include$lon[9] <- coords$Longitude[coords$Colony=="Small Island, Wadham Islands"]
colonies_to_include$lon[10] <- coords$Longitude[coords$Colony=="Penguin Island, South"]

mapview(st_as_sf(colonies_to_include,coords = c("lon", "lat"),crs=4326))

# order the colonies by latitude
# Convert Colony to factor ordered by latitude (north to south)
colonies_to_include$Colony <- factor(colonies_to_include$Colony, levels = colonies_to_include$Colony[order(-colonies_to_include$lat)])
# Check the factor levels
levels(colonies_to_include$Colony)
# Reorder the dataframe
colonies_to_include <- colonies_to_include %>%
  arrange(desc(lat), Colony)

# add province labels (make island of newfoundland NF so we can be consistent with ATPU where we have labrador that will be LAB)
colonies_to_include$Prov <- NA
colonies_to_include$Prov[1:8] <- "NF"
colonies_to_include$Prov[9] <- "NB"
colonies_to_include$Prov[10] <- "NS"
colonies_to_include$Colony.Prov <- paste0(colonies_to_include$Colony, ", ", colonies_to_include$Prov)

# Convert Colony to factor ordered by latitude (north to south)
colonies_to_include$Colony.Prov <- factor(colonies_to_include$Colony.Prov, levels = colonies_to_include$Colony.Prov[order(-colonies_to_include$lat)])
# Check the factor levels
levels(colonies_to_include$Colony.Prov)
# Reorder the dataframe
colonies_to_include <- colonies_to_include %>%
  arrange(desc(lat), Colony.Prov)

# Extract unique Colony levels from the first dataframe
ordered_levels <- levels(colonies_to_include$Colony.Prov)

# now add these to annual_summaries_colony and fit_samples_colony by matching to Colony
annual_summary_colony <- annual_summary_colony %>% left_join(colonies_to_include, by="Colony")
annual_summary_colony$Colony.Prov <- factor(annual_summary_colony$Colony.Prov, levels = ordered_levels)
levels(annual_summary_colony$Colony.Prov) # great
summary(annual_summary_colony$Colony.Prov)

fit_samples_colony <- fit_samples_colony %>% left_join(colonies_to_include, by="Colony")
fit_samples_colony$Colony.Prov <- factor(fit_samples_colony$Colony.Prov, levels = ordered_levels)
levels(fit_samples_colony$Colony.Prov) # great
summary(fit_samples_colony$Colony.Prov)
str(fit_samples_colony)

# Create plot of colony trajectories, overlaid with raw data
colony_trajectory_plot <- ggplot()+
  
  # Entire trajectory (shows the estimates outside the "reliable" window)
  #geom_ribbon(data = annual_summary_colony, 
              # aes(x = Year, ymin=N_q025, ymax = N_q975),
              # alpha = 0.2, fill = "gray80", col = "gray60", 
              # linetype = 2, linewidth = 0.5)+
  
 # geom_line(data = annual_summary_colony, 
            # aes(x = Year, y=N_med), linewidth = 1, col = "gray60")+
  
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
            aes(x = Year, y= N_med), linewidth = 1, col = "black")+
  
  # Thick dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_colony,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 0.5)+
  
  # Observed counts
  geom_point(data = spdat, aes(x = Year, y = Count), size = 2)+
  geom_errorbar(data = spdat, aes(x = Year, ymin = Count-1.96*SE_est, ymax = Count+1.96*SE_est), width = 0, linewidth = 1)+
  
  scale_y_continuous(labels = comma)+
  
  scale_color_manual(values=rep("grey50",length(unique(fit_samples_colony$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance") +
  facet_wrap(~Colony, scales = "free_y", ncol=3)

# Can't get facet_wrap to allow free y scale with the ordered Colony.Prov variable... frustrating.

colony_trajectory_plot

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trajectory_colony_all.years.png", plot=colony_trajectory_plot, 
       device="png", dpi=300, units="cm", width=30, height=25)


# Plot regional trajectory


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
# again, I don't think this does what we think?
samples_to_plot <- sample(unique(fit_samples_regional$samp),600)

# Choose years that are considered "the most reliable" for summarizing trends
# (i.e., the year range where surveys are available at the largest colonies)
t_start <- min(annual_summary_colony$Year)
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
  # geom_ribbon(data = annual_summary_regional, 
  #             aes(x = Year, ymin=N_q025, ymax = N_q975),
  #             alpha = 0.2, fill = "gray80", col = "gray60", 
  #             linetype = 2, linewidth = 0.5)+
  # 
  # geom_line(data = annual_summary_regional, 
  #           aes(x = Year, y=N_med), linewidth = 1, col = "gray60")+
  
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
            aes(x = Year, y= N_med), linewidth = 1, col = "black")+
  
  # Thick darkblue dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_regional,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 0.5)+
  
  coord_cartesian(ylim = c(ylim$min,ylim$max))+
  
  scale_y_continuous(labels = comma)+
  
  scale_color_manual(values=rep("grey50",length(unique(fit_samples_regional$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance")

regional_trajectory_plot

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trajectory_regional_all.years.png", plot=regional_trajectory_plot, 
       device="png", dpi=300, units="cm", width=20, height=20)



# Calculate regional trend ----

# 1973-2023 ----
# first for the whole monitoring period (1973-2023)
t_start <- 1973
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_all <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_all)
quantile(regional_trend_samples_all,c(0.025,0.975))

# Probability trend is negative
mean(regional_trend_samples_all <= 0)

# Probability trend is less than -1% per year
mean(regional_trend_samples_all <= -1)

# Probability trend is less than -2% per year
mean(regional_trend_samples_all <= -2)

# Probability trend is less than -3% per year
mean(regional_trend_samples_all <= -3)

# Histogram illustrating posterior trend estimate
regional_trend_samples_all <- as.data.frame(regional_trend_samples_all)
names(regional_trend_samples_all)<-c("trend")
trend_hist_all <- ggplot(regional_trend_samples_all, aes(x=trend)) + 
  geom_histogram() +
  geom_vline(aes(xintercept=median(trend)), linetype="solid", color = "blue", lwd=1) +
  geom_vline(aes(xintercept=quantile(trend,c(0.025))), linetype="dashed", color = "blue", lwd=1) +
  geom_vline(aes(xintercept=quantile(trend,c(0.975))), linetype="dashed", color = "blue", lwd=1) +
  geom_vline(aes(xintercept=0), linetype="solid", color = "black", lwd=1) +
  labs(x="Trend (% change per year)", y="Probability Density") 
trend_hist_all
ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trend_histogram_1973-2023.png", plot=trend_hist_all, 
       device="png", dpi=300, units="cm", width=20, height=20)

# Violin plot to visualize posterior trend estimate

trend_violin_plot_all <- ggplot(regional_trend_samples_all)+
  geom_hline(yintercept = 0, col = "gray80", linewidth = 2)+
  geom_violin(aes(y = trend, x = "LESP"),
              fill = "grey40",col = "black",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend (% change per year)")+
  #ggtitle("Posterior estimate of\nregional population trend")+
  coord_cartesian(ylim=c(-4,4))
trend_violin_plot_all

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trend_violinplot_1973-2023.png", plot=trend_violin_plot_all, 
       device="png", dpi=300, units="cm", width=20, height=20)

# 1984-2023 ----
# and again since 1984 when bacallieu monitoring started
t_start <- 1984
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_1984 <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_1984)
quantile(regional_trend_samples_1984,c(0.025,0.975))

# Probability trend is negative
mean(regional_trend_samples_1984 <= 0)

# Probability trend is less than -1% per year
mean(regional_trend_samples_1984 <= -1)

# Probability trend is less than -2% per year
mean(regional_trend_samples_1984 <= -2)

# Probability trend is less than -3% per year
mean(regional_trend_samples_1984 <= -3)

# Histogram illustrating posterior trend estimate
regional_trend_samples_1984 <- as.data.frame(regional_trend_samples_1984)
names(regional_trend_samples_1984)<-c("trend")
trend_hist_1984 <- ggplot(regional_trend_samples_1984, aes(x=trend)) + 
  geom_histogram() +
  geom_vline(aes(xintercept=median(trend)), linetype="solid", color = "blue", lwd=1) +
  geom_vline(aes(xintercept=quantile(trend,c(0.025))), linetype="dashed", color = "blue", lwd=1) +
  geom_vline(aes(xintercept=quantile(trend,c(0.975))), linetype="dashed", color = "blue", lwd=1) +
 geom_vline(aes(xintercept=0), linetype="solid", color = "black", lwd=1) +
  labs(x="Trend (% change per year)", y="Probability Density") 
trend_hist_1984
ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trend_histogram_1984-2023.png", plot=trend_hist_1984, 
       device="png", dpi=300, units="cm", width=20, height=20)

# Violin plot to visualize posterior trend estimate

trend_violin_plot_1984 <- ggplot(regional_trend_samples_1984)+
  geom_hline(yintercept = 0, col = "gray80", linewidth = 2)+
  geom_violin(aes(y = trend, x = "LESP"),
              fill = "grey40",col = "black",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend (% change per year)")+
  #ggtitle("Posterior estimate of\nregional population trend")+
  coord_cartesian(ylim=c(-4,4))
trend_violin_plot_1984

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trend_violinplot_1984-2023.png", plot=trend_violin_plot_1984, 
       device="png", dpi=300, units="cm", width=20, height=20)

# both periods ----
# we could combine them to show them together as well
regional_trend_samples_1984$period <- "1984-2023"
regional_trend_samples_all$period <- "1973-2023"

regional_trend_samples <- rbind(regional_trend_samples_1984, regional_trend_samples_all)

trend_violin_plot <- ggplot(regional_trend_samples, aes(group=period))+
  geom_hline(yintercept = 0, col = "gray80", linewidth = 2)+
  geom_violin(aes(y = trend, x = period),
              fill = "grey40",col = "black",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend (% change per year)")+
  #ggtitle("Posterior estimate of\nregional population trend")+
  coord_cartesian(ylim=c(-4,2))
trend_violin_plot

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trend_violinplot_2period.png", plot=trend_violin_plot, 
       device="png", dpi=300, units="cm", width=20, height=20)

# trajectory + violin trend ----

p1 <- regional_trajectory_plot + annotate("text", label = "A", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)
p2 <- trend_violin_plot + annotate("text", label = "B", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)

p <- p1|p2
p

ggsave(filename="output/figures/trajectory_and_trend_plots/LESP_trajectory.plus.trend.png", plot=p, 
       device="png", dpi=300, units="cm", width=30, height=15)

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
