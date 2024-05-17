# Trying out the 2a_ATPU on ATPU

# Load libraries ----


library(tidyverse) # for data formatting and plotting
library(readxl)    # for importing xlsx
library(jagsUI)    # for analysis
library(mgcv)      # for creating jagam object (bayesian gams)
library(scales)    # for plotting
library(here)      # for sanity
library(patchwork) # for stacking plots real nice
library(mapview) # for quick looks at spatial data
library(sf)

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

load(here("input", "Atlantic ATPU colonies clean.RData"))

# align variables with original code

names(spdat) <- tools::toTitleCase(names(spdat))

spdat <- subset(spdat, Year >= 1965) 

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

colonies_to_include = spdat %>%
  group_by(Colony) %>%
  summarize(mean_count = mean(Count),
            recent_count = Count[which.max(Year)],
            first_survey = min(Year),
            last_survey = max(Year),
            n_surveys = length(unique(Year)))

# keep them all for now

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

# Load fitted model ----
out <- readRDS(file = "output/ATPU_fitted.rds")

# How long it took to fit the model
out$mcmc.info$elapsed.mins

# ************************************************
# Model convergence ----
# ************************************************

# Which parameters have not converged?
unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)]
# and how bad are we off?
mean(unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)])#1.13
hist(unlist(out$Rhat)[which(unlist(out$Rhat)>1.1)], breaks=15)
# all between 1.1-1.5

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
# actually amazing at 0.45

# Plot results of Posterior predictive check: We want to see values clustered along the 1:1 line
lim <- range(c(out$sims.list$RMSE_actual,out$sims.list$RMSE_simulated))
RMSE_df = data.frame(actual = out$sims.list$RMSE_actual, simulated = out$sims.list$RMSE_simulated)
PostPredCheckPlot <- ggplot(data = RMSE_df,
                            aes(x = actual, y = simulated )) +
  geom_hex(binwidth = diff(lim)/50) +
  scale_fill_gradientn(colors = c("gray95","darkblue"), name = "Number of\nsimulated\ndatasets\n") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+
  coord_cartesian(xlim = lim,ylim=lim)+
  xlab("RMSE Actual Datasets")+
  ylab("RMSE Simulated Datasets")+
  theme_bw() + 
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
        legend.title = element_text(size=14))

PostPredCheckPlot

ggsave(filename="output/figures/goodness_of_fit/ATPU_PPC.png", plot=PostPredCheckPlot, 
       device="png", dpi=300, units="cm", width=20, height=15)

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
                          col = as.factor(Colony)))+
  geom_abline(slope = 1, intercept = 0, linetype="dashed")+
  geom_point()+
  geom_errorbar(width=0)+
  geom_errorbarh(height=0)+
  scale_colour_manual(values=viridis::plasma(length(unique(annual_summary_colony$Colony)))) +
  xlab("Observed Count (log-scale)")+
  ylab("Estimated Population Index (log-scale)")+
  scale_y_continuous(trans="log10", labels = comma)+
  scale_x_continuous(trans="log10", labels = comma)+
  theme(legend.position="none") +
  theme(axis.title.y = element_text(margin = margin(0,10,0,0),size=14),
        axis.title.x = element_text(margin = margin(10,0,0,0),size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12))

ObsPredPlot

ggsave(filename="output/figures/goodness_of_fit/ATPU_Obs_vs_Pred.png", plot=ObsPredPlot, 
       device="png", dpi=300, units="cm", width=20, height=20)

# results and plotting ----

# trajectories ----

# *colony level ----
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
# t_start <- min(annual_summary_colony$Year)
# t_end <- 2023

# we should re-order the colonies by latitude rather than alphabetical so that
# colonies geographically closer together are shown together

coords = read.csv("data/ATPU_colony_coordinates_for.maps.csv") %>%
  subset(!is.na(Lat)&!is.na(Lon)) %>%
  st_as_sf(coords = c("Lon", "Lat"),crs=4326, remove = FALSE) %>%
  subset(Country %in% c("Canada"))

# add the coordinates from coord to the colonies_to_include dataframe
colonies_to_include$lat <- NA
colonies_to_include$lat[1] <- coords$Lat[coords$Colony=="Bacalhao Island, LB"]
colonies_to_include$lat[2] <- coords$Lat[coords$Colony=="Baccalieu Island, NF"]
colonies_to_include$lat[3] <- coords$Lat[coords$Colony=="Coleman Island, NF"]
colonies_to_include$lat[4] <- coords$Lat[coords$Colony=="Gannet Clusters, GC 2"]
colonies_to_include$lat[5] <- coords$Lat[coords$Colony=="Gannet Clusters, GC 3"]
colonies_to_include$lat[6] <- coords$Lat[coords$Colony=="Gannet Clusters, GC 4"]
colonies_to_include$lat[7] <- coords$Lat[coords$Colony=="Gannet Clusters, GC 5"]
colonies_to_include$lat[8] <- coords$Lat[coords$Colony=="Gannet Clusters, GC 6"]
colonies_to_include$lat[9] <- coords$Lat[coords$Colony=="Great Island (Witless Bay), NF"]
colonies_to_include$lat[10] <- coords$Lat[coords$Colony=="Gull Island (Witless Bay), NF"]
colonies_to_include$lat[11] <- coords$Lat[coords$Colony=="Herring Islands 1 (Red Island)"]
colonies_to_include$lat[12] <- coords$Lat[coords$Colony=="Herring Islands 2 (Kipper Island)"]
colonies_to_include$lat[13] <- coords$Lat[coords$Colony=="Herring Islands 3 (Sardine Island)"]
colonies_to_include$lat[14] <- coords$Lat[coords$Colony=="Machias Seal Island, NB"]
colonies_to_include$lat[15] <- coords$Lat[coords$Colony=="North Bird Island, NF"]
colonies_to_include$lat[16] <- coords$Lat[coords$Colony=="North Green Island (Lab), LB"]
colonies_to_include$lat[17] <- 50.198 # entering manually for North Shore MBSs, as a central point
colonies_to_include$lat[18] <- coords$Lat[coords$Colony=="Pee Pee Island, NF"]
colonies_to_include$lat[19] <- coords$Lat[coords$Colony=="Puffin Islands (Lab), LB"]
colonies_to_include$lat[20] <- coords$Lat[coords$Colony=="Small Island, Wadham Is, NF"]
colonies_to_include$lat[21] <- coords$Lat[coords$Colony=="Penguin Island, South, NF"]
colonies_to_include$lat[22] <- coords$Lat[coords$Colony=="Tinker Island, nr. Indian Hr., LB"]
colonies_to_include$lon <- NA
colonies_to_include$lon[1] <- coords$Lon[coords$Colony=="Bacalhao Island, LB"]
colonies_to_include$lon[2] <- coords$Lon[coords$Colony=="Baccalieu Island, NF"]
colonies_to_include$lon[3] <- coords$Lon[coords$Colony=="Coleman Island, NF"]
colonies_to_include$lon[4] <- coords$Lon[coords$Colony=="Gannet Clusters, GC 2"]
colonies_to_include$lon[5] <- coords$Lon[coords$Colony=="Gannet Clusters, GC 3"]
colonies_to_include$lon[6] <- coords$Lon[coords$Colony=="Gannet Clusters, GC 4"]
colonies_to_include$lon[7] <- coords$Lon[coords$Colony=="Gannet Clusters, GC 5"]
colonies_to_include$lon[8] <- coords$Lon[coords$Colony=="Gannet Clusters, GC 6"]
colonies_to_include$lon[9] <- coords$Lon[coords$Colony=="Great Island (Witless Bay), NF"]
colonies_to_include$lon[10] <- coords$Lon[coords$Colony=="Gull Island (Witless Bay), NF"]
colonies_to_include$lon[11] <- coords$Lon[coords$Colony=="Herring Islands 1 (Red Island)"]
colonies_to_include$lon[12] <- coords$Lon[coords$Colony=="Herring Islands 2 (Kipper Island)"]
colonies_to_include$lon[13] <- coords$Lon[coords$Colony=="Herring Islands 3 (Sardine Island)"]
colonies_to_include$lon[14] <- coords$Lon[coords$Colony=="Machias Seal Island, NB"]
colonies_to_include$lon[15] <- coords$Lon[coords$Colony=="North Bird Island, NF"]
colonies_to_include$lon[16] <- coords$Lon[coords$Colony=="North Green Island (Lab), LB"]
colonies_to_include$lon[17] <- -60.20473 # entering manually for North Shore MBSs, as a central point
colonies_to_include$lon[18] <- coords$Lon[coords$Colony=="Pee Pee Island, NF"]
colonies_to_include$lon[19] <- coords$Lon[coords$Colony=="Puffin Islands (Lab), LB"]
colonies_to_include$lon[20] <- coords$Lon[coords$Colony=="Small Island, Wadham Is, NF"]
colonies_to_include$lon[21] <- coords$Lon[coords$Colony=="Penguin Island, South, NF"]
colonies_to_include$lon[22] <- coords$Lon[coords$Colony=="Tinker Island, nr. Indian Hr., LB"]

mapview(st_as_sf(colonies_to_include,coords = c("lon", "lat"),crs=4326))

# Convert Colony to factor ordered by latitude (north to south)
colonies_to_include$Colony <- factor(colonies_to_include$Colony, levels = colonies_to_include$Colony[order(-colonies_to_include$lat)])
# Check the factor levels
levels(colonies_to_include$Colony)
# Reorder the dataframe
colonies_to_include <- colonies_to_include %>%
  arrange(desc(lat), Colony)

# Extract unique Colony levels from the first dataframe
ordered_levels <- levels(colonies_to_include$Colony)

# now add these to annual_summaries_colony and fit_samples_colony and spdat by matching to Colony
annual_summary_colony <- annual_summary_colony %>% left_join(colonies_to_include, by="Colony")
annual_summary_colony$Colony <- factor(annual_summary_colony$Colony, levels = ordered_levels)
levels(annual_summary_colony$Colony) # great
summary(annual_summary_colony$Colony)

fit_samples_colony <- fit_samples_colony %>% left_join(colonies_to_include, by="Colony")
fit_samples_colony$Colony <- factor(fit_samples_colony$Colony, levels = ordered_levels)
levels(fit_samples_colony$Colony) # great
summary(fit_samples_colony$Colony)
str(fit_samples_colony)

spdat <- spdat %>% left_join(colonies_to_include, by="Colony")
spdat$Colony <- factor(spdat$Colony, levels = ordered_levels)
levels(spdat$Colony) # great
summary(spdat$Colony)
str(spdat)

samples_to_plot <- as.factor(samples_to_plot)
fit_samples_colony$samp <- as.factor(fit_samples_colony$samp)

# Create plot of colony trajectories, overlaid with raw data
colony_trajectory_plot <- ggplot()+
  
  # Plot 1000 trajectories from Bayesian posterior
geom_line(data = subset(fit_samples_colony, 
                        samp %in% samples_to_plot), 
          aes(x = Year, y = N_pred, col = factor(samp)),alpha = 0.1)+
  
  # Thick line for posterior median
  geom_line(data = annual_summary_colony, 
            aes(x = Year, y= N_med), linewidth = 1, col = "black")+
  
  # Thick dashed line for 95% CI
  geom_ribbon(data = annual_summary_colony, 
              aes(x = Year, ymin=N_q025, ymax = N_q975),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 0.5)+
  
  # Observed counts
  geom_point(data = spdat, aes(x = Year, y = Count), size = 2)+
  geom_errorbar(data = spdat, aes(x = Year, ymin = Count-1.96*SE_est, ymax = Count+1.96*SE_est), width = 0, linewidth = 1)+
  geom_point(data = subset(spdat, is.na(SE)), aes(x = Year, y = Count), size = 5, shape=1)+
  #geom_blank(data = limits_df, aes(x=x, y=y_max)) +
  scale_x_continuous(limits=c(1965,2023), expand=c(0,0)) +  
  scale_y_continuous(labels = comma) +
  scale_color_manual(values=rep("grey50",length(unique(fit_samples_colony$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance") +
  facet_wrap(~Colony, ncol=3, scales = "free_y") 

colony_trajectory_plot

ggsave(filename="output/figures/trajectory_and_trend_plots/ATPU_trajectory_colony_all.years.png", plot=colony_trajectory_plot, 
       device="png", dpi=300, units="cm", width=30, height=40)


# really want to get the y axes to shrink to have a max of just a fraction more than the upper CI
# but facet_wrap doesn't play nice with specifying different axes limits
# what if I plotted each colony separately in a loop, then stitched them all back together...
colonies <- unique(spdat$Colony)
# make a list to save the plots to
plot_list <- list()

for(i in 1:length(colonies)){
  a <- subset(fit_samples_colony, Colony==colonies[i])
  b <- subset(annual_summary_colony, Colony==colonies[i])
  c <- subset(spdat, Colony==colonies[i])
  
  p <- ggplot()+
    geom_line(data = subset(a,samp %in% samples_to_plot), 
              aes(x = Year, y = N_pred, group=factor(samp)),color="grey50", alpha = 0.1)+
    geom_line(data = b, 
              aes(x = Year, y= N_med), linewidth = 1, col = "black")+
    geom_ribbon(data = b, 
                aes(x = Year, ymin=N_q025, ymax = N_q975),
                fill = "transparent", col = "black", 
                linetype = 2, linewidth = 0.5)+
    geom_point(data = c, aes(x = Year, y = Count), size = 2)+
    geom_errorbar(data = c, aes(x = Year, ymin = Count-1.96*SE_est, ymax = Count+1.96*SE_est), width = 0, linewidth = 1)+
    geom_point(data = subset(c, is.na(SE)), aes(x = Year, y = Count), size = 5, shape=1) +
    scale_x_continuous(limits=c(1965,2023), expand=c(0,0)) +  
    scale_y_continuous(labels = comma, limits=c(0,max(b$N_q975)+0.1*max(b$N_q975))) +
    scale_color_manual(values=rep("grey50",length(unique(fit_samples_colony$samp))), 
                       guide = "none")+
    ylab("Index of Abundance") +
    xlab("Year") +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    
    facet_wrap(~Colony, ncol=1) 
  
  plot_list[[i]] <- p
}
plot_list[[3]] # Ok that's much better, not I just need to stack them all up!

# now add the axis titles and text just to the plots that need them
plot_list[[13]] <- plot_list[[13]] + ylab("Index of Abundance") + theme(axis.title.y=element_text())
plot_list[[22]] <- plot_list[[22]] + xlab("Year") + theme(axis.title.x=element_text(), axis.text.x=element_text())
plot_list[[21]] <- plot_list[[21]] + xlab("Year") + theme(axis.title.x=element_text(), axis.text.x=element_text())
plot_list[[20]] <- plot_list[[20]] + xlab("Year") + theme(axis.title.x=element_text(), axis.text.x=element_text())

p <- plot_grid(plotlist = plot_list, ncol = 3, align = "v")
p

# much better
ggsave(filename="output/figures/trajectory_and_trend_plots/ATPU_trajectory_colony_all.years.png", 
       plot=p, device="png", dpi=300, units="cm", width=30, height=40)


# *regional----


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

# turn off scientific notation
options(scipen=999)

# Create plot of regional trajectory
regional_trajectory_plot <- ggplot()+

  # Plot 1000 trajectories from Bayesian posterior
geom_line(data = subset(fit_samples_regional, 
                        samp %in% samples_to_plot & 
                          Year >= t_start & 
                          Year <= t_end), 
          aes(x = Year, y = N_pred/1000000, col = factor(samp)),alpha = 0.1)+
  
  # Thick darkblue line for posterior median
  geom_line(data = subset(annual_summary_regional, 
                          Year >= t_start & 
                            Year <= t_end), 
            aes(x = Year, y= N_med/1000000), linewidth = 1, col = "black")+
  
  # Thick darkblue dashed line for 95% CI
  geom_ribbon(data = subset(annual_summary_regional,Year >= t_start & 
                              Year <= t_end), 
              aes(x = Year, ymin=N_q025/1000000, ymax = N_q975/1000000),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 0.5)+
  
  
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  scale_x_continuous(limits=c(t_start,t_end), breaks = seq(1970, 2020, by = 10), expand = c(0, 0))+
  
  scale_color_manual(values=rep("grey50",length(unique(fit_samples_regional$samp))), 
                     guide = "none")+
  
  ylab("Index of Abundance\n(millions)")

regional_trajectory_plot

ggsave(filename="output/figures/trajectory_and_trend_plots/ATPU_trajectory_regional_all.years.png", plot=regional_trajectory_plot, 
       device="png", dpi=300, units="cm", width=20, height=20)


# Calculate regional trend ----

# first calculate trend in 5-day windows to find where the changepoints are (trend = 0)
trend.5.yr.windows <- fit_samples_regional %>%
  group_by(samp) %>%
  arrange(Year) %>%
  summarize(Year = Year, 
            N_pred_start = lag(N_pred, 2),  # 2 years before
         N_pred_end = lead(N_pred, 2),   # 2 years ahead
         trend = 100 * ((N_pred_end/N_pred_start)^(1/(5))-1)) %>% # 5-year window
  filter(!is.na(N_pred_start), !is.na(N_pred_end)) 

trend.5.yr.windows.summary <- trend.5.yr.windows %>% group_by(Year) %>%
  summarize(trend_med = median(trend),
            trend_q025 = quantile(trend,c(0.025)),
            trend_q975 = quantile(trend,c(0.975)))

trend_5.yr <- ggplot()+
  geom_line(data = trend.5.yr.windows, aes(x = Year, y= trend, col = factor(samp)),alpha = 0.1)+
  scale_color_manual(values=rep("grey50",length(unique(trend.5.yr.windows$samp))), 
                     guide = "none")+
  geom_line(data = trend.5.yr.windows.summary, aes(x = Year, y= trend_med), linewidth = 1, col = "black")+
  geom_ribbon(data = trend.5.yr.windows.summary,aes(x = Year, ymin=trend_q025, ymax = trend_q975),
              fill = "transparent", col = "black", 
              linetype = 2, linewidth = 0.5)+
  geom_hline(yintercept=0) +
  scale_x_continuous(limits=c(1965,2023), breaks = seq(1970, 2020, by = 10),expand = c(0, 0))+
  ylab("5-year Trend\n(% change per year)")
trend_5.yr
ggsave(filename="output/figures/trajectory_and_trend_plots/ATPU_trend_5.yr.windows.png", plot=trend_5.yr, 
       device="png", dpi=300, units="cm", width=20, height=20)

# so where does the median cross the 0 line?
subset(trend.5.yr.windows.summary, trend_med >-1 & trend_med <1)
# definitely around 1977-1978 and a switch again around 2005

# 1965-1978 ----
# first period 
t_start <- 1965
t_end <- 1978
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_1 <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_1)
quantile(regional_trend_samples_1,c(0.025,0.975))

# Probability trend is negative
mean(regional_trend_samples_1 < 0)

regional_trend_samples_1 <- as.data.frame(regional_trend_samples_1)
names(regional_trend_samples_1)<-c("trend")

# 1978-2005 ----
# second period
t_start <- 1978
t_end <- 2005
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_2 <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_2)
quantile(regional_trend_samples_2,c(0.025,0.975))

# Probability trend is positive
mean(regional_trend_samples_2 > 0)

regional_trend_samples_2 <- as.data.frame(regional_trend_samples_2)
names(regional_trend_samples_2)<-c("trend")
# 2005-2023 ----
# third period
t_start <- 2005
t_end <- 2023
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_3 <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_3)
quantile(regional_trend_samples_3,c(0.025,0.975))

# Probability trend is positive
mean(regional_trend_samples_3 > 0)

regional_trend_samples_3 <- as.data.frame(regional_trend_samples_3)
names(regional_trend_samples_3)<-c("trend")
# 1990-2023 ----
# and again since 1990 which is the last 3 generations, since generation time is 11 years according to Bird et al. 2020
t_start <- 1990
regional_indices_t_start <- subset(fit_samples_regional, Year == t_start)
regional_indices_t_end <- subset(fit_samples_regional, Year == t_end)
regional_trend_samples_4 <- 100 * ((regional_indices_t_end$N_pred/regional_indices_t_start$N_pred)^(1/(t_end-t_start))-1)

# Posterior median trend estimate (and 95% credible interval)
median(regional_trend_samples_4)
quantile(regional_trend_samples_4,c(0.025,0.975))

# Probability trend is positive
mean(regional_trend_samples_4 > 0)

regional_trend_samples_4 <- as.data.frame(regional_trend_samples_4)
names(regional_trend_samples_4)<-c("trend")
#multiple periods ----
# we could combine them to show them together as well
regional_trend_samples_4$period <- "1990-2023\n(last three generations)"
regional_trend_samples_3$period <- "2005-2023\n"
regional_trend_samples_2$period <- "1978-2005\n"
regional_trend_samples_1$period <- "1965-1978\n"

regional_trend_samples <- rbind(regional_trend_samples_1, regional_trend_samples_2, regional_trend_samples_3, regional_trend_samples_4)

# order the periods
regional_trend_samples$period <- factor(regional_trend_samples$period, levels=c("1965-1978\n","1978-2005\n","2005-2023\n","1990-2023\n(last three generations)"))

trend_violin_plot <- ggplot(regional_trend_samples, aes(group=period))+
  geom_hline(yintercept = 0)+
  geom_violin(aes(y = trend, x = period),
              fill = "grey40",col = "black",
              alpha = 0.5,
              draw_quantiles = c(0.025,0.5,0.975))+
  xlab("")+
  ylab("Trend\n(% change per year)")+
  #ggtitle("Posterior estimate of\nregional population trend")+
  coord_cartesian(ylim=c(-10,10))
trend_violin_plot

# trajectory + 5 yr trend + violin trend ----

p1 <- regional_trajectory_plot + annotate("text", label = "A", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)
p2 <- trend_5.yr + annotate("text", label = "B", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)
p3 <- trend_violin_plot + annotate("text", label = "C", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)

p <- p1/p2/p3
p

ggsave(filename="output/figures/trajectory_and_trend_plots/ATPU_trajectory.plus.trends.png", plot=p, 
       device="png", dpi=300, units="cm", width=20, height=30)

# saving ----

# save the workspace
save.image("input/ATPU_2a_workspace_04.30.2024")

# save the "colonies to include" dataframe for mapping along with puffins
write.csv(colonies_to_include, "input/ATPU_recent.counts.coords.csv", row.names = FALSE)

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
