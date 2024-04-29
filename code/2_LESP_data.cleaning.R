library(tidyverse)
library(readxl) 

# cleaning up the ATPU excel to match the LESP input file format

spdat = read_xlsx("data/LESP_trend_data.xlsx", sheet = 1) %>%
  dplyr::rename(Count = `Mature individuals`)

# Calculate SE for counts with upper and lower confidence limit but no SE
for(i in 1:nrow(spdat)){
  if(is.na(spdat$SE[i]) & !(is.na(spdat$`Lower CI`[i]))){
    spdat$SE[i] <- (spdat$`Upper CI`[i] - spdat$`Lower CI`[i])/3.92
  }
}

spdat <- spdat[,c(2:5)]

# save as RData

save(spdat, file="input/Atlantic LESP colonies clean.RData")
