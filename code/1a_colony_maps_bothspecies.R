

# map LESP and ATPU colonies together, either just the most recent counts for
# colonies included in models or all known colonies

library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(sf)
library(sp)
library(ggspatial) # for scalebars on ggplot maps
library(ggsflabel)
library(patchwork)
library(ggrepel)
library(mapview)
library(viridis)

options(scipen = 999)

# bring in monitored ATPU colonies ----
col.atpu.mon <- read.csv("input/ATPU_recent.counts.coords.csv")

sf_use_s2(FALSE)

col.atpu.mon <- st_as_sf(col.atpu.mon,coords = c("lon", "lat"),crs=4326) 
mapview(col.atpu.mon)
hist(col.atpu.mon$recent_count, breaks=50)
col.atpu.mon <- col.atpu.mon %>% mutate(count_f = case_when(
  recent_count %in% 0:1000 ~ "<1,000",
  recent_count %in% 1000:5000 ~ "1,000 to 5,000",
  recent_count %in% 5000:10000 ~ "5,000 to 10,000",
  recent_count %in% 10000:50000 ~ "10,000 to 50,000",
  recent_count %in% 50000:100000 ~ "50,000 to 100,000",
  recent_count %in% 100000:500000 ~ "100,000 to 500,000",
  recent_count %in% 500000:1000000 ~ "500,000 to 1,000,000",
  recent_count > 1000000 ~ "> 1,000,000"))
summary(as.factor(col.atpu.mon$count_f))
levels(as.factor(col.atpu.mon$count_f))
# make it an ordered factor
col.atpu.mon$count_f <- factor(col.atpu.mon$count_f, levels=c("<1,000",
                                                      "1,000 to 5,000",
                                                      "5,000 to 10,000",
                                                      "10,000 to 50,000",
                                                      "50,000 to 100,000",
                                                      "100,000 to 500,000",
                                                      "500,000 to 1,000,000",
                                                      "> 1,000,000"))
summary(col.atpu.mon$count_f)

# bring in monitored LESP colonies ----
col.lesp.mon <- read.csv("input/LESP_recent.counts.coords.csv")

col.lesp.mon <- st_as_sf(col.lesp.mon,coords = c("lon", "lat"),crs=4326) 
mapview(col.lesp.mon)

col.lesp.mon <- col.lesp.mon %>% mutate(count_f = case_when(
  recent_count %in% 0:1000 ~ "<1,000",
  recent_count %in% 1000:5000 ~ "1,000 to 5,000",
  recent_count %in% 5000:10000 ~ "5,000 to 10,000",
  recent_count %in% 10000:50000 ~ "10,000 to 50,000",
  recent_count %in% 50000:100000 ~ "50,000 to 100,000",
  recent_count %in% 100000:500000 ~ "100,000 to 500,000",
  recent_count %in% 500000:1000000 ~ "500,000 to 1,000,000",
  recent_count > 1000000 ~ "> 1,000,000"))
summary(as.factor(col.lesp.mon$count_f))

col.lesp.mon$count_f <- factor(col.lesp.mon$count_f, levels=c("<1,000",
                                                      "1,000 to 5,000",
                                                      "5,000 to 10,000",
                                                      "10,000 to 50,000",
                                                      "50,000 to 100,000",
                                                      "100,000 to 500,000",
                                                      "500,000 to 1,000,000",
                                                      "> 1,000,000"))
summary(col.lesp.mon$count_f)

# basemap 1 ----

bbox <- st_bbox(c(xmin = -68, xmax = -51, ymin = 42.5, ymax = 55), crs = 4326)

# downloaded from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
coastline <- st_read('data/spatial_data/Shoreline_Coastline_GSHHG/GSHHS_shp/f/GSHHS_f_L1.shp') %>% st_crop(bbox)
mapview(coastline)

basemap <- ggplot(data=coastline) + geom_sf(fill="gray95", size=0.01) + 
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"), 
        panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1), 
        panel.background = element_rect(fill = "aliceblue")) 
basemap 

basemap <- basemap +
  annotate(geom="text", x = -63.2, y=45.1, label = "Nova Scotia", size=4, colour="grey30") +
  annotate(geom="text", x = -66.3, y=46, label = "New Brunswick", size=4, colour="grey30") +
  annotate(geom="text", x = -56, y=48.5, label = "Newfoundland", size=4, colour="grey30") +
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) 
basemap

col_pal <- viridis::plasma(8)
col_pal <- rev(col_pal)

basemap + geom_point(data=col.lesp.mon, aes(geometry=geometry, color=count_f), size=3, stat="sf_coordinates") + 
  scale_color_manual(values=col_pal)

basemap + geom_point(data=col.atpu.mon, aes(geometry=geometry, color=count_f), size=3, stat="sf_coordinates") + 
  scale_color_manual(values=col_pal)

# map monitored colonies both species ----
# put them together in one dataframe

col.lesp.mon$sp <- "LESP"
col.atpu.mon$sp <- "ATPU"
col <- rbind(col.lesp.mon, col.atpu.mon)

basemap + geom_sf(data=col, aes(geometry=geometry, color=count_f, shape=sp), size=5, stroke=1.5) + 
  scale_color_manual(values=col_pal) + scale_shape_manual(values=c(1,2)) +
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) 

# nah that's too busy, we'll just stack them but keep the same basemap for easy comparison, but use this dataframe so we get the full legend

col_pal <- viridis::plasma(9)
col_pal <- rev(col_pal)
col_pal <- col_pal[c(2:9)]

p1 <- basemap + geom_sf(data=subset(col, sp=="LESP"), aes(geometry=geometry, color=count_f), size=5, stroke=1, shape=1) + 
  scale_color_manual(values=col_pal, drop = FALSE) +
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) +
  labs(color = "most recent count") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
        legend.box="vertical",
        legend.position="none", 
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=14))
p1

p2 <- basemap + geom_sf(data=subset(col, sp=="ATPU"), aes(geometry=geometry, color=count_f), size=5, stroke=1, shape=1) + 
  scale_color_manual(values=col_pal, drop = FALSE) +
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) +
  labs(color = "most recent count") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
        legend.box="vertical",
        #legend.position=c(0.8,0.17), 
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=14))
p2

p1 <- p1 + annotate("text", label = "A", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)
p2 <- p2 + annotate("text", label = "B", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)

p <- p1|p2
p

ggsave(filename="output/figures/maps/both.species.colonies.w.trends.png", plot=p, 
       device="png", dpi=300, units="cm", width=30, height=15)

# need to add province borders and labels but decent start

# bring in all ATPU colonies ----
col.atpu <- read.csv("input/ATPU_colony_coordinates_for.maps.csv")

col.atpu <- st_as_sf(col.atpu,coords = c("Lon", "Lat"),crs=4326) 
mapview(col.atpu)
head(col.atpu)
# remove the colonies with zero individuals
col.atpu <- subset(col.atpu, Individuals != 0)
# create a numeric recent_count from Individuals which will change "present" to NA
col.atpu$recent_count <- as.numeric(col.atpu$Individuals)

hist(col.atpu$recent_count, breaks=50)

# now we need to replace add all of the monitored colonies with the most recent count and more accurate position
col.atpu$sp <- "ATPU"
col.atpu <- col.atpu[,c("sp","Colony","recent_count","geometry")]
col.atpu.mon <- col.atpu.mon[,c("sp","Colony","recent_count","geometry")]
# unfortunately colony naming differs so we'll need to do add all of the monitored colonies then delete the duplicates
col.atpu <- rbind(col.atpu,col.atpu.mon)
col.atpu <- col.atpu[-c(11,16,17,19,20,21,25:29,55,56,58,63,66,69,71,72,84,111),] # also remove the generic combo NS MBS in QC and we'll just circle the general area

# now make a factor where we can give NAs a "present" factor
col.atpu <- col.atpu %>% mutate(count_f = case_when(
  is.na(recent_count) ~ "present",
  recent_count %in% 0:1000 ~ "<1,000",
  recent_count %in% 1000:5000 ~ "1,000 to 5,000",
  recent_count %in% 5000:10000 ~ "5,000 to 10,000",
  recent_count %in% 10000:50000 ~ "10,000 to 50,000",
  recent_count %in% 50000:100000 ~ "50,000 to 100,000",
  recent_count %in% 100000:500000 ~ "100,000 to 500,000",
  recent_count %in% 500000:1000000 ~ "500,000 to 1,000,000",
  recent_count > 1000000 ~ "> 1,000,000"))
summary(as.factor(col.atpu$count_f))
levels(as.factor(col.atpu$count_f))
# make it an ordered factor
col.atpu$count_f <- factor(col.atpu$count_f, levels=c("present","<1,000",
                                                              "1,000 to 5,000",
                                                              "5,000 to 10,000",
                                                              "10,000 to 50,000",
                                                              "50,000 to 100,000",
                                                              "100,000 to 500,000",
                                                              "500,000 to 1,000,000",
                                                              "> 1,000,000"))
summary(col.atpu$count_f)

mapview(col.atpu)

# bring in all LESP colonies ----
col.lesp <- read.csv("input/lesp_colony_coordinates_for.maps.csv")

col.lesp <- st_as_sf(col.lesp,coords = c("Longitude", "Latitude"),crs=4326) 
mapview(col.lesp) + mapview(col.lesp.mon, color="red")

head(col.lesp)
# remove the colonies with zero individuals
col.lesp <- subset(col.lesp, Individuals != 0)
# remove the commas from the counts
col.lesp$Individuals <- gsub(",", "", col.lesp$Individuals)

# create a numeric recent_count from Individuals which will change "present" to NA
col.lesp$recent_count <- as.numeric(col.lesp$Individuals)

hist(col.lesp$recent_count, breaks=100)

# now we need to replace add all of the monitored colonies with the most recent count and more accurate position
col.lesp$sp <- "LESP"
col.lesp <- col.lesp[,c("sp","Colony","recent_count","geometry")]
col.lesp.mon <- col.lesp.mon[,c("sp","Colony","recent_count","geometry")]
# unfortunately colony naming differs so we'll need to do add all of the monitored colonies then delete the duplicates
col.lesp <- rbind(col.lesp,col.lesp.mon)
row.names(col.lesp) <- c(1:nrow(col.lesp))
col.lesp <- col.lesp[-c(28,24,23,17,18,19,21,22,65,67,68),] 

# now make a factor where we can give NAs a "present" factor
col.lesp <- col.lesp %>% mutate(count_f = case_when(
  is.na(recent_count) ~ "present",
  recent_count %in% 0:1000 ~ "<1,000",
  recent_count %in% 1000:5000 ~ "1,000 to 5,000",
  recent_count %in% 5000:10000 ~ "5,000 to 10,000",
  recent_count %in% 10000:50000 ~ "10,000 to 50,000",
  recent_count %in% 50000:100000 ~ "50,000 to 100,000",
  recent_count %in% 100000:500000 ~ "100,000 to 500,000",
  recent_count %in% 500000:1000000 ~ "500,000 to 1,000,000",
  recent_count > 1000000 ~ "> 1,000,000"))
summary(as.factor(col.lesp$count_f))
levels(as.factor(col.lesp$count_f))
# make it an ordered factor
col.lesp$count_f <- factor(col.lesp$count_f, levels=c("present","<1,000",
                                                      "1,000 to 5,000",
                                                      "5,000 to 10,000",
                                                      "10,000 to 50,000",
                                                      "50,000 to 100,000",
                                                      "100,000 to 500,000",
                                                      "500,000 to 1,000,000",
                                                      "> 1,000,000"))
summary(col.lesp$count_f)

mapview(col.lesp)

# perfect

# basemap 2 ----

# larger extent

bbox <- st_bbox(c(xmin = -68, xmax = -51, ymin = 42.5, ymax = 58), crs = 4326)

# downloaded from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
coastline <- st_read('data/spatial_data/Shoreline_Coastline_GSHHG/GSHHS_shp/f/GSHHS_f_L1.shp') %>% st_crop(bbox)

basemap2 <- ggplot(data=coastline) + geom_sf(fill="gray95", size=0.001) + 
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 57), expand=FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"), 
        panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1), 
        panel.background = element_rect(fill = "aliceblue")) 
basemap2 

# basemap <- basemap +
#   annotate(geom="text", x = -63.2, y=45.1, label = "Nova Scotia", size=4, colour="grey30") +
#   annotate(geom="text", x = -66.3, y=46, label = "New Brunswick", size=4, colour="grey30") +
#   annotate(geom="text", x = -56, y=48.5, label = "Newfoundland", size=4, colour="grey30") +
#   coord_sf(xlim = c(-68, -51), ylim = c(42.5, 55), expand=FALSE) 
# basemap

# map all colonies both species ----
# put them together in one dataframe

col.all <- rbind(col.lesp[,c("Colony","count_f","sp","geometry")], col.atpu[,c("Colony","count_f","sp","geometry")])

col_pal <- viridis::plasma(10)
col_pal <- rev(col_pal)
col_pal <- col_pal[c(2:10)]

p1 <- basemap2 + geom_sf(data=subset(col.all, sp=="LESP"), aes(geometry=geometry, color=count_f, size=count_f), alpha=0.8, shape=20) + 
  scale_color_manual(values=col_pal, drop = FALSE) +
  scale_size_manual(values=c(3,3.5,4,4.5,5,5.5,6,6.5,7), guide="none") +
  geom_sf(data=subset(col, sp=="LESP"), aes(geometry=geometry), shape=1, size=6) + 
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 57), expand=FALSE) +
  labs(color = "Most Recent Count") +
  guides(color = guide_legend(override.aes = list(size = c(3,3.5,4,4.5,5,5.5,6,6.5,7)))) +
  theme(legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
        legend.box="vertical",
        legend.position="none", 
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
p1

p2 <- basemap2 + geom_sf(data=subset(col.all, sp=="ATPU"), aes(geometry=geometry, color=count_f, size=count_f), alpha=0.8, shape=20) + 
  scale_color_manual(values=col_pal, drop = FALSE) +
  scale_size_manual(values=c(3,3.5,4,4.5,5,5.5,6,6.5,7), guide="none") +
  geom_sf(data=subset(col, sp=="ATPU"), aes(geometry=geometry), shape=1, size=6) + 
  coord_sf(xlim = c(-68, -51), ylim = c(42.5, 57), expand=FALSE) +
  labs(color = "Most Recent Count") +
  guides(color = guide_legend(override.aes = list(size = c(3,3.5,4,4.5,5,5.5,6,6.5,7)))) +
  theme(legend.spacing.y = unit(1, "mm"), legend.direction="vertical",
        legend.box="vertical",
        #legend.position="none", 
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.box.margin = margin(0.4,0.4,0.4,0.4,"cm"),
        legend.background = element_rect(color = "white"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
p2

p1 <- p1 + annotate("text", label = "A", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)
p2 <- p2 + annotate("text", label = "B", x=-Inf, y=Inf, vjust=2, hjust=-1, size=8)

p <- p1|p2
p

ggsave(filename="output/figures/maps/both.species.all.colonies.w.trend.colonies.png", plot=p, 
       device="png", dpi=300, units="cm", width=30, height=15)

# maybe we should try adding fun clipart... 
library(png)
library(grid)
atpu <- readPNG(here("input/ATPU_clipart.png"), native=TRUE)
atpu <- rasterGrob(atpu, interpolate=TRUE)
lesp <- readPNG(here("input/LESP_clipart.png"), native=TRUE)
lesp <- rasterGrob(lesp, interpolate=TRUE)

p2.w.pic <- p2 + annotation_custom(atpu, xmin=-67.5, xmax=-61, ymin = 51, ymax=55) 
p1.w.pic <- p1 + annotation_custom(lesp, xmin=-67, xmax=-62, ymin = 51, ymax=56) 

p.w.pic <- p1.w.pic|p2.w.pic
p.w.pic

ggsave(filename="output/figures/maps/both.species.all.colonies.w.trend.colonies.clipart.png", plot=p.w.pic, 
       device="png", dpi=300, units="cm", width=30, height=15)
