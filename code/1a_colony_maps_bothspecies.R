

# map modelled LESP and ATPU col.lesponies together

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

# bring in ATPU ----
col.atpu <- read.csv("input/ATPU_recent.counts.coords.csv")

sf_use_s2(FALSE)

col.atpu <- st_as_sf(col.atpu,coords = c("lon", "lat"),crs=4326) 
mapview(col.atpu)
hist(col.atpu$recent_count, breaks=50)
col.atpu <- col.atpu %>% mutate(count_f = case_when(
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
col.atpu$count_f <- factor(col.atpu$count_f, levels=c("<1,000",
                                                      "1,000 to 5,000",
                                                      "5,000 to 10,000",
                                                      "10,000 to 50,000",
                                                      "50,000 to 100,000",
                                                      "100,000 to 500,000",
                                                      "500,000 to 1,000,000",
                                                      "> 1,000,000"))
summary(col.atpu$count_f)

# bring in LESP ----
col.lesp <- read.csv("input/LESP_recent.counts.coords.csv")

col.lesp <- st_as_sf(col.lesp,coords = c("lon", "lat"),crs=4326) 
mapview(col.lesp)

col.lesp <- col.lesp %>% mutate(count_f = case_when(
  recent_count %in% 0:1000 ~ "<1,000",
  recent_count %in% 1000:5000 ~ "1,000 to 5,000",
  recent_count %in% 5000:10000 ~ "5,000 to 10,000",
  recent_count %in% 10000:50000 ~ "10,000 to 50,000",
  recent_count %in% 50000:100000 ~ "50,000 to 100,000",
  recent_count %in% 100000:500000 ~ "100,000 to 500,000",
  recent_count %in% 500000:1000000 ~ "500,000 to 1,000,000",
  recent_count > 1000000 ~ "> 1,000,000"))
summary(as.factor(col.lesp$count_f))

col.lesp$count_f <- factor(col.lesp$count_f, levels=c("<1,000",
                                                      "1,000 to 5,000",
                                                      "5,000 to 10,000",
                                                      "10,000 to 50,000",
                                                      "50,000 to 100,000",
                                                      "100,000 to 500,000",
                                                      "500,000 to 1,000,000",
                                                      "> 1,000,000"))
summary(col.lesp$count_f)

# basemap ----

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

basemap + geom_point(data=col.lesp, aes(geometry=geometry, color=count_f), size=3, stat="sf_coordinates") + 
  scale_color_manual(values=col_pal)

basemap + geom_point(data=col.atpu, aes(geometry=geometry, color=count_f), size=3, stat="sf_coordinates") + 
  scale_color_manual(values=col_pal)

# put them together on one map

col.lesp$sp <- "LESP"
col.atpu$sp <- "ATPU"
col <- rbind(col.lesp, col.atpu)

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