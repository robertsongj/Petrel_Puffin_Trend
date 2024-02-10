# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse',
             'readxl',
             'RColorBrewer',
             'viridis',
             'jagsUI',
             'mgcv',
             'ggrepel',
             'scales',
             'sf',
             'ggspatial')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

theme_set(theme_bw())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

dirname <- thisPath()
setwd(dirname)
setwd("../")
`%!in%` <- Negate(`%in%`)


# ------------------------------------------------
# Load colony coordinates - prepare to create map
# ------------------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-60 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

coords = read_xlsx("data/LESP_colony_coordinates.xlsx", sheet = 1) %>%
  subset(!is.na(Latitude)&!is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),crs=4326, remove = FALSE) %>%
  subset(Country %in% c("Canada")) %>%
  st_transform(AEA_proj) %>%
  dplyr::rename(Count = `Estimated no. of mature individuals`)

ggplot(data = coords) + geom_sf()

bbox <- st_bbox(coords) %>% st_as_sfc() %>% st_buffer(10000)

coords$Count_numeric <- as.numeric(coords$Count)

# ------------------------------------------------
# PLOT COLONY LOCATIONS
# ------------------------------------------------

sf_use_s2(FALSE)

# downloaded from https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
coastline <- st_read('data/spatial_data/Shoreline_Coastline_GSHHG/GSHHS_shp/f/GSHHS_f_L1.shp')
coastline <- coastline %>%
  st_make_valid() %>%
  st_intersection(., st_transform(st_buffer(bbox,100000) , st_crs(.))) %>%
  st_transform(AEA_proj)


xlim <- range(as.data.frame(st_coordinates(bbox))$X)
xlim[2] <- xlim[2] + 200000
ylim <- range(as.data.frame(st_coordinates(bbox))$Y)
ylim[1] <- ylim[1]-100000

coords$count_levs <- cut(coords$Count_numeric,
                         breaks = c(-1,1000,10000,100000,1000000,1e+10),
                         labels = c("< 1,000",
                                    "1,000 to 10,000",
                                    "10,000 to 100,000",
                                    "100,000 to 1,000,000",
                                    "> 1,000,000"))



# https://icolorpalette.com/color/pale-blue
colony_map <- ggplot() + 
  geom_sf(data = coastline, col = "gray90", fill = "white") + 
  
  # Colonies with estimated abundances (plot these in layers, so that small colonies appear on top)
  geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[5]), aes(size = count_levs,fill= count_levs), col = "black", shape = 21, stroke = 0.2, alpha = 0.9)+
  #geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[5]), aes(size = count_levs,fill= count_levs), col = "black", shape = 10, stroke = 0.2)+
  
  geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[4]), aes(size = count_levs,fill= count_levs), col = "black", shape = 21, stroke = 0.2, alpha = 0.9)+
  #geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[4]), aes(size = count_levs,fill= count_levs),  col = "black", shape = 10, stroke = 0.2)+
  
  geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[3]), aes(size = count_levs,fill= count_levs), col = "black", shape = 21, stroke = 0.2, alpha = 0.9)+
  #geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[3]), aes(size = count_levs,fill= count_levs), col = "black", shape = 10, stroke = 0.2)+
  
  geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[2]), aes(size = count_levs,fill= count_levs), col = "black", shape = 21, stroke = 0.2, alpha = 0.9)+
  #geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[2]), aes(size = count_levs,fill= count_levs), col = "black", shape = 10, stroke = 0.2)+
  
  geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[1]), aes(size = count_levs,fill= count_levs), col = "black", shape = 21, stroke = 0.2, alpha = 0.9)+
  #geom_sf(data = subset(coords,count_levs == levels(coords$count_levs)[1]), aes(size = count_levs,fill= count_levs), col = "black", shape = 10, stroke = 0.2)+
  
  
  # Colonies with no estimated abundances
  geom_sf(data = subset(coords, Count == "Present"), col = "black",fill = "gray20", 
          size = 1, shape = 10, stroke = 0.2)+
  
  coord_sf(xlim = xlim, ylim = ylim)+
  scale_size_manual(name = "Colony Size",
                    values = c(1,2,3,5,8),
                    breaks = levels(coords$count_levs)
                    
  )+
  scale_fill_manual(name = "Colony Size",
                    #values = rev(magma(8)[2:6])
                    values = RColorBrewer::brewer.pal(5,"YlGnBu"),
                    breaks = levels(coords$count_levs)
  )+
  
  annotation_scale(style = "ticks")+
  annotation_north_arrow(which_north = "true",
                         location = "tr",
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering(text_col = 'gray20',
                                                                line_col = 'gray20',
                                                                fill = 'gray80'))+
  theme(panel.background = element_rect(fill = "#e3fefe",
                                        colour = "#e3fefe",
                                        size = 0.5, linetype = "solid"),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 0.2, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position=c(0.98,0.02),
        legend.justification = c("right","bottom"),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 6))+
  xlab("")+ylab("")


# Add map labels

map_labels <- rbind(data.frame(Latitude = 48,Longitude = -62,Label = "Gulf of St Lawrence", angle = 30, size = 3),
                    data.frame(Latitude = 50,Longitude = -51,Label = "North Atlantic\nOcean", angle = 10, size = 4)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),crs=4326, remove = FALSE)

for (i in 1:nrow(map_labels)){
  colony_map <- colony_map +
    geom_sf_text(data = map_labels[i,], aes(label = Label, angle = angle),
                 col = "gray30", alpha = 0.3, size = map_labels$size[i])
  
}

print(colony_map)

png(paste0("output/figures/maps/colony_map.png"), width=6, height=5.5, units="in", res=600, type="cairo")
print(colony_map)
dev.off()

# ------------------------------------------------
# Read in data, check out relationship between CV and colony size
# ------------------------------------------------

dat = read_xlsx("data/LESP_data_for_analysis.xlsx", sheet = 1) %>%
  dplyr::rename(Count = `Mature individuals`)

ggplot(data = subset(dat, Count > 0 & !is.na(CV)), aes(x = Count, y = CV))+
  geom_point(size=3)+
  scale_y_continuous(trans = "log10", labels = comma, name = "Coefficient of Variation")+
  scale_x_continuous(trans = "log10", labels = comma)+
  ggtitle("Relationship between colony abundance and survey error")+
  theme_bw()

# Reasonably linear relationship between log (CV) and log (count)
# - will hopefully allow error to be estimated for surveys that are missing that information