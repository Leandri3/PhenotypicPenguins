##---------------------------------------------------------------------------------------------------------
# This code: - Creates species-specific utilisations distribtuions (UDs)
# for:  Phenotypic plasticity chapter (1)
#       Gentoos
#       Each island separately
# In order to characterize an environmental space
# Calculates size of area used by species at both islands
#
# Plots the UDs against bathy background
#
# By Leandri de Kock, modified code from Chris Oosthuizen
# 07 August 2022
#--------------------------------------------------------------------------------------------------
# Utilization probabilities
# A utilization distribution is a surface of utilization probabilities. 
# https://www.r-bloggers.com/2016/05/adehabitathr-visualization/

library(dplyr) # load before adehabitarHR otherwise 'select' and 'id' get masked
library(adehabitatHR)
library(ggplot2)
library(sf)
library(viridis)
library(ggspatial)
library(raster)

source("./scripts/various/theme_rr.R")

# read in data
# gentoo_Nelson
gent1 <-list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
          pattern = '.*divestats_temp_gentoo_Nelson.*\\.rds$',
          full.names = TRUE)

print(gent1)

# make empty list
gentslist = list()

for (i in gent1) {             # to run on all individuals selected in divefiles
# Import divetemps data and save in a list
gentslist[[i]] = readRDS(i)
}

gent1 = as_tibble(bind_rows(gentslist))
str(gent1)

# make sure date.time is set to GMT (or UTC)
attr(gent1$begdesc, "tzone") # check time zone

gent1 = gent1 %>%
  mutate(species = "gentoo")%>%
  mutate(island = "nelson")

# create a new column for each breeding stage
gent1$stage = substr(gent1$ID , 8, nchar(gent1$ID)-6)
unique(gent1$stage)

# Change Gua to Bro
#--------------------------------------------------------------
# rename Gua = Bro
gent1$stage = if_else(gent1$stage == "Gua", "Bro", gent1$stage )
unique(gent1$stage)

# Look at sample size (uplinks) per stage # bro should be 96906
gent1 %>%
  group_by(stage) %>%
  dplyr::summarise(n= n()) # Bro = 90433
# works!
#----------------------------------------------------------------

plot(gent1$lon.x, gent1$lat.y, pch = ".", cex = 2)

# read in data
# gentoo_Kopaitic
gent2 <-list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
          pattern = '.*divestats_temp_gentoo_Kopaitic.*\\.rds$',
          full.names = TRUE)

print(gent2)

# make empty list
gentslist = list()

for (i in gent2) {             # to run on all individuals selected in divefiles
# Import divetemps data and save in a list
gentslist[[i]] = readRDS(i)
}

gent2 = as_tibble(bind_rows(gentslist))
str(gent2)

# make sure date.time is set to GMT (or UTC)
attr(gent2$begdesc, "tzone") # check time zone

gent2 = gent2 %>%
  mutate(species = "gentoo") %>%
  mutate(island = "kopaitic")


# create a new column for each breeding stage
gent2$stage = substr(gent2$ID , 10, nchar(gent2$ID)-6)
unique(gent2$stage)

plot(gent2$lon.x, gent2$lat.y, pch = ".", cex = 2)

dat = bind_rows(gent1,gent2)

dat = dat%>%
  dplyr::rename(id = ID) %>%
  dplyr::group_by(id) %>%
  arrange(begdesc) %>% # put the data in date order
  distinct(begdesc, .keep_all=TRUE) %>%   # remove duplicate date.time records for individuals
  dplyr::ungroup()
  
# create a new column for each breeding stage
#dat$stage = substr(dat$id , 8, nchar(dat$id)-6)
#unique(dat$stage)
# Look at sample size (uplinks) per stage
dat %>%
  group_by(stage) %>%
  dplyr::summarise(n= n()) # Bro = 111194, Cre = 34875, Inc = 81234

# Look at sample size (individuals) per stage
dat %>%
  group_by(stage) %>%
   dplyr::summarise(n= n_distinct(id))

# Look at sample size (individuals) per stage
dat %>%
  group_by(island) %>%
   dplyr::summarise(n= n_distinct(id))

str(dat)

# make a copy
divelocs.sp = dat %>%
  as_tibble(dat) %>%
  dplyr::select(begdesc, id, lon.x, lat.y, stage, species, island)%>%
  arrange(lat.y)

# remove rows with NA values in lat column 
divelocs.sp = divelocs.sp[complete.cases(divelocs.sp$lat.y),] # there were two rows with NaN values

str(divelocs.sp)

# Change Bro and Cre to chick-rearing phase:
#--------------------------------------------------------------
# rename 
divelocs.sp$phase = if_else(divelocs.sp$stage == "Bro", "Chick", divelocs.sp$stage)
divelocs.sp$phase = if_else(divelocs.sp$phase == "Cre", "Chick", divelocs.sp$phase)

unique(divelocs.sp$phase)
unique(divelocs.sp$stage)

# Look at sample size (uplinks) per stage # bro should be 96906
divelocs.sp %>%
  group_by(stage) %>%
  dplyr::summarise(n= n())

# chick should be:146067
divelocs.sp %>%
  group_by(phase) %>%
  dplyr::summarise(n= n()) # chick = 146067
# works!
#----------------------------------------------------------------


ggplot(data = divelocs.sp, aes(x = lon.x, y = lat.y, colour = island)) +
 geom_point()

# Quick visual checK:
sf_locs <- sf::st_as_sf(divelocs.sp, coords = c("lon.x","lat.y")) %>% 
  sf::st_set_crs(" +proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs ")

fig.NB = ggplot() + 
 annotation_map_tile(zoom = 1, cachedir = system.file("rosm.cache", package = "ggspatial")) +  # if using downloaded files
  layer_spatial(sf_locs, size = 0.75, shape = 16, alpha = 0.2, aes(color = phase)) +  # specified stage as ID above
  #scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  ggtitle("Breeding Gentoos", 
          subtitle = "Filtered tracks for analysis; colours = phase")

fig.NB 

#--------------------------------------------------------
## Utilization distributions
#--------------------------------------------------------

# Define the projections
#Set new projection to https://epsg.io/32721
utm.prj = " +proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs "   # Chris UTM King George
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#To assign a known CRS to spatial data:
utm.coord = SpatialPoints(cbind(divelocs.sp$lon.x, divelocs.sp$lat.y), proj4string=CRS(utm.prj))
utm.coord

## To transform from one CRS to another: from utm to lonlat
#wgs.coord <- st_transform(utm.coord, CRS(wgs84))
#wgs.coord
#
#sf::st_as_sf(coords=c("lon", "lat"), crs=4326) %>% # gpsdata (lat lon) is in wgs84 (4326)
#  sf::st_transform(32721) # transform lat lon data to UTM zone 21S https://epsg.io/32721 (crawl requires UTM)

# divelocs.sp$lon.x <- utm.coord$coords.x1
# divelocs.sp$lat.y <- utm.coord$coords.x2
head(divelocs.sp)

# Grid for UD
lms <- c(min(divelocs.sp$lon.x, na.rm = T) - 15000,
         max(divelocs.sp$lon.x, na.rm = T) + 15000,
         min(divelocs.sp$lat.y, na.rm = T) - 15000,
         max(divelocs.sp$lat.y, na.rm = T) + 15000)
rt <- raster(ext = extent(lms), crs = CRS(utm.prj), res = 250)
rt.sp <- as(rt, "SpatialPixelsDataFrame")

# group by breeding "phase" x "island"
divelocs.sp$group <-paste0(divelocs.sp$phase, " ", divelocs.sp$island)
head(divelocs.sp)
divelocs.sp <- divelocs.sp[ , c("lon.x", "lat.y", "group")]
unique(divelocs.sp$group)
#coordinates(divelocs.sp) <- c("lon.x", "lat.y")
#proj4string(divelocs.sp) <- utm.prj

# create a SpatialPointsDataFrame
coordinates(divelocs.sp) = ~lon.x+lat.y 	
proj4string(divelocs.sp) <-CRS(utm.prj)
head(coordinates(divelocs.sp)) 

plot(divelocs.sp)


kud <- kernelUD(divelocs.sp, h = 7000, grid = rt.sp) # mean h$scaleARS

#kud <- kernelUD(divelocs.sp, h = 5000, grid = rt.sp)
# from track2KBA FTP analysis h$scaleARS = 5 km: 
                    # h$mag = 2.34 km, 
                    #h$med.max.dist = 10.41 km, 
                    #h$ref = 1.95 (underestimates)


image(kud)

# Create rasters of each breeding stage
vud <- getvolumeUD(kud)
head(vud)

# From RR_GP
# Get contours
#all.contour <- getverticeshr(kud, percent = 95)
#all.contour.50 <- getverticeshr(kud, percent = 50)

vud.gt.nel.inc.raster <- raster(as(vud$`Inc nelson`, "SpatialPixelsDataFrame"))
vud.gt.kop.inc.raster <- raster(as(vud$`Inc kopaitic`, "SpatialPixelsDataFrame"))

vud.gt.nel.chick.raster <- raster(as(vud$`Chick nelson`, "SpatialPixelsDataFrame"))
vud.gt.kop.chick.raster <- raster(as(vud$`Chick kopaitic`, "SpatialPixelsDataFrame"))

#vud.gt.nel.cre.raster <- raster(as(vud$`Cre nelson`, "SpatialPixelsDataFrame"))
#vud.gt.kop.cre.raster <- raster(as(vud$`Cre kopaitic`, "SpatialPixelsDataFrame"))

#plot(vud.kop.gt.inc.raster)
#plot(vud.kop.gt.bro.raster)
#plot(vud.kop.gt.cre.raster)

# Write the rasters
writeRaster(vud.gt.nel.inc.raster, "./outputs/utilisation distributions/phenotypic/vud.gt.nel.inc.grd", format = "raster", overwrite = T)
writeRaster(vud.gt.kop.inc.raster, "./outputs/utilisation distributions/phenotypic/vud.gt.kop.inc.grd", format = "raster", overwrite = T)

writeRaster(vud.gt.nel.chick.raster, "./outputs/utilisation distributions/phenotypic/vud.gt.nel.chick.grd", format = "raster", overwrite = T)
writeRaster(vud.gt.kop.chick.raster, "./outputs/utilisation distributions/phenotypic/vud.gt.kop.chick.grd", format = "raster", overwrite = T)

#writeRaster(vud.gt.nel.cre.raster, "vud.gt.nel.cre.grd", format = "raster", overwrite = T)
#writeRaster(vud.gt.kop.cre.raster, "vud.gt.kop.cre.grd", format = "raster", overwrite = T)

# Import rasters
# vud.incubation <- raster("./output/vud.kop.inc.grd")
# vud.brood <- raster("./output/vud.kop.bro.grd")
# vud.creche <- raster("./output/vud.kop.cre.grd")


# calculate the size of the areas 
 kernel.area(kud, percent = c(50, 95),
             unin = c("m"),
             unout = c("km2"), standardize = FALSE)

#plot(vud.gt.kop.inc.raster, col = terrain.colors(100))
#plot(vud.kop.gt.bro.raster, col = terrain.colors(100))
#plot(vud.kop.gt.cre.raster, col = terrain.colors(100))

vud.raster <-stack(vud.gt.nel.inc.raster,
                   vud.gt.kop.inc.raster,
                   vud.gt.nel.chick.raster,
                   vud.gt.kop.chick.raster)#,
                   #vud.gt.nel.cre.raster,
                   #vud.gt.kop.cre.raster)
names(vud.raster) <-c("nel incubation", "kop incubation",  
                      "nel chick", "kop chick")#,
                      #"nel creche","kop creche")

plot(vud.raster, col = terrain.colors(100))

#### Maps of utilization distributions # can do by breeding stage
# nelson
# chick
chick.nel <- rasterToPoints(vud.gt.nel.chick.raster)
chick.nel <- data.frame(chick.nel)
colnames(chick.nel) <- c("lon", "lat", "val")
chick.nel$ph <- "chick"
chick.nel$isl <- "nelson"

# incubation
incubation.nel <- rasterToPoints(vud.gt.nel.inc.raster)
incubation.nel <- data.frame(incubation.nel)
colnames(incubation.nel) <- c("lon", "lat", "val")
incubation.nel$ph <- "incubation"
incubation.nel$isl <- "nelson"

# creche
#creche.nel <- rasterToPoints(vud.gt.nel.cre.raster)
#creche.nel <- data.frame(creche.nel)
#colnames(creche.nel) <- c("lon", "lat", "val")
#creche.nel$ph <- "creche"
#creche.nel$isl <- "nelson"

# kopaitic
# brood
chick.kop <- rasterToPoints(vud.gt.kop.chick.raster)
chick.kop <- data.frame(chick.kop)
colnames(chick.kop) <- c("lon", "lat", "val")
chick.kop$ph <- "chick"
chick.kop$isl <- "kopaitic"

# incubation
incubation.kop <- rasterToPoints(vud.gt.kop.inc.raster)
incubation.kop <- data.frame(incubation.kop)
colnames(incubation.kop) <- c("lon", "lat", "val")
incubation.kop$ph <- "incubation"
incubation.kop$isl <- "kopaitic"

# creche
#creche.kop <- rasterToPoints(vud.gt.kop.cre.raster)
#creche.kop <- data.frame(creche.kop)
#colnames(creche.kop) <- c("lon", "lat", "val")
#creche.kop$ph <- "creche"
#creche.kop$isl <- "kopaitic"

kern<-rbind(incubation.nel, incubation.kop,
            chick.nel, chick.kop)#,
            #creche.nel, creche.kop)
kern <- kern[kern$val <= 95, ]



