# Session info: R version 4.0.1 (2020-06-06) See Things Now
#### 1. Build urban data ####

library(sf) # v0.9.2
library(rgeos) # v0.5.3
library(raster) # v3.1.5
library(rgdal) # v1.5.10
library(tidyverse)
library(viridis)

#### Load data ####
gdata.amph <- read.csv("amph_data.csv", head = TRUE) # (final dataset)

# Spatial dataframe (point feature)
sites <- st_as_sf(gdata.amph, coords = c("lon", "lat"), crs = 4326)

# Human population density
# https://neo.sci.gsfc.nasa.gov/view.php?datasetId=SEDAC_POP (year 2000)
popden <- read.asciigrid("PopDen1.asc")
popden <- raster(popden)

# Human Footprint Index
# https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic (1995-2004)
hfi <- raster("hfp_N_America_grid/hfp_n_amer/w001001.adf")

# % impervious surface
## 1 km resolution
# https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1 (2010)
# 1km resolution
gmis <- raster("gmis_impervious_surface/1km_geographic_polygon/gmis_impervious_surface_percentage_geographic_1000m.tif")
NAvalue(gmis) <- 255 ## 255 is NoData
gmis[gmis == 200] <- 0 ## '200' = not in Global Human Built-up And Settlement Extent; can be set to 0

## 30 m resolution
usa30 <- raster("cusa_30m.tif")
can30 <- raster("can_30m.tif")
alaska30 <- raster('USAW3_gmis_impervious_surface_percentage_geographic_30m.tif')

# Roads
# https://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1 (1980-2010)
roads1 <- read_sf("groads/groads-v1-americas-shp/gROADS-v1-americas.shp")
roads <- roads1 %>% 
  filter(EXS == 1) ## Filter EXS = 1 for "existence = definite" 2 = doubtful, 0 = unspecified

# Urban areas
# https://hifld-geoplatform.opendata.arcgis.com/datasets/d7f5ddbc05324e058500327a13abbbfb
usa_uas <- read_sf("tl_2019_us_uac10/tl_2019_us_uac10.shp") # USA
# https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2011-eng.cfm
popctrs <- read_sf("lpc_000a16a_e.shp") # Canada

## Merge urban areas:
usajoin <- usa_uas %>% 
          select(c(GEOID10, NAME10, geometry)) %>% # remove and rename columns to join
          rename(ID = GEOID10, name = NAME10)

canjoin <- popctrs %>% 
          select(c(PCUID, PCNAME, geometry)) %>% # remove and rename columns to join
          rename(ID = PCUID, name = PCNAME)

canjoin <- st_transform(canjoin, crs = st_crs(usajoin)) # make projections the same

uas <- rbind(usajoin, canjoin) # merge shapefiles

#### Taxonomy ####
amph$order <- ifelse(amph$species %in% c("Ambystoma_barbouri", "Ambystoma_maculatum",
                                         "Desmognathus_fuscus", "Dicamptodon_aterrimus",
                                         "Dicamptodon_copei", "Ensatina_eschscholtzii",
                                         "Hydromantes_brunus", "Hydromantes_platycephalus",
                                         "Plethodon_albagula", "Plethodon_cinereus", 
                                         "Taricha_granulosa"), "Caudata", "Anura")

#### Urbanization variables ####
# Extract values per site within 1km, 5km, 10km, 15km
buffers <- list(1000, 5000, 10000, 15000)

## 1) population density
popdenmulti <- lapply(buffers, function(x) raster::extract(popden, sites, fun=mean, buffer=x, na.rm=TRUE, df=TRUE))
popdendf <- as.data.frame(popdenmulti)
popdendf <- popdendf[,-c(1,3,5,7)] # remove 'ID' columns; keep data only
names(popdendf) <- c("pd_1_km", "pd_5_km", "pd_10_km", "pd_15_km")

## 2) Human Footprint Index
hfimulti <- lapply(buffers, function(x) raster::extract(hfi, sites, fun=mean, buffer=x, na.rm=TRUE, df=TRUE))
hfidf <- as.data.frame(hfimulti)
hfidf <- hfidf[,-c(1,3,5,7)]
names(hfidf) <- c("hfi_1_km", "hfi_5_km", "hfi_10_km", "hfi_15_km")

## 3) % Impervious surface
impmulti <- lapply(buffers, function(x) raster::extract(gmis, sites, fun=mean, buffer=x, na.rm=TRUE, df=TRUE))
impdf <- as.data.frame(impmulti)
impdf <- impdf[,-c(1,3,5,7)]
names(impdf) <- c("imp_1_km", "imp_5_km", "imp_10_km", "imp_15_km")

amph <- Reduce(function(x,y) cbind(x,y), list(gdata.amph, popdendf, hfidf, impdf))

## 4) % Impervious surface (30 m)
ex.usa <- raster::extract(usa30, sites, fun=mean, buffer=0, na.rm=TRUE, df=TRUE)
ex.can <- raster::extract(can30, sites, fun=mean, buffer=0, na.rm=TRUE, df=TRUE)
ex.ak <-  raster::extract(alaska30, sites, fun=mean, buffer=0, na.rm=TRUE, df=TRUE)

sites$imp_30m <- ifelse(is.na(ex.usa$cusa_30m), ex.can$can_30m, ex.usa$cusa_30m) # replace NA Canada sites
sites$imp_30m <- ifelse(is.na(sites$imp_30m), ex.ak$USAW3_gmis_impervious_surface_percentage_geographic_30m, 
                        sites$imp_30m) # replace NA Alaska sites

## % impervious surface is from 0-100
## values that are 200 are *not* 'HBASE' (Global Human Built-up And Settlement Extent)
## Thus we can make them 0% impervious
sites$imp_30m <- ifelse(sites$imp_30m==200, 0, sites$imp_30m)
amph$imp_30m <- sites$imp_30m

# 5) Road length
## Reproject into planar projection
AEA <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
sites_AEA <- st_transform(sites, AEA)
roads_AEA <- st_transform(roads, AEA)

roadlist <- lapply(buffers, function(x) st_buffer(sites_AEA, dist = x)) 
roadmulti <- lapply(roadlist, function(x) st_intersection(x, roads_AEA)) # clip roads within sites
roadmulti_len <- mapply(cbind, roadmulti, len = lapply(roadmulti, st_length), SIMPLIFY = FALSE) # 1 row per road (multiple rows per pop)


roadmultidf <- lapply(roadmulti_len, function(x) as.data.frame(x) %>% 
                        group_by(pop) %>%
                        summarise(road_length = (sum(len))/1000) %>% # sum total road length within each pop; convert to km
                        mutate(road_length = as.numeric(road_length))
)

roadbufs <- reduce(roadmultidf, full_join, by = "pop")
roadbufs[is.na(roadbufs)] <- 0
names(roadbufs)[c(2:5)] <- c("road_1_km", "road_5_km", "road_10_km", "road_15_km")

amphrb <- merge(amph, roadbufs, by = "pop", all = TRUE)

# 6) Urban areas
## sites within 10km of urban area
uas_AEA <- st_transform(uas, AEA)
s10km <- st_is_within_distance(sites_AEA, uas_AEA, dist = 10000)
urban <- ifelse(lengths(s10km)>0, 1, 0)

amphrb$urban <- urban

#### Final dataset ####
#write.csv(amphrb, "amph_data.csv", row.names = F, quote = F)
