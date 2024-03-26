## Luke Ozsanlav-Harris

## Created: 28/03/2023

## Create various maps for publication to show Islay location, habitat and roads
## Create Figure 1: map of Islay with RAMSAR and GBG roost sites marked and entire UK as inlay
## Create sup figure: map of Islay with all habitat shown
## Create sup figure: map of Islay with only agricultural land, saltmarsh and roads shown

## Load required packages
pacman::p_load(tidyverse, data.table, sf, raster, ggspatial, cowplot, MetBrewer, ggnewscale, terra)




##-------------------------------------##
#### 1. Read in all Spatial Data set ####
##-------------------------------------##

## Ramsar outline on Islay
Ram1 <- st_read("Landcover Data/Ramsar Outlines/Gruinart_Layer.shp")
Ram2 <- st_read("Landcover Data/Ramsar Outlines/The_Oa.shp")
st_crs(Ram1); st_crs(Ram2)

## Roads of the UK
Road <- st_read("Landcover Data/Roads/roads.shp")
st_crs(Road)

## Coastline of the UK
UK <- st_read("Landcover Data/High-Res Coastline/UK_Coatline.shp")
st_crs(UK)
## Coastline of Islay
Islay <- st_read("Landcover Data/High-Res Coastline/Islay_Coatline.shp")
st_crs(Islay)

## Read in bounding boxes of UK and Islay
UKBBox <- st_read("Landcover Data/High-Res Coastline/BoundingBoxes/UK_BoundingBox.shp")
IslayBBox <- st_read("Landcover Data/High-Res Coastline/BoundingBoxes/Islay_BoundingBox.shp")
plot(IslayBBox)

## Landcover data for a middle year
Land <- raster("Landcover Data/Islay landcover data/LandCov2018_Islay.grd")
st_crs(Land)


## Create the three GBG roosts sites in lat long
Roost <- data.frame(Loc = c("Gruinart", "Bridgend", "Lagan"), 
                    Lat = c(55.829241, 55.775400, 55.688561), 
                    Long = c(-6.332374, -6.262273, -6.277231))
Roost <- st_as_sf(Roost, coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant")




##-----------------------------------------##
#### 2. Reproject/Crop/Join Spatial Data ####
##-----------------------------------------##

##-- Joining --##

## Join the two RAMSAR areas
Ram <- rbind(Ram1, Ram2)
Ram$id <- c("Gruinart", "Oa")


##-- Crop --##

## crop the road to just the Islay
Road_cr <- st_crop(Road, xmin = -6.6, xmax = -6, ymin = 55.55, ymax = 56)
Road_cr <- st_intersection(IslayBBox, Road_cr)
ggplot() +geom_sf(data = Road_cr, aes(geometry = geometry)) + theme_light()

## Re-project the raster to lat/long
LandR <- projectRaster(Land, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method = "ngb")
IslayBBoxR <- st_transform(IslayBBox, crs = crs(Land))
## crop the raster to exlclude the  bit of Jura
LandI <- mask(Land, IslayBBoxR)
## check this worked
plot(LandI)



##---------------------------------##
#### 3. Reclassify habitat codes ####
##---------------------------------##

## convert the raster to a data frame for plotting
test_spdf <- as(LandI, "SpatialPixelsDataFrame")
Habitat <- as.data.frame(test_spdf)
colnames(Habitat) <- c("value", "x", "y")

## check this has worked
a <- ggplot() +  
  geom_raster(data=Habitat, aes(x=x, y=y, fill=value), alpha=1) + theme_light()
## check that some weird interpolation did not happen
unique(Habitat$value)

## classify the habitats
Habitat <- Habitat %>% mutate(value = as.character(value),
                              Hab = fct_collapse(value,
                                                 #marine = c("13"),
                                                 #woodland = c("1", "2"), # move it in "other" as some woodland classified sites are mostly fields with areas of woodland inside/around (small fields)
                                                 freshwater = c("14"),
                                                 #urban = c("20", "21") # place in other
                                                 #saltmarsh = c("19"),
                                                 coastal_saltmarsh = c("15", "16", "17", "18", "19"),
                                                 other_grassland = c("7","10"),
                                                 improved_grassland = c("4"),
                                                 arable_hort  = c("3"),
                                                 bog = c("11"),
                                                 other = c("1","2","9", "20", "21", "13")))

## Assign meaningful names
Habitat$Hab <- dplyr::recode(Habitat$Hab, arable_hort = "Arable", bog = "Upland Bog", improved_grassland = "Improved grassland", 
                             other_grassland = "Other grassland", coastal_saltmarsh = "Saltmarsh & Coastal", freshwater = "Freshwater", other = "Other")

## check the habitats have the correct dates
unique(Habitat$Hab)



##-----------------------------------##
#### 4. Create Figure 1: Islay Map ####
##-----------------------------------##

## get a three colour palette
x <- met.brewer("Egypt", n = 4)
x[1];x[2];x[3];x[4]

## Create a layer of just improved grassland
## Read in the habitat raster again
Land <- rast("Landcover Data/Islay landcover data/LandCov2018_Islay.grd")
## Re-project the raster and the bounding box to lat/long
LandR <- project(Land, crs(Islay), method = "near")
IslayBBoxR <- st_transform(IslayBBox, crs = crs(Land))
## crop the raster to exclude the  bit of Jura
LandI <- mask(Land, IslayBBoxR)
## check this worked
plot(LandI)
## Now turn all values that are not improved grassland (4) to nA
values(LandI) <- ifelse(values(LandI)==4, 4, NA)
## Now convert the raster to polygons for plotting and then convert to a sf object to plot
LandIPoly <- as.polygons(LandI)
plot(LandIPoly)
LandIPoly <- st_as_sf(LandIPoly)



## create the map of Islay
IslayMap <- ggplot() +
          # Add the coastline
          geom_sf(data = Islay, aes(geometry = geometry)) + 
  

          
          # Add the Habitat data
          geom_sf(data = LandIPoly, aes(geometry = geometry, fill = "#43b284"), colour = "#43b284", alpha = 0.8) +
          scale_fill_manual(name = "Improved Grassland",
                    values =c("#43b284"="#43b284"),
                    labels = c("")) +
          new_scale_fill() +
          
          # add other spatial objects
          geom_sf(data = Ram, aes(geometry = geometry, fill = "#0f7ba2"), alpha = 0.6) +
          scale_fill_manual(name = "Nature Reserves",
                            values =c("#0f7ba2"="#0f7ba2"),
                            labels = c("")) +
          
          geom_sf(data = Road_cr, aes(geometry = geometry, colour = "#fab255"), fill = NA, alpha = 0.9, size = 0.9) +
          scale_colour_manual(name = "Public Roads",
                            values =c("#fab255"="#fab255"),
                            labels = c("")) +
          new_scale_color() +
  
          geom_sf(data = Roost, aes(geometry = geometry, colour = "#dd5129"), size = 4, alpha = 1, shape = 16) +
          scale_colour_manual(name = "GBG Roost Sites",
                              values =c("#dd5129"="#dd5129"),
                              labels = c("")) +

          # set the plot limits to add margin for inset
          coord_sf(xlim = c(-6.89, -6), 
                   ylim = c(55.55, 55.95), crs = 4326, expand = F) +
          
          # add map extras
          annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.2, "in")) +
          annotation_north_arrow(location = "bl", which_north = "true",
                                 pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                                 style = north_arrow_orienteering,
                                 height = unit(1, "cm"), width = unit(0.7, "cm"),) +
  
          # Set the labels
          xlab("Longitude") + ylab("Latitude") +
          # set the theme
          theme_light() +
          theme(legend.position = "bottom", panel.grid.major = element_blank())


## create the map of UK
UKMap <- ggplot() +
      # Add the coastline
      geom_sf(data = UK, aes(geometry = geometry) , alpha = 0.6, size = 0.25) +
      geom_sf(data = Islay, aes(geometry = geometry), fill = "#DC3220", colour = "#DC3220") + 
  
      # set the theme
      theme_light() +
      theme(axis.text = element_blank(), panel.grid.major = element_blank(),
            axis.ticks = element_blank())


## inset the UK map onto the Islay map
plot.with.inset <-
  ggdraw() +
  draw_plot(IslayMap) +
  draw_plot(UKMap, x = 0.1, y = .57, width = .4, height = .4)


## save the plot
ggsave(filename = "Plots/Maps/Figure1_IslayMap.png", 
       plot = plot.with.inset,
       width = 21, 
       height = 15,
       units = "cm",
       dpi = 300)



##--------------------------##
#### 5. Islay Habitat Map ####
##--------------------------##

## transform the roads and Islay coat to the same crs as the habitat map
Road_tr <- st_transform(Road_cr, crs = crs(Land))
Islay_tr <- st_transform(Islay, crs = crs(Land))
  
# Add the coastline
# geom_sf(data = Islay, aes(geometry = geometry)) + 

## Get a colour palette
MetPal <- met.brewer("Archambault", length(unique(Habitat$Hab)))
c(MetPal)
## Changes the order 
MetPalR <-c("#ab3329", "#381a61", "#7c4b73", "#ed968c", "#88a0dc", "#e78429", "#f9d14a")

## Create Islay Habitat Map
IsayHab <- ggplot() +  
            
            ## add roads underneath so that scale bvar is correct
            geom_sf(data = Road_tr, aes(geometry = geometry), fill = NA, alpha = 0.8, size = 1.5, linetype = "longdash") +
  
            # Add the Habitat data
            geom_raster(data=Habitat, aes(x=x, y=y, fill=Hab), alpha=1) + 
            scale_fill_manual(values=MetPalR) +
  
            # add map extras
            annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.2, "in")) +
            annotation_north_arrow(location = "bl", which_north = "true",
                                   pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                                   style = north_arrow_orienteering,
                                   height = unit(1, "cm"), width = unit(0.7, "cm"),) +
            
            # Set the labels
            xlab("Longitude") + ylab("Latitude") + labs(fill = "Habitat") + 
      
            # set the plot limits to add margin for inset
            # coord_sf(xlim = c(-6.6, -5.95), 
            #          ylim = c(55.55, 55.95), crs = 4326, expand = F) +
  
            # set the theme
            theme_light() +
            theme(legend.position = "right", panel.grid.major = element_blank(),
                  legend.box.background=element_rect(colour = "#BDC3C7"),legend.box.margin=margin(1,1,1,1),
                  legend.spacing.y = unit(0.5, 'cm'), legend.title = element_text(size =14),
                  axis.title = element_text(size =14), axis.text = element_text(size =12),
                  panel.grid.minor = element_blank()) +
            guides(fill = guide_legend(byrow = TRUE))


## save Islay habitaty plot
ggsave(filename = "Plots/Maps/SupFigure_IslayHabitats.png", 
       plot = IsayHab,
       width = 21, 
       height = 19,
       units = "cm",
       dpi = 300)
          
          
          
##---------------------------------##
#### 6. Islay Habitat & Road Map ####
##---------------------------------##

## remove all the other habitat besides the two were selection changes in response to shooting and distance to road
Habitat2 <- Habitat
Habitat2 <- mutate(Habitat2, Hab = ifelse(Hab== "Improved grassland" | Hab == "Saltmarsh & Coastal", paste0(Hab), NA))

## Get the hex codes for the colours
MatVals <- met.brewer("Archambault", length(unique(Habitat$Hab)))
c(MatVals)


## make the plot
RoadHabs <- ggplot() +  
  
  # Add the Habitat data
  geom_raster(data=Habitat2, aes(x=x, y=y, fill=Hab), alpha=1) + 
  scale_fill_manual(name = "Habitat",
                    values =c("Improved grassland" = "#f9d14a",
                              "Saltmarsh & Coastal" = "#88a0dc"),
                    labels = c("Improved grassland",
                               "Saltmarsh & Coastal"), na.value = "transparent") +
  
  ## Add the coastline
  geom_sf(data = Islay_tr, aes(geometry = geometry), fill = NA) +
  geom_sf(data = Road_tr, aes(geometry = geometry, colour = "black"), fill = NA, alpha = 0.6, size = 0.7) +
  scale_colour_manual(name = "Roads",
                      values =c("black" = "black"),
                      labels = c("")) +
  
  # add map extras
  annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.2, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_orienteering,
                         height = unit(1, "cm"), width = unit(0.7, "cm"),) +
  
  # Set the labels
  xlab("Longitude") + ylab("Latitude") + labs(fill = "Habitat") + 
  
  # set the theme
  theme_light() +
  theme(legend.position = "right", panel.grid.major = element_blank(),
        legend.box.background=element_rect(colour = "#BDC3C7"),legend.box.margin=margin(1,1,1,1),
        legend.spacing.y = unit(0.5, 'cm'), legend.title = element_text(size =14),
        axis.title = element_text(size =14), axis.text = element_text(size =12),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))


## save Islay habitaty plot
ggsave(filename = "Plots/Maps/SupFigure_Roads&Habs.png", 
       plot = RoadHabs,
       width = 21, 
       height = 19,
       units = "cm",
       dpi = 300)