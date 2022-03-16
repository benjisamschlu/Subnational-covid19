


#######################################################################################
##      PROJECT: Sub-national estimation of covid impact in BEL at district level
##      -- by Benjamin-Samuel Schlüter --
##      UCLouvain
##      09/04/2021
#######################################################################################
#
# CODE AIM: Create maps for outputs of models
#
#######################################################################################
#
#
# Notes:
# 1) 

rm(list=ls())

############################################################################################################################################################################################



#######################################################################################################################################################################################################################
# ----- Load pkg ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################################################################################################################################################################################################################


# packages 

library(mapview)
library(ggmap)
library(zipcodeR)
library(rgeos)
library(sp)
library(maptools)
library(BelgiumMaps.StatBel)
library(leaflet)
library(readr)
library(mapview)
library(leafsync)
library(viridis)
library(manipulateWidget)
library(htmlwidgets)
library(psych)
library(StatMeasures)
library(tidyr)
library(broman)

data(BE_ADMIN_DISTRICT, package = "BelgiumMaps.StatBel")
data(BE_ADMIN_MUNTY, package = "BelgiumMaps.StatBel")


####################################################################################################################
# ----- MAPS FOR LLE0 --------------------------------------------------------------------------------------------------
####################################################################################################################


# data coming from estimates of lle0
dta = readRDS("./data/estimates/lle0_for_maps.rda")
dta <- dta %>% 
        mutate(TX_ADM_DSTR_DESCR_FR = dist,
               TX_ADM_DSTR_DESCR_FR = str_replace(TX_ADM_DSTR_DESCR_FR, "'", "’"))
# Edit districts according to new admin from 1st Jan 2019
# Create La Louvière district 
# and rename Tournai and Mouscron to Tournai-Mouscron
BE_ADMIN_MUNTY@data <- BE_ADMIN_MUNTY@data %>% 
        mutate(TX_ADM_DSTR_DESCR_FR = ifelse(TX_MUNTY_DESCR_FR %in% c("La Louvière", "Binche", "Estinnes", "Morlanwelz"), "Arrondissement de La Louvière",
                                             TX_ADM_DSTR_DESCR_FR),
               TX_ADM_DSTR_DESCR_FR = case_when(TX_ADM_DSTR_DESCR_FR == "Arrondissement de Mouscron" ~ "Arrondissement de Tournai-Mouscron",
                                                TX_ADM_DSTR_DESCR_FR == "Arrondissement de Tournai" ~ "Arrondissement de Tournai-Mouscron",
                                                TRUE ~ TX_ADM_DSTR_DESCR_FR)
               )
BE_ADMIN_DISTRICT@data <- BE_ADMIN_DISTRICT@data %>%
        mutate(TX_ADM_DSTR_DESCR_FR = case_when(TX_ADM_DSTR_DESCR_FR == "Arrondissement de Mouscron" ~ "Arrondissement de Tournai-Mouscron",
                                                TX_ADM_DSTR_DESCR_FR == "Arrondissement de Tournai" ~ "Arrondissement de Tournai-Mouscron",
                                                TRUE ~ TX_ADM_DSTR_DESCR_FR))
# Merge polygons of municipality into districts to obtain correct lines
CORRECT_DIST = aggregate(BE_ADMIN_MUNTY, by = list(ID = BE_ADMIN_MUNTY@data$TX_ADM_DSTR_DESCR_FR), FUN = "mean")

# Exctract info for the labels of districts
centroids <- getSpPPolygonsLabptSlots(BE_ADMIN_DISTRICT)

labels_df <- tibble(dist = BE_ADMIN_DISTRICT$TX_ADM_DSTR_DESCR_FR,
                    lat = centroids[,2],
                    lng = centroids[,1]) %>% 
        mutate(dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d’" ~ substr(dist, 18, nchar(dist))),
               # Slight modifications to avoid overlap of labels
               lat = ifelse(dist == "Hal-Vilvorde", 50.96, lat),
               lat = ifelse(dist == "Ath", 50.66, lat)) %>% 
        rbind(tibble(dist = "La Louvière",
                     lat = 50.4,
                     lng = 4.2))

# merge with the file mymaps 
mymap <- merge(BE_ADMIN_MUNTY, dta, by = "TX_ADM_DSTR_DESCR_FR", all.x=F, all.y=F) #merge your data file with the shapefile

pal <- colorNumeric(palette = heat.colors(50),  domain = mymap$median, na.color = "#808080",
                    reverse=F) # palette definition
m1 <-                                          # <br> stands for enter; %s stands for number
        leaflet(mymap,
                options = leafletOptions(zoomControl = FALSE)) %>%
        #addTiles() %>% # avoid ploting surrounding maps
        addLegend(title = "Posterior median (in years)", 
                  pal = pal,
                  values = ~ median, 
                  position = "bottomleft", na.label = "Missing") %>%
        # polygons of color are added at the municipality level so they are correct
        addPolygons(color = ~pal(median), 
                    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.85) %>%  
        # add only lables without a symbol (ie an indicator)
        addLabelOnlyMarkers(lng = labels_df$lng,
                            lat = labels_df$lat,
                            label = labels_df$dist,
                            labelOptions = labelOptions(noHide = T, textOnly = TRUE,
                                                        direction = "center",
                                                        textsize = "12px",
                                                        style = list("font-weight" = "bold")))
m1 <- addPolylines(m1, data = subset(CORRECT_DIST), weight = 1, color = "black")  # add your required borders 
m1


####################################################################################################################
# ----- MAPS FOR LLE0 UNCERTAINTY --------------------------------------------------------------------------------------------------
####################################################################################################################

pal <- colorNumeric(palette = viridis(20, option = "D"),  domain = mymap$width_u, na.color = "#808080",
                    reverse=T) # palette definition
m2 <-                                          # <br> stands for enter; %s stands for number
        leaflet(mymap) %>%
        #addTiles() %>% # avoid ploting surrounding maps
        addLegend(title = "Width of 95% CI (in years)", 
                  pal = pal,
                  values = ~ width_u, 
                  position = "bottomleft", na.label = "Missing") %>%
        # polygons of color are added at the municipality level so they are correct
        addPolygons(color = ~pal(width_u), 
                    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.85) %>%  
        # add only lables without a symbol (ie an indicator)
        addLabelOnlyMarkers(lng = labels_df$lng,
                            lat = labels_df$lat,
                            label = labels_df$dist,
                            labelOptions = labelOptions(noHide = T, textOnly = TRUE,
                                                        direction = "center",
                                                        textsize = "12px",
                                                        style = list("color" = "#black",
                                                                     "font-weight" = "bold")))
m2 <- addPolylines(m2, data = subset(CORRECT_DIST), weight = 1, color = "black")  # add your required borders 
m2



####################################################################################################################
# ----- MAPS FOR SMR --------------------------------------------------------------------------------------------------
####################################################################################################################


dta = readRDS("./data/estimates/smr_for_maps.rda")
dta <- dta %>% 
        mutate(TX_ADM_DSTR_DESCR_FR = dist,
               TX_ADM_DSTR_DESCR_FR = str_replace(TX_ADM_DSTR_DESCR_FR, "'", "’")) %>% 
        filter(year == "2020") %>% 
        mutate(width_u = upper - lower)

# merge with the file mymaps 
mymap <- merge(BE_ADMIN_MUNTY, dta, by = "TX_ADM_DSTR_DESCR_FR", all.x=F, all.y=F) #merge your data file with the shapefile

pal <- colorNumeric(palette = heat.colors(50),  domain = mymap$median, na.color = "#808080",
                    reverse=T) # palette definition
m3 <-                                          # <br> stands for enter; %s stands for number
        leaflet(mymap,
                options = leafletOptions(zoomControl = FALSE)) %>%
        #addTiles() %>% # avoid ploting surrounding maps
        addLegend(title = "Posterior median", 
                  pal = pal,
                  values = ~ median, 
                  position = "bottomleft", na.label = "Missing") %>%
        # polygons of color are added at the municipality level so they are correct
        addPolygons(color = ~pal(median), 
                    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.85) %>%  
        # add only lables without a symbol (ie an indicator)
        addLabelOnlyMarkers(lng = labels_df$lng,
                            lat = labels_df$lat,
                            label = labels_df$dist,
                            labelOptions = labelOptions(noHide = T, textOnly = TRUE,
                                                        direction = "center",
                                                        textsize = "12px",
                                                        style = list("font-weight" = "bold")))
m3 <- addPolylines(m3, data = subset(CORRECT_DIST), weight = 1, color = "black")  # add your required borders 
m3


####################################################################################################################
# ----- MAPS FOR SMR UNCERTAINTY --------------------------------------------------------------------------------------------------
####################################################################################################################

pal <- colorNumeric(palette = viridis(20, option = "D"),  domain = mymap$width_u, na.color = "#808080",
                    reverse=T) # palette definition
m4 <-                                          # <br> stands for enter; %s stands for number
        leaflet(mymap,
                options = leafletOptions(zoomControl = FALSE)) %>%
        #addTiles() %>% # avoid ploting surrounding maps
        addLegend(title = "Width of 95% CI", 
                  pal = pal,
                  values = ~ width_u, 
                  position = "bottomleft", na.label = "Missing") %>%
        # polygons of color are added at the municipality level so they are correct
        addPolygons(color = ~pal(width_u), 
                    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.85) %>%  
        # add only lables without a symbol (ie an indicator)
        addLabelOnlyMarkers(lng = labels_df$lng,
                            lat = labels_df$lat,
                            label = labels_df$dist,
                            labelOptions = labelOptions(noHide = T, textOnly = TRUE,
                                                        direction = "center",
                                                        textsize = "12px",
                                                        style = list("font-weight" = "bold")))
m4 <- addPolylines(m4, data = subset(CORRECT_DIST), weight = 1, color = "black")  # add your required borders 
m4
