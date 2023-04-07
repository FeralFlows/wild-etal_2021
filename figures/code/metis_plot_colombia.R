##################################################################
# Plot Demeter Output (land allocation) to spatial mapping       #
# Author: Mengqi Zhao                                            #
# Email: mengqiz@umd.edu                                         #
# Last Update: 2020-10-21                                        #  
##################################################################



rm(list = ls())

if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("metis" %in% rownames(installed.packages()) == F){devtools::install_github(repo="JGCRI/metis")}
library(metis)
if("rmap" %in% rownames(installed.packages()) == F){devtools::install_github(repo="JGCRI/rmap")}
library(rmap)
if('gcamdata' %in% rownames(installed.packages()) == F){devtools::install_github(repo='JGCRI/gcamdata')}
library(gcamdata)
if('rgdal' %in% rownames(installed.packages()) == F){install.packages('rgdal')}
library(rgdal)
if('raster' %in% rownames(installed.packages()) == F){install.packages('raster')}
library(raster)
if('dplyr' %in% rownames(installed.packages()) == F){install.packages('dplyr')}
library(dplyr)
if('tidyr' %in% rownames(installed.packages()) == F){install.packages('tidyr')}
library(tidyr)
if('data.table' %in% rownames(installed.packages()) == F){install.packages('data.table')}
library(data.table)
if('snow' %in% rownames(installed.packages()) == F){install.packages('snow')}
library(snow)
if('sp' %in% rownames(installed.packages()) == F){install.packages('sp')}
library(sp)
if('stringr' %in% rownames(installed.packages()) == F){install.packages('stringr')}
library(stringr)
if('purrr' %in% rownames(installed.packages()) == F){install.packages('purrr')}
library(purrr)

# tidyverse_conflicts()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work.dir <- getwd()
demeter_output <- 'E:/NEXO-UA/Demeter/example/outputs'
# GCAM v5.3
select_folders <- c('gcam5p3-stash_GFDL-ESM2M_rcp2p6_Reference_2021-05-10_16h55m37s',
                    'gcam5p3-stash_GFDL-ESM2M_rcp2p6_Impacts_2021-05-10_17h42m47s',
                    'gcam5p3-stash_GFDL-ESM2M_rcp2p6_Policy_2021-05-10_18h41m52s')

# select_folders <- c('gcam5p3_IPSL-CM5A-LR_rcp8p5_Reference_2021-04-05_09h08m49s',
#                     'gcam5p3_IPSL-CM5A-LR_rcp8p5_Impacts_2021-04-05_09h56m10s',
#                     'gcam5p3_IPSL-CM5A-LR_rcp8p5_Policy_2021-04-05_10h43m36s')
# GCAM v5.1
# select_folders <- c('reference_2020-10-28_20h29m05s',
#                     'HadGEM2-ES_rcp8p5_Impacts_2020-12-04_08h58m54s',
#                     'HadGEM2-ES_rcp8p5_Policy_2020-12-04_13h07m37s')
scenarios <- c('reference', 'impacts', 'policy')
select_output <- paste(demeter_output, select_folders, sep = '/')
output_files <- list.files(select_output, pattern = '0p5deg', recursive = TRUE, full.names = TRUE)




# Create Basins and hydrologic watersheds within Colombia (better solution)
m_basin <- metis::mapGCAMBasins
m_hydroshed <- metis::mapHydroShed3
m_country <- metis::mapCountries
m_states <- metis::mapStates

m_colombia <- m_country[m_country@data$subRegion %in% c("Colombia"),]; sp::plot(m_colombia)
m_colombia_ext <- m_country[m_country@data$subRegion %in% c("Colombia", 'Ecuador', 'Peru', 'Brazil', 'Venezuela', 'Panama'),];
sp::plot(m_colombia_ext)
m_colombia_states <- m_states[m_states@data$region %in% c('Colombia'), ]; sp::plot(m_colombia_states)

# # Departments within Colombia
# m_colombia_states <- sp::spTransform(m_colombia, raster::crs(m_states))
# m_colombia_states <- raster::crop(m_states, m_colombia)
# m_colombia_states@data <- droplevels(m_colombia_states@data)
# m_colombia_states <- m_colombia_states[m_colombia_states@data$region %in% c('Colombia'), ]
# sp::plot(m_colombia_states)

# Basins within Colombia
m_colombia_basin <- sp::spTransform(m_colombia, raster::crs(m_basin))
m_colombia_basin <- raster::crop(m_basin, m_colombia_basin)
m_colombia_basin@data <-  droplevels(m_colombia_basin@data)
sp::plot(m_colombia_basin)

# Basins within extended Colombia Region
m_colombia_ext_basin <- sp::spTransform(m_colombia_ext, raster::crs(m_basin))
m_colombia_ext_basin <- raster::crop(m_basin, m_colombia_ext_basin)
m_colombia_ext_basin@data <- droplevels(m_colombia_ext_basin@data)
sp::plot(m_colombia_ext_basin)

# Countrys within extended Colombia Region
m_colombia_ext_country <- sp::spTransform(m_colombia_ext, raster::crs(m_country))
m_colombia_ext_country <- raster::crop(m_country, m_colombia_ext_country)
m_colombia_ext_country@data <- droplevels(m_colombia_ext_country@data)
sp::plot(m_colombia_ext_country)

# Hydrosheds within Colombia
m_colombia_hydroshed <- sp::spTransform(m_colombia, raster::crs(m_hydroshed))
m_colombia_hydroshed <- raster::crop(m_hydroshed, m_colombia_hydroshed)
m_colombia_hydroshed@data <- droplevels(m_colombia_hydroshed@data)
sp::plot(m_colombia_hydroshed)

# Hydrosheds within extended Colombia
m_colombia_ext_hydroshed <- sp::spTransform(m_colombia, raster::crs(m_hydroshed))
m_colombia_ext_hydroshed <- raster::crop(m_hydroshed, m_colombia_ext_hydroshed)
m_colombia_ext_hydroshed <- raster::bind(m_colombia_ext_country, m_colombia_ext_hydroshed)
m_colombia_ext_hydroshed@data <- droplevels(m_colombia_ext_hydroshed@data)
sp::plot(m_colombia_ext_hydroshed)


# Plot Colombia Maps at basin and hydroshed level
metis.map(m_colombia, fileName = 'colombia', labels = F, labelsSize = 0.6, pdfpng = 'pdf')
metis.map(m_colombia_states, fileName = 'colombia_states', labels = T, labelsSize = 0.6, pdfpng = 'pdf')
metis.map(m_colombia_basin, fileName = 'colombia_basin', labels = F, labelsSize = 0.6, pdfpng = 'pdf')
metis.map(m_colombia_ext_basin, fileName = 'colombia_ext_basin', labels = T, labelsSize = 0.6)
metis.map(m_colombia_ext_country, fileName = 'colombia_ext_country', labels = T, labelsSize = 0.6)
metis.map(m_colombia_hydroshed, fileName = 'colombia_hydroshed', labels = F, labelsSize = 0.2, pdfpng = 'pdf')
metis.map(m_colombia_ext_hydroshed, fileName = 'colombia_ext_hydroshed', labels = T, labelsSize = 0.2)

# Shapefiles for selected region
# shape <- m_colombia_hydroshed
shape <- m_colombia_states
# shape_ext <- m_colombia_ext_hydroshed # extended region for plotting
# boundary <- raster::extent(bbox(m_colombia_hydroshed))
boundary <- raster::extent(bbox(m_colombia_states))

# Create bounday shapefile in SpatialPolygonsDataFrame format
x1 <- boundary@xmin - 1.0
x2 <- boundary@xmax + 1.0
y1 <- boundary@ymin - 1.0
y2 <- boundary@ymax + 1.0
shape_boundary <- sp::Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))
shape_boundary <- sp::Polygons(list(shape_boundary), ID = "A")
shape_boundary <- sp::SpatialPolygons(list(shape_boundary))
sp::proj4string(shape_boundary) <- m_colombia@proj4string
df_boundary <- matrix(data = NA)
rownames(df_boundary) <- 'A'
shape_boundary <- sp::SpatialPolygonsDataFrame(shape_boundary, data = as.data.frame(df_boundary))

region <- 'Colombia'
subRegType_i <- 'states'
# gcm <- 'IPSL-CM5A-LR' # 'GFDL-ESM2M'  # 'HadGEM2-ES' # 'IPSL-CM5A-LR' # 'MIROC-ESM-CHEM'
# rcp <- 'rcp8p5' # 'rcp2p6'
gcm <- 'GFDL-ESM2M'
rcp <- 'rcp2p6'

# GCAM 5.3
gcam_version <- 'gcam5p3-stash'
main_folder <- paste(region, gcam_version, gcm, rcp, sep = '_')
# GCAM 5.1
# main_folder <- paste(region, gcm, rcp, sep = '_')


#------------------------------- For Demeter -------------------------------#
#Read Demeter data in paralell####
cl <- makeSOCKcluster(4)

#function to read before run####
read_files <- function(file, header = TRUE, sep = 'auto', ...){
  data <- data.table::fread(file, header = header, sep = sep, ...)
  filename <- basename(file[1])
  scenario <- tolower(strsplit(basename(filename), '\\_|\\.')[[1]][2])
  year <- strsplit(basename(filename), '\\_|\\.')[[1]][4]
  data$year <- year
  data$scenario <- scenario
  # data$file_type <- filename
  return(data)
}

# basename is a function to choose file by name
all_data <- parLapply(cl, output_files, read_files)
stopCluster(cl)

col_gather <- c('water', 'forest', 'shrub', 'grass',
                'urban', 'snow', 'sparse', 'corn_irr',
                'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
                'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
                'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
                'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
                'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
                'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
                'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')

df_data <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::rename(gridcode = pkey_0p5_deg, lon = longitude, lat = latitude) %>%
  tidyr::gather(key = 'crop', value = 'value', col_gather)

# Filter data for Colombia region
df_data_colombia <- df_data %>% 
  dplyr::filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)

#------------------------------- DEMETER MAP PLOTTING -------------------------------#
# Plot all the crops for 3 scenarios: reference, impacts, and policy
# The output figures include crop land allocation fraction under individual scenario,
# and the absolute difference and percent difference between (1) impacts and reference;
# and (2) policy and reference
# Note: this plotting may take hours. You may select less crop types in order to reduce run time
if(T){
  crop_i <-c('water', 'forest', 'shrub', 'grass',
             'urban', 'snow', 'sparse', 'corn_irr',
             'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
             'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
             'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
             'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
             'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
             'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
             'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  # crop_i <- c('water', 'forest', 'shrub', 'grass', 'urban', 'snow', 'sparse')
  # crop_i <- c('root_tuber_irr', 'root_tuber_rfd', 'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  
  # crop_i <-c('corn_irr',
  #            'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #            'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
  #            'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
  #            'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #            'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
  #            'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd',
  #            'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  # crop_i <- c('biomass_grass_irr', 'biomass_tree_irr')
  # # 
  # crop_i <- c('fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #              'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #              'wheat_rfd',
  #              'biomass_grass_irr', 'biomass_grass_rfd')

  # crop_i <- c('shrub', 'sparse')
  # crop_i <- c('misccrop_irr', 'fodderherb_irr')
  
  nameAppend_i <- c('_LandUse')
  df <- df_data_colombia
  scenario_i <- c('reference', 'impacts', 'policy')
  year_i <- c(2020, 2030, 2040, 2050)

  
  data_2 <- df %>% 
    filter(scenario %in% scenario_i & year %in% year_i & crop %in% crop_i) %>% 
    tidyr::unite(class, c(crop, scenario, year), sep = '-', remove = TRUE) %>% 
    dplyr::select(-gridcode); str(data_2)
  
  poly_table <- metis.grid2poly(gridFiles = data_2,
                                subRegShape = shape,
                                aggType = 'depth',
                                subRegCol = 'subRegion',
                                subRegType = subRegType_i,
                                nameAppend = '_Demeter',
                                folderName = main_folder,
                                saveFiles = T)
  
  # Find and add missing polygons
  poly_table_temp <- poly_table[0,]
  for (crop in crop_i){
    for (scenario in scenario_i){
      for (year in year_i){
        class_temp <- paste(crop, scenario, year, sep = '-')
        temp <- poly_table %>% 
          filter(class == class_temp) %>% 
          as.data.frame()
        missing <- dplyr::setdiff(shape$subRegion, temp$subRegion) %>% 
          as.data.frame() %>% 
          dplyr::rename(subRegion = '.')
        n_miss <- length(missing$subRegion)
        print(paste0('Missing number of subregions for ', crop, ' | ', scenario, ' | ', year, ' are: ', n_miss))
        if(n_miss > 0){
          df_missing <- data.frame(subRegion = missing$subRegion,
                                   value = rep(1e-10, times = n_miss),
                                   class = rep(class_temp, times = n_miss),
                                   scenario = rep('scenario', times = n_miss),
                                   x = rep('x', times = n_miss),
                                   param = rep('param', times = n_miss),
                                   aggType = rep('depth', times = n_miss),
                                   subRegType = rep(subRegType_i, times = n_miss))
          temp <- rbind(temp, df_missing)
        }
        poly_table_temp <- rbind(poly_table_temp, temp)
      }
    }
  }
  
  # regular plot
  poly_table_2 <- poly_table_temp %>%
    tidyr::separate(class, sep = '-', into = c('class', 'scenario', 'x')) %>% 
    dplyr::mutate(param = class,
                  classPalette = 'pal_green',
                  value = if_else(value < 1e-6, 1e-6, value))

  
  metis.mapsProcess(polygonTable = poly_table_2,
                    subRegShape = shape,
                    subRegCol = 'subRegion',
                    subRegType = subRegType_i,
                    scenRef = 'reference',
                    nameAppend = nameAppend_i,
                    folderName = main_folder,
                    xRange = year_i,
                    scaleRange = c(0,1),
                    mapTitleOn = F,
                    legendFixedBreaks = 8,
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extdendedLabelSize = 0.7,
                    cropToBoundary = F,
                    extension = T,
                    expandPercent = 25,
                    # classPalette = 'pal_green',
                    classPaletteDiff = 'pal_div_BrGn',
                    facetCols = 4,
                    animateOn = F,
                    pdfpng = 'pdf')
  
  # rmap::map(data = poly_table_2,
  #           shape = shape,
  #           subRegCol = 'subRegion',
  #           subRegType = subRegType_i,
  #           scenRef = 'reference',
  #           # nameAppend = nameAppend_i,
  #           folder = paste(getwd(), 'outputs/rmap', main_folder, 'LandUse', sep = '/'),
  #           xRange = year_i,
  #           scaleRange = c(0, 1),
  #           mapTitleOn = F,
  #           # legendFixedBreaks = 8,
  #           boundaryRegShape = shape_boundary,
  #           background = T,
  #           extendedFillColor = "Greys",
  #           extendedLabels = T,
  #           # extdendedLabelSize = 0.7,
  #           cropToBoundary = T,
  #           expandPercent = 25,
  #           # classPalette = 'pal_green',
  #           classPaletteDiff = 'pal_div_BrGn',
  #           facetCols = 4,
  #           animateOn = F)
  
  if(T){
    # aggregate biomass_grass_irr and biomass_tree_irr together
    poly_table_biomass_irr <- poly_table_2 %>% 
      filter(class %in% c('biomass_grass_irr', 'biomass_tree_irr')) %>% 
      mutate(class = 'biomass_irr',
             param = 'biomass_irr',) %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    poly_table_biomass_rfd <- poly_table_2 %>% 
      filter(class %in% c('biomass_grass_rfd', 'biomass_tree_rfd')) %>% 
      mutate(class = 'biomass_rfd',
             param = 'biomass_rfd') %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    
    # plot by crop_irr, crop_rfd
    poly_table_crop <- poly_table_2
    poly_table_crop$param[grepl('irr', poly_table_crop$param)] <- 'crop_irr'
    poly_table_crop$param[grepl('rfd', poly_table_crop$param)] <- 'crop_rfd'
    poly_table_crop$class <- poly_table_crop$param
    poly_table_crop <- poly_table_crop %>% 
      dplyr::filter(class %in% c('crop_irr', 'crop_rfd')) %>% 
      dplyr::group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      dplyr::ungroup()
    
    # plot by crop total
    poly_table_crop_all <- poly_table_crop
    poly_table_crop_all[grepl('crop_irr|crop_rfd', poly_table_crop_all)] <- 'crop_all'
    poly_table_crop_all <- poly_table_crop_all %>% 
      dplyr::group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      dplyr::ungroup()
    
    # plot by land allocation type
    # c('water', 'forest', 'shrub', 'grass', 'urban', 'snow', 'sparse')
    poly_table_landAlloc <- poly_table_2 %>% 
      dplyr::filter(class %in% c('forest', 'shrub', 'grass')) %>% 
      dplyr::bind_rows(poly_table_crop_all %>% dplyr::mutate(class = 'crop')) %>% 
      dplyr::mutate(param = 'landAlloc')
    
  }

  
  poly_table_list <- list(poly_table_biomass_irr, poly_table_biomass_rfd, poly_table_crop, poly_table_crop_all, poly_table_landAlloc)
  # poly_table_list <- list(poly_table_crop, poly_table_crop_all)
  

  
  # Plot all the crop land allocations by crop types
  for(poly_table_i in poly_table_list){
    metis.mapsProcess(polygonTable = poly_table_i,
                      subRegShape = shape,
                      subRegCol = 'subRegion',
                      subRegType = subRegType_i, # 'hydroshed',
                      scenRef = 'reference',
                      nameAppend = nameAppend_i,
                      folderName = main_folder,
                      xRange = year_i,
                      mapTitleOn = F,
                      legendFixedBreaks = 6,
                      boundaryRegShape = shape_boundary,
                      extendedLabels = T,
                      extdendedLabelSize = 0.7,
                      cropToBoundary = F,
                      extension = T,
                      expandPercent = 25,
                      # classPalette = 'pal_green',
                      classPaletteDiff = 'pal_div_BrGn',
                      facetCols = 4,
                      animateOn = F,
                      pdfpng = 'pdf')
  }
  
}



#============================= Spatial Mapping =================================
# Debugging parameters
# gridFiles <- grid_data
# subRegShape <- shape
# subRegShapeExt <- shape
# boundaryRegShape <- shape_boundary
# aggType <- 'depth'
# subRegCol <- 'subRegion'
# subRegType <- subRegType_i
# nameAppend <- '_Tethys_Reference_Dept'
# folderName <- main_folder
# years <- years

find_missing <- function(poly_table, subRegShape, aggType, subRegType, years){
  # Find missing polygons
  shape_area <- data.frame(subRegion = subRegShape@data$subRegion, area = raster::area(subRegShape)) # area is in m2
  
  poly_table_temp <- poly_table[0,]
  for(scenario_t in unique(poly_table$scenario)){
    for(year_t in unique(years)){
      for(param_t in unique(poly_table$param)){
        for(class_t in unique(poly_table$class)){
          temp <- poly_table %>% 
            filter(class == class_t & param == param_t & scenario == scenario_t & x == year_t) %>% 
            as.data.frame()
          missing <- dplyr::setdiff(subRegShape$subRegion, temp$subRegion) %>% 
            as.data.frame() %>% 
            dplyr::rename(subRegion = '.')
          n_miss <- length(missing$subRegion)
          if(n_miss > 0){
            sprintf('-------------------- Adding %s missing polygons --------------------', n_miss)
            missing <- data.frame(subRegion = missing$subRegion,
                                  value = rep(1e-10, times = n_miss),
                                  year = paste('X', rep(year_t, times = n_miss), sep = ''),
                                  scenario = rep(scenario_t, times = n_miss),
                                  x = rep(year_t, times = n_miss),
                                  class = rep(class_t, times = n_miss),
                                  param = rep(param_t, times = n_miss),
                                  aggType = rep(aggType, times = n_miss),
                                  subRegType = rep(subRegType, times = n_miss)) %>% 
              dplyr::select(names(temp))
            
            # Add missing polygons and convert unit from volume to depth
            temp <- rbind(temp, missing)
          }
          poly_table_temp <- rbind(poly_table_temp, temp)
        }
      }
    }
  }
  
  return(poly_table_temp)
}


# Spatial mapping function for Tethys and Xanthos output in km3
sp_mapping <- function(gridFiles, subRegShape, subRegShapeExt, boundaryRegShape,
                       aggType, subRegCol, subRegType, nameAppend, folderName, years, mapPlot=T, ...){
  # Aggregate Tethys 0.5 degree cells to selected polygons (e.g. hydroshed)
  poly_table <- metis.grid2poly(gridFiles = gridFiles,
                                subRegShape = subRegShape,
                                aggType = aggType, # Tethys output in km3, so use 'vol' to aggregate
                                subRegCol = subRegCol,
                                subRegType = subRegType,
                                nameAppend = nameAppend,
                                folderName = folderName,
                                saveFiles = T)
  poly_table_temp <- find_missing(poly_table = poly_table,
                                  subRegShape = subRegShape,
                                  aggType = aggType,
                                  subRegType = subRegType,
                                  years = years)
  
  if(aggType == 'vol'){
    poly_table_2 <- poly_table_temp %>% 
      left_join(shape_area, by = subRegCol) %>% 
      mutate(value = (value*1000^3/area)*1000) # convert from vol (km3) to depth (mm)
  }else{
    poly_table_2 <- poly_table_temp
  }
  
  # facetCols <- ifelse(length(years) <= 4, length(years), ceiling(length(years)/2))
  facetCols <- if (length(years) <= 4) {
    length(years)
  }
  else if (length(years) <= 10) {
    ceiling(length(years)/2)
  }
  else {
    ceiling(sqrt(length(years)))
  }
  # Uncomment if plotting by class
  if(unique(poly_table_2$class) != 'class'){
    facetCols <- ceiling(length(unique(poly_table_2$class))/2)
  }

  if(mapPlot == T){
    # Plot maps
    metis.mapsProcess(polygonTable = poly_table_2,
                      subRegShape = subRegShapeExt,
                      subRegCol = subRegCol,
                      subRegType = subRegType,
                      nameAppend = nameAppend,
                      folderName = folderName,
                      # scenRef = 'Reference',
                      mapTitleOn = F,
                      legendFixedBreaks = 7,
                      # scaleRange = c(0,1),
                      # xRange = years,
                      boundaryRegShape = boundaryRegShape,
                      extendedLabels = T,
                      extdendedLabelSize = 0.7,
                      extension = T,
                      # expandPercent = 25,
                      cropToBoundary = F,
                      classPalette = 'pal_wet', # pal_wet for water, pal_green for cropland
                      animateOn = F,
                      facetCols = facetCols,
                      pdfpng = 'pdf')
  }

  
  return(poly_table_2)
}

# Calculate grid cell area and convert from water volumn to depth
grid_vol_to_depth <- function(grid_data, shape){
  grid <- sp::SpatialPointsDataFrame(sp::SpatialPoints(
    coords = (cbind(grid_data$lon, grid_data$lat))),
    data = grid_data %>% dplyr::select(lat, lon))
  
  sp::gridded(grid) <- TRUE
  grid_raster <- raster::stack(grid); grid_raster
  grid_poly <- raster::rasterToPolygons(grid_raster)
  sp::proj4string(grid_poly) <- sp::proj4string(shape)
  
  # Calculate grid cell area and left join based on lat and lon
  poly_area <- grid_poly@data %>%
    mutate(area=raster::area(grid_poly)) # area in m2
  grid_data <- grid_data %>%
    left_join(poly_area, by=c('lat', 'lon')) %>% 
    mutate(value = value*1000^3/area*1000) %>% 
    dplyr::select(-area)# depth in mm
  
  return(grid_data)
}

#================================ Tethys Output ================================
# wdirr_data_reference <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Reference/wdirr_km3peryr.csv', header = TRUE)
# wdtotal_data_reference <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Reference/wdtotal_km3peryr.csv', header = TRUE)
# wdtotal_data_impacts <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Impacts/wdtotal_km3peryr.csv', header = TRUE)
# wdtotal_data_policy <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Policy/wdtotal_km3peryr.csv', header = TRUE)


years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))
# Total Water Demand (irrigation, domestic, electricity, livestock, manufacturing, mining, non-ag)
tethys_mapping <- function(wd_sector, years, boundary, scenario, gcm, rcp, gcam_version){
  tethys_output <- 'E:/NEXO-UA/Tethys/example/Output'
  gcam_v <- gcam_version
  time_scale <- '2005-2050'
  folder_name <- paste(gcam_v, gcm, rcp, time_scale, scenario, sep = '_')
  file_name <- paste0(wd_sector, '_km3peryr', '.csv')
  file_path <- paste(tethys_output, folder_name, file_name, sep = '/')
  wd_data <- data.table::fread(file_path, header = TRUE)
  
  wd <- wd_data %>% 
    as.data.frame() %>% 
    tidyr::gather(key = 'year', value = 'value', years) %>% 
    dplyr::select(-ID, -ilon, -ilat)
  
  # Filter data to Colombia area
  wd_colombia <- wd %>% 
    dplyr::filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)
  
  # Convert volume to depth based on each grid cell
  wd_colombia_depth <- grid_vol_to_depth(wd_colombia, shape)
  # wd_colombia_depth <- wd_colombia_depth %>%
  #   dplyr::mutate(units = 'Water Demand (mm)')
  
  nameAppend <- paste('_Tethys', scenario, sep = '_')
  
  demand <- sp_mapping(wd_colombia_depth, shape, shape, shape_boundary, 'depth',
                      'subRegion', subRegType_i, nameAppend, main_folder, years, mapPlot = TRUE)
  
  
  return(demand)
}
demand_reference <- tethys_mapping('wdtotal', years, boundary, 'Reference', 'GFDL-ESM2M', 'rcp2p6', 'gcam5p3-stash')
demand_reference <- tethys_mapping('wdtotal', years, boundary, 'Reference', gcm, rcp, 'gcam5p3-stash')
demand_impacts <- tethys_mapping('wdtotal', years, boundary, 'Impacts', gcm, rcp, 'gcam5p3-stash')
demand_policy <- tethys_mapping('wdtotal', years, boundary, 'Policy', gcm, rcp, 'gcam5p3-stash')


demand_reference <- demand_reference %>%
  dplyr::mutate(units = 'Water Demand (mm)')
# Plot GIF Maps
metis.mapsProcess(polygonTable = demand_reference,
                  subRegShape = shape,
                  subRegCol = 'subRegion',
                  subRegType = subRegType_i,
                  nameAppend = '_Tethys_Reference',
                  folderName = main_folder,
                  mapTitleOn = F,
                  legendFixedBreaks = 8,
                  # scaleRange = c(0,1),
                  # xRange = years,
                  boundaryRegShape = shape_boundary,
                  extendedLabels = T,
                  extdendedLabelSize = 0.7,
                  extension = T,
                  # expandPercent = 25,
                  cropToBoundary = F,
                  classPalette = 'pal_wet', # pal_wet for water, pal_green for cropland
                  animateOn = F, # FALSE
                  # fps = 2, # comment if not using animation
                  # legendOutsideSingle = F, # comment if not using animation
                  # legendPosition = c('left', 'bottom'), # comment if not using animation
                  facetCols = 5,
                  pdfpng = 'png')

# Plot each water withdrawal sector in one frame
tethys_indv_mapping <- function(gcm, rcp, scenario, gcam_version, boundary, shape, scen_all=F){
  tethys_output <- 'E:/NEXO-UA/Tethys/example/Output'
  gcam_v <- gcam_version
  time_scale <- '2005-2050'
  folder_name <- paste(gcam_v, gcm, rcp, time_scale, scenario, sep = '_')
  folder_path <- paste(tethys_output, folder_name, sep = '/')
  file_list <- list.files(folder_path, pattern = 'km3peryr', recursive = TRUE, full.names = TRUE)
  file_list <- file_list[!grepl('wdnonag|wdtotal', file_list)]
  
  # Read Tethys files
  # Read Demeter data in paralell####
  cl <- makeSOCKcluster(4)
  clusterExport(cl, list('%>%', 'gather'), envir = .GlobalEnv)
  
  # function to read before run####
  read_files_tethys <- function(file, header = TRUE, sep = 'auto', ...){
    data <- data.table::fread(file, header = header, sep = sep, ...)
    filename <- basename(file[1])
    class <- tolower(strsplit(basename(filename), '\\_|\\.')[[1]][1])
    patterns <- c('wddom', 'wdelec', 'wdirr', 'wdliv', 'wdmfg', 'wdmin')
    replacement <- c('Domestic', 'Electric', 'Irrigation', 'Livestock', 'Manufacturing', 'Mining')
    for (i in seq_along(patterns)){
      class <- gsub(patterns[i], replacement[i], class, perl = TRUE)
    }
    if(scen_all == F){
      data$class <- class
    }else{
      data$class <- paste(class, strsplit(dirname(file[1]), '_')[[1]][5], sep = '-')
    }
    
    years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))
    data <- data %>% 
      as.data.frame() %>% 
      tidyr::gather(key = 'year', value = 'value', years) %>% 
      dplyr::select(-ID, -ilon, -ilat)
    
    return(data)
  }
  
  # basename is a function to choose file by name
  wd_indv <- parLapply(cl, file_list, read_files_tethys)
  stopCluster(cl)
  
  wd_indv_df <- rbindlist(wd_indv, idcol = FALSE)
  
  # Filter data to boundary area
  wd_indv_boundary <- wd_indv_df %>% 
    dplyr::filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)
  
  # Convert volume to depth based on each grid cell
  wd_indv_boundary_depth <- grid_vol_to_depth(wd_indv_boundary, shape)
  
  if(scen_all == F){
    nameAppend <- paste('_Tethys', scenario, 'Indv', sep = '_')
    demand <- sp_mapping(wd_indv_boundary_depth, shape, shape, shape_boundary, 'depth',
                         'subRegion', subRegType_i, nameAppend, main_folder, years)
  }else{
    nameAppend <- paste('_Tethys', 'Indv_ALL', sep = '_')
    poly_table <- metis.grid2poly(gridFiles = wd_indv_boundary_depth,
                                  subRegShape = shape,
                                  aggType = 'depth', # Tethys output in km3, so use 'vol' to aggregate
                                  subRegCol = 'subRegion',
                                  subRegType = subRegType_i,
                                  nameAppend = nameAppend,
                                  folderName = main_folder,
                                  saveFiles = T)
    poly_table_temp <- find_missing(poly_table = poly_table,
                                    subRegShape = shape, 
                                    aggType = 'depth',
                                    subRegType = subRegType_i,
                                    years = years)
    
    poly_table_2 <- poly_table_temp %>%
      tidyr::separate(class, sep = '-', into = c('class', 'scenario')) %>% 
      mutate(param = 'WaterDemands',
             classPalette = 'pal_wet')
    
    metis.mapsProcess(polygonTable = poly_table_2,
                      subRegShape = shape,
                      subRegCol = 'subRegion',
                      subRegType = subRegType_i,
                      scenRef = 'Reference',
                      nameAppend = nameAppend,
                      folderName = main_folder,
                      xRange = year_i,
                      mapTitleOn = F,
                      legendFixedBreaks = 5,
                      boundaryRegShape = shape_boundary,
                      extendedLabels = T,
                      extdendedLabelSize = 0.7,
                      cropToBoundary = F,
                      extension = T,
                      expandPercent = 25,
                      classPaletteDiff = 'pal_div_BrGn',
                      facetCols = 6,
                      animateOn = F,
                      pdfpng = 'pdf')
    
  }
  
}

tethys_indv_mapping(gcm, rcp, c('Reference', 'Impacts', 'Policy'), 'gcam5p3-stash', boundary, shape, scen_all=T)
tethys_indv_mapping(gcm, rcp, 'Reference', 'gcam_5p3', boundary, shape, scen_all=F)
tethys_indv_mapping(gcm, rcp, 'Impacts', 'gcam_5p3', boundary, shape, scen_all=F)
tethys_indv_mapping(gcm, rcp, 'Policy', 'gcam_5p3', boundary, shape, scen_all=F)


#================================ Xanthos Output ===============================
# for original Xanthos output
run <- 'clim_impacts'
runoff_var_name <- 'q_km3peryear'
# gcm <- 'MIROC-ESM-CHEM'
# rcp <- 'rcp6p0'
time_scale <- '1950_2099'

xanthos_folder <- paste(run, gcm, rcp, sep = '_')
xanthos_output <- paste('E:/NEXO-UA/Xanthos/example/output', xanthos_folder, sep = '/')
runoff_file <- paste0(paste(runoff_var_name, gcm, rcp, time_scale, sep = '_'), '.csv')


# Read Xanthos runoff output in km3
runoff <- data.table::fread(file = paste(xanthos_output, runoff_file, sep = '/'), header = TRUE)
coord_0p5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
  dplyr::rename(gridcode = V1,
         longitude = V2,
         latitude = V3,
         x = V4,
         y = V5)

runoff_years <- sprintf('%s', seq(from = 1950, to = 2099, by = 1))
years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))
# years <- sprintf('%s', seq(from = 1985, to = 2050, by = 5))

## Calculate historical mean from GCMs and RCPs for the Reference scenario
runoff_hist <- data.frame(value = rowMeans(runoff[, which(colnames(runoff)=='1950') : which(colnames(runoff)=='2005')]))
runoff_hist <- runoff_hist %>% 
  dplyr::mutate(lat = coord_0p5deg$latitude,
                lon = coord_0p5deg$longitude,
                year = 'historical') %>% 
  dplyr::filter(lat <= boundary@ymax, lat >= boundary@ymin, 
                lon <= boundary@xmax, lon >= boundary@xmin)
runoff_hist_depth <- grid_vol_to_depth(runoff_hist, shape)

supply_reference <- sp_mapping(runoff_hist_depth, shape, shape, shape_boundary, 'depth',
                     'subRegion', subRegType_i, '_Xanthos_HistoricalMean', main_folder, years='historical', mapPlot = T)

## Calculate runoff for GCMs and RCPs projections
runoff_colombia <- runoff %>% 
  dplyr::mutate(lat = coord_0p5deg$latitude,
         lon = coord_0p5deg$longitude) %>% 
  tidyr::gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
  dplyr::filter(lat <= boundary@ymax, lat >= boundary@ymin, 
                lon <= boundary@xmax, lon >= boundary@xmin)
  # dplyr::filter(year %in% years) %>% 
  # dplyr::select(-id)


# Use this to smooth the runoff
runoff_colombia_smooth <- runoff_colombia %>% 
  # tidyr::gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
  dplyr::filter(year >= 2000) %>%
  dplyr::group_by(id) %>% 
  tidyr::nest() %>%
  dplyr::mutate(smooth = data %>% purrr::map(~loess(value ~ year, data = .))) %>% 
  dplyr::mutate(pred = purrr::map2(smooth, data, predict)) %>%
  tidyr::unnest(c(pred, data)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(lat, lon, year, pred) %>% 
  dplyr::rename(value = pred)

runoff_colombia_depth <- grid_vol_to_depth(runoff_colombia_smooth, shape)

supply <- sp_mapping(runoff_colombia_depth, shape, shape, shape_boundary, 'depth',
                     'subRegion', subRegType_i, '_Xanthos_Smooth', main_folder, years, mapPlot = T)

supply <- supply %>%
  dplyr::mutate(units = 'Runoff (mm)')
# Plot GIF Maps
metis.mapsProcess(polygonTable = supply,
                  subRegShape = shape,
                  subRegCol = 'subRegion',
                  subRegType = subRegType_i,
                  nameAppend = '_Xanthos',
                  folderName = main_folder,
                  mapTitleOn = F,
                  legendFixedBreaks = 7,
                  # scaleRange = c(0,1),
                  # xRange = years,
                  boundaryRegShape = shape_boundary,
                  extendedLabels = T,
                  extdendedLabelSize = 0.7,
                  extension = T,
                  # expandPercent = 25,
                  cropToBoundary = F,
                  classPalette = 'pal_wet', # pal_wet for water, pal_green for cropland
                  animateOn = F, # FALSE
                  # fps = 1, # comment if not using animation
                  # legendOutsideSingle = F, # comment if not using animation
                  # legendPosition = c('left', 'bottom'), # comment if not using animation
                  facetCols = 5,
                  pdfpng = 'png')


if(F){
  # for historical water supply
  gcm <- 'watch+wfdei'
  rcp <- '1970_2010'
  time_scale <- '1970_2010'
  xanthos_folder <- paste(run, gcm, time_scale, sep = '_')
  xanthos_output <- paste('E:/NEXO-UA/Xanthos/example/output', xanthos_folder, sep = '/')
  runoff_file <- paste0(paste(runoff_var_name, gcm, time_scale, sep = '_'), '.csv')

  # Read Xanthos runoff output in km3
  runoff <- data.table::fread(file = paste(xanthos_output, runoff_file, sep = '/'), header = TRUE)
  coord_0p5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
    dplyr::rename(gridcode = V1,
           longitude = V2,
           latitude = V3,
           x = V4,
           y = V5)
  
  years_hist <- c('2005', '2010')
  runoff_years <- sprintf('%s', seq(from = 1970, to = 2010, by = 1))
  runoff_colombia_hist <- runoff %>% 
    mutate(lat = coord_0p5deg$latitude,
           lon = coord_0p5deg$longitude) %>% 
    gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
    filter(year %in% runoff_years, lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin) %>% 
    dplyr::select(-id)
  
  runoff_colombia_hist_depth <- grid_vol_to_depth(runoff_colombia_hist, shape)
  
  supply_hist <- sp_mapping(runoff_colombia_hist_depth, shape, shape, shape_boundary, 'depth',
                       'subRegion', subRegType_i, '_XanthosHist', 'Colombia_Hist', runoff_years)
}






#------------------------------- Water Scarcity -------------------------------#
# demand_clean <- demand$shapeTbl %>% 
#   dplyr::select(-units, -region, -classPalette, -multiFacetCol, -multiFacetRow)
# supply_clean <- supply$shapeTbl %>% 
#   dplyr::select(-units, -region, -classPalette, -multiFacetCol, -multiFacetRow)

scarcity_mapping <- function(demand, supply, scenario, ext_append){
  if(scenario == 'Reference' & ext_append == 'Historical'){
    poly_table_scarcity <- demand %>% 
      left_join(supply, by = c('subRegion',
                               'class',
                               'scenario',
                               'param',
                               'aggType',
                               'subRegType')) %>% 
      dplyr::mutate(value = value.x/value.y,
                    units = 'Water Scarcity') %>% 
      dplyr::rename(year = year.x,
                    x = x.x) %>% 
      dplyr::select(-value.x, -value.y, -year.y, -x.y) 
  } else {
    poly_table_scarcity <- demand %>% 
      left_join(supply, by = c('subRegion',
                               'year',
                               'x',
                               'class',
                               'scenario',
                               'param',
                               'aggType',
                               'subRegType')) %>% 
      dplyr::mutate(value = value.x/value.y,
                    units = 'Water Scarcity') %>% 
      dplyr::select(-value.x, -value.y) 
  }
  
  
  numeric2Cat_param <- list("param")
  numeric2Cat_breaks <- list(c(0,0.1,0.2,0.4,Inf))
  numeric2Cat_labels <- list(c("None","Low","Moderate","Severe"))
  numeric2Cat_palette <- list(c("#3288BD", "#ABDDA4", "#FDAE61", "#9E0142"))
  numeric2Cat_legendTextSize <- list(c(1))
  numeric2Cat_list <-list(numeric2Cat_param=numeric2Cat_param,
                          numeric2Cat_breaks=numeric2Cat_breaks,
                          numeric2Cat_labels=numeric2Cat_labels,
                          numeric2Cat_palette=numeric2Cat_palette,
                          numeric2Cat_legendTextSize=numeric2Cat_legendTextSize)
  nameAppend <- paste('_Scarcity', scenario, ext_append, sep = '_')
  
  # Plot maps
  metis.mapsProcess(polygonTable = poly_table_scarcity,
                    subRegShape = shape,
                    subRegCol = 'subRegion',
                    subRegType = subRegType_i,
                    nameAppend = nameAppend,
                    folderName = main_folder,
                    mapTitleOn = F,
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extdendedLabelSize = 0.7,
                    extension = T,
                    cropToBoundary = F,
                    facetCols = 5,
                    animateOn = F,
                    # fps = 2, # comment if not using animation
                    # legendOutsideSingle = F, # comment if not using animation
                    # legendPosition = c('left', 'bottom'), # comment if not using animation
                    numeric2Cat_list = numeric2Cat_list,
                    pdfpng = 'pdf')
}

scarcity_mapping(demand_reference, supply_reference, 'Reference', 'Historical')
scarcity_mapping(demand_reference, supply, 'Reference', 'Smooth')
scarcity_mapping(demand_impacts, supply, 'Impacts', 'Smooth')
scarcity_mapping(demand_policy, supply, 'Policy', 'Smooth')



# sp_scarcity <- sp::SpatialPointsDataFrame(sp::SpatialPoints(coords = (cbind(datax$lon, 
#                                                              datax$lat))), data = datax)
# metis.map(dataPolygon = poly_table_scarcity,
#           fillColumn = 'subRegion',
#           shpFile = shape_ext,
#           fillPalette = 'pal_ScarcityCat',
#           legendBreaks = c(0,0.1,0.2,0.4),
#           fileName = 'Colombia_Scarcity')



#------------------------------- Demeter Historical Land Allocation -------------------------------#
modis <- data.table::fread(file = 'E:/NEXO-UA/Demeter/example/inputs/observed/gcam_reg32_basin235_modis_mirca_2000_0p5deg_MZ.csv',
                           header = TRUE)
df_modis <- modis %>% 
  rename(lat = latitude,
         lon = longitude) %>% 
  dplyr::select(-pkey_0p5_deg) %>% 
  gather(key = 'class', value = 'value', append(col_gather, 'crops')) %>% 
  filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)

df_modis_select <- df_modis %>% 
  filter(class %in% c('grass', 'forest', 'shrub', 'urban', 'crops'))


poly_table_modis <- sp_mapping(df_modis_select, shape, shape, shape_boundary, 'vol',
           'subRegion', 'hydroshed', '_DemeterModis', 'Colombia', '2000')





#------------------------------- NOT IN USE -------------------------------#
if(F){
  
  basins <- unique(m_colombia_basin@data$subRegion)
  value_col <- df_data_by_basin %>% 
    filter(basin_name %in% basins & year %in% 2005 & scenario %in% 'reference') %>% 
    select(value=corn_irr)
  
  
  data <- data.frame(subRegion = value_col$basin_name,
                     year = rep(2005, 12),
                     value = value_col$value)
  
  
  
  # metis.prepGrid (Fail)
  # This function is designed to be used with specific open-source downscaling models (Xanthos [18],
  # Demeter [19], and Tethys [20]) that downscale GCAM data to the grid level. The function takes
  # outputs from these various models and processes them into the format required for providing
  # input to the metis.mapsProcess.R function
  metis.prepGrid(demeterFolders = paste(demeter_output, '/reference_2020-10-05_15h02m50s',sep = ''),
                 demeterScenarios = 'reference',
                 demeterTimesteps = seq(from = 2005, to = 2050, by = 5),
                 demeterUnits = 'fraction',
                 dirOutputs = paste(getwd(), '/demeter', sep = ''))
  
  # metis.grid2poly (Works)
  # Function used to crop and aggregate gridded data by a given polygon shape file. If no grid is
  # provided, the function can still be used to produce regional and subregional maps.
  
  metis.grid2poly(gridFiles = 'E:/NEXO-UA/Results/metis/outputs/Colombia/Grid2Poly/demeter_test.csv',
                  subRegShape = m_colombia_hydroshed,
                  subRegCol = 'subRegion',
                  subRegType = 'hydrosheds',
                  aggType = 'depth',
                  folderName = 'Colombia',
                  saveFiles = T)
  
  data_hydroshed <- data.table::fread('E:/NEXO-UA/Results/metis/outputs/Colombia/Grid2Poly/poly_scenario_hydrosheds_param.csv',
                                      header = TRUE)
  
  str(data_hydroshed)
  col.names <- c('value', 'forest',
                 'shrub', 'grass', 'urban', 'snow', 'sparse', 'corn_irr',
                 'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
                 'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
                 'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
                 'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
                 'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
                 'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
                 'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  
  data_hydroshed[, (col.names) := lapply(.SD, as.numeric), .SDcols = col.names]
  
  
  data <- data_hydroshed %>% 
    group_by(subRegion) %>% 
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup()
  
  data <- data.frame(subRegion = as.character(m_colombia_hydroshed$subRegion),
                     year = rep(2005, length(m_colombia_hydroshed$subRegion)),
                     value = m_colombia_hydroshed$SUB_AREA)
  
  metis.mapsProcess(polygonTable = data_hydroshed,
                    subRegShape = m_colombia_hydroshed,
                    folderName = 'Colombia',
                    subRegCol = 'subRegion',
                    # subRegType = 'hydrosheds',
                    mapTitleOn = F,
                    legendFixedBreaks = 6,
                    cropToBoundary = F,
                    extension = T,
                    expandPercent = 15,
                    classPalette = 'pal_green')
  
  
  # read basin codes and coordinates
  basin <- data.table::fread(file = 'basin.csv', header = TRUE); str(basin)
  basin_names <- data.table::fread(file = 'gcam_basin_lookup.csv', header = TRUE) %>% 
    rename(basin = basin_id); str(basin_names)
  coord_5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
    rename(gridcode = V1,
           longitude = V2,
           latitude = V3,
           x = V4,
           y = V5) %>% 
    append(basin) %>% 
    as.data.frame() %>% 
    left_join(basin_names, by = 'basin') %>% 
    select(-x, -y); str(coord_5deg)
  
  # assign the basin codes according to the coordinates
  df_data_2 <- df_data %>% 
    as_tibble() %>% 
    gcamdata::left_join_error_no_match(coord_5deg, by = c('gridcode', 'longitude', 'latitude'))
  
  df_data_by_basin <- df_data_2 %>% 
    group_by(scenario, year, basin, basin_name, glu_name) %>% 
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% 
    select(-longitude, -latitude, -gridcode)
  str(df_data_by_basin)
  
  # Read in interception of World Basins and 32 region shape file
  map_inter_basin_32reg <- metis::mapIntersectGCAMBasin32Reg; head(map_inter_basin_32reg@data) 
  map_inter_basin_32reg <- map_inter_basin_32reg[map_inter_basin_32reg@data$subRegion_GCAMReg32 %in% c("Colombia"),]
  map_inter_basin_32reg@data <- droplevels(map_inter_basin_32reg@data)
  metis.map(map_inter_basin_32reg,labels=F) # View custom shape
  map_inter_basin_32reg@data <- map_inter_basin_32reg@data %>% 
    dplyr::rename(basins = subRegion_GCAMBasin)
  str(map_inter_basin_32reg@data)
  map_inter_basin_32reg@data <- map_inter_basin_32reg@data %>% 
    dplyr::mutate(subRegion=basins);  map_inter_basin_32reg@data
}

if(F){
  # Plot crop land allocation fraction in group of four crops for years and scenarios
  # of your choice
  sp_alloc <- function(args_ls, df, scenario_i, year_i, shape){
    # browser()
    class_i <- args_ls[1:4]
    nameAppend_i <- args_ls[5]
    
    data_2 <- df %>% 
      filter(scenario %in% scenario_i & year %in% year_i & class %in% class_i) %>% 
      tidyr::unite(class, c(class, scenario, year), sep = '-', remove = TRUE) %>% 
      dplyr::select(-gridcode); str(data_2)
    
    poly_table <- metis.grid2poly(gridFiles = data_2,
                                  subRegShape = shape,
                                  aggType = 'depth',
                                  subRegCol = 'subRegion',
                                  subRegType = 'hydroshed',
                                  nameAppend = '_ColombiaFrac',
                                  folderName = 'Colombia',
                                  saveFiles = T)
    
    # poly_table_2 <- poly_table %>%
    #   tidyr::separate(class, sep = '-', into = c('class', 'scenario', 'x')) %>% 
    #   mutate(param = class)
    
    metis.mapsProcess(polygonTable = poly_table,
                      subRegShape = shape,
                      subRegCol = 'subRegion',
                      subRegType = 'hydroshed',
                      nameAppend = nameAppend_i,
                      folderName = 'Colombia',
                      mapTitleOn = F,
                      legendFixedBreaks = 9,
                      cropToBoundary = F,
                      extension = T,
                      expandPercent = 25,
                      classPalette = 'pal_green',
                      facetCols = 8,
                      animateOn = F)
  }
  
  
  scenario_i <- 'reference'
  year_i <- c(2020, 2030, 2040, 2050)
  # class_i <- c('biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  
  args_ls <- list(C1 = c('water', 'forest', 'shrub', 'grass', '_Group1'),
                  C2 = c('urban', 'snow', 'sparse', 'corn_irr', '_Group2'),
                  C3 = c('fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr', '_Group3'),
                  C4 = c('oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr', '_Group4'),
                  C5 = c('root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd', '_Group5'),
                  C6 = c('fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd', '_Group6'),
                  C7 = c('oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd', '_Group7'),
                  C8 = c('root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland', '_Group8'),
                  C9 = c('biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd', '_Group9'))
  
  lapply(args_ls, sp_alloc, df_data_colombia, scenario_i, year_i, m_colombia_hydroshed)
}

if(F){
  # Aggreagte gridcells to selected polygons
  poly_table_runoff <- metis.grid2poly(gridFiles = runoff_colombia,
                                       subRegShape = shape,
                                       aggType = 'vol', # Xanthos output in km3, so use 'vol' to aggregate
                                       subRegCol = 'subRegion',
                                       subRegType = 'hydroshed',
                                       nameAppend = '_ColombiaXanthos',
                                       folderName = 'Colombia',
                                       saveFiles = T)
  
  # Find missing polygons
  missing <- dplyr::setdiff(shape$subRegion, poly_table$subRegion) %>% 
    as.data.frame() %>% 
    rename(subRegion = '.') %>% 
    left_join(shape_area, by = 'subRegion')
  if(length(missing$subRegion) > 0){
    rep <- length(years) * length(missing$subRegion)
    n_miss <- length(missing$subRegion)
    missing <- data.frame(subRegion = rep(missing$subRegion, times = length(years)),
                          value = rep(0, times = rep),
                          year = rep(years, each = n_miss),
                          class = rep('class', times = rep),
                          scenario = rep('scenario', times = rep),
                          x = paste('X', rep(years, each = n_miss), sep = ''),
                          param = rep('param', times = rep),
                          aggType = rep('vol', times = rep),
                          subRegType = rep('hyroshed', times = rep))
  }
  
  
  # Add missing polygons and convert unit from volume to depth
  poly_table_supply <- poly_table_runoff %>% 
    # dplyr::union(missing) %>% 
    left_join(shape_area, by = 'subRegion') %>% 
    # filter(subRegion %in% m_colombia_hydroshed@data$subRegion & !subRegion %in% m_chile_hydroshed$subRegion) %>% 
    # dplyr::select(subRegion %in% m_chile_hydroshed$subRegion) %>% 
    mutate(value = (value/area)*1000*1000) # convert from km to mm
  
  # Plot maps
  metis.mapsProcess(polygonTable = poly_table_supply,
                    subRegShape = shape_ext,
                    subRegCol = 'subRegion',
                    subRegType = 'hydroshed',
                    nameAppend = '_Xanthos',
                    folderName = 'Colombia',
                    mapTitleOn = F,
                    legendFixedBreaks = 7,
                    # scaleRange = c(0,254),
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extdendedLabelSize = 1,
                    extension = T,
                    expandPercent = 25,
                    cropToBoundary = F,
                    classPalette = 'pal_wet',
                    facetCols = 5
  )
}