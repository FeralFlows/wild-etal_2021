
# metis.master.R
# Script to run different parts of the metis package.

rm(list=ls())

#----------------------------
# Install necessary packages
#----------------------------
if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("metis" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/metis")}
library(metis)
if("rgcam" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/rgcam")}
library(rgcam)
if("tibble" %in% rownames(installed.packages()) == F){install.packages("tibble")}
library(tibble)
if("rgdal" %in% rownames(installed.packages()) == F){install.packages("rgdal")}
library(rgdal)
if("tmap" %in% rownames(installed.packages()) == F){install.packages("tmap")}
library(tmap)
if("zoo" %in% rownames(installed.packages()) == F){install.packages("zoo")}
library(zoo)
if("RSQLite" %in% rownames(installed.packages()) == F){install.packages("RSQLite")}
library(RSQLite)
if("ggplot2" %in% rownames(installed.packages()) == F){install.packages("ggplot2")}
library(ggplot2)
if("ggalluvial" %in% rownames(installed.packages()) == F){install.packages("ggalluvial")}
library(ggalluvial)

if("dplyr" %in% rownames(installed.packages()) == F){install.packages("dplyr")}
library(dplyr)
if("dbplyr" %in% rownames(installed.packages()) == F){install.packages("dbplyr")}
library(dbplyr)
if("tidyr" %in% rownames(installed.packages()) == F){install.packages("tidyr")}
library(tidyr)
library(plutus)

#----------------------------
# Read GCAM Data (metis.readgcam.R)
#---------------------------

setwd('E:/NEXO-UA/Results/metis/gcam_database')

# Connect to gcam database or project
# gcamdatabasePath_i <- 'E:/NEXO-UA/GCAM-Workspace/gcam-core_LAC_v02_5Nov2019/output/FinalRuns' # Use if gcamdatabase is needed
# gcamdatabaseName_i <- 'IDBNexus_HadGEM2-ES_rcp8p5' # "Reference_originalSW" Use if gcamdatabse is needed
gcamdatabasePath_i <- 'E:/NEXO-UA/GCAM-Workspace/gcam-core-gcam-v5.3/output/FinalRuns' # Use if gcamdatabase is needed
gcamdatabaseName_i <- 'gcam5p3-stash_GFDL-ESM2M_rcp2p6' # "Reference_originalSW" Use if gcamdatabse is needed
dataProjPath_i <- getwd() # Path to dataProj file.
dataProj_i <- "dataProj_gcam5p3-stash_GFDL-ESM2M_rcp2p6.proj" # "dataProj_gcam5p1_HadGEM2-ES_rcp8p5.proj"  # Use if gcamdata has been saved as .proj file

# Get list of scenarios and rename if desired.
# conn <- rgcam::localDBConn(gcamdatabasePath_i,gcamdatabaseName_i) # if connecting directly to gcam database
# prj <- rgcam::addScenario(conn, proj=paste(dataProjPath_i, "IDBNexus.proj", sep='/'), scenario=c('Reference', 'Impacts', 'Policy'))
# dataProjLoaded <- loadProject(paste(dataProjPath_i, "/",dataProj_i , sep = ""))
#  listScenarios(dataProjLoaded)  # List of Scenarios in GCAM database

scenOrigNames_i = c("Reference", "Impacts", "Policy") #  c(, "DeepDecarb1Mkt_2DS")  #'GCAMOriginal',
scenNewNames_i = c("Reference", "Climate Impacts", "Climate Policy") # "Original",   c(, "DeepDecarb1Mkt_2DS"")  # Names to replace the original names for final figures.

# Choose Parameters or set to "All" for all params. For complete list see ?metis.readgcam
paramsList <-  metis.mappings()$mapParamQuery
# paramsSelect_i <- c('aggLandAlloc')
# paramsSelect_i <- c('watWithdrawBySec')
# paramsSelect_i <- c('watWithdrawBySec',
#                     'watWithdrawByCrop',
#                     'watIrrWithdrawBasin',
#                     'landAlloc',
#                     'landAllocByCrop',
#                     'elecByTechTWh',
#                     'elecNewCapCost',
#                     "elecCapByFuel",
#                     "elecFinalBySecTWh",
#                     "elecFinalByFuelTWh",
#                     "elecNewCapGW",
#                     "elecAnnualRetPrematureCost",
#                     "elecAnnualRetPrematureGW",
#                     "elecCumCapCost",
#                     "elecCumCapGW",
#                     "elecCumRetPrematureCost",
#                     "elecCumRetPrematureGW",
#                     'energyPrimaryByFuelEJ',
#                     "energyPrimaryRefLiqProdEJ",
#                     'energyFinalConsumBySecEJ',
#                     "energyFinalByFuelBySectorEJ",
#                     "energyFinalSubsecByFuelTranspEJ",
#                     "energyFinalSubsecByFuelBuildEJ",
#                     "energyFinalSubsecByFuelIndusEJ",
#                     "energyFinalSubsecBySectorBuildEJ",
#                     'agProdByCrop'
#                     )


if(T){
  ## If selecting all parameters, switch to if(T)
  paramesSelect_energy <- c("energyPrimaryByFuelEJ","energyPrimaryRefLiqProdEJ", "energyFinalConsumBySecEJ","energyFinalByFuelBySectorEJ","energyFinalSubsecByFuelTranspEJ", "energyFinalSubsecByFuelBuildEJ", "energyFinalSubsecByFuelIndusEJ","energyFinalSubsecBySectorBuildEJ", "energyPrimaryByFuelMTOE","energyPrimaryRefLiqProdMTOE", "energyFinalConsumBySecMTOE","energyFinalbyFuelMTOE","energyFinalSubsecByFuelTranspMTOE", "energyFinalSubsecByFuelBuildMTOE", "energyFinalSubsecByFuelIndusMTOE","energyFinalSubsecBySectorBuildMTOE", "energyPrimaryByFuelTWh","energyPrimaryRefLiqProdTWh", "energyFinalConsumBySecTWh","energyFinalbyFuelTWh","energyFinalSubsecByFuelTranspTWh", "energyFinalSubsecByFuelBuildTWh", "energyFinalSubsecByFuelIndusTWh","energyFinalSubsecBySectorBuildTWh")
  
  paramsSelect_electricity <- c("elecByTechTWh","elecCapByFuel","elecFinalBySecTWh","elecFinalByFuelTWh", "elecNewCapCost","elecNewCapGW","elecAnnualRetPrematureCost","elecAnnualRetPrematureGW","elecCumCapCost","elecCumCapGW","elecCumRetPrematureCost","elecCumRetPrematureGW")
  
  paramsSelect_transport <- c( "transportPassengerVMTByMode", "transportFreightVMTByMode", "transportPassengerVMTByFuel", "transportFreightVMTByFuel")
  
  paramsSelect_water <- c("watConsumBySec", "watWithdrawBySec", "watWithdrawByCrop", "watBioPhysCons", "watIrrWithdrawBasin","watIrrConsBasin")
  
  paramsSelect_socioecon <- c("gdpPerCapita", "gdp", "gdpGrowthRate", "pop")
  
  paramsSelect_ag <- c("agProdbyIrrRfd", "agProdBiomass", "agProdForest","agProdByCrop")
  
  paramsSelect_livestock <- c("livestock_MeatDairybyTechMixed","livestock_MeatDairybyTechPastoral","livestock_MeatDairybyTechImports", "livestock_MeatDairybySubsector")
  
  paramsSelect_land <- c("landIrrRfd", "landIrrCrop","landRfdCrop", "landAlloc","landAllocByCrop")
  
  paramsSelect_emissions <- c("emissLUC", "emissNonCO2BySectorGWPAR5","emissNonCO2BySectorGTPAR5","emissNonCO2BySectorOrigUnits", "emissNonCO2ByResProdGWPAR5", "emissBySectorGWPAR5FFI","emissMethaneBySourceGWPAR5", "emissByGasGWPAR5FFI", "emissByGasGWPAR5LUC", "emissBySectorGWPAR5LUC", "emissNonCO2ByResProdGTPAR5", "emissBySectorGTPAR5FFI","emissMethaneBySourceGTPAR5", "emissByGasGTPAR5FFI", "emissByGasGTPAR5LUC","emissBySectorGTPAR5LUC", "emissCO2BySectorNoBio")
  
  paramsSelect_all <- c(paramesSelect_energy, paramsSelect_electricity, paramsSelect_transport, paramsSelect_water, paramsSelect_socioecon, paramsSelect_ag, paramsSelect_livestock, paramsSelect_land, paramsSelect_emissions)
  
  paramsSelect_i <- paramsSelect_all
}


queriesSelect_i <- c("All")

# Select regions from the 32 GCAM regions.
# regionsSelect_i <- c("Colombia", "Argentina", "Uruguay")
# regionsSelect_i <- c("Argentina")
regionsSelect_i <- c("Colombia")


# Reading in the no bio query so it works with Rgcam
# dataGCAM <- metis.readgcam(#dataProjFile = paste(dataProjPath_i, dataProj_i, sep = "/"),
#                            gcamdatabase = paste(gcamdatabasePath_i, gcamdatabaseName_i, sep='/'),
#                            reReadData = TRUE,
#                            scenOrigNames = scenOrigNames_i,
#                            scenNewNames = scenNewNames_i,
#                            regionsSelect = regionsSelect_i,
#                            paramsSelect = paramsSelect_i,
#                            queryFile = NULL)

if(file.exists(paste(dataProjPath_i, 'outputs', dataProj_i, sep = "/"))){
  dataGCAM <- metis.readgcam(dataProjFile = paste(dataProjPath_i, 'outputs', dataProj_i, sep = "/"),
                             scenOrigNames = scenOrigNames_i,
                             scenNewNames = scenNewNames_i,
                             regionsSelect = regionsSelect_i,
                             paramsSelect = paramsSelect_i)
}else{
  dataGCAM <- metis.readgcam(gcamdatabase = paste(gcamdatabasePath_i, gcamdatabaseName_i, sep = '/'),
                             reReadData = TRUE,
                             scenOrigNames = scenOrigNames_i,
                             scenNewNames = scenNewNames_i,
                             regionsSelect = regionsSelect_i,
                             paramsSelect = paramsSelect_i,
                             queryFile = NULL,
                             nameAppend = gcamdatabaseName_i)
  file.rename(paste(getwd(),'outputs', 'dataProj.proj', sep = '/'),
              paste(getwd(), 'outputs', dataProj_i, sep = '/'))
  
}

if(file.exists(paste(dataProjPath_i, 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = "/"))){
  invest <- plutus::gcamInvest(dataProjFile = paste(dataProjPath_i, 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = "/"),
                               reReadData = T,
                               scenOrigNames = scenOrigNames_i,
                               scenNewNames = scenNewNames_i,
                               regionsSelect = regionsSelect_i,
                               saveData = F)
}else{
  invest <- plutus::gcamInvest(gcamdatabase = paste(gcamdatabasePath_i, gcamdatabaseName_i, sep = '/'),
                               reReadData = T,
                               scenOrigNames = scenOrigNames_i,
                               scenNewNames = scenNewNames_i,
                               regionsSelect = regionsSelect_i,
                               saveData = F)
  file.rename(paste(getwd(),'outputs', 'dataProj.proj', sep = '/'),
              paste(getwd(), 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = '/'))
}

dataGCAM$data # To view the data read that was read.
invest$data
#------------------------------------------------------------------------------------------
# Charts Process (metis.chartsProcess.R)
#------------------------------------------------------------------------------------------

# Can also add data .csv outputs from metis.readgcam.R which are autmatically saved in
# ./metis/outputs/readGCAMTables/Tables_gcam
# for each of the regions selected.
# gcamDataTable_Argentina.csv, gcamDataTable_China.csv, gcamDataTable_Pakistan.csv
# This would be added to dataTables_i as:
# dataTables_i = c(paste(getwd(), "/outputs/readGCAMTables/Tables_local/local_Regional_Colombia.csv", sep = "")
#                  #paste(getwd(), "/outputs/readGCAMTables/Tables_gcam/gcamDataTable_Colombia.csv", sep = "")
# )

# Read in the data from the function metis.readgcam.
rTable_i <- rbind(dataGCAM$data, invest$data)
rTable_i <- rTable_i[complete.cases(rTable_i$class1),] # remove rows with NA values

# Get classes that have color codes in metis color palette pal_metis
class_metis <- unlist(attributes(metis.colors()$'pal_metis'), use.names=FALSE)
# Find parameters no have all classes in calss_metis, but the assigned palette is pal_metis
param_pal_16 <- unique(rTable_i$param[!(rTable_i$class1 %in% class_metis) & rTable_i$classPalette1 == 'pal_metis'])
param_pal_16_select <- c('livestock_MeatDairybyTechMixed', 'livestock_MeatDairybySubsector', 'livestock_MeatDairybyTechPastoral', 'watConsumBySec', 'transportPassengerVMTByMode', 'transportFreightVMTByMode', 'watConsumBySec')
rTable_i$classPalette1[rTable_i$param %in% param_pal_16_select] <- 'pal_16'

rTable_i$classPalette1[rTable_i$param %in% c('watWithdrawBySec', 'watBioPhysCons')] <- 'pal_metis'

# For GCAM 5.1
# rTable_i$class1[grepl('biomass_tree|biomass_grass', rTable_i$class1)] <- 'biomass'

# For GCAM 5.3
rTable_i$class1[grepl('biomassTree|biomassGrass', rTable_i$class1)] <- 'biomass'
rTable_i$class1[grepl('naturalOtherTree|naturalOtherGrass', rTable_i$class1)] <- 'biomass'
rTable_i$class1[grepl('^landAlloc$', rTable_i$param) & grepl('RootTuber|biomass', rTable_i$class1)] <- 'crops'
# rTable_i$class1[grepl('RootTuber', rTable_i$class1[rTable_i$param == 'landAlloc'])] <- 'crops'
  
# Choose Parameters or set to "All" for all params. For complete list see ?metis.chartsProcess

# paramsSelect_i <- c("finalNrgbySec", "TranspFinalNrgByFuel", "BuildFinalNrgByFuel",
#                     "IndFinalNrgByFuel", "primNrgConsumByFuel", "elecByTech", "watWithdrawBySec",
#                     "aggLandAlloc", "LUCemiss", "nonco2emissionBySectorGWPAR5",
#                     "finalNrgbyFuel","finalElecbySec","finalElecbyFuel",
#                     "NonCo2EmissionsByResProdGWPAR5",
#                     "TotalFFIEmissBySec", "CO2BySector_NonCO2Gases_GWPAR5", "CO2BySector_NonCO2Gases_GWPAR5_LUC",
#                     "TotalEmissBySec", "LandAllocByCrop", "MethaneBySource", "PassengerVMTByMode", "FreightVMTByMode",
#                     "BuildFinalNrgBySector",
#                     "co2emissionBySectorNoBio", "PassengerVMTByFuel", "FreightVMTByFuel", "RefiningByLiq")
# paramsSelect_i <- "All"
# paramsSelect_i <- "elecByTechTWh"  # "elecInvest"
# paramsSelect_i <- c('landAlloc', 'landAllocByCrop', 'landIrrCrop', 'landRfdCrop')
# paramsSelect_i <- param_pal_16_select
# paramsSelect_i <- c('elecAnnualRetPrematureCost', 'elecAnnualRetPrematureGW')

# Select regions from the 32 GCAM regions.
# paramsSelect_i <- c('BuildFinalNrgBySector')
# Charts Process
#regionsSelect_i <- c('Colombia')
# paramsSelect_i <- c('elecNewCapCost')

charts <- metis.chartsProcess(rTable = rTable_i, # Default is NULL
                              #dataTables=dataTables_i, # Default is NULL
                              paramsSelect = paramsSelect_i, # Default is "All"
                              regionsSelect = regionsSelect_i, # Default is "All"
                              xCompare = c("2020", "2030", "2040", "2050"), # Default is c("2015","2030","2050","2100")
                              scenRef = "Reference", # Default is NULL
                              dirOutputs = paste(getwd(), "/outputs", sep = ""), # Default is paste(getwd(),"/outputs",sep="")
                              regionCompareOnly = 0, # Default 0. If set to 1, will only run comparison plots and not individual
                              scenarioCompareOnly = 0,
                              regionCompare = 0,
                              useNewLabels = 0,
                              folderName = "IDBNexus",
                              xRange = c(2020, 2030, 2040, 2050), # c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050), #
                              colOrder1 = c("Reference", "Climate Impacts", "Climate Policy"), #"Original",
                              colOrderName1 = "scenario",
                              pdfpng = 'pdf',
                              # sizeBarLines = 1,
                              facetLabelSize = 36,
                              # plotBGcolor = 'white',
                              multiPlotFigsOnly = T) # Default 0. If set to 1, will only run comparison plots and not individual


