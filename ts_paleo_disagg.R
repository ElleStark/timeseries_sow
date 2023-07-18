# Space-time disagreggation for Paleo flows within time series selections
# Elle Stark July 2023

# see https://github.com/rabutler-usbr/knnstdisagg.git

# Install disagreggation package and CO River Natural flow package if needed: 
#library(devtools)
#install_github('rabutler-usbr/knnstdisagg')
#install_github('rabutler-usbr/CoRiverNF')


library(knnstdisagg)
library(CoRiverNF)
library(tidyverse)

############ OBTAIN AND SET UP DATA ########################

# Annual flow data
flow_data <- readRDS('Data/inputs/annual_flow_long.rds')

# List of folder names for finding data and creating CRSS input files
selected_sow_list <- c('badHydro_badSCD', 'badHydro_goodSCD')

# MAYBE START LOOP HERE TO GO THROUGH EACH FOLDER AUTOMATICALLY
# setwd(paste0('Data/outputs/paleo_disagg/', selected_sow_list[m]))

# sow as selected in 'ts_clustering.R'
all_selected_sow <- readRDS('Data/outputs/opt_sow_badH_badIC.rds')

# Find sow index for paleo flows (already have disaggregated data for other sources)
paleo_idx <- all_selected_sow %>%
  filter(Scenario=='PCNF' | Scenario=='DPNF') 
paleo_idx <- paleo_idx$sow_idx

# Select annual flow data for desired Paleo traces
ann_flow_all <- flow_data %>%
  filter(sow_idx %in% paleo_idx) %>%
  dplyr::select(sow_idx, Year, LF_Annual)

# setup annual data (see knnstdisagg documentation)
annual_index <- CoRiverNF::wyAnnTot$LeesFerry
yrs <- as.numeric(format(index(wyAnnTot$LeesFerry), "%Y"))
annual_index <- as.matrix(annual_index)
annual_index <- cbind(yrs, annual_index)

# setup monthly data (see knnstdisagg documentation)
last_month <- paste0("/", max(yrs), "-09")
monthly_data <- CoRiverNF::monthlyInt[last_month]  


##### Test by selecting one trace before setting up loop
ann_flow <- ann_flow_all[which(ann_flow_all$sow_idx==paleo_idx[1]),]
ann_flow$Year <- as.numeric(ann_flow$Year)
ann_flow$LF_Annual <- ann_flow$LF_Annual*1000000
ann_flow <- as.matrix(dplyr::select(ann_flow, -sow_idx))

# create disaggregation object
disagg <- knn_space_time_disagg(
  ann_flow = ann_flow,
  ann_index_flow = annual_index,
  mon_flow = monthly_data,
  start_month = 10,
  nsim = 1,
  scale_sites = 1:20,
  k_weights = knn_params_default(nrow(annual_index))
)

tail(knnst_get_disagg_data(disagg)[,5:10])
i=1
write_knnst(disagg, 'Data/outputs/paleo_disagg/badHydro_badIC')

disagg_df <- read.csv('Data/outputs/paleo_disagg/badHydro_badIC/disagg_1.csv')
disagg_df <- round(disagg_df[,-1], digits=4)

#### Create files for CRSS input
# Write header lines into txt file (CHANGE SIMULATION INFO IF NEEDED)
writeLines(c('start_date: 2024-1-31 24:00', 'units: acre-ft/month'), 
           'Data/outputs/paleo_disagg/badHydro_badSCD/headertext.txt')

for(i in 1:ncol(disagg_df)){
  # Get filename
  inflow_pt <- colnames(disagg_df[i])
  inflow_file <- paste0('Data/outputs/paleo_disagg/badHydro_badSCD/', inflow_pt, 'NF.inflow')
  
  # Save flows for one column as txt 
  f_temp <- paste0('Data/outputs/paleo_disagg/badHydro_badSCD/', inflow_pt, 'flow_only.txt')
  write.table(disagg_df[i], file = f_temp, row.names = FALSE, col.names = FALSE)
  
  # Append header info 
  file.copy('Data/outputs/paleo_disagg/badHydro_badSCD/headertext.txt', 
            inflow_file)
  file.append(inflow_file, f_temp)
  
  # Delete temporary file
  file.remove(f_temp)
}



  
  