# Space-time disagreggation for Paleo flows within time series selections
# Elle Stark July 2023

# see https://github.com/rabutler-usbr/knnstdisagg.git

# Install disagreggation package and CO River Natural flow package if needed: 
#library(devtools)
#install_github('rabutler-usbr/knnstdisagg')
install_github('rabutler-usbr/CoRiverNF')


library(knnstdisagg)
library(CoRiverNF)
library(tidyverse)

############ OBTAIN AND SET UP DATA ########################
# sow as selected in 'ts_clustering.R'
all_selected_sow <- readRDS('Data/outputs/opt_sow_badH_badIC.rds')
flow_data <- readRDS('Data/inputs/annual_flow_long.rds')

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

head(knnst_get_disagg_data(disagg)[,5:10])

write_knnst(disagg, 'Data/outputs/paleo_disagg/badHydro')
  
  