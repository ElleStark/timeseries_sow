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
setwd('C:/Users/elles/Documents/CU_Boulder/Research/ROAM_code/timeseries_sow/Data')
flow_data <- readRDS('inputs/annual_flow_long.rds')

# List of folder names for finding data and creating CRSS input files
sow_folder_list <- c('badHydro_badSCD', 'badHydro_goodSCD', 
                     'goodHydro_badSCD', 'goodHydro_goodSCD')

for(m in 1:length(sow_folder_list)){
  # LOOP THROUGH EACH FOLDER AUTOMATICALLY
  folder <- paste0('C:/Users/elles/Documents/CU_Boulder/Research/ROAM_code/timeseries_sow/Data/outputs/paleo_disagg/', 
                   sow_folder_list[m])
  
  
  # sow as selected in 'ts_clustering.R'
  all_selected_sow <- readRDS(paste0(folder,'/opt_sow.rds'))
  
  # Find sow index for paleo flows (already have disaggregated data for other sources)
  paleo_sow <- all_selected_sow %>%
    filter(Scenario=='PCNF' | Scenario=='DPNF') 
  paleo_idx <- paleo_sow$sow_idx
  scen_list <- paleo_sow$Scenario
  trace_list <- paleo_sow$TraceNumber
  
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
  
  ##### Loop through each paleo trace in the folder
  for(t in 1:length(paleo_idx)){
    # Get scenario, trace number, sow index
    scen <- scen_list[t]
    trace <- trace_list[t]
    sow <- paleo_idx[t]
    # Create folder to store CRSS trace input files
    dir.create(paste0(folder, '/', scen, trace, '_sow', sow, '/'))
    
    # Set up annual flow data based on built-in LF NF data from knnstdisagg 
    ann_flow <- ann_flow_all[which(ann_flow_all$sow_idx==sow),]
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
    
    #tail(knnst_get_disagg_data(disagg)[,5:10]) # Check data
    write_knnst(disagg, folder)
    
    disagg_df <- read.csv(paste0(folder,'/disagg_1.csv'))
    # Start file in January - uses water year by default
    disagg_df <- tail(disagg_df, -3)
    disagg_df <- round(disagg_df[,-1], digits=4) # exclude first column ('Year')
    
    #### Create files for CRSS input
    # Write header lines into txt file (CHANGE SIMULATION INFO IF NEEDED)
    writeLines(c('start_date: 2024-1-31 24:00', 'units: acre-ft/month'), 
               'headertext.txt')
    
    for(i in 1:ncol(disagg_df)){
      # Get filename
      inflow_pt <- colnames(disagg_df[i])
      inflow_file <- paste0(folder, '/', scen, trace, '_sow', sow, '/', inflow_pt, 'NF.inflow')
      
      # Save flows for one column as txt 
      f_temp <- paste0(folder, '/', inflow_pt, 'flow_only.txt')
      write.table(disagg_df[i], file = f_temp, row.names = FALSE, col.names = FALSE)
      
      # Append header info 
      file.copy('headertext.txt', 
                inflow_file)
      file.append(inflow_file, f_temp)
      
      # Delete temporary file
      file.remove(f_temp)
    }
  }
}



  
  