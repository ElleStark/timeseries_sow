# Given initial SOW set, cluster based on cumulative annual LF flow time series
# Elle Stark July 2023

library(tidyverse)
library(ggplot2)
library(GGally)

################# OBTAIN DATA #####################
initial_sow_set <- readRDS('Data/inputs/initial_sow_set_uclhs_500.rds')

# Read in Lee Ferry annual flow data
flow_data_alltraces <- read.csv('Data/inputs/hydrology_all_annual.csv')
# Filter to 2024-2056 (simulation period)
flow_data_alltraces <- subset(flow_data_alltraces, Year>=2024 & Year<=2056) %>%
  dplyr::select(-c('X', 'ID'))

wide_flow_data_alltraces <- pivot_wider(flow_data_alltraces, names_from = Year, values_from = CNF_LF)

# Look up annual flow data for traces included in initial sow set
flow_data_wide <- data.frame()

for(i in 1:nrow(initial_sow_set)){
  scen <- initial_sow_set[i, 'Scenario']
  trace <- initial_sow_set[i, 'TraceNumber']
  tracedata <- filter(wide_flow_data_alltraces, Scenario==scen & TraceNumber==trace)
  flow_data_wide <- rbind(flow_data_wide, tracedata)
}

# Save annual flow data for initial SOW set
saveRDS(flow_data_wide, 'Data/outputs/annual_flow_wide.rds')

################ TIME SERIES CLUSTERING #############



