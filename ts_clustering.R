# Given initial SOW set, cluster based on cumulative annual LF flow time series
# Elle Stark July 2023

library(tidyverse)
library(ggplot2)
library(GGally)
library(TSclust)
library(cluster)
library(purrr)

################## FUNCTIONS #########################

### function to normalize a single column (or vector) ###
# such that values are between 0 and 1
normalize_col <- function(x){
  (x - min(x))/(max(x)-min(x))
}

### Function to convert distance object to dataframe with each row as one pairwise distance 
# Copied from: https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
dist_2_df <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    distance = as.vector(inDist))
}

################# OBTAIN DATA #####################
initial_sow_set <- readRDS('Data/inputs/initial_sow_set_uclhs_500.rds')

# Read in Lee Ferry annual flow data
flow_data_alltraces <- read.csv('Data/inputs/hydrology_all_annual.csv')
# Filter to 2024-2056 (simulation period)
sim_start = 2024
sim_end = 2056
flow_data_alltraces <- subset(flow_data_alltraces, Year>=sim_start & Year<=(sim_end+1)) %>%
  dplyr::select(-c('X', 'ID'))

wide_flow_data_alltraces <- pivot_wider(flow_data_alltraces, 
                                        names_from = Year, values_from = CNF_LF)

# Look up annual flow data for traces included in initial sow set
# Loop through manually to include duplicated traces in SOW set
flow_data_wide <- data.frame()

for(i in 1:nrow(initial_sow_set)){
  scen <- initial_sow_set[i, 'Scenario']
  trace <- initial_sow_set[i, 'TraceNumber']
  tracedata <- filter(wide_flow_data_alltraces, 
                      Scenario==scen & TraceNumber==trace)
  flow_data_wide <- rbind(flow_data_wide, tracedata)
}

flow_data_wide <- mutate(flow_data_wide, sow_idx=seq(1, 500, by=1))

flow_data_long <- pivot_longer(flow_data_wide, 
                               cols = -c('Scenario', 'TraceNumber', 'sow_idx'), 
                               names_to = 'Year', values_to = 'LF_Annual')

# Save annual flow data for initial SOW set
saveRDS(flow_data_wide, 'Data/inputs/annual_flow_wide.rds')
saveRDS(flow_data_long, 'Data/inputs/annual_flow_long.rds')

################ TIME SERIES CLUSTERING #############

# Calculate cumulative flow at LF for each year, wide and long format
cumulative_flow_df_long <- flow_data_long %>%
  group_by(sow_idx) %>%
  mutate(LF_cumulative=cumsum(LF_Annual)) %>%
  ungroup() %>%
  select(Scenario, TraceNumber, sow_idx, Year, LF_cumulative) 

cumulative_flow_df_wide <- pivot_wider(cumulative_flow_df_long, names_from = Year, 
                                       values_from = LF_cumulative)

saveRDS(cumulative_flow_df_wide, 'Data/outputs/cumulative_flow_wide.rds')
saveRDS(cumulative_flow_df_long, 'Data/outputs/cumulative_flow_long.rds')

# Transpose data frame for use in CID calculation. cnf = cumulative natural flow
timeseries_cnf <- as.data.frame(t(select(cumulative_flow_df_wide, 
                                         -Scenario, -TraceNumber, -sow_idx))) 
# Test clustering on annual flow data
timeseries_annual <- as.data.frame(t(select(flow_data_wide, 
                                            -Scenario, -TraceNumber, -sow_idx)))


# Create dissimilarity matrix using TSclust package - 
# CID (complexity-invariant distance) method with PAM (partitioning around medoids)
# See Steinmann et al., 2020, for use of CID for Scenario Discovery
flow_dist_cid <- diss(timeseries_cnf, 'CID')

####### SELECT NUMBER OF CLUSTERS #############

# Test a range of number of clusters and analyze using elbow and silhouette methods
# get metrics of interest for plotting
avg_diss <- c()
max_diss <- c()
diameter <- c()
separation <- c()

for(i in 2:10){
  clusters <- pam(flow_dist_cid, i, diss = TRUE)
  
  cluster_info <- clusters$clusinfo
  avg_diss[i-1] <- mean(cluster_info[,'av_diss'])
  max_diss[i-1] <- mean(cluster_info[,'max_diss'])
  diameter[i-1] <- mean(cluster_info[,'diameter'])
  separation[i-1] <- mean(cluster_info[,'separation'])

  # Cumulative flow data, clustered
  flow_clustered_wide <- mutate(cumulative_flow_df_wide,
                                Cluster = as.character(clusters$clustering))

  flow_clustered_long <- pivot_longer(flow_clustered_wide,
                                 cols = -c('Scenario', 'TraceNumber', 'sow_idx', 'Cluster'),
                                 names_to = 'Year', values_to = 'LF_cumulative')

  # Annual flow data, clustered
  cluster_idx <- select(flow_clustered_wide, sow_idx, Cluster)
  annual_clustered_long <- left_join(flow_data_long, cluster_idx, by='sow_idx')

  # Plot cumulative flow clustered
  cumulative_plot <- ggplot(flow_clustered_long, aes(x=Year, y=LF_cumulative, col=Cluster,
                                                group=sow_idx)) +
    geom_path() +
    ylab('Cumulative Flow at Lee Ferry (MAF)')
  f_name <- paste0('Data/temp/ts_clust_size/k.', i, '/cumulative_ts.png')
  ggsave(f_name, plot = cumulative_plot, width = 12, height = 7)

  # Plot individual clusters of annual LF flow
  for(j in 1:i){
    annual_plot <- ggplot(annual_clustered_long[annual_clustered_long$Cluster==as.character(j),],
                          aes(x=Year, y=LF_Annual, group=sow_idx)) +
      geom_path(alpha=0.4) +
      ylab('Annual Flow at Lee Ferry (MAF)') +
      ylim(0, 45)
    f_name <- paste0('Data/temp/ts_clust_size/k.', i, '/annual_ts_clust', j, '.png')
    ggsave(f_name, plot = annual_plot, width = 12, height = 7)
  }

  # Plot silhouette plot of clusters.
  # Need to save manually due to RStudio issues with silhouette plotting
  windows()
  plot(clusters)
}

clust_metrics_df <- data.frame(avg_diss, max_diss, diameter, separation)
clust_metrics_df <- clust_metrics_df %>%
  mutate(across(everything(), normalize_col)) 

n_clust <- 2:10
elbow_df <- as.data.frame(cbind("n_clust" = n_clust, clust_metrics_df))
elbow_df_long <- pivot_longer(elbow_df, !n_clust, names_to = 'metric', values_to = 'value')

elbow_plot <- ggplot(elbow_df_long, aes(x=n_clust, y=value, col=metric)) +
  geom_line() +
  ylab('normalized metric value') +
  xlab('number of clusters') +
  scale_x_continuous(breaks = seq(1,10,1)) +
  theme(panel.grid.minor = element_blank())

############### End number of clusters analysis #############

# select 4 clusters based on above results
k=4

# final clustering
clusters <- pam(flow_dist_cid, k, diss = TRUE)

# Cumulative flow data, clustered
flow_clustered_wide <- mutate(cumulative_flow_df_wide,
                              Cluster = as.character(clusters$clustering))

flow_clustered_long <- pivot_longer(flow_clustered_wide,
                                    cols = -c('Scenario', 'TraceNumber', 'sow_idx', 'Cluster'),
                                    names_to = 'Year', values_to = 'LF_cumulative')
saveRDS(flow_clustered_wide, 'Data/outputs/cumulative_flow_wide_clustered.RDS')

# Cluster indexes
cluster_idx <- select(flow_clustered_wide, sow_idx, Cluster)

# Save cluster number for each SOW index 
saveRDS(cluster_idx, 'Data/outputs/ts_clusters_idx.rds')

n_clust <- 4


### Create optimization groups that take into account demand and initial conditions 

# Label SOW uncertainty dataframe with clusters
clustered_sow <- mutate(initial_sow_set, 'Cluster'=cluster_idx$Cluster, 'sow_idx'=seq(1,500, by=1))

# Pairwise scatter of clusters mapped onto uncertainty metrics
cluster_pairs <- ggpairs(dplyr::select(clustered_sow, -Scenario, -TraceNumber, -sow_idx), 
                         ggplot2::aes(col=Cluster, alpha=0.5))


### Create 4 optimization SOW sets for Borg Run: 
# 1. Good hydrology, demand, and initial conditions
# 2. Good hydrology, bad demand, and bad initial conditions
# 3. Bad hydrology, good demand, and good initial conditions
# 4. Bad hydrology, bad demand, and bad initial conditions

# Boxplots of Demand, Mead, Powell within each cluster
# May help determine Good/Bad thresholds
mead_boxplot <- ggplot(data=clustered_sow) +
  geom_boxplot(mapping = aes(y=mead_pe, fill=Cluster)) +
  theme_bw() +
  ylab('Mead Pool Elevation') +
  theme(axis.text.x = element_blank())

powell_boxplot <- ggplot(data=clustered_sow) +
  geom_boxplot(mapping = aes(y=powell_pe, fill=Cluster)) +
  theme_bw() +
  ylab('Powell Pool Elevation') +
  theme(axis.text.x = element_blank())

demand_boxplot <- ggplot(data=clustered_sow) +
  geom_boxplot(mapping = aes(y=demand, fill=Cluster)) +
  theme_bw() +
  ylab('Demand (MAF)') +
  theme(axis.text.x = element_blank())

# Test values for Good/Bad Threshold
mead_bad <- c(1000, 1025, 1050, 1060, 1075, 1090)
powell_bad <- c(3490, 3525, 3535, 3550, 3575, 3590, 3600)
demand_bad <- c(4.5, 4.75, 5.0, 5.1, 5.25, 5.5)

# Create dataframe with number of 'bad' SOW for each combo of thresholds
bad_sow_list_all <- list() # list of bad sow cluster lists
bad_sow_list_c <- list() # bad sow list for each cluster
good_sow_list_c <- list()
good_sow_list_all <- list()
thresholds <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(thresholds) <- c('Mead', 'Powell', 'Demand')

for (c in 1:n_clust){
  cluster_sow <- filter(clustered_sow, Cluster==c)
  for(i in 1:length(mead_bad)){
    bad_sow_m <- filter(cluster_sow, mead_pe<=mead_bad[i])
    good_sow_m <- filter(cluster_sow, mead_pe>mead_bad[i])
    for(j in 1:length(powell_bad)){
      bad_sow_p <- filter(bad_sow_m, powell_pe<=powell_bad[j])
      good_sow_p <- filter(good_sow_m, powell_pe>powell_bad[j])
      for(k in 1:length(demand_bad)){
        bad_sow <- filter(bad_sow_p, demand>=demand_bad[k])
        good_sow <- filter(good_sow_p, demand<demand_bad[k])
        bad_sow_list_c <- append(bad_sow_list_c, list(bad_sow))
        good_sow_list_c <- append(good_sow_list_c, list(good_sow))
        thresholds[nrow(thresholds)+1,] <- c(mead_bad[i], powell_bad[j], demand_bad[k])
      }
    }
  }
  bad_sow_list_all <- append(bad_sow_list_all, list(bad_sow_list_c))
  good_sow_list_all <- append(good_sow_list_all, list(good_sow_list_c))
  bad_sow_list_c <- list()
  good_sow_list_c <- list()
}

bad_idx_list <- list()
good_idx_list <- list()
for (i in 1:n_clust){
  # find index of SOW lists in each cluster with at least 8 rows
  # see if other clusters have any matching indexes
  b_sow_list <- bad_sow_list_all[[i]]
  bad_idx_list[[i]] <- which(sapply(b_sow_list, function(x) nrow(x)>=8))
  g_sow_list <- good_sow_list_all[[i]]
  good_idx_list[[i]] <- which(sapply(g_sow_list, function(x) nrow(x)>=8))
}

common_idx_b <- Reduce(intersect, bad_idx_list) # indexes for which all clusters have at least 8 'bad' SOW
common_idx_g <- Reduce(intersect, good_idx_list) # indexes for which all clusters have at least 8 'good' SOW
common_idx <- intersect(common_idx_g, common_idx_b) # overlap in indexes
thresh_opt <- thresholds[common_idx,] # Options for threshold values for enough 'good' and 'bad' in each cluster

# Choose approx. median values: 1060, 3590, 5.1 (index 160)
selected_good_sow = list()
selected_bad_sow = list()
thresh_idx = 160

for(i in 1:n_clust){
  selected_good_sow[[i]] <- good_sow_list_all[[i]][thresh_idx]
  selected_bad_sow[[i]] <- bad_sow_list_all[[i]][thresh_idx]
}

sow_bad_df <- rbind(selected_bad_sow[[1]][[1]], selected_bad_sow[[2]][[1]], 
                    selected_bad_sow[[3]][[1]], selected_bad_sow[[4]][[1]])
sow_bad_df['ic_cat'] <- 'bad' # ic_cat stands for initial conditions & demand category

sow_good_df <- rbind(selected_good_sow[[1]][[1]], selected_good_sow[[2]][[1]], 
                    selected_good_sow[[3]][[1]], selected_good_sow[[4]][[1]])
sow_good_df['ic_cat'] <- 'good' # ic_cat stands for initial conditions & demand category

categorized_sow <- rbind(sow_bad_df, sow_good_df)
ic_cat_idx <- dplyr::select(categorized_sow, sow_idx, ic_cat)
saveRDS(categorized_sow, 'Data/outputs/SOW_metrics_ICcat.rds')
saveRDS(ic_cat_idx, 'Data/outputs/IC_cat_sow_idx.rds')

# Get long format LF cumulative flow data for above sow categories
cat_flow_long <- inner_join(flow_clustered_long, ic_cat_idx, by='sow_idx')
cat_flow_wide <- inner_join(flow_clustered_wide, ic_cat_idx, by='sow_idx')
#cat_flow_long <- left_join(flow_clustered_long, ic_cat_idx, by='sow_idx') %>% 
#  mutate(ic_cat=replace_na(ic_cat, 'a_mix'))

# Plot each cluster showing 'good' and 'bad' demand/i.c.
cluster_plot_IC <- ggplot(cat_flow_long[cat_flow_long$Cluster=='1',], aes(x=Year, y=LF_cumulative, col=ic_cat,
                                                   group=sow_idx)) +
  geom_path() +
  ylab('Cumulative Flow at Lee Ferry (MAF)') +
  ylim(0, 600)




############# SUBSET GROUPS FOR OPTIMIZATION ###################

# use 8 sow to match number in baseline optimization
# also allows for 18-core machine to run 2 simultaneous Borg runs 
n_sow = 8 

### Need to go through each group and select traces
# 1. Good hydrology, demand, and initial conditions: 
      # Cluster 4, good ic_cat
# 2. Good hydrology, bad demand, and bad initial conditions:
      # Cluster 4, bad ic_cat
# 3. Bad hydrology, good demand, and good initial conditions
      # Cluster 1, good ic_cat
# 4. Bad hydrology, bad demand, and bad initial conditions:
      # Cluster 1, bad ic_cat
badH_badIC_df <- categorized_sow[which(categorized_sow$Cluster=='1' & categorized_sow$ic_cat=='bad'),]
badH_goodIC_df <- categorized_sow[which(categorized_sow$Cluster=='1' & categorized_sow$ic_cat=='good'),]
goodH_badIC_df <- categorized_sow[which(categorized_sow$Cluster=='4' & categorized_sow$ic_cat=='bad'),]
goodH_goodIC_df <- categorized_sow[which(categorized_sow$Cluster=='4' & categorized_sow$ic_cat=='good'),]

df_groups_list <- list('badH_badIC_df'=badH_badIC_df, 'badH_goodIC_df'=badH_goodIC_df, 
                       'goodH_badIC_df'=goodH_badIC_df, 'goodH_goodIC_df'=goodH_goodIC_df)
# Labels used for saving results later
group_labels <- c('badH_badIC', 'badH_goodIC', 'goodH_badIC', 'goodH_goodIC')

bad_hydr_idx <- c(1,2)
good_hydr_idx <- c(3,4)

# Loop through each group to select based on hydrology and plot results
for(i in 1:length(df_groups_list)){
  df = df_groups_list[[i]]
  
  ggpairs(dplyr::select(badH_badIC_df, 
                        -Scenario, -TraceNumber, -sow_idx, -ic_cat, -Cluster))
  
  flow_df <- filter(flow_clustered_wide, sow_idx %in% df$sow_idx)
  

  flow_for_dist <- as.data.frame(t(select(flow_df, -Cluster,
                                          -Scenario, -TraceNumber, -sow_idx)))
  
  # Calculate pairwise distances in group based on time series of LF cumulative flow
  cid_matrix <- diss(flow_for_dist, 'CID')
  # Convert dissimilartiy matrix to dataframe to process - function at top of script
  cid_df <- dist_2_df(cid_matrix)
  
  # If in quadrant with 'good' hydrology, choose traces closest to trace with max cumulative LF flow
  if(i %in% good_hydr_idx){
    flow_idx <- which.max(as.matrix(flow_df[,(ncol(flow_df)-1)]))
  }
  # If in quadrant with 'bad' hydrology, choose traces closest to trace with min cumulative LF flow
  if(i %in% bad_hydr_idx){
    flow_idx <- which.min(as.matrix(flow_df[,(ncol(flow_df)-1)]))
  }
  
  # Filter dissimilarity matrix to only comparisons with selected extreme SOW
  compare_df <- filter(cid_df, (row==paste0('V',flow_idx)) | col==paste0('V', flow_idx))
  # Select nearest traces to get desired number of SOW
  nearest <- compare_df %>% slice_min(distance, n=n_sow-1)
  # Process list to get SOW index for referencing other dataframes
  nearest_list <- as.vector(as.matrix(nearest[,c('row', 'col')])) 
  nearest_list <- nearest_list[!nearest_list==paste0('V',flow_idx)]
  nearest_list <- as.integer(sub('V', '', nearest_list))
  selected_temp_idx <- c(flow_idx, nearest_list)
  selected_flows <- flow_df[selected_temp_idx,]
  selected_sow_idx <- selected_flows$sow_idx
  
  # FINAL DATA FRAME of SOW for that group, labeling those selected for optimization
  df <- mutate(df, 'selected'=(df$sow_idx %in% selected_sow_idx))
  
  # Save selected SOW as RDS and csv
  selected_df <- df[which(df$selected==TRUE),]
  f_name <- paste0('Data/outputs/opt_sow_', group_labels[i], '.rds')
  saveRDS(selected_df, f_name)
  f_name <- paste0('Data/outputs/opt_sow_', group_labels[i], '.csv')
  write.csv(selected_df, f_name)
  
  # plot selected flows - pairwise
  pairs_selected <- ggpairs(dplyr::select(df, -Scenario, -TraceNumber, -sow_idx, 
                        -ic_cat, -Cluster), ggplot2::aes(col=selected))
  
  f_name <- paste0('Data/outputs/pairs_group', i, '.png')
  ggsave(f_name, pairs_selected, width = 12, height = 7)
  
  # plot selected flows - cumulative
  cluster_plot_selected <- ggplot() +
    geom_path(cat_flow_long[cat_flow_long$sow_idx %in% df$sow_idx,], mapping=aes(x=Year, y=LF_cumulative,
                                                                         group=sow_idx), col='darkgray') +
    geom_path(data = cat_flow_long[cat_flow_long$sow_idx %in% selected_sow_idx,], mapping=aes(x=Year, y=LF_cumulative,
                                                                                group=sow_idx), col='red') +
    ylab('Cumulative Flow at Lee Ferry (MAF)') +
    ylim(0, 600)
  
  f_name <- paste0('Data/outputs/flow_group', i, '.png')
  ggsave(f_name, cluster_plot_selected, width = 11, height = 7)
}






