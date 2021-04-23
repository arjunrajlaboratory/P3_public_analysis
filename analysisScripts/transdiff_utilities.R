bootstrap_high_cells <- function(per_cell_counts,             # data.frame with counts per cell in all wells
                                 threshold,                   # count threshold at or above which calling a cell "positive"
                                 negative_condition,          # name of condition corresponding to negative control wells, e.g., 'PO'
                                 base_condition,              # name of condition corresponding to reference wells, e.g., '7F only'
                                 n_iter = 1000,               # number of bootstrap iterations, default = 1000
                                 seed_val = 74004003,         # seed value for pseudorandom reproducibility, default = 74004003
                                 all_results_flag = F,        # whether to return full bootstrap resample value set (will be large), default FALSE
                                 CIwidth = 0.9){              # width of confidence intervals (e.g., 0.9 for a 90% CI), default = 0.9
  
  # check structure of dataset
  # needs to be a data.frame with columns: cellID, target, counts, wellID, condition
  if(sum(colnames(per_cell_counts) %in% c('cellID', 'target', 'wellID', 'condition', 'counts')) < 5) {
    cat(paste0('Error: per_cell_counts must contain columns:cellID, target, wellID, condition, counts\n'))
    break
  }
  
  if(sum(per_cell_counts$condition == negative_condition) == 0){
    cat(paste0('Error: no entries in condition column are the specified negative_condition, ', negative_condition, '\n'))
    break
  }
  
  if(sum(per_cell_counts$condition == base_condition) == 0){
    cat(paste0('Error: no entries in condition column are the specified base_condition, ', base_condition, '\n'))
    break
  }
 
  set.seed(seed_val)
  
  # results will be a list containing results
  # if all_results = T, includes full bootstrap simulation in final list entry
  # if all_results = F, includes only entries with:
  # - p-value for threshold-exceeding fraction larger than negative_condition average
  # - p-value for threshold-exceeding fraction larger than base_condition average
  # - resampled average and 90% CI on threshold-exceeding fraction
  # - fold-change in threshold-exceeding fraction over negative_condition, resampled average fold-change, and 90% CI on FC
  # - fold-change in threshold-exceeding fraction over base_condition, resampled average fold-change, and 90% CI on FC
  results <- list()
  
  cells_per_well <- per_cell_counts %>%
    group_by(wellID, condition) %>%
    summarise(Nc = length(wellID),
              obs_nGEthresh = sum(counts >= threshold),
              obs_fracGEthresh = obs_nGEthresh/Nc)
  
  cells_per_well$wellNum <- 1:nrow(cells_per_well)
    
  # main loop
  # for each iteration:
  # 1. resample each well's per-cell counts with replacement Nc times, where Nc is the number of cells in that well
  # 2. calculate fraction above specified threshold for each well
  all_resamps <- list()
  for (well in cells_per_well$wellNum) {
    
    temp_Nc <- cells_per_well %>% filter(wellNum == well)
    
    tempDat <- per_cell_counts %>% 
      filter(wellID == temp_Nc$wellID) 
    
    cat(paste0('Working on well ', as.character(temp_Nc$wellID), '...\n'))
    cat(paste0('Resampling ', as.character(temp_Nc$Nc), ' cells with replacement ', as.character(n_iter), ' times.\n'))
    temp_resamps <- matrix(0, temp_Nc$Nc, n_iter)
    colnames(temp_resamps) <- paste0('iter_', as.character(1:n_iter))
    
    for (i in 1:n_iter){
      temp_resamps[,i] <- sample(tempDat$counts, temp_Nc$Nc, replace = T)
    }
    
    temp_resamps <- as.data.frame(temp_resamps)
    temp_resamps$wellID <- temp_Nc$wellID
    
    if(is.null(dim(all_resamps))) {
      all_resamps <- temp_resamps
    } else {
      all_resamps %<>% bind_rows(temp_resamps)
    }
    
  }
  
  if(all_results_flag) {
    results[['all_resamples']] <- all_resamps
  }
  
  cat('Calculating statistics...\n')
  
  boot_thresh <- all_resamps %>%
    gather('iteration', 'resampled_count', 1:n_iter) %>%
    group_by(iteration, wellID) %>%
    summarise(nGEthresh = sum(resampled_count > threshold))
  
  lo_ci = (1 - CIwidth)/2
  hi_ci = 1 - lo_ci
  
  sum_boot_basicstats_perWell <- boot_thresh %>%
    inner_join(cells_per_well, by = 'wellID') %>%
    group_by(wellID, Nc, obs_nGEthresh, obs_fracGEthresh) %>%
    summarise(mean_resampled_nGEthresh = mean(nGEthresh),
              CI_lo_resampled_nGEthresh = quantile(nGEthresh, probs = lo_ci),
              CI_hi_resampled_nGEthresh = quantile(nGEthresh, probs = hi_ci)) %>%
    mutate(mean_resampled_fracGEthresh = mean_resampled_nGEthresh/Nc,
           CI_lo_resampled_fracGEthresh = CI_lo_resampled_nGEthresh/Nc,
           CI_hi_resampled_fracGEthresh = CI_hi_resampled_nGEthresh/Nc,
           CI_width = CIwidth)
  
  sum_boot_basicstats_perCondition <- boot_thresh %>%
    inner_join(cells_per_well, by = 'wellID') %>%
    mutate(fracGEthresh = nGEthresh/Nc) %>%
    group_by(iteration, condition) %>%
    summarise(mean_resampled_fracGEthresh = mean(fracGEthresh),
              mean_obs_fracGEthresh = mean(obs_fracGEthresh)) %>%
    group_by(condition, mean_obs_fracGEthresh) %>%
    summarise(mean_resampled_mean_fracGEthresh = mean(mean_resampled_fracGEthresh),
              CI_lo_resampled_mean_fracGEthresh = quantile(mean_resampled_fracGEthresh, probs = lo_ci),
              CI_hi_resampled_mean_fracGEthresh = quantile(mean_resampled_fracGEthresh, probs = hi_ci)) %>%
    mutate(CI_width = CIwidth)
  
  boot_thresh_negative <- boot_thresh %>%
    inner_join(cells_per_well, by = 'wellID') %>%
    filter(condition == negative_condition) %>%
    mutate(resampled_fracGEthresh = nGEthresh/Nc) %>%
    group_by(iteration) %>%
    summarise(mean_resampled_negative_fracGEthresh = mean(resampled_fracGEthresh))
 
  boot_thresh_base <- boot_thresh %>%
    inner_join(cells_per_well, by = 'wellID') %>%
    filter(condition == base_condition) %>%
    mutate(resampled_fracGEthresh = nGEthresh/Nc) %>%
    group_by(iteration) %>%
    summarise(mean_resampled_base_fracGEthresh = mean(resampled_fracGEthresh))
  
  sum_boot_comparisons_perCondition <- boot_thresh %>%
    inner_join(cells_per_well, by = 'wellID') %>%
    mutate(fracGEthresh = nGEthresh/Nc) %>%
    group_by(iteration, condition) %>%
    summarise(mean_resampled_fracGEthresh = mean(fracGEthresh),
              mean_resampled_nGEthresh = mean(nGEthresh)) %>%
    inner_join(boot_thresh_negative, by = 'iteration') %>%
    inner_join(boot_thresh_base, by = 'iteration') %>%
    mutate(fc_vs_negative = ifelse(mean_resampled_negative_fracGEthresh > 0, mean_resampled_fracGEthresh / mean_resampled_negative_fracGEthresh, mean_resampled_nGEthresh),
           fc_vs_base = ifelse(mean_resampled_base_fracGEthresh > 0, mean_resampled_fracGEthresh / mean_resampled_base_fracGEthresh, mean_resampled_nGEthresh)) %>%
    group_by(condition) %>%
    summarise(fracGEthresh_vs_negative_fractionGnegative = sum(mean_resampled_fracGEthresh > mean_resampled_negative_fracGEthresh)/n_iter,
              fracGEthresh_vs_base_fractionGbase = sum(mean_resampled_fracGEthresh > mean_resampled_base_fracGEthresh)/n_iter,
              mean_foldChange_vs_negative = mean(fc_vs_negative),
              ci_lo_foldChange_vs_negative = quantile(fc_vs_negative, probs = lo_ci),
              ci_hi_foldChange_vs_negative = quantile(fc_vs_negative, probs = hi_ci),
              mean_foldChange_vs_base = mean(fc_vs_base),
              ci_lo_foldChange_vs_base = quantile(fc_vs_base, probs = lo_ci),
              ci_hi_foldChange_vs_base = quantile(fc_vs_base, probs = hi_ci)) %>%
    mutate(p_val_vs_negative = 1 - fracGEthresh_vs_negative_fractionGnegative,
           p_val_vs_base = 1 - fracGEthresh_vs_base_fractionGbase,
           n_iterations = n_iter,
           seed_value = seed_val,
           count_threshold = threshold)
  
  cat('Done.\n')
  
  results[['sum_boot_basicstats_perWell']] <- sum_boot_basicstats_perWell
  results[['sum_boot_basicstats_perCondition']] <- sum_boot_basicstats_perCondition
  results[['sum_boot_comparisons_perCondition']] <- sum_boot_comparisons_perCondition
  
  return(results)
}
