

save_tuning_sigma_ar <- function(ar_ParmNvaried, sigma_alltune_ParmNvaried, path, rank_in_tune, n_ranks){
  fname <- file.path(path, paste0("sigma_ar_intune.rank", rank_in_tune, ".dat"))
  
  con <- file(fname, open = "a")
  
  cat("################################################\n", file = con)
  cat("ALLTUNE sigma variations:\n", file = con)
  
  for (j in 1:N_parm) {
   
    cat(sigma_alltune_ParmNvaried[j, 1:n_ranks], " ", file = con)
    
    cat("\n", file = con)
  }
  
  cat("ARs of them:\n", file = con)
  for (j in 1:N_parm) {
    
    cat(ar_ParmNvaried[j, 1:n_ranks], " ", file = con)
    
    cat("\n", file = con)
  }
  
  close(con)
  
  return(0)
}

#将后验参数链保存到tunex.dat.中
save_the_batch_tune <- function(chain_IterParm, n_iter_a_batch, path, i_rank, rank_in_tune, logpost){
  
  fname <- file.path(path, paste0("tune.", FoutPre, rank_in_tune, ".", i_rank, FoutSuf))
  
  con <- file(fname, open = "a")
  
  for (i in 1:n_iter_a_batch) {
    cat(sprintf("%.12e", chain_IterParm[i, 1:N_parm]), file = con)
    cat("", sprintf("%f\n", logpost[i]), "", file = con)
  }
  
  close(con)
  
  return(0)
  
}

#计算每个参数的在所有rank的接受率均值标准差，可以评估每个参数的接受率的变异情况，从而决定哪些参数需要进一步调优
calc_SD_allParm <- function(std_Parm, ar_ParmNvaried, n_ranks){
  
  mean_oneparm <- rowMeans(ar_ParmNvaried)
  
  variance_oneparm <- rowSums((ar_ParmNvaried - mean_oneparm)^2) / n_ranks
  
  std_Parm <- sqrt(variance_oneparm)
  
  return(std_Parm)
}

argmax <- function(array) {
  return(which.max(array))
}

argmin <- function(array) {
  return(which.min(array))
}

#目的是在所有参数和所有rank的接受率（ar）中，找到最接近目标接受率 ar_best 的值，并返回这个最接近的值及其对应的参数索引和rank索引
calc_closest_ar <- function(ar_ParmNvaried, iparm_closest, irank_closest) {
  
  abs_diff <- abs(ar_ParmNvaried - ar_best)
  
  closest_value <- min(abs_diff)
  closest_indices <- which(abs_diff == closest_value, arr.ind = TRUE)
  
  iparm_closest <- closest_indices[1, "row"]
  irank_closest <- closest_indices[1, "col"]
      
  
  return(list(closest_value = closest_value, iparm_closest = iparm_closest, irank_closest = irank_closest))
}

#目的是找到给定参数（由 iparm_stdmax 指定）在所有rank上的接受率（ar）中，最接近目标接受率 ar_best 的值，并返回这个最接近的差值，及其对应的rank索引。
calc_closest_ar_oneparm <- function(ar_ParmNvaried, irank_closest_oneparm, iparm_stdmax) {
  
  ar_row <- ar_ParmNvaried[iparm_stdmax, ]
  
  abs_diff <- abs(ar_row - ar_best)
  
  closest_value_oneparm <- min(abs_diff)
  
  irank_closest_oneparm <- which.min(abs_diff)

  
  return(list(closest_value_oneparm = closest_value_oneparm, irank_closest_oneparm = irank_closest_oneparm))
}

#目的是根据接受率（ar）矩阵 ar_ParmNvaried 决定需要改变的参数和rank。具体来说，这个函数通过计算标准差和找到最接近目标接受率的参数和rank来决定需要改变的参数iparm和对应的irank。
decide_sigma_to_change <- function(ar_ParmNvaried, n_ranks, iparm_change, irank_change) {
  
  std_Parm <- double(N_parm)
  #计算每个parm的ar标准差，可以评估每个参数的ar在allrank上的偏差
  std_Parm <- calc_SD_allParm(std_Parm, ar_ParmNvaried, n_ranks)
  
  irank_closest_oneparm <- 1
  #找到arstd最大的parm索引
  iparm_stdmax <- argmax(std_Parm)
  #第一种方案下，arstd最大的parm
  iparm_change1 <- iparm_stdmax
  
  #目的是先找到arstd最大的parm（由 iparm_stdmax 指定），再找到其在allrank上的ar中，最接近目标接受率 ar_best 的值，并返回这个最接近的差值，及其对应的rank索引
  results_oneparm <- calc_closest_ar_oneparm(ar_ParmNvaried, irank_closest_oneparm, iparm_stdmax)
  v_oneparm_closest <- results_oneparm$closest_value_oneparm
  irank_change1 <- results_oneparm$irank_closest_oneparm
  
  #方案二：直接在allparm和allrank的ar中，找到最接近目标接受率 ar_best 的值，并返回这个最接近的值及其对应的parm索引和rank索引（而无需先算出每个参数接受率标准差std）
  iparm_closest <- 1
  irank_closest <- 1
  results_2d <- calc_closest_ar(ar_ParmNvaried, iparm_closest, irank_closest)
  v_2d_closest <- results_2d$closest_value
  iparm_change2 <- results_2d$iparm_closest
  irank_change2 <- results_2d$irank_closest
  
  #第一个方案含义：在arstd最大的parm中，找到使得ar最接近ar_best的rank，那么可以说，这个parm这个rank的sigma可以使得ar最接近ar_best。
  #第二个方案含义：在所有parm的所有rank中，找到使得ar最接近ar_best的parm以及rank，那么可以说，这个parm这个rank的sigma可以使得ar最接近ar_best。
  
  #ar_accept_diff在input.ini中给出=0.1,如果第一个方案这个接受率差值<0.1
  if (v_oneparm_closest < ar_accept_diff) {
    #第一个方案and第二个方案的找到的接受率差值也<0.1，随机选择第一种方案或者第二种方案，来得到iparm,irank下的sigma去调整原来的矩阵
    if (v_2d_closest < ar_accept_diff) {
      rand_unif <- runif(1)
      if (rand_unif > 0.5) {
        iparm_change <- iparm_change1
        irank_change <- irank_change1
      } else {
        iparm_change <- iparm_change2
        irank_change <- irank_change2
      }
    } else {
      
      #只有第一个方案的接受率差值<0.1，选择第一个方案的iparm，irank的sigma去调整原来的矩阵
      iparm_change <- iparm_change1
      irank_change <- irank_change1
    }
  } else {
    
    #第一个方案>0.1，第二个方案<0.1，选择第二个方案的iparm，irank的sigma去调整原来的矩阵
    if (v_2d_closest < ar_accept_diff) {
      iparm_change <- iparm_change2
      irank_change <- irank_change2
    } else {
      
      #第一个方案和第二个方案>0.1
      rand_unif <- runif(1)
      if (rand_unif > 0.5) {
        iparm_change <- -1 
        irank_change <- 1 
      } else {
        iparm_change <- -1 
        irank_change <- -1 
      }
    }
  }
  
  return(list(iparm_change = iparm_change, irank_change = irank_change))
}

#为每个参数找到最接近目标接受率 ar_best 的sigma值
race_all_parm <- function(sigma_prop_rankintune, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks) {
  
  sigma_prop_rankintune <- sapply(1:N_parm, function(i){
    j_min <- which.min(abs(ar_ParmNvaried[i, ] - ar_best))
    return(sigma_alltune_ParmNvaried[i, j_min])
  })
  
  return(sigma_prop_rankintune)
}


modify_sigma_prop_rankintune <- function(rank_in_tune, iparm_change, irank_change, sigma_RanksParm, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks) {
  
  if (iparm_change >= 1) {
    sigma_RanksParm[rank_in_tune, iparm_change] <- sigma_alltune_ParmNvaried[iparm_change, irank_change]
  } else {
    sigma_prop_rankintune <- numeric(N_parm)
    if (irank_change >= 1) {
      sigma_prop_rankintune <- init_gaussian_proposal(sigma_prop_rankintune, init_gp_ratio)
    } else {
      #为每个参数找到最接近目标接受率 ar_best 的sigma值
      sigma_prop_rankintune <- race_all_parm(sigma_prop_rankintune, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks)
    }
    
    sigma_RanksParm[rank_in_tune, ] <- sigma_prop_rankintune
    
  }
  
  return(sigma_RanksParm)
}



create_logspace_array <- function(base_value, half_scale, n_points, scaled_arr) {
  
  log_min <- log(base_value / half_scale)
  log_max <- log(base_value * half_scale)
  log_values <- seq(from = log_min, to = log_max, length.out = n_points)
  
  return(exp(log_values))
}

check_bounceInside_sigma_boundary <- function(scaled_arr, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, n_ranks, i_parm){
  half_scale <- sqrt(sigma_jumpin_ratio)
  if (scaled_arr[1] < sigma_parm_min[i_parm]) {
    scaled_arr <- create_logspace_array(sigma_parm_min[i_parm] * half_scale, half_scale, n_ranks, scaled_arr)
  }
  if (scaled_arr[n_ranks] > sigma_parm_max[i_parm]) {
    scaled_arr <- create_logspace_array(sigma_parm_max[i_parm] / half_scale, half_scale, n_ranks, scaled_arr)
    
  }
  return(scaled_arr)
}

#gen sigma all rank
gen_sigma_alltune <- function(sigma_alltune_ParmNvaried, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, sigma_RanksParm, rank_in_tune, n_ranks, N_parm){
  
  for (i in 1:N_parm) {
    scaled_arr <- create_logspace_array(sigma_RanksParm[rank_in_tune, i], sigma_scale_half_ratio, n_ranks, scaled_arr)
    scaled_arr <- check_bounceInside_sigma_boundary(scaled_arr, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, n_ranks, i)
    sigma_alltune_ParmNvaried[i, ] <- scaled_arr
  }
  scaled_arr <- NULL
  return(sigma_alltune_ParmNvaried)
}



iter_batch_mh_tune <- function(chain_IterParm, sigma_prop, n_iter_in_tune, nline_data, data_NlineNdim, rank_in_tune, i_accumul_accept_tune, logpost_old, i_rank) {
  
  save_debug <- 1
  save_allch_ll <- 1
  
  one_chain_old <- chain_IterParm[1, ]
  one_chain_new <- double(N_parm)
  
  logpost <- double(n_iter_in_tune)
  logpost[1] <- logpost_old
  
  for (i in 1:n_iter_in_tune) {
    copy_old <- 1
    if (i > 1) {
      logpost_old <- logpost[i - 1]
      one_chain_old <- chain_IterParm[i - 1, ]
    } 
    
    one_chain_new <- do_gaussian_propose(one_chain_old, one_chain_new, sigma_prop)
    one_chain_new <- para_boundary(one_chain_new)
    q_factor <- 1
    
    logll_tempered_new <- logll_beta(one_chain_new, nline_data, data_NlineNdim, rank_in_tune)
    logprior_new <- log_prior(one_chain_new)
    logpost_new <- logll_tempered_new + logprior_new
    
    if ((logpost_new - logpost_old) > log(1 / q_factor)) {
      chain_IterParm[i, ] <- one_chain_new
      logpost[i] <- logpost_new
      i_accumul_accept_tune <- i_accumul_accept_tune + 1
      copy_old <- 0
    } else {
      H <- exp(logpost_new - logpost_old) * q_factor
      rand_unif <- runif(1)
      if (rand_unif < H) {
        chain_IterParm[i, ] <- one_chain_new
        logpost[i] <- logpost_new
        i_accumul_accept_tune <- i_accumul_accept_tune + 1
        copy_old <- 0
      }
    }
    
    if (copy_old && i > 1) {
      chain_IterParm[i, ] <- chain_IterParm[i - 1, ]
      logpost[i] <- logpost_old
    }
    
    if (save_allch_ll) {
      save_log_posterior(one_chain_new, logll_tempered_new, logprior_new, results_dir, i_rank)
    }
  }
  
  if (save_debug) {
    save_the_batch_tune(chain_IterParm, n_iter_in_tune, results_dir, i_rank, rank_in_tune, logpost)
  }
  
  ar <- i_accumul_accept_tune / n_iter_in_tune
  
  return(ar)
}
