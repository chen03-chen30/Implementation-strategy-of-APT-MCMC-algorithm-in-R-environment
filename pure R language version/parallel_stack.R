
#将do_swap保存到swap_decision.dat文件中
save_debug_stack_doswap <- function(so_swap){
  fname <- file.path(results_dir, "swap_decision.dat")
  
  con <- file(fname, open = "a")
  
  cat(do_swap, "\n", file = con)
  
  close(con)
  
  return(0)
}

#将ptr_i_accumul和i_swap保存到swap_sequence.dat文件中
save_debug_stack_sequence <- function(i_accumul, i_swap){
  fname <- file.path(results_dir, "swap_sequence.dat")
  
  con <- file(fname, open = "a")
  
  cat(i_accumul, " ", i_swap, "\n", file = con)
  
  close(con)
  
  return(0)
}

#Swap two parm chains
# swap_two_chains <- function(ptr_ichain, ptr_jchain, N_parm){
#   if (length(ptr_ichain) != length(ptr_jchain) || length(ptr_ichain) != N_parm) {
#     stop("The lengths of the two chains must be the same and equal to N_parm.")
#     
#     temp <- ptr_ichain
#     ptr_ichain <- ptr_jchain
#     ptr_jchain <- temp
#     
#     return(list(ptr_ichain = ptr_ichain, ptr_jchain = ptr_jchain))
#   }
# }

#// min() Find minimum between two numbers.
#// max() Find maximum between two numbers.

#判断两个链是否可以交换，并进行交换
judge_and_swap <- function(transit_BetaParm, nline_data, data_NlineNdim, logpost_all_ranks, i_swap, j_swap, Swapmode, n_ranks){
  
  debug_save_swap <- 0
  do_swap <- 0
  H <- 0
  rand_unif <- 0
  
  #input.ini中给出swapmode=0,即相邻链交换
  if (Swapmode == 0) {
    j_swap <- i_swap + 1
  } 
  
  #随机交换
  if (Swapmode == 1) {
    j_swap <- i4_unif_ab(1, n_ranks - 1)
    if (j_swap >= i_swap) {
      j_swap <- j_swap + 1
    }
  }
    
  #i链使用ibeta的logpost和j链使用jbeta的logpost
  logpost_ichain_ibeta <- logpost_all_ranks[i_swap]
  logpost_jchain_jbeta <- logpost_all_ranks[j_swap]
    
    #复制i_swap链的参数到chain_Parm_ichain
  chain_Parm_ichain <- double(N_parm)
  chain_Parm_ichain <- transit_BetaParm[i_swap, ]
    
    #计算似然函数，但用的是j_swap的beta值了
  logll_tempered_mix_ij <- logll_beta(chain_Parm_ichain, nline_data, data_NlineNdim, j_swap)
  logprior_i <- log_prior(chain_Parm_ichain)
    #计算出i链使用jbeta的对数后验概率
  logpost_ichain_jbeta <- logll_tempered_mix_ij + logprior_i
    
  chain_Parm_jchain <- double(N_parm)
  chain_Parm_jchain <- transit_BetaParm[j_swap, ]
  logll_tempered_mix_ji <- logll_beta(chain_Parm_jchain, nline_data, data_NlineNdim, i_swap)
  logprior_j <- log_prior(chain_Parm_jchain)
  logpost_jchain_ibeta <- logll_tempered_mix_ji + logprior_j
    
    #据metropolis-Hastings准则判断是否进行交换
  if ((logpost_ichain_jbeta + logpost_jchain_ibeta - logpost_ichain_ibeta - logpost_jchain_jbeta) > 0) {
    do_swap <- 1
  } else {
    H <- exp(logpost_ichain_jbeta + logpost_jchain_ibeta - logpost_ichain_ibeta - logpost_jchain_jbeta)
    rand_unif <- runif(1)
    if (rand_unif < H) {
      do_swap <- 1
    }
  }
    
    #如果决定交换
  if (do_swap) {
    temp <- transit_BetaParm[i_swap, ]
    transit_BetaParm[i_swap, ] <- transit_BetaParm[j_swap, ]
    transit_BetaParm[j_swap, ] <- temp
    logpost_all_ranks[i_swap] <- logpost_jchain_ibeta
    logpost_all_ranks[j_swap] <- logpost_ichain_jbeta
  }
    
    #如果debug_save_swap=1,则记录交换信息do_swap到swap_decision.dat文件中
  if (debug_save_swap) {
    save_debug_stack_doswap(do_swap)
  }
    
  return(list(transit_BetaParm = transit_BetaParm, logpost_all_ranks = logpost_all_ranks))
}
