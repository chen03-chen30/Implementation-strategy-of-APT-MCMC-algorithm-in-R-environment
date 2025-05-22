
save_the_batch <- function(chain_IterParm, n_iter_a_batch, path, i_rank, logpost, accumul, accumul_accept){
  
  fname <- file.path(path, paste0(FoutPre, i_rank, FoutSuf))
  
  con <- file(fname, open = "a")
  
  for (i in 1:n_iter_a_batch) {
    cat(sprintf("%.12e", chain_IterParm[i,1:N_parm]), file = con, sep = " ")
    
    cat("", logpost[i], "", file = con, sep = " ")
    cat("", accumul[i], "", file = con, sep = " ")
    cat("", accumul_accept[i], "", file = con, sep = " ")

    cat("\n", file = con)
  }
  
  close(con)
  
  return(0)
}

#将参数链，新的对数先验，新的后验保存到chainx.dat.all.ll中
save_log_posterior <- function(one_chain_new, logll_temper_new, logprior_new, path, i_rank){
  
  fname <- file.path(path, paste0(FoutPre, i_rank, FoutSuf, ".all.ll"))
  
  con <- file(fname, open = "a")
  
  
  cat(sprintf("%.12e", one_chain_new[1:N_parm]), file = con, sep = " ")
  
  cat("", logll_temper_new, "", file = con, sep = " ")
  cat("", logprior_new, "", file = con, sep = " ")
  cat("", logprior_new + logll_temper_new, "", file = con, sep = "\n")
  
  close(con)
  
  return(0)
}

#将当前参数链和对应的标准差保存到debug_gaussian_prop中
save_debug_gaussian_proposal <- function(one_chain, sigma_prop, results_dir, N_parm){
  
  fname <- file.path(results_dir, "debug_gaussian_prop")
  
  con <- file(fname, open = "a")
  
  cat(sprintf("%.12e", one_chain[1:N_parm]), file = con, sep = " ")
  cat(sigma_prop[1:N_parm], file = con, sep = " ")

  cat("\n", file = con)
  
  close(con)
  
  return(0)
}

#实现了一个高斯提议分布，用于生成新的参数链
do_gaussian_propose <- function(one_chain_old, one_chain_new, sigma_prop){
  #Gaussian proposal distribution.
  #This proposal is nearly as simple as it gets, mathematically it is:
  #  q(x∗∣xi)=Normal(xi,σ2),
  #that is, a Gaussian centered on the current position xi with variance given by a standard deviation parameter σ.
  #  Input:
  #  ptr_one_chain: poiter to the part of Parameter array in 2d chain_iterparm.
  #ptr_one_chain_new: poiter to a Parameter array of parms.
  #ptr_sigma_prop: poiter to the array of the Standard deviation of Gaussian distribution.
  #  Main Work: propose parameter set of in ptr_one_chain_new.
  #  Return: 
  #  q_factor: ratio of proposal densities.
  
  for(i in 1:N_parm){
    x <- one_chain_old[i]
    sigma_prop_x <- sigma_prop[i]
    dx <- r8_normal_01()*sigma_prop_x
    one_chain_new[i] <- x + dx
  }
  
  debug <- 0
  if (debug) {
    save_debug_gaussian_proposal(one_chain_new, sigma_prop, results_dir = "./output", N_parm)
  }
  
  return(one_chain_new)
}

#实现一个batch的MH算法，就是将parm的后验参数链进行一个batch的迭代，最后返回一个batch最后的logpost
iter_batch_mh <- function(chain_IterParm, sigma_prop, n_iter_a_batch, nline_data, data_NlineNdim, i_rank, i_accumul, i_accumul_accept, i_save_begin, logpost_old){
  accumul_accept <- integer(n_iter_a_batch)
  accumul <- integer(n_iter_a_batch)
  logpost <- double(n_iter_a_batch)
  
  logpost_final <- 0
  logll_tempered_new <- 0
  logprior_new <- 0
  logpost_new <- 0
  q_factor <- 0
  H <- 0
  rand_unif <- 0
  
  #当前迭代的参数链输入
  one_chain_old <- chain_IterParm[1, ]
  
  #debug int used to turn on int save_log_posterior
  save_allch_ll <- 0
  
  #set the initial first logpost 
  logpost[1] <- logpost_old
  
  #开始迭代循环
  for (i in 1:n_iter_a_batch) {
    copy_old <- 1
    
    #如果不是第一次迭代，更新logpost_old和one_chain_old
    if (i > 1) {
      logpost_old <- logpost[i - 1]
      one_chain_old <- chain_IterParm[i - 1, ]
    }
    
    #新的参数链
    one_chain_new <- double(N_parm)
    
    #实现了一个高斯提议分布，通过旧参数链，来生成新的参数链，并且给出高斯提议比率q=1
    one_chain_new <- do_gaussian_propose(one_chain_old, one_chain_new, sigma_prop)
    q_factor <- 1
    
    #检查参数是否在边界内
    one_chain_new <- para_boundary(one_chain_new)
    
    #用new_chain计算新的对数似然，对数先验，对数后验
    logll_tempered_new <- logll_beta(one_chain_new, nline_data, data_NlineNdim, i_rank)
    logprior_new <- log_prior(one_chain_new)
    logpost_new <- logll_tempered_new + logprior_new
    
    #MH-algorithm ///////////MH算法
    #NOTE: not compute Hastings ratio in the 1st_criterion to avoid overflow of exp function
    #save_chain(one_chain_new, logpost_new, sigma_prop , results_dir, i_rank)
    
    #如果新的对数后验概率>旧的对数后验概率，则直接接受新的参数链
    if ((logpost_new - logpost_old) > log(1 / q_factor)) {
      
      #1st setting, total 3//更新当前设置的参数链为new_chain的参数，更新当前的logpost为new_logpost，接受迭代次数+1
      chain_IterParm[i, ] <- one_chain_new
      logpost[i] <- logpost_new
      i_accumul_accept[i_rank] <- i_accumul_accept[i_rank] + 1
      copy_old <- 0
    }else{
      
      ##如果新的对数后验概率<旧的对数后验概率，计算MH比率H判断
      H <- exp(logpost_new - logpost_old) * q_factor
      rand_unif <- runif(1)
      
      #如果生成的随机数小于H，则接受新的参数链
      if (rand_unif < H) {
        
        ##2nd是 setting, total 3
        #更新当前设置的参数链为new_chain的参数，更新当前的logpost为new_logpost，接受迭代次数+1
        chain_IterParm[i, ] <- one_chain_new
        logpost[i] <- logpost_new
        i_accumul_accept[i_rank] <- i_accumul_accept[i_rank] + 1
        copy_old <- 0
      }
    }
    
    #如果没有接受新的参数链，则复制旧的参数链和对数后验概率
    if (copy_old && i > 1 ) {
      chain_IterParm[i, ] <- chain_IterParm[i - 1, ]
      logpost[i] <- logpost_old
    }
    
    if (save_allch_ll) {
      save_log_posterior(one_chain_new, logll_tempered_new, logprior_new, results_dir, i_rank)
    }
    
    #更新ptr_i_accumul，ptr_i_accumul_accept
    i_accumul[i_rank] <- i_accumul[i_rank] + 1
    accumul[i] <- i_accumul[i_rank]
    accumul_accept[i] <- i_accumul_accept[i_rank]
  }
  
  #如果累计迭代次数>=500000，则保存当前batch的结果
  if (i_accumul[i_rank] >= i_save_begin) {
    save_the_batch(chain_IterParm, n_iter_a_batch, results_dir, i_rank, logpost, accumul, accumul_accept)
  }
  #在chain.dat文件中，会发现有的行parm、logpost变了，但是i_accumul_accept没有增加的情况，这是因为这一行只是与别的chain将parm进行swap了。
  
  #设置最终的对数后验概率为最后一个迭代的对数后验概率
  logpost_final <- logpost[n_iter_a_batch]
  
  #返回一个batch的最终的对数后验概率，累计接受迭代次数
  return(list(i_rank = i_rank, 
              final_parm = chain_IterParm[n_iter_a_batch, ],
              logpost_final = logpost_final, 
              rank_i_accumul = i_accumul[i_rank],
              rank_i_accumul_accept = i_accumul_accept[i_rank]))
}
