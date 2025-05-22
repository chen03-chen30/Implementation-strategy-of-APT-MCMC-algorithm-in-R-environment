library(foreach)
library(doParallel)
source("readin.R")
source("data_loader.R")
source("parallel_init.R")
source("rand_func.R")
source("user_logll.R")
source("user_prior.R")
source("parallel_flow.R")
source("parallel_stack.R")
source("parallel_batch.R")
source("parallel_tune.R")

#Print debug
debug <- 1

#Initialize the MPI environment
n_ranks <- 4
cl <- makeCluster(32)
registerDoParallel(cl)

#Runtime monitering start point
runtime_start <- Sys.time()

#Read the input parameters (every process should do that to load global parameters)
read_input_ini("input.ini")

#Check_and_make results_dir.
make_dir(results_dir)

#Set a rolling random seed
if (init_rand_seed <= 0) {
  rolling_seed <- as.integer(Sys.time())
} else {
  rolling_seed <- init_rand_seed
}

#get the value of nline_data
nline_data <- get_nlines_of_file(Data_file)

#read the data using root=0 and then broadcast to other ranks.
data_NlineNdim <- read_2d_data(Data_file, nline_data, ndim_data, Delimiter)

#read the parm range
read_parm_range("input.ini")

#cal the sigma_parm_min & max
sigma_parm_min <- numeric(N_parm)
sigma_parm_max <- numeric(N_parm)
sigma_parm <- calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max)
sigma_parm_min <- sigma_parm$sigma_parm_min
sigma_parm_max <- sigma_parm$sigma_parm_max

#Generate a initial random N_beta*N_parm parm array ：transit_BetaParm
transit_BetaParm <- matrix(0, nrow = n_ranks, ncol = N_parm)
transit_BetaParm <- foreach(i_rank = 1:n_ranks, .combine = rbind) %dopar% {
  
  #rolling seed
  rolling_seed <- rolling_seed + i_rank
  
  #save seed of rand
  if (debug) {
    save_the_seed(rolling_seed, results_dir, i_rank)
  }
  
  #gen rank's init_parm
  init_parm <- init_parm_set(rolling_seed, transit_BetaParm[i_rank, ])
  init_parm
}

#save all rank's init_parm
if (debug) {
  save_init_parm(transit_BetaParm, results_dir, n_ranks, N_parm)
}


#calc the initial logpost_all_ranks
logpost_all_ranks <- double(n_ranks)
logpost_all_ranks <- foreach(i_rank = 1:n_ranks, .combine = "c") %dopar% {
  debug_arr_pointer <- 1
  
  #copy init_parm to transit_betaparm
  init_Parm <- double(N_parm)
  init_Parm <- transit_BetaParm[i_rank, ]
  
  #验证chain_Parm_root数组和transit_BetaParm二维数组中根进程行的第一个元素的值是否相等。
  if (debug_arr_pointer) {
    cat("init_Parm, transit_BP:", init_Parm, transit_BetaParm[i_rank, ], "\n")
  }
  
  #计算对数似然、先验，后验
  logll_tempered <- logll_beta(init_Parm, nline_data, data_NlineNdim, i_rank)
  logprior <- log_prior(init_Parm)
  logpost <- logll_tempered + logprior
  
  save_first_chain(init_Parm, results_dir, i_rank, logpost, N_parm)
  
  logpost
}


#init_gaussian_prop 
sigma_gaussian_prop <- double(N_parm)
sigma_gaussian_prop <- init_gaussian_proposal(sigma_gaussian_prop, init_gp_ratio)
sigma_RanksParm <- matrix(0, nrow = n_ranks, ncol = N_parm)
for (i in 1:n_ranks) {
  save_sigma_gauss_prop(sigma_gaussian_prop, i)
  sigma_RanksParm[i, ] <- sigma_gaussian_prop
}

#accumulate i for the entire Markov chain
i_accumul <- numeric(n_ranks)

#accumulate N_accept for the entire Markov chain
i_accumul_accept <- numeric(n_ranks)

#Run a flow of stacks of batches. Where we tune sigma_prop between stacks and swap chains between batches

#The_Entire_flow
save_debug <- 1
processing_debug <- 1
n_accept_old <- numeric(n_ranks)
n_accept_now <- numeric(n_ranks)
accept_rate_a_stack <- numeric(n_ranks)

tune_ranks <- integer(n_ranks)

i_tmp_stack <- 0
i_next_stack <- 0

#Flow begin////////////////////////////////////////////////////

while (i_next_stack < N_iter) {
  
  #Cal i_next_stack
  i_next_stack <- i_tmp_stack + n_iter_a_stack
  if (i_next_stack > N_iter) {
    n_iter_a_stack <- N_iter - i_tmp_stack
    i_next_stack <- N_iter
  }
  
  n_accept_old <- i_accumul_accept
  
  for (i in 1:n_ranks) {
    sigma_prop <- sigma_RanksParm[i, ]
    save_sigma_gauss_prop(sigma_prop, i)
  }
  
  #Stack begin////////////////////////////////////////////////
  save_debug <- 1
  i_tmp <- 0
  i_next <- 0
  n_iter_a_batch <- 0
  i_swap <- 1
  j_swap <- 1
  
  while (i_next < n_iter_a_stack) {
    
    #Cal i_next
    n_iter_a_batch_adjust <- i4_unif_ab(-n_iter_a_batch_rand, n_iter_a_batch_rand)
    n_iter_a_batch <- n_iter_a_batch_base + n_iter_a_batch_adjust
    i_next <- i_tmp + n_iter_a_batch
    
    if (i_next > n_iter_a_stack) {
      n_iter_a_batch <- n_iter_a_stack - i_tmp
      i_next <- n_iter_a_stack
    }
    
    #Batch begin/////////////////////////////////////////////
    batch_results <- foreach(i_rank = 1:n_ranks, .combine = c) %dopar% {
      
      chain_IterParm <- matrix(0, n_iter_a_batch, N_parm)
      chain_IterParm[1, ] <- transit_BetaParm[i_rank, ]
      
      logpost_old <- logpost_all_ranks[i_rank]
      
      sigma_prop <- sigma_RanksParm[i_rank, ]
      
      iter_batch_results <- iter_batch_mh(chain_IterParm, sigma_prop, n_iter_a_batch, nline_data, data_NlineNdim, i_rank, i_accumul, i_accumul_accept, i_save_begin, logpost_old)
      
      iter_batch_results
    }
    
    transit_BetaParm <- do.call(rbind, batch_results[names(batch_results) == "final_parm"])
    logpost_all_ranks <- sapply(batch_results[names(batch_results) == "logpost_final"], function(x) x[[1]])
    i_accumul <- sapply(batch_results[names(batch_results) == "rank_i_accumul"], function(x) x[[1]])
    i_accumul_accept <- sapply(batch_results[names(batch_results) == "rank_i_accumul_accept"], function(x) x[[1]])
    
    #Swap begin////////////////////////////////////////////////
    i_swap <- i4_unif_ab(1, n_ranks - 1)
    if (save_debug) {
      save_debug_stack_sequence(i_accumul[i_swap], i_swap)
    }
    
    swap_results <- judge_and_swap(transit_BetaParm, nline_data, data_NlineNdim, logpost_all_ranks, i_swap, j_swap, Swapmode, n_ranks)
    transit_BetaParm <- swap_results$transit_BetaParm
    logpost_all_ranks <- swap_results$logpost_all_ranks
    
    i_tmp <- i_next
  }
  
  
  if (processing_debug) {
    cat(sprintf("%9d out of total %d iterations have completed.\n", i_next_stack, N_iter))
  }
  
  n_accept_now <- i_accumul_accept
  accept_rate_a_stack <- (n_accept_now - n_accept_old) / n_iter_a_stack
  
  if (save_debug) {
    save_ar_stack(i_next_stack, n_iter_a_stack, accept_rate_a_stack, n_ranks)
  }
  
  #Tune begin//////////////////////////////////////////////////////
  if ((i_next_stack < N_stoptune) && (i_next_stack < N_iter)) {
    
    #Find rank to tune
    debug <- 1
    for (i in 1:n_ranks) {
      if ((accept_rate_a_stack[i] < ar_ok_upper) && (accept_rate_a_stack[i] > ar_ok_lower)) {
        tune_ranks[i] <- 0
      } else {
        tune_ranks[i] <- 1
      }
    }
    
    for (rank_in_tune in 1:n_ranks) {
      if (tune_ranks[rank_in_tune] == 1) {
        if (processing_debug) {
          cat(sprintf("     tuning rank %d...\n", rank_in_tune))
        }
        
        #gen sigma all rank
        #sigma_alltune_ParmNvaried: gen sigma for every parm in every rank
        save_tuning_debug <- 1
        sigma_alltune_ParmNvaried <- matrix(0, nrow = N_parm, ncol = n_ranks)
        sigma_alltune_ParmNvaried <- gen_sigma_alltune(sigma_alltune_ParmNvaried, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, sigma_RanksParm, rank_in_tune, n_ranks, N_parm)
        
        
        logpost_single_tune <- 0
        logpost_single_tune <- logpost_all_ranks[rank_in_tune]
        chain_single_tune <- double(N_parm)
        chain_single_tune <- transit_BetaParm[rank_in_tune, ]
        ar_ParmNvaried <- matrix(0, nrow = N_parm, ncol = n_ranks)
        
        sigma_tune1parm_NvariedParm <- matrix(0, nrow = n_ranks, ncol = N_parm)
        sigma_tune1parm_NvariedParm <- matrix(rep(sigma_RanksParm[rank_in_tune, ], n_ranks), nrow = n_ranks, ncol = N_parm, byrow = TRUE)
        
        for(j_parm in 1:N_parm) {
          sigma_tune1parm_NvariedParm[ , j_parm] <- sigma_alltune_ParmNvaried[j_parm, ]
          
          #cal ar for different parm with different sigma(from sigma_tune1parm_NvariedParm)
          chain_IterParm <- matrix(0, nrow = n_iter_in_tune, ncol = N_parm)
          chain_IterParm[1, ] <- chain_single_tune
          logpost_old <- logpost_single_tune
          i_accumul_accept_tune <- 0
          
          ar_jparm_ranks <- foreach(i_rank = 1:n_ranks, .combine = c) %dopar% {
            sigma_prop <- sigma_tune1parm_NvariedParm[i_rank, ]
            ar <- iter_batch_mh_tune(chain_IterParm, sigma_prop, n_iter_in_tune, nline_data, data_NlineNdim, rank_in_tune, i_accumul_accept_tune, logpost_old, i_rank)
            ar
          }
          
          ar_ParmNvaried[j_parm, ] <- ar_jparm_ranks
        }
        
        #choose the iparm,irank's sigma(if it is used would close to ar_best) by ar_ParmNvaried
        iparm_change <- 1
        irank_change <- 1
        sigma_decision <- decide_sigma_to_change(ar_ParmNvaried, n_ranks, iparm_change, irank_change)
        iparm_change <- sigma_decision$iparm_change
        irank_change <- sigma_decision$irank_change
        
        #modify_sigma
        sigma_RanksParm <- modify_sigma_prop_rankintune(rank_in_tune, iparm_change, irank_change, sigma_RanksParm, sigma_alltune_ParmNvaried, ar_ParmNvaried, n_ranks)
        
        if (save_tuning_debug) {
          save_tuning_sigma_ar(ar_ParmNvaried, sigma_alltune_ParmNvaried, results_dir, rank_in_tune, n_ranks)
        }
      }
    }
  }
  
  i_tmp_stack <- i_next_stack
}

if (debug) {
  cat(sprintf("Total: ranks,i_accumul,i_accumul_accept: %d %d %d\n", n_ranks, i_accumul, i_accumul_accept));
}

runtime_end <- Sys.time()
time_spent <- as.numeric(runtime_end - runtime_start, units = "secs")

if (debug) {
  cat(sprintf("Time spent: %f \n", time_spent));
}

stopCluster(cl)

