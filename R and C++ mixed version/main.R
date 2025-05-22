#部分使用了Cpp，因为需要将整个Cycle转换成Cpp以提高效率
library(Rcpp)
source("readin.R")
source("data_loader.R")
source("user_prior.R")
source("user_logll.R")
library(cppfiles)

#//////解决了cpp文件之间，cpp与r文件的链接问题：将cpp文件全部打包成一个R包

#//////由于R-foreach-cpp遇到了困难，R-Rcppparallel-cpp也遇到了困难：
#//////主要是R与C在并行时对于数组的传递与识别有区别而导致的。
#//////故尝试R-OpenMP-cpp：1.OpenMP并行无法启动，实际上还是串行，2.拟合结果很不好。
#使用了OPENMP后总是崩溃，找到原因是：核心问题在于 OpenMP 并行循环内部对 Rcpp 对象（特别是 Rcpp::NumericVector 的创建和赋值，以及可能通过 Rcpp::checkUserInterrupt 对 R 内部状态的访问）的并发操作存在不稳定性或未预期的副作用（如 stack imbalance 或内存损坏），这些副作用在迭代次数增多时最终导致崩溃。【这可能也是其他并行也无法运行的根本原因：并行&Rcpp对象操作不稳定】
#拟解决方法：1.避免了 Rcpp 对象交互。2.减少了 R API 调用。3.最小化了栈压力。
#范围缩小：iter_mh_batch和iter_mh_tune部分会导致崩溃
#最终解决方法：并行部分必须全部使用Cpp语句，否则会报错。

#Print debug
debug <- 1

#Runtime monitering start point
runtime_start <- Sys.time()

##///////////////////////////////
##////////////READ///////////////
##///////////////////////////////

#Initialize the MPI environment
n_ranks <- 8

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

para_min <- c(para0_min, para1_min, para2_min, para3_min, para4_min,
              para5_min, para6_min, para7_min, para8_min, para9_min,
              para10_min, para11_min, para12_min, para13_min, para14_min)
para_max <- c(para0_max, para1_max, para2_max, para3_max, para4_max,
              para5_max, para6_max, para7_max, para8_max, para9_max,
              para10_max, para11_max, para12_max, para13_max, para14_max)

#cal the sigma_parm_min & max
sigma_parm_min <- numeric(N_parm)
sigma_parm_max <- numeric(N_parm)
sigma_parm <- calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max)
sigma_parm_min <- sigma_parm$sigma_parm_min
sigma_parm_max <- sigma_parm$sigma_parm_max


##///////////////////////////////////////
##////////////Init_setting///////////////
##///////////////////////////////////////


#Generate a initial random N_beta*N_parm parm array ：transit_BetaParm
transit_BetaParm <- matrix(0, nrow = n_ranks, ncol = N_parm)
for (i_rank in 0:(n_ranks - 1)) {
  # rolling seed
  local_rolling_seed <- rolling_seed + i_rank

  # save seed of rand
  #if (debug) {
    #save_the_seed_cpp(local_rolling_seed, results_dir, i_rank)
  #}

  # gen rank's init_parm
  init_parm_rank <- init_parm_set(local_rolling_seed, transit_BetaParm[i_rank + 1, ])
  transit_BetaParm[i_rank + 1, ] <- init_parm_rank # 将结果存储在预先分配的矩阵中
}

#save all rank's init_parm
#if (debug) {
#  save_init_parm_cpp(transit_BetaParm, results_dir, n_ranks, N_parm)
#}


#calc the initial logpost_all_ranks
logpost_all_ranks <- double(n_ranks)
for (i_rank in 0:(n_ranks - 1)) {
  debug_arr_pointer <- 0

  #copy init_parm to transit_betaparm
  init_Parm <- double(N_parm)
  init_Parm <- transit_BetaParm[i_rank + 1, ]

  #验证chain_Parm_root数组和transit_BetaParm二维数组中根进程行的第一个元素的值是否相等。
  if (debug_arr_pointer) {
    cat("init_Parm, transit_BP:", init_Parm, transit_BetaParm[i_rank + 1, ], "\n")
  }

  #计算对数似然、先验，后验
  logll_tempered <- logll_beta(init_Parm, nline_data, data_NlineNdim, i_rank)
  logprior <- log_prior(init_Parm)
  logpost <- logll_tempered + logprior

  #save_first_chain_cpp(init_Parm, results_dir, i_rank, logpost, N_parm, FoutPre, FoutSuf)

  logpost_all_ranks[i_rank] <- logpost # 将结果存储在预先分配的向量中
}


#init_gaussian_prop
sigma_gaussian_prop <- double(N_parm)
sigma_gaussian_prop <- init_gaussian_proposal_cpp(para_min, para_max, init_gp_ratio)
sigma_RanksParm <- matrix(0, nrow = n_ranks, ncol = N_parm)
for (i in 0:(n_ranks - 1)) {
  #save_sigma_gauss_prop_cpp(sigma_gaussian_prop, i, results_dir, N_parm)
  sigma_RanksParm[i + 1, ] <- sigma_gaussian_prop
}

##///////////////////////////////////////
##////////////Flow///////////////////////
##///////////////////////////////////////

#accumulate i for the entire Markov chain
i_accumul <- rep(0L, n_ranks) # 使用 rep(0L, n_ranks) 初始化为整数零向量

#accumulate N_accept for the entire Markov chain
i_accumul_accept <- rep(0L, n_ranks) # 使用 rep(0L, n_ranks) 初始化为整数零向量


#The_Entire_flow
cat("Starting MCMC using run_mcmc_cpp...\n")

my_seed <- 520
set.seed(my_seed)
cpp_seed <- as.integer(my_seed)

mcmc_results <- run_mcmc_cpp(
  N_iter = N_iter,                   # 总迭代次数 (来自 read_input_ini)
  n_ranks = n_ranks,                 # rank 数量 (来自 read_input_ini 或脚本开始处设置)
  N_parm = N_parm,                   # 参数数量 (来自 read_input_ini)
  nline_data = nline_data,           # 数据行数 (来自 data_loader.R)
  data_NlineNdim = data_NlineNdim,   # 数据矩阵 (来自 data_loader.R)
  sigma_RanksParm = sigma_RanksParm, # 初始 sigma 提议矩阵 (来自上面的初始化)
  transit_BetaParm = transit_BetaParm,# 初始参数矩阵 (来自上面的初始化)
  logpost_all_ranks = logpost_all_ranks, # 初始对数后验概率向量 (来自上面的初始化)
  i_accumul = i_accumul,             # 初始累积迭代计数 (初始化为 0) - C++ 函数需要 IntegerVector
  i_accumul_accept = i_accumul_accept, # 初始累积接受计数 (初始化为 0) - C++ 函数需要 IntegerVector
  n_iter_a_stack = n_iter_a_stack,   # 每个堆栈的基础迭代次数 (来自 read_input_ini)
  n_iter_a_batch_base = n_iter_a_batch_base, # 每个批次的基础迭代次数 (来自 read_input_ini)
  n_iter_a_batch_rand = n_iter_a_batch_rand, # 批次迭代次数随机范围 (来自 read_input_ini)
  Swapmode = Swapmode,               # 交换模式 (来自 read_input_ini)
  i_save_begin = i_save_begin,       # 开始保存批次的迭代数 (来自 read_input_ini)
  N_stoptune = N_stoptune,           # 停止调优的迭代数 (来自 read_input_ini)
  ar_ok_upper = ar_ok_upper,         # 接受率上界 (来自 read_input_ini)
  ar_ok_lower = ar_ok_lower,         # 接受率下界 (来自 read_input_ini)
  sigma_parm_min = sigma_parm_min,   # Sigma 最小值向量 (来自 calc_sigma_scale_boundary)
  sigma_parm_max = sigma_parm_max,   # Sigma 最大值向量 (来自 calc_sigma_scale_boundary)
  sigma_jumpin_ratio = sigma_jumpin_ratio, # Sigma 跳跃比例 (来自 read_input_ini)
  sigma_scale_half_ratio = sigma_scale_half_ratio, # Sigma 缩放半比率 (来自 read_input_ini)
  n_iter_in_tune = n_iter_in_tune,   # 调优内部迭代次数 (来自 read_input_ini)
  ar_best = ar_best,                 # 最佳接受率 (来自 read_input_ini)
  ar_accept_diff = ar_accept_diff,   # 可接受的接受率差异 (来自 read_input_ini)
  init_gp_ratio = init_gp_ratio,     # 初始高斯提议比例 (来自 read_input_ini)
  para_min = para_min,               # 参数下界向量 (来自 read_parm_range)
  para_max = para_max,               # 参数上界向量 (来自 read_parm_range)
  Beta_Values = Beta_Values,         # Beta 温度值向量 (来自 read_input_ini)
  results_dir = results_dir          # 结果目录路径 (来自 read_input_ini)
  
)

final_sigma_RanksParm <- mcmc_results$sigma_RanksParm
final_transit_BetaParm <- mcmc_results$transit_BetaParm
final_logpost_all_ranks <- mcmc_results$logpost_all_ranks
final_i_accumul <- mcmc_results$i_accumul
final_i_accumul_accept <- mcmc_results$i_accumul_accept

cat("\n--- Final MCMC State ---\n")
cat("Final total iterations per rank:\n")
print(final_i_accumul)
cat("\nFinal total accepted iterations per rank:\n")
print(final_i_accumul_accept)
cat("\nFinal overall acceptance rate per rank:\n")
print(final_i_accumul_accept / final_i_accumul) # 避免除零错误 (如果 final_i_accumul 可能为 0)
cat("\nFinal log posterior per rank:\n")
print(final_logpost_all_ranks)
cat("\nFinal parameters (first few rows/cols):\n")
print(head(final_transit_BetaParm))
cat("\nFinal sigma parameters (first few rows/cols):\n")
print(head(final_sigma_RanksParm))

if (debug) {
  cat(sprintf("Total: ranks, i_accumul, i_accumul_accept: %d %s %s\n",
              n_ranks,
              paste(final_i_accumul, collapse=","), # 打印向量
              paste(final_i_accumul_accept, collapse=",") # 打印向量
  ));
}

runtime_end <- Sys.time()
time_spent <- as.numeric(runtime_end - runtime_start, units = "secs")

if (debug) {
  cat(sprintf("Time spent: %f \n", time_spent));
}


