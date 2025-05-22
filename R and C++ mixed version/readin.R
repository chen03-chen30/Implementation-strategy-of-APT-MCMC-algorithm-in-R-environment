#CORRECT

N_parm <- NULL
N_iter <- NULL
N_stoptune <- NULL
N_beta <- NULL
Fout_Len <- NULL
FoutPre <- NULL
FoutSuf <- NULL
results_dir <- NULL
Data_file <- NULL
ndim_data <- NULL
Delimiter <- NULL
Beta_Values <- NULL
init_gp_ratio <- NULL
init_rand_seed <- NULL
i_save_begin <- NULL
n_iter_a_stack <- NULL
n_iter_a_batch_base <- NULL
n_iter_a_batch_rand <- NULL
n_iter_in_tune <- NULL
ar_ok_lower <- NULL
ar_ok_upper <- NULL
ar_best <- NULL
ar_accept_diff <- NULL
sigma_scale_half_ratio <- NULL
sigma_scale_min <- NULL
sigma_scale_max <- NULL
sigma_jumpin_ratio <- NULL
N_swap <- NULL
Swapmode <- NULL




#返回input.ini中特定的para_name的那一行
read_onepara <- function(path, para_name){
  lines <- readLines(path)
  
  para_line <- NULL
  found <- 0
  
  for (line in lines) {
    if (grepl(paste0("^", para_name, ":\\s*"), line, perl = TRUE)) {
      para_line <- line
      found <- found + 1
    }
  }
  
  if (found > 1) {
    stop(paste("Error:duplicate definition of '", para_name, "'!", sep = ""))
  } else if (found == 0) {
    stop(paste("Error: '", para_name, "' not found!", sep = ""))
  }
  
  return(para_line)
}

#将input.ini参数的值读取到那个参数中
read_chains_grid <- function(path){
  
  para_name <- "N_parm"
  para_line <- read_onepara(path, para_name)
  N_parm <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "N_iter"
  para_line <- read_onepara(path, para_name)
  N_iter <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "N_stoptune"
  para_line <- read_onepara(path, para_name)
  N_stoptune <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "N_beta"
  para_line <- read_onepara(path, para_name)
  N_beta <<- as.integer(strsplit(para_line, ":")[[1]][2])
}

#使用read_onepara将input.ini中file_name读取到变量中
read_chain_out_name <- function(path){
  
  para_name <- "Fout_Len"
  para_line <- read_onepara(path, para_name)
  Fout_Len <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "FoutPre"
  para_line <- read_onepara(path, para_name)
  FoutPre <<- trimws(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "FoutSuf"
  para_line <- read_onepara(path, para_name)
  FoutSuf <<- trimws(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "results_dir"
  para_line <- read_onepara(path, para_name)
  results_dir <<- trimws(strsplit(para_line, ":")[[1]][2])
}

#使用read_onepara将input.ini中datafile参量读取到变量中
# 定义读取数据描述的函数
read_data_desc <- function(path) {
  
  para_name <- "Data_file"
  para_line <- read_onepara(path, para_name)
  Data_file <<- trimws(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "ndim_data"
  para_line <- read_onepara(path, para_name)
  ndim_data <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "Delimiter"
  para_line <- read_onepara(path, para_name)
  Delimiter <<- trimws(strsplit(para_line, ":")[[1]][2])
  if (Delimiter == "blank") {
    Delimiter <<- " "
  }
}

#读取input.ini中的beta_values值并存进去
read_beta_values <- function(path){
  
  para_name <- "Beta_Values"
  para_line <- read_onepara(path, para_name)
  beta_values_str <- strsplit(para_line, ":")[[1]][2]
  Beta_Values <<- as.numeric(unlist(strsplit(beta_values_str, ",")))
  
  para_name <- "N_beta"
  para_line <- read_onepara(path, para_name)
  N_beta <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  if  (length(Beta_Values) != N_beta) {
    stop("Error: number of Beta values does not match N_beta")
  }
}

#读取input.ini中采样参数存储到变量中
# 定义读取采样参数的函数
read_sampling_para <- function(path) {
  
  para_name <- "init_gp_ratio"
  para_line <- read_onepara(path, para_name)
  init_gp_ratio <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "init_rand_seed"
  para_line <- read_onepara(path, para_name)
  init_rand_seed <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "i_save_begin"
  para_line <- read_onepara(path, para_name)
  i_save_begin <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "n_iter_a_stack"
  para_line <- read_onepara(path, para_name)
  n_iter_a_stack <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "n_iter_a_batch_base"
  para_line <- read_onepara(path, para_name)
  n_iter_a_batch_base <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "n_iter_a_batch_rand"
  para_line <- read_onepara(path, para_name)
  n_iter_a_batch_rand <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "n_iter_in_tune"
  para_line <- read_onepara(path, para_name)
  n_iter_in_tune <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "ar_ok_lower"
  para_line <- read_onepara(path, para_name)
  ar_ok_lower <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "ar_ok_upper"
  para_line <- read_onepara(path, para_name)
  ar_ok_upper <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "ar_best"
  para_line <- read_onepara(path, para_name)
  ar_best <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "ar_accept_diff"
  para_line <- read_onepara(path, para_name)
  ar_accept_diff <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "sigma_scale_half_ratio"
  para_line <- read_onepara(path, para_name)
  sigma_scale_half_ratio <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "sigma_scale_min"
  para_line <- read_onepara(path, para_name)
  sigma_scale_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "sigma_scale_max"
  para_line <- read_onepara(path, para_name)
  sigma_scale_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "sigma_jumpin_ratio"
  para_line <- read_onepara(path, para_name)
  sigma_jumpin_ratio <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "N_swap"
  para_line <- read_onepara(path, para_name)
  N_swap <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "Swapmode"
  para_line <- read_onepara(path, para_name)
  Swapmode <<- as.integer(strsplit(para_line, ":")[[1]][2])
  
}

make_dir <- function(path) {
  # 尝试创建目录（默认不递归创建父目录，权限受系统限制）
  if (!dir.create(path, showWarnings = FALSE)) {
    # 当目录已存在或创建失败时执行
    message("\nWARNING: dir: ", path, " exists!")
    message("-------  may overwrite previous results!")
    message("-------  MSG from R function.")
    message()
    return(1)  # 返回1表示目录已存在或创建失败
  }
  return(0)  # 返回0表示创建成功
}

#【集成以上】读取input.ini中所有的参数,所有参数都是全局变量
read_input_ini <- function(path){
  read_chains_grid(path)
  read_chain_out_name(path)
  read_beta_values(path)
  read_data_desc(path)
  read_sampling_para(path)
}


