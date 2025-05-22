#CORRECT

#Part 0
#General part, no need to change in most cases
#source("readin.R")
#source("rand_func.R")

#Part 1
#define how many parameters in the model
#their values are read from input files by read_parm_range(char *path)
para0_min <- NULL
para1_min <- NULL
para2_min <- NULL
para3_min <- NULL
para4_min <- NULL
para5_min <- NULL
para6_min <- NULL
para7_min <- NULL
para8_min <- NULL
para9_min <- NULL
para10_min <- NULL
para11_min <- NULL
para12_min <- NULL
para13_min <- NULL
para14_min <- NULL

para0_max <- NULL
para1_max <- NULL
para2_max <- NULL
para3_max <- NULL
para4_max <- NULL
para5_max <- NULL
para6_max <- NULL
para7_max <- NULL
para8_max <- NULL
para9_max <- NULL
para10_max <- NULL
para11_max <- NULL
para12_max <- NULL
para13_max <- NULL
para14_max <- NULL

#Part 2: read prior range
#NOTE: Input Min and Max in default.
#NOTE: Add your own control para (mean, std, etc) in case of needed.
read_parm_range <- function(path) {
  ################### para0
  para_name <- "para0_max"
  para_line <- read_onepara(path, para_name)
  para0_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para0_min"
  para_line <- read_onepara(path, para_name)
  para0_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para1
  para_name <- "para1_max"
  para_line <- read_onepara(path, para_name)
  para1_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para1_min"
  para_line <- read_onepara(path, para_name)
  para1_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para2
  para_name <- "para2_max"
  para_line <- read_onepara(path, para_name)
  para2_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para2_min"
  para_line <- read_onepara(path, para_name)
  para2_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para3
  para_name <- "para3_max"
  para_line <- read_onepara(path, para_name)
  para3_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para3_min"
  para_line <- read_onepara(path, para_name)
  para3_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para4
  para_name <- "para4_max"
  para_line <- read_onepara(path, para_name)
  para4_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para4_min"
  para_line <- read_onepara(path, para_name)
  para4_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para5
  para_name <- "para5_max"
  para_line <- read_onepara(path, para_name)
  para5_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para5_min"
  para_line <- read_onepara(path, para_name)
  para5_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para6
  para_name <- "para6_max"
  para_line <- read_onepara(path, para_name)
  para6_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para6_min"
  para_line <- read_onepara(path, para_name)
  para6_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para7
  para_name <- "para7_max"
  para_line <- read_onepara(path, para_name)
  para7_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para7_min"
  para_line <- read_onepara(path, para_name)
  para7_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para8
  para_name <- "para8_max"
  para_line <- read_onepara(path, para_name)
  para8_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para8_min"
  para_line <- read_onepara(path, para_name)
  para8_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para9
  para_name <- "para9_max"
  para_line <- read_onepara(path, para_name)
  para9_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para9_min"
  para_line <- read_onepara(path, para_name)
  para9_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para10
  para_name <- "para10_max"
  para_line <- read_onepara(path, para_name)
  para10_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para10_min"
  para_line <- read_onepara(path, para_name)
  para10_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para11
  para_name <- "para11_max"
  para_line <- read_onepara(path, para_name)
  para11_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para11_min"
  para_line <- read_onepara(path, para_name)
  para11_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para12
  para_name <- "para12_max"
  para_line <- read_onepara(path, para_name)
  para12_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para12_min"
  para_line <- read_onepara(path, para_name)
  para12_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para13
  para_name <- "para13_max"
  para_line <- read_onepara(path, para_name)
  para13_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para13_min"
  para_line <- read_onepara(path, para_name)
  para13_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  ################### para14
  para_name <- "para14_max"
  para_line <- read_onepara(path, para_name)
  para14_max <<- as.numeric(strsplit(para_line, ":")[[1]][2])
  
  para_name <- "para14_min"
  para_line <- read_onepara(path, para_name)
  para14_min <<- as.numeric(strsplit(para_line, ":")[[1]][2])
}

#Part 3: initialize prior, part 1 
#NOTE: uniform initialization in default.
#NOTE: Specify your own init function in case of needed.

para0_init <- function(para0_min, para0_max) {
  para0 <- r8_unif_ab(para0_min, para0_max)
  return(para0)
}

para1_init <- function(para1_min, para1_max) {
  para1 <- r8_unif_ab(para1_min, para1_max)
  return(para1)
}

para2_init <- function(para2_min, para2_max) {
  para2 <- r8_unif_ab(para2_min, para2_max)
  return(para2)
}

para3_init <- function(para3_min, para3_max) {
  para3 <- r8_unif_ab(para3_min, para3_max)
  return(para3)
}

para4_init <- function(para4_min, para4_max) {
  para4 <- r8_unif_ab(para4_min, para4_max)
  return(para4)
}

para5_init <- function(para5_min, para5_max) {
  para5 <- r8_unif_ab(para5_min, para5_max)
  return(para5)
}

para6_init <- function(para6_min, para6_max) {
  para6 <- r8_unif_ab(para6_min, para6_max)
  return(para6)
}

para7_init <- function(para7_min, para7_max) {
  para7 <- r8_unif_ab(para7_min, para7_max)
  return(para7)
}

para8_init <- function(para8_min, para8_max) {
  para8 <- r8_unif_ab(para8_min, para8_max)
  return(para8)
}

para9_init <- function(para9_min, para9_max) {
  para9 <- r8_unif_ab(para9_min, para9_max)
  return(para9)
}

para10_init <- function(para10_min, para10_max) {
  para10 <- r8_unif_ab(para10_min, para10_max)
  return(para10)
}

para11_init <- function(para11_min, para11_max) {
  para11 <- r8_unif_ab(para11_min, para11_max)
  return(para11)
}

para12_init <- function(para12_min, para12_max) {
  para12 <- r8_unif_ab(para12_min, para12_max)
  return(para12)
}

para13_init <- function(para13_min, para13_max) {
  para13 <- r8_unif_ab(para13_min, para13_max)
  return(para13)
}

para14_init <- function(para14_min, para14_max) {
  para14 <- r8_unif_ab(para14_min, para14_max)
  return(para14)
}

#Part 4 initialize, part 2
#NOTE: uniform initialization in default.
#NOTE: Specify your own init function in case of needed.

#Part 5: calc the scaling boundary of each para

calc_sigma_scale_boundary <- function(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max){
  
  sigma_parm_min[1] <- sigma_scale_min * (para0_max - para0_min)
  sigma_parm_min[2] <- sigma_scale_min * (para1_max - para1_min)
  sigma_parm_min[3] <- sigma_scale_min * (para2_max - para2_min)
  sigma_parm_min[4] <- sigma_scale_min * (para3_max - para3_min)
  sigma_parm_min[5] <- sigma_scale_min * (para4_max - para4_min)
  sigma_parm_min[6] <- sigma_scale_min * (para5_max - para5_min)
  sigma_parm_min[7] <- sigma_scale_min * (para6_max - para6_min)
  sigma_parm_min[8] <- sigma_scale_min * (para7_max - para7_min)
  sigma_parm_min[9] <- sigma_scale_min * (para8_max - para8_min)
  sigma_parm_min[10] <- sigma_scale_min * (para9_max - para9_min)
  sigma_parm_min[11] <- sigma_scale_min * (para10_max - para10_min)
  sigma_parm_min[12] <- sigma_scale_min * (para11_max - para11_min)
  sigma_parm_min[13] <- sigma_scale_min * (para12_max - para12_min)
  sigma_parm_min[14] <- sigma_scale_min * (para13_max - para13_min)
  sigma_parm_min[15] <- sigma_scale_min * (para14_max - para14_min)
  
  sigma_parm_max[1] <- sigma_scale_max * (para0_max - para0_min)
  sigma_parm_max[2] <- sigma_scale_max * (para1_max - para1_min)
  sigma_parm_max[3] <- sigma_scale_max * (para2_max - para2_min)
  sigma_parm_max[4] <- sigma_scale_max * (para3_max - para3_min)
  sigma_parm_max[5] <- sigma_scale_max * (para4_max - para4_min)
  sigma_parm_max[6] <- sigma_scale_max * (para5_max - para5_min)
  sigma_parm_max[7] <- sigma_scale_max * (para6_max - para6_min)
  sigma_parm_max[8] <- sigma_scale_max * (para7_max - para7_min)
  sigma_parm_max[9] <- sigma_scale_max * (para8_max - para8_min)
  sigma_parm_max[10] <- sigma_scale_max * (para9_max - para9_min)
  sigma_parm_max[11] <- sigma_scale_max * (para10_max - para10_min)
  sigma_parm_max[12] <- sigma_scale_max * (para11_max - para11_min)
  sigma_parm_max[13] <- sigma_scale_max * (para12_max - para12_min)
  sigma_parm_max[14] <- sigma_scale_max * (para13_max - para13_min)
  sigma_parm_max[15] <- sigma_scale_max * (para14_max - para14_min)
  
  sigma_bound_ok <- 0
  for (i in 1:N_parm){
    if (sigma_parm_min[i] < 0) {
      sigma_bound_ok <- sigma_bound_ok + 1
    }
    if (sigma_parm_max[i] < 0) {
      sigma_bound_ok <- sigma_bound_ok + 1
    }
    if ((sigma_parm_max[i] - sigma_parm_min[i]) < 0) {
      sigma_bound_ok <- sigma_bound_ok + 1
    }
  }
  
  if (sigma_bound_ok > 0) {
    cat("ERR: bad sigma_parm_min/max", sigma_bound_ok, "TIMES!\n")
    return(1)
  }else{
    return(list(sigma_parm_min = sigma_parm_min, sigma_parm_max = sigma_parm_max))
  }
}

#set init parameter at N_ITER = 0
init_parm_set <- function(seed, chain_parm){
  set.seed(seed)
  
  # input_file <- "input.ini"
  # read_parm_range(input_file)
  
  chain_parm[1] <- para0_init(para0_min, para0_max)
  chain_parm[2] <- para1_init(para1_min, para1_max)
  chain_parm[3] <- para2_init(para2_min, para2_max)
  chain_parm[4] <- para3_init(para3_min, para3_max)
  chain_parm[5] <- para4_init(para4_min, para4_max)
  chain_parm[6] <- para5_init(para5_min, para5_max)
  chain_parm[7] <- para6_init(para6_min, para6_max)
  chain_parm[8] <- para7_init(para7_min, para7_max)
  chain_parm[9] <- para8_init(para8_min, para8_max)
  chain_parm[10] <- para9_init(para9_min, para9_max)
  chain_parm[11] <- para10_init(para10_min, para10_max)
  chain_parm[12] <- para11_init(para11_min, para11_max)
  chain_parm[13] <- para12_init(para12_min, para12_max)
  chain_parm[14] <- para13_init(para13_min, para13_max)
  chain_parm[15] <- para14_init(para14_min, para14_max)
  
  # sigma_parm_min <- numeric(N_parm)
  # sigma_parm_max <- numeric(N_parm)
  # result <- calc_sigma_scale_boundary(sigma_scale_min, sigma_scale_max, sigma_parm_min, sigma_parm_max)
  # sigma_parm_min <- result$sigma_parm_min
  # sigma_parm_max <- result$sigma_parm_max
  return(chain_parm)
}







#Part 6: init gaussian proposal
#set init gaussian proposal for the sampling
init_gaussian_proposal <- function(ptr_sigma_prop, init_gp_ratio){
  ptr_sigma_prop[1] <- (para0_max - para0_min) * init_gp_ratio
  ptr_sigma_prop[2] <- (para1_max - para1_min) * init_gp_ratio
  ptr_sigma_prop[3] <- (para2_max - para2_min) * init_gp_ratio
  ptr_sigma_prop[4] <- (para3_max - para3_min) * init_gp_ratio
  ptr_sigma_prop[5] <- (para4_max - para4_min) * init_gp_ratio
  ptr_sigma_prop[6] <- (para5_max - para5_min) * init_gp_ratio
  ptr_sigma_prop[7] <- (para6_max - para6_min) * init_gp_ratio
  ptr_sigma_prop[8] <- (para7_max - para7_min) * init_gp_ratio
  ptr_sigma_prop[9] <- (para8_max - para8_min) * init_gp_ratio
  ptr_sigma_prop[10] <- (para9_max - para9_min) * init_gp_ratio
  ptr_sigma_prop[11] <- (para10_max - para10_min) * init_gp_ratio
  ptr_sigma_prop[12] <- (para11_max - para11_min) * init_gp_ratio
  ptr_sigma_prop[13] <- (para12_max - para12_min) * init_gp_ratio
  ptr_sigma_prop[14] <- (para13_max - para13_min) * init_gp_ratio
  ptr_sigma_prop[15] <- (para14_max - para14_min) * init_gp_ratio
  
  # 返回 ptr_sigma_prop 向量
  return(ptr_sigma_prop)
}

#Part 6: check if a proposed point is within the range
#TODO: add a new function for this ugly long boring function.
save_debug_para_boundary <- function(p_s, p, p_min, p_max){
  fname <- file.path(result_dir, paste0(p_s, ".debug_para_boundary"))
  con <- file(fname, open = "a")
  
  cat(sprintf("%.12e", p), file = con, sep = " ")
  cat(sprintf("%.12e", p_min), file = con, sep = " ")
  cat(sprintf("%.12e", p_max), file = con, sep = "\n")
  
  close(con)
  
  return(0)
}

bounce_inside <- function(para, para_min, para_max){
  if (para > para_max) {
    para <- para_max - (para - para_max)
    para <- max(para, para_min)
  }
  if (para < para_min) {
    para <- para_min + (para_min - para)
    para <- min(para, para_max)
  }
  return(para)
}

para_boundary <- function(ptr_one_chain_new){
  para0 <- ptr_one_chain_new[1]
  para1 <- ptr_one_chain_new[2]
  para2 <- ptr_one_chain_new[3]
  para3 <- ptr_one_chain_new[4]
  para4 <- ptr_one_chain_new[5]
  para5 <- ptr_one_chain_new[6]
  para6 <- ptr_one_chain_new[7]
  para7 <- ptr_one_chain_new[8]
  para8 <- ptr_one_chain_new[9]
  para9 <- ptr_one_chain_new[10]
  para10 <- ptr_one_chain_new[11]
  para11 <- ptr_one_chain_new[12]
  para12 <- ptr_one_chain_new[13]
  para13 <- ptr_one_chain_new[14]
  para14 <- ptr_one_chain_new[15]
  
  debug <- 0
  
  ptr_one_chain_new[1] <- bounce_inside(para0, para0_min, para0_max)
  ptr_one_chain_new[2] <- bounce_inside(para1, para1_min, para1_max)
  ptr_one_chain_new[3] <- bounce_inside(para2, para2_min, para2_max)
  ptr_one_chain_new[4] <- bounce_inside(para3, para3_min, para3_max)
  ptr_one_chain_new[5] <- bounce_inside(para4, para4_min, para4_max)
  ptr_one_chain_new[6] <- bounce_inside(para5, para5_min, para5_max)
  ptr_one_chain_new[7] <- bounce_inside(para6, para6_min, para6_max)
  ptr_one_chain_new[8] <- bounce_inside(para7, para7_min, para7_max)
  ptr_one_chain_new[9] <- bounce_inside(para8, para8_min, para8_max)
  ptr_one_chain_new[10] <- bounce_inside(para9, para9_min, para9_max)
  ptr_one_chain_new[11] <- bounce_inside(para10, para10_min, para10_max)
  ptr_one_chain_new[12] <- bounce_inside(para11, para11_min, para11_max)
  ptr_one_chain_new[13] <- bounce_inside(para12, para12_min, para12_max)
  ptr_one_chain_new[14] <- bounce_inside(para13, para13_min, para13_max)
  ptr_one_chain_new[15] <- bounce_inside(para14, para14_min, para14_max)
  
  if (debug) {
    save_debug_para_boundary("para0", ptr_one_chain_new[1], para0_min, para0_max)
    save_debug_para_boundary("para1", ptr_one_chain_new[2], para1_min, para1_max)
    save_debug_para_boundary("para2", ptr_one_chain_new[3], para2_min, para2_max)
    save_debug_para_boundary("para3", ptr_one_chain_new[4], para3_min, para3_max)
    save_debug_para_boundary("para4", ptr_one_chain_new[5], para4_min, para4_max)
    save_debug_para_boundary("para5", ptr_one_chain_new[6], para5_min, para5_max)
    save_debug_para_boundary("para6", ptr_one_chain_new[7], para6_min, para6_max)
    save_debug_para_boundary("para7", ptr_one_chain_new[8], para7_min, para7_max)
    save_debug_para_boundary("para8", ptr_one_chain_new[9], para8_min, para8_max)
    save_debug_para_boundary("para9", ptr_one_chain_new[10], para9_min, para9_max)
    save_debug_para_boundary("para10", ptr_one_chain_new[11], para10_min, para10_max)
    save_debug_para_boundary("para11", ptr_one_chain_new[12], para11_min, para11_max)
    save_debug_para_boundary("para12", ptr_one_chain_new[13], para12_min, para12_max)
    save_debug_para_boundary("para13", ptr_one_chain_new[14], para13_min, para13_max)
    save_debug_para_boundary("para14", ptr_one_chain_new[15], para14_min, para14_max)
  }
  return(ptr_one_chain_new)
}

#Part 7: log prior
#NOTE: Write your own prior function of all parameters
#NOTE: para_boundary ensures: ((a<=a_max) && (a>=a_min))
#in two-planet orbit retrieval:
#para0:  cos_i   p1
#para1:  ecc     p1
#para2:  anO     p1
#para3:  po      p1
#para4:  M0      p1
#para5:  mp      p1
#para6:  var_uke 
#para7:  period  p1
#para8:  cos_i   p2
#para9:  ecc     p2
#para10: anO     p2
#para11: po      p2
#para12: M0      p2
#para13: mp      p2
#para14: period  p2
prior_para0 <- function(para0_min, para0_max) {
  return(1 / (para0_max - para0_min))
}

prior_para1 <- function(para1_min, para1_max) {
  return(1 / (para1_max - para1_min))
}

prior_para2 <- function(para2_min, para2_max) {
  return(1 / (para2_max - para2_min))
}

prior_para3 <- function(para3_min, para3_max) {
  return(1 / (para3_max - para3_min))
}

prior_para4 <- function(para4_min, para4_max) {
  return(1 / (para4_max - para4_min))
}

prior_para5 <- function(para5, para5_min, para5_max) {
  return(1 / (para5 * log(para5_max / para5_min)))
}

prior_para6 <- function(para6, para6_min, para6_max) {
  return(1 / (para6 * log(para6_max / para6_min)))
}

prior_para7 <- function(para7, para7_min, para7_max) {
  return(1 / (para7 * log(para7_max / para7_min)))
}

prior_para8 <- function(para8_min, para8_max) {
  return(1 / (para8_max - para8_min))
}

prior_para9 <- function(para9_min, para9_max) {
  return(1 / (para9_max - para9_min))
}

prior_para10 <- function(para10_min, para10_max) {
  return(1 / (para10_max - para10_min))
}

prior_para11 <- function(para11_min, para11_max) {
  return(1 / (para11_max - para11_min))
}

prior_para12 <- function(para12_min, para12_max) {
  return(1 / (para12_max - para12_min))
}

prior_para13 <- function(para13, para13_min, para13_max) {
  return(1 / (para13 * log(para13_max / para13_min)))
}

prior_para14 <- function(para14, para14_min, para14_max) {
  return(1 / (para14 * log(para14_max / para14_min)))
}

#Part 8: Combine all prior distributions 
log_prior <- function(ptr_one_chain){
  para0 <- ptr_one_chain[1]
  para1 <- ptr_one_chain[2]
  para2 <- ptr_one_chain[3]
  para3 <- ptr_one_chain[4]
  para4 <- ptr_one_chain[5]
  para5 <- ptr_one_chain[6]
  para6 <- ptr_one_chain[7]
  para7 <- ptr_one_chain[8]
  para8 <- ptr_one_chain[9]
  para9 <- ptr_one_chain[10]
  para10 <- ptr_one_chain[11]
  para11 <- ptr_one_chain[12]
  para12 <- ptr_one_chain[13]
  para13 <- ptr_one_chain[14]
  para14 <- ptr_one_chain[15]
  
  log_prior <- log(prior_para0(para0_min, para0_max)) +
    log(prior_para1(para1_min, para1_max)) +
    log(prior_para2(para2_min, para2_max)) +
    log(prior_para3(para3_min, para3_max)) +
    log(prior_para4(para4_min, para4_max)) +
    log(prior_para5(para5, para5_min, para5_max)) +
    log(prior_para6(para6, para6_min, para6_max)) +
    log(prior_para7(para7, para7_min, para7_max)) +
    log(prior_para8(para8_min, para8_max)) +
    log(prior_para9(para9_min, para9_max)) +
    log(prior_para10(para10_min, para10_max)) +
    log(prior_para11(para11_min, para11_max)) +
    log(prior_para12(para12_min, para12_max)) +
    log(prior_para13(para13, para13_min, para13_max)) +
    log(prior_para14(para14, para14_min, para14_max))
  
  debug <- 0
  
  if (debug) {
    cat("para0,para1,para2,para3,para4,para5,para6,para7,para8,para9,para10,para11,para12,para13,para14: ", 
        para0, para1, para2, para3, para4, para5, para6, para7, para8, para9, para10, para11, para12, para13, para14, "\n")
    cat("their prior and sum: ", 
        log(prior_para0(para0_min, para0_max)), 
        log(prior_para1(para1_min, para1_max)), 
        log(prior_para2(para2_min, para2_max)), 
        log(prior_para3(para3_min, para3_max)), 
        log(prior_para4(para4_min, para4_max)), 
        log(prior_para5(para5, para5_min, para5_max)), 
        log(prior_para6(para6, para6_min, para6_max)), 
        log(prior_para7(para7, para7_min, para7_max)), 
        log(prior_para8(para8_min, para8_max)), 
        log(prior_para9(para9_min, para9_max)), 
        log(prior_para10(para10_min, para10_max)), 
        log(prior_para11(para11_min, para11_max)), 
        log(prior_para12(para12_min, para12_max)), 
        log(prior_para13(para13, para13_min, para13_max)), 
        log(prior_para14(para14, para14_min, para14_max)), 
        log_prior, "\n")
  }
  
  return(log_prior)
}
