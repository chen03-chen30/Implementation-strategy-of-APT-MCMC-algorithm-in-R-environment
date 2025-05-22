#CORRECT

period_to_au <- function(period_days, ms_Msun){
  PI <- as.double(3.14159265358)
  GG <- as.double(6.67259e-8)
  Msun <- as.double(1.9891e33)
  AU2cm <- as.double(1.4959787e13)
  
  period_second <- as.double(period_days) * 24.0 * 3600.0
  
  a <- as.double(period_second/(2.0 * PI))
  b <- as.double(GG * Msun * ms_Msun)
  r <- as.double((a ^ 2.0) * b)
  a_AU <- (r ^ (1/3))/AU2cm
  
  return(a_AU)
}
calc_osi <- function(mp_Mearth, ms_Msun, a_AU, d_pc){
  osi_alpha <- 3.0 * (mp_Mearth) * (ms_Msun^-1.0) * (a_AU) * (d_pc^-1.0)
  return(osi_alpha)
}

newton_solver <- function(EE, e){
  tolerance <- 1e-8
  eps <- 1
  M <- as.double(EE)
  
  while (abs(eps) > tolerance) {
    E1 <- EE - (EE - e * sin(EE) - M) / (1 - e * cos(EE))
    eps <- as.double(E1 - EE)
    EE <- as.double(E1)
  }
  
  return(EE)
}

func_as <- function(time_con, Nline_time, ecc, osi, cosi, OmegaO, M0, omega, per, Ntime_radec){
  deg2rad <- 0.0174532925
  PI <- 3.14159265358
  
  omega <- omega * deg2rad
  M0 <- M0 * deg2rad
  OmegaO <- OmegaO * deg2rad
  
  coso <- cos(omega)
  sino <- sin(omega)
  cosOg <- cos(OmegaO)
  sinOg <- sin(OmegaO)
  
  A <- osi * (coso * cosOg - sino * sinOg * cosi)
  B <- osi * (coso * sinOg + sino * cosOg * cosi)
  F <- osi * (-sino * cosOg - coso * sinOg * cosi)
  G <- osi * (-sino * sinOg + coso * cosOg * cosi)
  
  for (i in 1:Nline_time) {
    EE <- newton_solver(((2 * PI) / per * time_con[i] - M0), ecc)
    X <- cos(EE) - ecc
    Y <- sqrt(1 - ecc^2) * sin(EE)
    Ntime_radec[i, 1] <- B * X + G * Y
    Ntime_radec[i, 2] <- A * X + F * Y
  }
  
  return(Ntime_radec)
}

log_likelihood <- function(Nline_time, Ntime_radec, data_radec, var_uk){
  PI2 <- 6.28318530716
  log_llhd <- 0
  sig_power <- var_uk * var_uk
  AC_twice <- 2.0 * log(PI2^(-0.5)) + 2.0 * log(sig_power^(-0.5))
  
  for (i in 1:Nline_time) {
    ra_once <- -((data_radec[i, 1] - Ntime_radec[i, 1])^2) / (2.0 * sig_power)
    dec_once <- -((data_radec[i, 2] - Ntime_radec[i, 2])^2) / (2.0 * sig_power)
    log_llhd <- log_llhd + ra_once + dec_once + AC_twice
  }
  
  return(log_llhd)
}

logll_beta <- function(ptr_one_chain, nline_data, data_NlineNdim, i_rank){
  ms_Msun <- 1.0
  d_pc <- 3.0
  radec_int2 <- 2
  ndim_data <- 2
  
  cos_inc1 <- ptr_one_chain[1]
  ecc1 <- ptr_one_chain[2]
  an_Omg1 <- ptr_one_chain[3]
  p_omg1 <- ptr_one_chain[4]
  M0_1 <- ptr_one_chain[5]
  pl_m1 <- ptr_one_chain[6]
  var_uk <- ptr_one_chain[7]
  period1 <- ptr_one_chain[8]
  cos_inc2 <- ptr_one_chain[9]
  ecc2 <- ptr_one_chain[10]
  an_Omg2 <- ptr_one_chain[11]
  p_omg2 <- ptr_one_chain[12]
  M0_2 <- ptr_one_chain[13]
  pl_m2 <- ptr_one_chain[14]
  period2 <- ptr_one_chain[15]
  
  logll <- 0
  Nline_time <- nline_data / 2
  time_con <- double(Nline_time)
  data_radec <- matrix(0, nrow = Nline_time, ncol = radec_int2)
  
  for (i in 1:Nline_time) {
    time_con[i] <- data_NlineNdim[(i - 1) * ndim_data + 1, 1]
    data_radec[i, 1] <- data_NlineNdim[(i - 1) * ndim_data + 1, 2]
    data_radec[i, 2] <- data_NlineNdim[(i - 1) * ndim_data + 2, 2]
  }
  
  a_AU1 <- period_to_au(period1, ms_Msun)
  osi1 <- calc_osi(pl_m1, ms_Msun, a_AU1, d_pc)
  
  a_AU2 <- period_to_au(period2, ms_Msun)
  osi2 <- calc_osi(pl_m2, ms_Msun, a_AU2, d_pc)
  
  Ntime_radec1 <- matrix(0, nrow = Nline_time, ncol = radec_int2)
  Ntime_radec2 <- matrix(0, nrow = Nline_time, ncol = radec_int2)
  
  Ntime_radec1 <- func_as(time_con, Nline_time, ecc1, osi1, cos_inc1, an_Omg1, M0_1, p_omg1, period1, Ntime_radec1)
  Ntime_radec2 <- func_as(time_con, Nline_time, ecc2, osi2, cos_inc2, an_Omg2, M0_2, p_omg2, period2, Ntime_radec2)
  
  Ntime_radec <- matrix(0, nrow = Nline_time, ncol = radec_int2)
  
  for (i in 1:Nline_time) {
    Ntime_radec[i, 1] <- Ntime_radec1[i, 1] + Ntime_radec2[i, 1]
    Ntime_radec[i, 2] <- Ntime_radec1[i, 2] + Ntime_radec2[i, 2]
  }
  
  logll <- log_likelihood(Nline_time, Ntime_radec, data_radec, var_uk)
  logll <- logll * Beta_Values[i_rank]
  
  Ntime_radec <- NULL
  Ntime_radec1 <- NULL
  Ntime_radec2 <- NULL
  data_radec <- NULL
  time_con <- NULL
  
  return(logll)
}