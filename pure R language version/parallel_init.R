
#保存至chainx.dat。将每个rank的first_chain保存到指定路径的文件中
save_first_chain <- function(chain_parm, path, i_rank, logpost_first, N_parm){
  fname <- file.path(path, paste0(FoutPre, i_rank, FoutSuf))
  
  con <- file(fname, open = "a")
  
  cat(sprintf("%.12e", chain_parm[1:N_parm]), file = con, sep = " ")
  cat("", sprintf("%.12e", logpost_first), "", file = con, sep = " ")
  cat(0, "", file = con, sep = " ")
  cat(0, file = con)
  
  cat("\n", file = con)
  
  close(con)
  
  return(0)
}

#将每个rank的init.parm保存到init.parm
save_init_parm <- function(transit_BetaParm, path, n_ranks, N_parm){
  fname <- file.path(path, "init.parm")
  
  con <- file(fname, open = "a")
  
  for (i in 1:n_ranks) {
    cat(paste0("init para set", i, ": "), file = con)
    
    cat(sprintf("%f", transit_BetaParm[i, 1:N_parm]), file = con, sep = " ")
    
    cat("\n", file = con)
  }
  close(con)
  
  return(0)
}

#将每个rank的seed保存到init.randseed
save_the_seed <- function(seed, path, i_rank){
  fname <- file.path(path, "init.randseed")
  
  con <- file(fname, open = "a")
  
  cat(sprintf("seed for rank %d is %d\n", i_rank, seed), file = con)
  
  close(con)
  
  return(0)
}
