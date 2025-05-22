
#保存至gaussian_prop.chain。将每个rank的sigma_gauss保存到指定路径的文件中。
save_sigma_gauss_prop <- function(ptr_sigma_prop, i_rank){
  fname <- file.path(results_dir, paste0("gaussian_prop.chain", i_rank))
  con <- file(fname, open = "a")
  
  cat(sprintf("%.12e", ptr_sigma_prop[1:N_parm]), file = con, sep = " ")
  cat("\n", file = con)
  
  close(con)
  
  return(0)
}

#将当前stack的接收率信息保存到accept_rate_stacks.chain中。
save_ar_stack <- function(i_next_stack, n_iter_a_stack, accept_rate_a_stack, n_ranks){
  fname <- file.path(results_dir, paste0("accept_rate_stacks.chain"))
  
  con <- file(fname, open = "a")
  
  cat("", i_next_stack, "", file = con, sep = " ")
  cat("", n_iter_a_stack, "", file = con, sep = " ")
  cat(sprintf("%.12e", accept_rate_a_stack[1:n_ranks]), file = con, sep = " ")
  cat("\n", file = con)
  
  close(con)
  
  return(0)
}
