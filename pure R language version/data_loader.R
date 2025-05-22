#CORRECT

#library(Rmpi)

# 获取datafile的行数
get_nlines_of_file <- function(path) {
  con <- file(path, "r") 
  on.exit(close(con)) 
  line_number <- length(readLines(con)) 
  return(line_number) 
}

# 读取datafile数据，并存储在data_NlineNdim中
read_2d_data <- function(path, nline_data, ndim_data, delimiter){
  con <- file(path,"r")
  on.exit(close(con))
  lines <- readLines(con)
  data_NlineNdim <- matrix(0, nrow = nline_data, ncol = ndim_data)
  
  for (iline_local in 1:nline_data) {
    tokens <- strsplit(lines[iline_local], delimiter)[[1]]
    data_NlineNdim[iline_local, ] <- as.numeric(tokens)
  }
  
  return(data_NlineNdim)
}
