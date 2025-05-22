c_program_dir <- "/home/R External Interface version"

n_processes <- 8
mpi_executable_name <- "main"
mpi_executable_path <- file.path(c_program_dir, mpi_executable_name)

command <- "mpirun"
args <- c("-np", as.character(n_processes), "-oversubscribe", mpi_executable_path)

setwd(c_program_dir)

system2(command, args = args, stdout = "", stderr = "")