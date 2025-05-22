# 使用 C++ 编译器
CXX = mpicxx

# C++ 编译选项
CXXFLAGS = -Wall -Werror 
CXXFLAGS += -Wextra
CXXFLAGS += -O2 	#-O3
CXXFLAGS += -lm 

# 目标文件
mtmc: main.cpp mpi_init.cpp mpi_batch.cpp mpi_stack.cpp mpi_flow.cpp mpi_tune.cpp data_loader.cpp alloc.cpp rand_func.cpp readin.cpp user_logll.cpp user_prior.cpp 
    $(CXX) main.cpp mpi_init.cpp mpi_batch.cpp mpi_stack.cpp mpi_flow.cpp mpi_tune.cpp data_loader.cpp alloc.cpp rand_func.cpp readin.cpp user_logll.cpp user_prior.cpp $(CXXFLAGS) -o a.out

# 清理规则
clean:
    rm -rf *.o
    rm a.out

