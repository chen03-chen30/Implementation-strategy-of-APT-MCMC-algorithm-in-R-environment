#include <Rcpp.h>
#include <Rmath.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

using namespace Rcpp;


/////////////////init/////////////////////
// 保存至chainx.dat。将每个rank的first_chain保存到指定路径的文件中
/// 
// [[Rcpp::export]]
int save_first_chain_cpp(NumericVector chain_parm, std::string path, int i_rank, double logpost_first, int N_parm, std::string FoutPre, std::string FoutSuf) {
  std::string fname = path + "/" + FoutPre + std::to_string(i_rank) + FoutSuf;
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < N_parm; ++i) {
      con << std::fixed << std::setprecision(12) << chain_parm[i] << (i == N_parm - 1 ? "" : " ");
    }
    con << " " << std::fixed << std::setprecision(12) << logpost_first << " 0 0" << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

// 将每个rank的init.parm保存到init.parm
/// 
// [[Rcpp::export]]
int save_init_parm_cpp(NumericMatrix transit_BetaParm, std::string path, int n_ranks, int N_parm) {
  std::string fname = path + "/init.parm";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < n_ranks; ++i) {
      con << "init para set" << (i + 1) << ": ";
      for (int j = 0; j < N_parm; ++j) {
        con << std::fixed << std::setprecision(6) << transit_BetaParm(i, j) << (j == N_parm - 1 ? "" : " ");
      }
      con << std::endl;
    }
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

// 将每个rank的seed保存到init.randseed
// [[Rcpp::export]]
int save_the_seed_cpp(int seed, std::string path, int i_rank) {
  std::string fname = path + "/init.randseed";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    con << "seed for rank " << i_rank << " is " << seed << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

/////////////////init/////////////////////
/////////////////flow/////////////////////
// 保存至gaussian_prop.chain。将每个rank的sigma_gauss保存到指定路径的文件中。
/// 
// [[Rcpp::export]]
int save_sigma_gauss_prop_cpp(NumericVector ptr_sigma_prop, int i_rank, std::string results_dir, int N_parm) {
  std::string fname = results_dir + "/gaussian_prop.chain" + std::to_string(i_rank);
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < N_parm; ++i) {
      con << std::fixed << std::setprecision(12) << ptr_sigma_prop[i] << (i == N_parm - 1 ? "" : " ");
    }
    con << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

// 将当前stack的接收率信息保存到accept_rate_stacks.chain中。
/// 
// [[Rcpp::export]]
int save_ar_stack_cpp(int i_next_stack, int n_iter_a_stack, NumericVector accept_rate_a_stack, int n_ranks, std::string results_dir) {
  std::string fname = results_dir + "/accept_rate_stacks.chain";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    con << " " << i_next_stack << " " << n_iter_a_stack << " ";
    for (int i = 0; i < n_ranks; ++i) {
      con << std::fixed << std::setprecision(12) << accept_rate_a_stack[i] << (i == n_ranks - 1 ? "" : " ");
    }
    con << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

/////////////////flow/////////////////////
/////////////////batch/////////////////////
/// 
// [[Rcpp::export]]
int save_the_batch_cpp(const std::vector<std::vector<double>>& chain_IterParm_cpp, int n_iter_a_batch, std::string path, int i_rank, const std::vector<double>& logpost, const std::vector<int>& accumul, const std::vector<int>& accumul_accept, int N_parm, std::string FoutPre, std::string FoutSuf) {

  std::string fname = path + "/" + FoutPre + std::to_string(i_rank) + FoutSuf;
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < n_iter_a_batch; ++i) {
      for (int	 j = 0; j < N_parm; ++j) {
        con << std::fixed << std::setprecision(12) << chain_IterParm_cpp[i][j] << (j == N_parm - 1 ? "" : " ");
      }
      con << " " << logpost[i] << " " << accumul[i] << " " << accumul_accept[i] << std::endl;
    }
    con.close();
    return 0;
  } else {
    std::cerr << "ERROR: Unable to open file: " << fname << std::endl;
    return 1; // Indicate an error
  }
}

// 将参数链，新的对数先验，新的后验保存到chainx.dat.all.ll中
/// 
// [[Rcpp::export]]
int save_log_posterior_cpp(const std::vector<double>& one_chain_new_cpp, double logll_temper_new, double logprior_new, std::string path, int i_rank, int N_parm, std::string FoutPre, std::string FoutSuf) {
  std::string fname = path + "/" + FoutPre + std::to_string(i_rank) + FoutSuf + ".all.ll";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < N_parm; ++i) {
      con << std::fixed << std::setprecision(12) << one_chain_new_cpp[i] << (i == N_parm - 1 ? "" : " ");
    }
    con << " " << logll_temper_new << " " << logprior_new << " " << (logprior_new + logll_temper_new) << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

// 将当前参数链和对应的标准差保存到debug_gaussian_prop中
/// 
// [[Rcpp::export]]
int save_debug_gaussian_proposal_cpp(NumericVector one_chain, NumericVector sigma_prop, std::string results_dir, int N_parm) {
  std::string fname = results_dir + "/debug_gaussian_prop";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < N_parm; ++i) {
      con << std::fixed << std::setprecision(12) << one_chain[i] << (i == N_parm - 1 ? "" : " ");
    }
    for (int i = 0; i < N_parm; ++i) {
      con << " " << sigma_prop[i] << (i == N_parm - 1 ? "" : " ");
    }
    con << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}
/////////////////batch/////////////////////
/////////////////stack/////////////////////
// 将do_swap保存到swap_decision.dat文件中
/// 
// [[Rcpp::export]]
int save_debug_stack_doswap_cpp(bool so_swap, std::string results_dir) {
  std::string fname = results_dir + "/swap_decision.dat";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    con << so_swap << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

// 将ptr_i_accumul和i_swap保存到swap_sequence.dat文件中
/// 
// [[Rcpp::export]]
int save_debug_stack_sequence_cpp(int i_accumul, int i_swap, std::string results_dir) {
  std::string fname = results_dir + "/swap_sequence.dat";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    con << i_accumul << " " << i_swap << std::endl;
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

/////////////////stack/////////////////////
/////////////////tune/////////////////////
/// 
// [[Rcpp::export]]
int save_tuning_sigma_ar_cpp(NumericMatrix ar_ParmNvaried, NumericMatrix sigma_alltune_ParmNvaried, std::string path, int rank_in_tune, int n_ranks, int N_parm) {
  std::string fname = path + "/sigma_ar_intune.rank" + std::to_string(rank_in_tune) + ".dat";
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    con << "################################################" << std::endl;
    con << "ALLTUNE sigma variations:" << std::endl;

    for (int j = 0; j < N_parm; ++j) {
      for (int k = 0; k < n_ranks; ++k) {
        con << sigma_alltune_ParmNvaried(j, k) << (k == n_ranks - 1 ? "" : " ");
      }
      con << std::endl;
    }

    con << "ARs of them:" << std::endl;
    for (int j = 0; j < N_parm; ++j) {
      for (int k = 0; k < n_ranks; ++k) {
        con << ar_ParmNvaried(j, k) << (k == n_ranks - 1 ? "" : " ");
      }
      con << std::endl;
    }

    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

/// 
// [[Rcpp::export]]
int save_the_batch_tune_cpp(const std::vector<std::vector<double>>& chain_IterParm_cpp, int n_iter_a_batch, std::string path, int i_rank, int rank_in_tune, NumericVector logpost, int N_parm, std::string FoutPre, std::string FoutSuf) {
  std::string fname = path + "/tune." + FoutPre + std::to_string(rank_in_tune) + "." + std::to_string(i_rank) + FoutSuf;
  std::ofstream con(fname, std::ios::app);

  if (con.is_open()) {
    for (int i = 0; i < n_iter_a_batch; ++i) {
      for (int j = 0; j < N_parm; ++j) {
        con << std::fixed << std::setprecision(12) << chain_IterParm_cpp[i][j] << (j == N_parm - 1 ? "" : " ");
      }
      con << " " << std::fixed << std::setprecision(6) << logpost[i] << std::endl;
    }
    con.close();
    return 0;
  } else {
    Rcpp::warning("Unable to open file: %s", fname.c_str());
    return 1; // Indicate an error
  }
}

/////////////////tune/////////////////////
