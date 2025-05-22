#include <Rcpp.h>
#include <omp.h>      // 包含 OpenMP 头文件
#include <vector>
#include <cmath>      
#include <numeric>    
#include <limits>     
#include <algorithm>  
#include <string>
#include <stdexcept>  
#include <iomanip>    
#include <fstream>    



// From save.cpp (假设签名如你所述)
int save_sigma_gauss_prop_cpp(Rcpp::NumericVector ptr_sigma_prop, int i_rank, std::string results_dir, int N_parm);
int save_debug_stack_sequence_cpp(int i_accumul, int i_swap, std::string results_dir);
int save_ar_stack_cpp(int i_next_stack, int n_iter_a_stack, Rcpp::NumericVector accept_rate_a_stack, int n_ranks, std::string results_dir);
int save_debug_stack_doswap_cpp(bool do_swap, std::string results_dir); // 假设 do_swap 是 bool
int save_log_posterior_cpp(const std::vector<double>& one_chain_new_cpp, double logll_temper_new, double logprior_new, std::string path, int i_rank, int N_parm, std::string FoutPre = "chain", std::string FoutSuf = ".dat");
int save_the_batch_cpp(const std::vector<std::vector<double>>& chain_IterParm_cpp, int n_iter_a_batch, std::string path, int i_rank, const std::vector<double>& logpost, const std::vector<int>& accumul, const std::vector<int>& accumul_accept, int N_parm, std::string FoutPre = "chain", std::string FoutSuf = ".dat");
int save_the_batch_tune_cpp(const std::vector<std::vector<double>>& chain_IterParm_cpp, int n_iter_a_batch, std::string path, int i_rank, int rank_in_tune, Rcpp::NumericVector logpost, int N_parm, std::string FoutPre = "tune_chain", std::string FoutSuf = ".dat");
int save_tuning_sigma_ar_cpp(Rcpp::NumericMatrix ar_ParmNvaried, Rcpp::NumericMatrix sigma_alltune_ParmNvaried, std::string results_dir, int rank_in_tune, int n_ranks, int N_parm);
int save_debug_gaussian_proposal_cpp(const std::vector<double>& one_chain_cpp, const std::vector<double>& sigma_prop_cpp, std::string results_dir, int N_parm);
  
// From rand_func.cpp (假设是线程安全的)
inline int i4_unif_ab(int a, int b) {
    if (a > b) throw std::range_error("i4_unif_ab: a > b");
    if (a == b) return a;
    double range = static_cast<double>(b) - a + 1.0;
    int offset = static_cast<int>(std::floor(R::runif(0.0, 1.0) * range));
    int result = a + offset;
    return std::min(result, b); // Clamp for safety
}
inline double r8_normal_01() { return R::norm_rand(); }
inline double r8_unif_ab(double a, double b) { return R::runif(a, b); }

// From other .cpp files (根据 R 代码函数名和你的描述)
void do_gaussian_propose_cpp(const std::vector<double>& one_chain_old_cpp,
                             std::vector<double>& one_chain_new_cpp, // Modified in place
                             const std::vector<double>& sigma_prop_cpp,
                             int N_parm, std::string results_dir);
                             
std::vector<double> para_boundary_cpp(std::vector<double> one_chain_new_cpp, // Pass by value to modify copy
                                     const std::vector<double>& para_min_cpp, // Keep Rcpp for bounds - read only
                                     const std::vector<double>& para_max_cpp
                                     ); // Added N_parm
                                     
double logll_beta_cpp(const std::vector<double>& ptr_one_chain_cpp, // Use C++ vector
                     int nline_data,
                     const std::vector<std::vector<double>>& data_cpp, // Keep Rcpp - read only
                     int i_rank,
                     const std::vector<double>& beta_values_cpp);    // Keep Rcpp - read only

double log_prior_cpp(const std::vector<double>& ptr_one_chain_cpp, // Use C++ vector
                    const std::vector<double>& para_min_cpp,         // Keep Rcpp - read only
                    const std::vector<double>& para_max_cpp);         // Keep Rcpp - read only

Rcpp::NumericVector init_gaussian_proposal_cpp(const Rcpp::NumericVector& para_min,
                                              const Rcpp::NumericVector& para_max,
                                              double init_gp_ratio); // This likely stays Rcpp


// --- Helper Function Translations (R -> C++) ---

// 实现了一个高斯提议分布，用于生成新的参数链
// 注意：在 C++ 中，通常修改传入的 one_chain_new，而不是返回新向量，以避免复制
void do_gaussian_propose_cpp(const std::vector<double>& one_chain_old_cpp,
                             std::vector<double>& one_chain_new_cpp, // Modified in place
                             const std::vector<double>& sigma_prop_cpp,
                             int N_parm, std::string results_dir) {
    
    for(int i = 0; i < N_parm; ++i) {
        double x = one_chain_old_cpp[i];
        double sigma_prop_x = sigma_prop_cpp[i];
        double dx = r8_normal_01() * sigma_prop_x; // Assumes thread-safe RNG
        one_chain_new_cpp[i] = x + dx;
    }

    bool debug = false;
    //if (debug) {
        // Conversion needed if save function expects Rcpp type
         // Rcpp::NumericVector one_chain_rcpp = Rcpp::wrap(one_chain_new_cpp);
         // Rcpp::NumericVector sigma_prop_cpp = Rcpp::wrap(sigma_prop_cpp);
         // save_debug_gaussian_proposal_cpp(one_chain_rcpp, sigma_prop_cpp, results_dir, N_parm);
        // OR modify the save function signature (as done in forward decl)
        //save_debug_gaussian_proposal_cpp(one_chain_new_cpp, sigma_prop_cpp, results_dir, N_parm);
    //}
}

// 定义一个结构体来存储 iter_batch_mh 的并行结果 (每个rank一个)
struct IterBatchResult {
    int i_rank; // 0-based rank index
    std::vector<double> final_parm_cpp; // 存储参数结果
    double logpost_final;
    int rank_i_accumul; // 使用int可能更好，如果迭代次数很大用long long
    int rank_i_accumul_accept;
};


// 实现一个batch的MH算法 (修改为返回结构体，并在并行循环外处理)
// 注意：这个函数现在是*串行*的，将在 OpenMP 循环内部为每个 rank 调用
IterBatchResult iter_batch_mh_cpp(
    const std::vector<double>& initial_parm_cpp, // 只需传入当前rank的起始参数
    const std::vector<double>& sigma_prop_cpp,
    int n_iter_a_batch,
    int nline_data,
    const std::vector<std::vector<double>>& data_cpp,
    int i_rank, // 0-based index
    int i_accumul_start,      // 传入当前rank的起始累积计数
    int i_accumul_accept_start, // 传入当前rank的起始接受计数
    int i_save_begin,
    double logpost_start, // 传入当前rank的起始logpost
    const std::vector<double>& para_min_cpp, // 需要 para_min/max for boundary/prior
    const std::vector<double>& para_max_cpp,
    const std::vector<double>& beta_values_cpp, // 需要 Beta_Values for logll
    std::string results_dir,
    int N_parm
) {
    
    
    std::vector<std::vector<double>> chain_IterParm_cpp(n_iter_a_batch, std::vector<double>(N_parm));
    std::vector<double> logpost_cpp(n_iter_a_batch);
    std::vector<int> accumul_cpp(n_iter_a_batch); // 使用 Rcpp IntegerVector
    std::vector<int> accumul_accept_cpp(n_iter_a_batch);

    
    double logll_tempered_new = 0.0;
    double logprior_new = 0.0;
    double logpost_new = 0.0;
    double H = 0.0;
    double rand_unif = 0.0;
    double logpost_old = 0.0;

    std::vector<double> one_chain_old_cpp = initial_parm_cpp;
    std::vector<double> one_chain_new_cpp(N_parm);
    

    bool save_allch_ll = false; // R code had 0

    logpost_cpp[0] = logpost_start;
    
    chain_IterParm_cpp[0] = one_chain_old_cpp;
    
    int current_i_accumul = i_accumul_start;
    int current_i_accumul_accept = i_accumul_accept_start;

    // 开始迭代循环 (0-based loop for n_iter_a_batch iterations)
    for (int i = 0; i < n_iter_a_batch; ++i) {
        bool copy_old = true;

        // one_chain_old 和 logpost_old 已经在每次循环开始时/或上一次迭代结束时设置好
        if (i > 0) {
             one_chain_old_cpp = chain_IterParm_cpp[i - 1];
             logpost_old = logpost_cpp[i - 1];
        } else {
             logpost_old = logpost_start; // For the very first iteration (i=0)
             // one_chain_old is already set correctly before loop
        }

        // 实现了一个高斯提议分布
        do_gaussian_propose_cpp(one_chain_old_cpp, one_chain_new_cpp, sigma_prop_cpp, N_parm, results_dir);
        // q_factor is 1 for symmetric Gaussian proposal
        double q_factor = 1.0;
        
        // 检查参数是否在边界内
        one_chain_new_cpp = para_boundary_cpp(one_chain_new_cpp, para_min_cpp, para_max_cpp);

        // 用new_chain计算新的对数似然，对数先验，对数后验
        // 注意：logll_beta_cpp 需要 i_rank (0-based) 和 Beta_Values
        logll_tempered_new = logll_beta_cpp(one_chain_new_cpp, nline_data, data_cpp, i_rank, beta_values_cpp);
        logprior_new = log_prior_cpp(one_chain_new_cpp, para_min_cpp, para_max_cpp);
        logpost_new = logll_tempered_new + logprior_new;

        // MH-algorithm
        // (log(1/q_factor) is 0 since q_factor is 1)
        if ((logpost_new - logpost_old) > log(1 / q_factor)) {
            chain_IterParm_cpp[i] = one_chain_new_cpp; // Store accepted C++ vector
            logpost_cpp[i] = logpost_new;
            current_i_accumul_accept++;
            copy_old = false;
        } else {
            H = std::exp(logpost_new - logpost_old) * q_factor; // * q_factor (which is 1)
            rand_unif = r8_unif_ab(0.0, 1.0); 

            if (rand_unif < H) {
                chain_IterParm_cpp[i] = one_chain_new_cpp; // Store accepted C++ vector
                logpost_cpp[i] = logpost_new;
                current_i_accumul_accept++;
                copy_old = false;
            }
        }

        if (copy_old) { // If not accepted
             if (i > 0) {
                  chain_IterParm_cpp[i] = chain_IterParm_cpp[i - 1]; // Copy previous C++ vector
                  logpost_cpp[i] = logpost_cpp[i - 1]; // Use the previous logpost
             } else {
                  // For i=0, if not accepted, keep initial state
                  chain_IterParm_cpp[i] = one_chain_old_cpp;
                  logpost_cpp[i] = logpost_start;           // Keep initial logpost
             }
        }


        if (save_allch_ll) {
            save_log_posterior_cpp(one_chain_new_cpp, logll_tempered_new, logprior_new, results_dir, i_rank, N_parm); // R i_rank 是 1-based, C++ 是 0-based. 确认 save 函数期望哪个
        }

        current_i_accumul++;
        accumul_cpp[i] = current_i_accumul;
        accumul_accept_cpp[i] = current_i_accumul_accept;
    }

    // 如果累计迭代次数 >= i_save_begin，则保存当前 batch 的结果
    // 注意：R 的 i_accumul 是全局的，这里 current_i_accumul 是这个 rank 的
    if (current_i_accumul >= i_save_begin) {
         save_the_batch_cpp(chain_IterParm_cpp, n_iter_a_batch, results_dir, i_rank, logpost_cpp, accumul_cpp, accumul_accept_cpp, N_parm); // 再次检查 i_rank 基准
    }


    IterBatchResult result;
    result.i_rank = i_rank; // Store 0-based index
    result.final_parm_cpp = chain_IterParm_cpp[n_iter_a_batch - 1]; // Store final C++ vector
    result.logpost_final = logpost_cpp[n_iter_a_batch - 1];
    result.rank_i_accumul = current_i_accumul;
    result.rank_i_accumul_accept = current_i_accumul_accept;

    return result;
}

// 定义 judge_and_swap 的返回结构体
struct SwapResult {
    Rcpp::NumericMatrix transit_BetaParm;
    Rcpp::NumericVector logpost_all_ranks;
};

// R code for judge_and_swap -> C++
SwapResult judge_and_swap_cpp(
    Rcpp::NumericMatrix transit_BetaParm, // Pass copies or handle mutation carefully
    Rcpp::NumericVector logpost_all_ranks,
    int nline_data,
    const std::vector<std::vector<double>>& data_cpp,
    int i_swap, // 0-based index
    int j_swap, // 0-based index (calculated below or passed in)
    int Swapmode,
    int n_ranks,
    const std::vector<double>& para_min_cpp,    // <<< C++ input
    const std::vector<double>& para_max_cpp,    // <<< C++ input
    const std::vector<double>& beta_values_cpp, // <<< C++ input
    std::string results_dir,
    int N_parm
) {
    bool debug_save_swap = false; // R code had 0
    bool do_swap = false;
    double H = 0.0;
    double rand_unif = 0.0;
    
    
    // Note: i_swap was 1-based in R, now should be 0-based.
    // The R code generated j_swap based on 1-based i_swap. Adjust logic.
    if (Swapmode == 0) { // Adjacent swap
        // i_swap is 0 to n_ranks-2. j_swap should be i_swap + 1.
        j_swap = i_swap + 1;
    } else if (Swapmode == 1) { // Random swap
        // Choose j_swap between 0 and n_ranks-1, ensuring j_swap != i_swap
        j_swap = i4_unif_ab(0, n_ranks - 2); // Generate in [0, n_ranks-2]
        if (j_swap >= i_swap) {
            j_swap = j_swap + 1; // Map to [0, i_swap-1] U [i_swap+1, n_ranks-1]
        }
         // Double check bounds: if i_swap is n_ranks-1, i4_unif_ab(0, n_ranks-2) is fine.
         // If i_swap is 0, i4_unif_ab(0, n_ranks-2) gives [0, n_ranks-2]. If result is 0, j_swap becomes 1. Correct.
    } else {
         Rcpp::stop("Unsupported Swapmode.");
    }
    
    
    double logpost_ichain_ibeta = logpost_all_ranks[i_swap];
    double logpost_jchain_jbeta = logpost_all_ranks[j_swap];
    
    
    Rcpp::NumericVector chain_Parm_ichain = transit_BetaParm(i_swap, Rcpp::_);
    std::vector<double> chain_Parm_ichain_cpp = Rcpp::as<std::vector<double>>(chain_Parm_ichain); // Convert
    double logll_tempered_mix_ij = logll_beta_cpp(chain_Parm_ichain_cpp, nline_data, data_cpp, j_swap, beta_values_cpp);
    double logprior_i = log_prior_cpp(chain_Parm_ichain_cpp, para_min_cpp, para_max_cpp);
    double logpost_ichain_jbeta = logll_tempered_mix_ij + logprior_i;

    Rcpp::NumericVector chain_Parm_jchain = transit_BetaParm(j_swap, Rcpp::_);
    std::vector<double> chain_Parm_jchain_cpp = Rcpp::as<std::vector<double>>(chain_Parm_jchain); // Convert
    double logll_tempered_mix_ji = logll_beta_cpp(chain_Parm_jchain_cpp, nline_data, data_cpp, i_swap, beta_values_cpp);
    double logprior_j = log_prior_cpp(chain_Parm_jchain_cpp, para_min_cpp, para_max_cpp);
    double logpost_jchain_ibeta = logll_tempered_mix_ji + logprior_j;

    double log_accept_ratio = logpost_ichain_jbeta + logpost_jchain_ibeta - logpost_ichain_ibeta - logpost_jchain_jbeta;

    if (log_accept_ratio > 0.0) { // Corresponds to H > 1
        do_swap = true;
    } else {
        H = std::exp(log_accept_ratio);
        rand_unif = r8_unif_ab(0.0, 1.0); // Again, assumes thread safety if called in parallel context (but swap is likely serial)
        if (rand_unif < H) {
            do_swap = true;
        }
    }

    // Create copies to return modified versions
    Rcpp::NumericMatrix next_transit_BetaParm = Rcpp::clone(transit_BetaParm);
    Rcpp::NumericVector next_logpost_all_ranks = Rcpp::clone(logpost_all_ranks);


    if (do_swap) {
        // Swap rows in the matrix copy
        Rcpp::NumericVector temp = next_transit_BetaParm(i_swap, Rcpp::_);
        next_transit_BetaParm(i_swap, Rcpp::_) = next_transit_BetaParm(j_swap, Rcpp::_);
        next_transit_BetaParm(j_swap, Rcpp::_) = temp;

        // Update logposts in the vector copy
        next_logpost_all_ranks[i_swap] = logpost_jchain_ibeta;
        next_logpost_all_ranks[j_swap] = logpost_ichain_jbeta;

    }

    if (debug_save_swap) {
        save_debug_stack_doswap_cpp(do_swap, results_dir);
    }

    SwapResult result;
    result.transit_BetaParm = next_transit_BetaParm;
    result.logpost_all_ranks = next_logpost_all_ranks;
    return result;
}


// Helper for calc_SD_allParm
Rcpp::NumericVector colMeans(const Rcpp::NumericMatrix& x) {
    int nc = x.ncol();
    int nr = x.nrow();
    Rcpp::NumericVector means(nc);
    for (int j = 0; j < nc; ++j) {
        double sum = 0;
        for (int i = 0; i < nr; ++i) {
            sum += x(i, j);
        }
        means[j] = sum / nr;
    }
    return means;
}

// R code for calc_SD_allParm -> C++
// Calculates *column* standard deviations (assuming rows=parm, cols=rank in ar_ParmNvaried)
// But R code calculated *row* standard deviations. Let's match R.
Rcpp::NumericVector calc_SD_allParm_cpp(const Rcpp::NumericMatrix& ar_ParmNvaried) { // n_ranks is ar_ParmNvaried.ncol()
    int N_parm = ar_ParmNvaried.nrow();
    int n_ranks = ar_ParmNvaried.ncol();
    Rcpp::NumericVector std_Parm(N_parm);

    if (n_ranks == 0) return std_Parm; // Avoid division by zero

    for (int i = 0; i < N_parm; ++i) {
        Rcpp::NumericVector row = ar_ParmNvaried(i, Rcpp::_);
        double mean_oneparm = Rcpp::mean(row); // Rcpp sugar for mean
        double variance_oneparm = 0.0;
        for (int j = 0; j < n_ranks; ++j) {
            variance_oneparm += std::pow(ar_ParmNvaried(i, j) - mean_oneparm, 2);
        }
        std_Parm[i] = std::sqrt(variance_oneparm / n_ranks); // Use population variance formula as in R code
    }
    return std_Parm;
}

// R argmax/argmin -> C++ indices (0-based)
int argmax_cpp(const Rcpp::NumericVector& array) {
    if (array.size() == 0) return -1; // Or throw error
    return std::distance(array.begin(), std::max_element(array.begin(), array.end()));
}

int argmin_cpp(const Rcpp::NumericVector& array) {
     if (array.size() == 0) return -1; // Or throw error
     return std::distance(array.begin(), std::min_element(array.begin(), array.end()));
}


struct ClosestArResult {
    double closest_value; // The minimum absolute difference
    int iparm_closest; // 0-based row index
    int irank_closest; // 0-based column index
};

// R calc_closest_ar -> C++
ClosestArResult calc_closest_ar_cpp(const Rcpp::NumericMatrix& ar_ParmNvaried, double ar_best) {
    ClosestArResult result = {std::numeric_limits<double>::max(), -1, -1};

    int N_parm = ar_ParmNvaried.nrow();
    int n_ranks = ar_ParmNvaried.ncol();

    for (int i = 0; i < N_parm; ++i) {
        for (int j = 0; j < n_ranks; ++j) {
            double abs_diff = std::abs(ar_ParmNvaried(i, j) - ar_best);
            if (abs_diff < result.closest_value) {
                result.closest_value = abs_diff;
                result.iparm_closest = i;
                result.irank_closest = j;
            }
        }
    }
    return result;
}


struct ClosestArOneParmResult {
    double closest_value_oneparm;
    int irank_closest_oneparm; // 0-based index
};

// R calc_closest_ar_oneparm -> C++
ClosestArOneParmResult calc_closest_ar_oneparm_cpp(const Rcpp::NumericMatrix& ar_ParmNvaried, int iparm_stdmax, double ar_best) { // iparm_stdmax is 0-based
    ClosestArOneParmResult result = {std::numeric_limits<double>::max(), -1};
    int n_ranks = ar_ParmNvaried.ncol();

    Rcpp::NumericVector ar_row = ar_ParmNvaried(iparm_stdmax, Rcpp::_);
    Rcpp::NumericVector abs_diff = Rcpp::abs(ar_row - ar_best); // Rcpp sugar

    if (abs_diff.size() > 0) {
         result.irank_closest_oneparm = argmin_cpp(abs_diff); // Find index of minimum difference
         result.closest_value_oneparm = abs_diff[result.irank_closest_oneparm]; // Get the minimum difference itself
    }

    return result;
}

struct DecideSigmaResult {
    int iparm_change; // 0-based index, or -1
    int irank_change; // 0-based index, or -1
};

// R decide_sigma_to_change -> C++
DecideSigmaResult decide_sigma_to_change_cpp(
    const Rcpp::NumericMatrix& ar_ParmNvaried,
    double ar_best,
    double ar_accept_diff,
    int N_parm
    
) {
    

    
    DecideSigmaResult result = {-1, -1}; // Default to special values
    int n_ranks = ar_ParmNvaried.ncol();

    Rcpp::NumericVector std_Parm = calc_SD_allParm_cpp(ar_ParmNvaried);

    int iparm_stdmax = argmax_cpp(std_Parm);

    // --- Scheme 1 ---
    int iparm_change1 = iparm_stdmax;
    ClosestArOneParmResult results_oneparm = calc_closest_ar_oneparm_cpp(ar_ParmNvaried, iparm_stdmax, ar_best);
    double v_oneparm_closest = results_oneparm.closest_value_oneparm;
    int irank_change1 = results_oneparm.irank_closest_oneparm;


    // --- Scheme 2 ---
    ClosestArResult results_2d = calc_closest_ar_cpp(ar_ParmNvaried, ar_best);
    double v_2d_closest = results_2d.closest_value;
    int iparm_change2 = results_2d.iparm_closest;
    int irank_change2 = results_2d.irank_closest;

    // --- Decision Logic ---
    if (v_oneparm_closest < ar_accept_diff) {
        if (v_2d_closest < ar_accept_diff) {
            // Both good, choose randomly
            double rand_unif = r8_unif_ab(0.0, 1.0);
            if (rand_unif > 0.5) {
                result.iparm_change = iparm_change1;
                result.irank_change = irank_change1;
            } else {
                result.iparm_change = iparm_change2;
                result.irank_change = irank_change2;
            }
        } else {
            // Only scheme 1 is good
            result.iparm_change = iparm_change1;
            result.irank_change = irank_change1;
        }
    } else {
        if (v_2d_closest < ar_accept_diff) {
            // Only scheme 2 is good
            result.iparm_change = iparm_change2;
            result.irank_change = irank_change2;
        } else {
            // Neither is good, random special values
            double rand_unif = r8_unif_ab(0.0, 1.0);
            if (rand_unif > 0.5) {
                result.iparm_change = -1;
                result.irank_change = 0; // Use 1/-1 as in R code
            } else {
                result.iparm_change = -1;
                result.irank_change = -1; // Use 1/-1 as in R code
            }
        }
    }

    return result;
}

// R race_all_parm -> C++
Rcpp::NumericVector race_all_parm_cpp(const Rcpp::NumericMatrix& sigma_alltune_ParmNvaried, const Rcpp::NumericMatrix& ar_ParmNvaried, double ar_best, int N_parm) {
    Rcpp::NumericVector sigma_prop_rankintune(N_parm);
    int n_ranks = ar_ParmNvaried.ncol();
    if (n_ranks == 0) return sigma_prop_rankintune; // Return empty/zero vector

    for (int i = 0; i < N_parm; ++i) {
         Rcpp::NumericVector ar_row = ar_ParmNvaried(i, Rcpp::_);
         Rcpp::NumericVector diff_row = Rcpp::abs(ar_row - ar_best);
         int j_min = argmin_cpp(diff_row); // 0-based index of best rank for this param
         if (j_min != -1) {
             sigma_prop_rankintune[i] = sigma_alltune_ParmNvaried(i, j_min);
         } else {
             // Handle error: couldn't find min? Maybe default or warn.
             sigma_prop_rankintune[i] = 0.0; // Or some default
             Rcpp::warning("Could not find best sigma for parameter %d in race_all_parm_cpp", i + 1);
         }
    }
    return sigma_prop_rankintune;
}

// R modify_sigma_prop_rankintune -> C++
// Modifies sigma_RanksParm in place
void modify_sigma_prop_rankintune_cpp(
    int rank_in_tune, // 0-based index
    int iparm_change, // 0-based index, or -1
    int irank_change, // 0-based index, or -1 (from decide_sigma)
    Rcpp::NumericMatrix& sigma_RanksParm, // Modify in place
    const Rcpp::NumericMatrix& sigma_alltune_ParmNvaried,
    const Rcpp::NumericMatrix& ar_ParmNvaried,
    double ar_best, // Needed for race_all_parm
    double init_gp_ratio, // Needed for init_gaussian_proposal
    const Rcpp::NumericVector& para_min, // Needed for init_gaussian_proposal
    const Rcpp::NumericVector& para_max, // Needed for init_gaussian_proposal
    int N_parm
) {

    if (iparm_change >= 0) { // Corresponds to R iparm_change >= 1
        // Ensure irank_change is valid before using
        if(irank_change >= 0)
        {
            sigma_RanksParm(rank_in_tune, iparm_change) = sigma_alltune_ParmNvaried(iparm_change, irank_change);
        }

    } else { // iparm_change < 0 (was -1 in R)
        Rcpp::NumericVector sigma_prop_rankintune(N_parm);
        if (irank_change >= 0) { // Corresponds to R irank_change >= 1 (was 1 in R code) - check if this mapping is correct
            // R code called init_gaussian_proposal with init_gp_ratio
            // Assuming this initializes based on parameter ranges
            sigma_prop_rankintune = init_gaussian_proposal_cpp(para_min, para_max, init_gp_ratio);
        } else { // irank_change < 0 (was -1 in R code)
            // R code called race_all_parm
            sigma_prop_rankintune = race_all_parm_cpp(sigma_alltune_ParmNvaried, ar_ParmNvaried, ar_best, N_parm);
        }

        sigma_RanksParm(rank_in_tune, Rcpp::_) = sigma_prop_rankintune;
        
    }
    // No return needed, sigma_RanksParm modified by reference
}

// R create_logspace_array -> C++
Rcpp::NumericVector create_logspace_array_cpp(double base_value, double half_scale, int n_points) {
    Rcpp::NumericVector result(n_points);

    double log_min = std::log(base_value / half_scale);
    double log_max = std::log(base_value * half_scale);
    double log_step = (log_max - log_min) / (n_points - 1);

    for (int i = 0; i < n_points; ++i) {
        result[i] = std::exp(log_min + i * log_step);
    }
    return result;
}

// R check_bounceInside_sigma_boundary -> C++
Rcpp::NumericVector check_bounceInside_sigma_boundary_cpp(
    Rcpp::NumericVector scaled_arr, // Pass by value (copy) as it might be replaced
    const Rcpp::NumericVector& sigma_parm_min,
    const Rcpp::NumericVector& sigma_parm_max,
    double sigma_jumpin_ratio, // Note: R code uses this name but value seems like 'scale'
    int i_parm // 0-based index for the parameter being checked
) {
    int n_ranks = scaled_arr.size();

    double half_scale = std::sqrt(sigma_jumpin_ratio); // As in R code
    double current_min = sigma_parm_min[i_parm];
    double current_max = sigma_parm_max[i_parm];


    // Check lower bound
    if (scaled_arr[0] < current_min) {
         // Rcpp::Rcout << "Sigma boundary bounce (low) for param " << i_parm << std::endl;
        scaled_arr = create_logspace_array_cpp(current_min * half_scale, half_scale, n_ranks);
        // Need to re-check the upper bound in case the bounce pushed it too high
    }

     // Check upper bound (potentially after bouncing from low)
    if (scaled_arr[n_ranks - 1] > current_max) {
         // Rcpp::Rcout << "Sigma boundary bounce (high) for param " << i_parm << std::endl;
        scaled_arr = create_logspace_array_cpp(current_max / half_scale, half_scale, n_ranks);
         // Re-check lower bound? Could potentially oscillate if bounds are very tight.
         // The R code didn't re-check the lower bound here. Let's stick to that for now.
          // if (n_ranks > 0 && scaled_arr[0] < current_min) {
          //      Rcpp::warning("Sigma boundary oscillation detected for param %d", i_parm + 1);
          //      // Maybe just clamp?
          //      scaled_arr = create_logspace_array_cpp(current_min * half_scale, half_scale, n_ranks);
          // }
    }


    return scaled_arr;
}


// R gen_sigma_alltune -> C++
Rcpp::NumericMatrix gen_sigma_alltune_cpp(
    const Rcpp::NumericVector& sigma_parm_min,
    const Rcpp::NumericVector& sigma_parm_max,
    double sigma_jumpin_ratio,
    double sigma_scale_half_ratio, // Need this value from R code context
    const Rcpp::NumericMatrix& sigma_RanksParm,
    int rank_in_tune, // 0-based index
    int n_ranks,
    int N_parm
) {
    Rcpp::NumericMatrix sigma_alltune_ParmNvaried(N_parm, n_ranks);

    for (int i = 0; i < N_parm; ++i) {
        Rcpp::NumericVector scaled_arr = create_logspace_array_cpp(sigma_RanksParm(rank_in_tune, i), sigma_scale_half_ratio, n_ranks);
        scaled_arr = check_bounceInside_sigma_boundary_cpp(scaled_arr, sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, i);
        sigma_alltune_ParmNvaried(i, Rcpp::_) = scaled_arr;
    }
    // scaled_arr goes out of scope, no need to NULLify
    return sigma_alltune_ParmNvaried;
}


// R iter_batch_mh_tune -> C++
// Returns the acceptance rate (double)
double iter_batch_mh_tune_cpp(
    const std::vector<double>& chain_single_tune_cpp, // Initial chain state
    const std::vector<double>& sigma_prop_cpp, // Sigma for this specific tuning iteration
    int n_iter_in_tune,
    int nline_data,
    const std::vector<std::vector<double>>& data_cpp,
    int rank_in_tune, // 0-based index: the rank being tuned
    double logpost_old_start, // Initial log posterior
    int i_rank, // 0-based index: used only for saving output? Matches the loop variable in R.
    const std::vector<double>& para_min_cpp,
    const std::vector<double>& para_max_cpp,
    const std::vector<double>& beta_values_cpp, // Needed for logll
    std::string results_dir,
    int N_parm
    ) {
    

    
    bool save_debug = false; // R code had 1
    bool save_allch_ll = false; // R code had 1

    std::vector<std::vector<double>> chain_IterParm_cpp(n_iter_in_tune, std::vector<double>(N_parm));
    std::vector<double> logpost(n_iter_in_tune);

    std::vector<double> one_chain_old_cpp = chain_single_tune_cpp;
    std::vector<double> one_chain_new_cpp(N_parm);
    


    double logpost_old = logpost_old_start;
    logpost[0] = logpost_old;
    chain_IterParm_cpp[0] = one_chain_old_cpp; // Store initial C++ vector


    int i_accumul_accept_tune = 0; // Local acceptance counter for this tune run

    for (int i = 0; i < n_iter_in_tune; ++i) {
        bool copy_old = true;

        if (i > 0) {
            logpost_old = logpost[i - 1];
            one_chain_old_cpp = chain_IterParm_cpp[i - 1];
        } else {
             logpost_old = logpost_old_start;
             
        }

        do_gaussian_propose_cpp(one_chain_old_cpp, one_chain_new_cpp, sigma_prop_cpp, N_parm, results_dir); // Generate proposal
        one_chain_new_cpp = para_boundary_cpp(one_chain_new_cpp, para_min_cpp, para_max_cpp); // Check boundaries
        double q_factor = 1.0; // Symmetric proposal

        // Calculate log likelihood using the beta of the rank being tuned
        double logll_tempered_new = logll_beta_cpp(one_chain_new_cpp, nline_data, data_cpp, rank_in_tune, beta_values_cpp);
        double logprior_new = log_prior_cpp(one_chain_new_cpp, para_min_cpp, para_max_cpp);
        double logpost_new = logll_tempered_new + logprior_new;


        if ((logpost_new - logpost_old) > log(1 / q_factor)) { // log(1/q_factor) = 0
          // Bounds check for assignment
            chain_IterParm_cpp[i] = one_chain_new_cpp;
            logpost[i] = logpost_new;
            i_accumul_accept_tune++;
            copy_old = false;
        } else {
            double H = std::exp(logpost_new - logpost_old) * q_factor;
            double rand_unif = r8_unif_ab(0.0, 1.0); // Assumes thread safety
            if (rand_unif < H) {
                chain_IterParm_cpp[i] = one_chain_new_cpp;
                logpost[i] = logpost_new;
                i_accumul_accept_tune++;
                copy_old = false;
            }
        }

        if (copy_old) { // Not accepted
            if (i > 0) {
                chain_IterParm_cpp[i] = chain_IterParm_cpp[i - 1];
                logpost[i] = logpost[i - 1];
            } else {
                // If first iteration is rejected, keep initial state
                chain_IterParm_cpp[i] = one_chain_old_cpp;
                logpost[i] = logpost_old_start;
            }
        }


        if (save_allch_ll) {
            // i_rank here seems to correspond to the sigma being tested (from the outer loop in R)
            save_log_posterior_cpp(one_chain_new_cpp, logll_tempered_new, logprior_new, results_dir, i_rank, N_parm); // Check i_rank base
        }
    }

    //if (save_debug) {
        // i_rank again seems to be the sigma index, rank_in_tune is the actual chain rank
      //  save_the_batch_tune_cpp(chain_IterParm_cpp, n_iter_in_tune, results_dir, i_rank, rank_in_tune, logpost, N_parm); // Check indices bases
    //}

    double ar = (double)i_accumul_accept_tune / n_iter_in_tune;
    return ar;
}


// --- Main Exported Function ---

/// 
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(
    int N_iter,
    int n_ranks,
    int N_parm,
    int nline_data,
    const Rcpp::NumericMatrix& data_NlineNdim,
    Rcpp::NumericMatrix sigma_RanksParm, // Modifiable copy
    Rcpp::NumericMatrix transit_BetaParm, // Modifiable copy
    Rcpp::NumericVector logpost_all_ranks, // Modifiable copy
    Rcpp::IntegerVector i_accumul,        // Modifiable copy (use IntegerVector)
    Rcpp::IntegerVector i_accumul_accept, // Modifiable copy (use IntegerVector)
    int n_iter_a_stack_base, // Renamed from n_iter_a_stack to avoid confusion with loop variable
    int n_iter_a_batch_base,
    int n_iter_a_batch_rand,
    int Swapmode,
    int i_save_begin,
    int N_stoptune,
    double ar_ok_upper,
    double ar_ok_lower,
    const Rcpp::NumericVector& sigma_parm_min,
    const Rcpp::NumericVector& sigma_parm_max,
    double sigma_jumpin_ratio, // Assume this is the scale factor
    double sigma_scale_half_ratio, // Need this value
    int n_iter_in_tune,
    double ar_best,
    double ar_accept_diff,
    double init_gp_ratio,
    const Rcpp::NumericVector& para_min,
    const Rcpp::NumericVector& para_max,
    const Rcpp::NumericVector& Beta_Values, // Temp scaling for loglik
    std::string results_dir
    
    
) {
    
    
    std::vector<std::vector<double>> data_cpp(data_NlineNdim.nrow(), std::vector<double>(data_NlineNdim.ncol()));
    for(int i=0; i<data_NlineNdim.nrow(); ++i) { for (int j=0; j<data_NlineNdim.ncol(); ++j) { data_cpp[i][j] = data_NlineNdim(i,j); }}
    std::vector<double> para_min_cpp = Rcpp::as<std::vector<double>>(para_min);
    std::vector<double> para_max_cpp = Rcpp::as<std::vector<double>>(para_max);
    std::vector<double> beta_values_cpp = Rcpp::as<std::vector<double>>(Beta_Values);

    // --- Initialization ---
    int i_tmp_stack = 0;
    int i_next_stack = 0;
    // Rcpp::IntegerVector i_accumul_accept_int = Rcpp::as<Rcpp::IntegerVector>(i_accumul_accept); // Work with integers

    // Variables for inner loop
    int i_tmp = 0;
    int i_next = 0;
    int n_iter_a_batch = 0;

    // Tuning related
    Rcpp::IntegerVector tune_ranks(n_ranks, 0); // Initialize with 0

    bool save_debug = false; // From R code context
    bool processing_debug = true; // Use bool
    

    // --- Outer Stack Loop ---
    while (i_next_stack < N_iter) {

        // Calculate iterations for this stack
        // Use a local variable for n_iter_a_stack within this loop iteration
        int n_iter_a_stack = n_iter_a_stack_base; // Assuming n_iter_a_stack was constant per stack in R? Or should it be calculated? Let's assume base value.
                                                 // If it depends on previous state, logic needs adjustment.

        i_next_stack = i_tmp_stack + n_iter_a_stack;
        if (i_next_stack > N_iter) {
            n_iter_a_stack = N_iter - i_tmp_stack;
            i_next_stack = N_iter;
        }
        
        // Store acceptance count at start of stack for rate calculation
         Rcpp::IntegerVector n_accept_old = Rcpp::clone(i_accumul_accept);


        // Save initial sigma props for this stack (seems to happen once per stack in R)
        for (int i = 0; i < n_ranks; ++i) {
            Rcpp::NumericVector sigma_prop = sigma_RanksParm(i, Rcpp::_);
            // Assuming save_sigma_gauss_prop_cpp expects 0-based i_rank
            save_sigma_gauss_prop_cpp(sigma_prop, i, results_dir, N_parm);
        }

        // --- Inner Batch Loop (within a stack) ---
        i_tmp = 0; // Reset inner loop counter for each stack
        i_next = 0; // Reset inner loop counter

        while (i_next < n_iter_a_stack) {

            // Calculate iterations for this batch
            int n_iter_a_batch_adjust = i4_unif_ab(-n_iter_a_batch_rand, n_iter_a_batch_rand); // Assumes i4_unif_ab is inclusive [a, b]
            n_iter_a_batch = n_iter_a_batch_base + n_iter_a_batch_adjust;
            n_iter_a_batch = std::max(1, n_iter_a_batch); // Ensure at least 1 iteration
            i_next = i_tmp + n_iter_a_batch;

            if (i_next > n_iter_a_stack) {
                n_iter_a_batch = n_iter_a_stack - i_tmp;
                i_next = n_iter_a_stack;
            }

            // --- Batch MH Step (Parallel) ---
            
            std::vector<std::vector<double>> initial_parm_cpp(n_ranks, std::vector<double>(N_parm));
	    std::vector<std::vector<double>> sigma_prop_cpp(n_ranks, std::vector<double>(N_parm));
	    
	    for (int r = 0; r < n_ranks; ++r) {
	         Rcpp::NumericVector initial_parm_rcpp = transit_BetaParm(r, Rcpp::_);
	         Rcpp::NumericVector sigma_prop_rcpp = sigma_RanksParm(r, Rcpp::_);
	    // Perform conversion SERIALLY here
	         initial_parm_cpp[r] = Rcpp::as<std::vector<double>>(initial_parm_rcpp);
	         sigma_prop_cpp[r] = Rcpp::as<std::vector<double>>(sigma_prop_rcpp);
	     }
	     
	     std::vector<IterBatchResult> batch_results(n_ranks); // Store results temporarily
	     
	     

	     #pragma omp parallel for schedule(static) // Use static schedule for potentially balanced load
	     for (int i_rank = 0; i_rank < n_ranks; ++i_rank) {
		        
		    
		        // Get initial state for this rank for this batch
		         
		double logpost_start = logpost_all_ranks[i_rank];
		int i_accumul_start = i_accumul[i_rank];
		int i_accumul_accept_start = i_accumul_accept[i_rank];
		         

		        // Call the *serial* MH batch function for this rank
		        
		        
		batch_results[i_rank] = iter_batch_mh_cpp(
		   initial_parm_cpp[i_rank], sigma_prop_cpp[i_rank], n_iter_a_batch,
		   nline_data, data_cpp,
		   i_rank, // Pass 0-based index
		   i_accumul_start, i_accumul_accept_start,
		   i_save_begin, logpost_start,
		   para_min_cpp, para_max_cpp, beta_values_cpp,
		   results_dir, N_parm
		);
		        
	}
	     // End of parallel for loop

            // --- Combine results from parallel execution ---
            // This part is serial again
            
            for (int i_rank = 0; i_rank < n_ranks; ++i_rank) {
                const IterBatchResult& res = batch_results[i_rank];
                // Ensure results are for the correct rank (optional check)
                if (res.i_rank == i_rank) {
                    Rcpp::NumericVector final_parm_rcpp = Rcpp::wrap(res.final_parm_cpp);
                    transit_BetaParm(i_rank, Rcpp::_) = final_parm_rcpp;
                    logpost_all_ranks[i_rank] = res.logpost_final;
                    i_accumul[i_rank] = res.rank_i_accumul;
                    i_accumul_accept[i_rank] = res.rank_i_accumul_accept;
                } else {
                    Rcpp::warning("Mismatch in batch result rank index! Expected %d, got %d", i_rank, res.i_rank);
                }
            }
            

            // --- Swap Step (Serial) ---
            // Swap only makes sense if there's more than one rank
            // Choose i_swap (0-based index from 0 to n_ranks-2 for adjacent swap mode)
            int i_swap = i4_unif_ab(0, n_ranks - 2); // R code used 1 to n_ranks-1, so 0 to n_ranks-2 for 0-based adjacent swap
            int j_swap = -1; // Will be determined by judge_and_swap_cpp based on Swapmode
            
            if (save_debug) {
              // Save sequence before potential swap. Use 0-based i_swap.
              save_debug_stack_sequence_cpp(i_accumul[i_swap], i_swap, results_dir);
            }
            
            // Call swap function (passing copies/getting results back)
            SwapResult swap_results = judge_and_swap_cpp(
              transit_BetaParm, logpost_all_ranks,
              nline_data, data_cpp,
              i_swap, j_swap, Swapmode, n_ranks,
              para_min_cpp, para_max_cpp, beta_values_cpp,
              results_dir, N_parm
            );
            
            // Update state with swap results
            transit_BetaParm = swap_results.transit_BetaParm;
            logpost_all_ranks = swap_results.logpost_all_ranks;
            
            // Update inner loop counter
            i_tmp = i_next;

        } // End of inner batch while loop

        // --- Post-Stack Processing ---
        if (processing_debug) {
            Rprintf("%9d out of total %d iterations have completed.\n", i_next_stack, N_iter);
        }

        // Calculate acceptance rate for the completed stack
        Rcpp::NumericVector accept_rate_a_stack(n_ranks);
        Rcpp::IntegerVector n_accept_now = i_accumul_accept; // Current accept counts

        for(int i=0; i < n_ranks; ++i) {
          accept_rate_a_stack[i] = static_cast<double>(n_accept_now[i] - n_accept_old[i]) / n_iter_a_stack;
        }

        if (save_debug) {
            save_ar_stack_cpp(i_next_stack, n_iter_a_stack, accept_rate_a_stack, n_ranks, results_dir);
        }

        // --- Tuning Step (Serial overall, but parallel inside) ---
        if ((i_next_stack < N_stoptune) && (i_next_stack < N_iter)) {

            // Find ranks to tune
             std::fill(tune_ranks.begin(), tune_ranks.end(), 0); // Reset tune_ranks
            for (int i = 0; i < n_ranks; ++i) {
                if ((accept_rate_a_stack[i] < ar_ok_lower) || (accept_rate_a_stack[i] > ar_ok_upper)) {
                    tune_ranks[i] = 1;
                }
            }

            // Iterate through ranks needing tuning
            for (int rank_in_tune = 0; rank_in_tune < n_ranks; ++rank_in_tune) {
                if (tune_ranks[rank_in_tune] == 1) {

                    if (processing_debug) {
                        Rprintf("     tuning rank %d...\n", rank_in_tune + 1); // Print 1-based rank
                    }

                    // Generate candidate sigmas for this rank
                    Rcpp::NumericMatrix sigma_alltune_ParmNvaried = gen_sigma_alltune_cpp(
                        sigma_parm_min, sigma_parm_max, sigma_jumpin_ratio, sigma_scale_half_ratio,
                        sigma_RanksParm, rank_in_tune, n_ranks, N_parm
                    );

                    // Get current state for the rank being tuned
                    double logpost_single_tune = logpost_all_ranks[rank_in_tune];
                    Rcpp::NumericVector chain_single_tune_rcpp = transit_BetaParm(rank_in_tune, Rcpp::_);
                    std::vector<double> chain_single_tune_cpp = Rcpp::as<std::vector<double>>(chain_single_tune_rcpp);

                    // Matrix to store AR results (Rows: Param, Cols: Sigma candidate index/Rank)
                    Rcpp::NumericMatrix ar_ParmNvaried(N_parm, n_ranks);

                    // --- Parallel Tuning Simulation ---
                    // Iterate through each parameter j_parm
                    for(int j_parm = 0; j_parm < N_parm; ++j_parm) {

                         // Create sigma matrix for this parameter test:
                         // Start with current sigma for the rank being tuned
                         Rcpp::NumericMatrix sigma_tune1parm_NvariedParm(n_ranks, N_parm);
                         Rcpp::NumericVector current_sigma_for_rank = sigma_RanksParm(rank_in_tune, Rcpp::_);
                         for(int k=0; k<n_ranks; ++k) {
                              sigma_tune1parm_NvariedParm(k, Rcpp::_) = current_sigma_for_rank;
                         }
                         // Now replace the column for j_parm with the candidate sigmas
                         Rcpp::NumericVector candidate_sigmas_for_jparm = sigma_alltune_ParmNvaried(j_parm, Rcpp::_);
                         for(int k=0; k<n_ranks; ++k) {
                              sigma_tune1parm_NvariedParm(k, j_parm) = candidate_sigmas_for_jparm[k];
                         }

                        // Vector to store AR for this j_parm across different sigma candidates (ranks)
                        
                        std::vector<std::vector<double>> sigma_tune1parm_NvariedParm_cpp(n_ranks, std::vector<double>(N_parm));
			for(int k=0; k<n_ranks; ++k) {
			    Rcpp::NumericVector sigma_row_rcpp = sigma_tune1parm_NvariedParm(k, Rcpp::_);
			    sigma_tune1parm_NvariedParm_cpp[k] = Rcpp::as<std::vector<double>>(sigma_row_rcpp);
			}
			
			std::vector<double> ar_jparm_ranks_cpp(n_ranks);
			
			


		        #pragma omp parallel for schedule(static)
		            for (int i_rank_sigma = 0; i_rank_sigma < n_ranks; ++i_rank_sigma) {
		                    // sigma_prop is the i_rank_sigma'th row from sigma_tune1parm_NvariedParm
		                    
		                    
		                const std::vector<double>& sigma_prop_cpp = sigma_tune1parm_NvariedParm_cpp[i_rank_sigma];


		                    // Call the tuning MH function
		                double result_ar = iter_batch_mh_tune_cpp(
		                    chain_single_tune_cpp, sigma_prop_cpp, n_iter_in_tune,
		                    nline_data, data_cpp,
		                    rank_in_tune, // The actual rank being tuned
		                    logpost_single_tune,
		                    i_rank_sigma, // Index for the sigma being tested (used for saving)
		                    para_min_cpp, para_max_cpp, beta_values_cpp,
		                    results_dir, N_parm
		                );
		                    
		                ar_jparm_ranks_cpp[i_rank_sigma] = result_ar;
		                    
		            } // End parallel loop for sigma candidates
                        
                        Rcpp::NumericVector ar_jparm_ranks(n_ranks);
                        std::copy(ar_jparm_ranks_cpp.begin(), ar_jparm_ranks_cpp.end(), ar_jparm_ranks.begin());

                        // Store results for this j_parm
                        ar_ParmNvaried(j_parm, Rcpp::_) = ar_jparm_ranks;

                    } // End loop over j_parm


                    // --- Decide which sigma to change (Serial) ---
                    DecideSigmaResult sigma_decision = decide_sigma_to_change_cpp(
                        ar_ParmNvaried, ar_best, ar_accept_diff, N_parm
                    );
                    int iparm_change = sigma_decision.iparm_change; // 0-based or -1
                    int irank_change = sigma_decision.irank_change; // 0-based or -1/1


                    // --- Modify Sigma (Serial) ---
                    modify_sigma_prop_rankintune_cpp(
                        rank_in_tune, iparm_change, irank_change,
                        sigma_RanksParm, // Modify in place
                        sigma_alltune_ParmNvaried, ar_ParmNvaried, ar_best, init_gp_ratio,
                        para_min, para_max, N_parm
                    );

                    // --- Save Tuning Debug Info (Serial) ---
                    bool save_tuning_debug = false; // From R code context
                    if (save_tuning_debug) {
                         // Need to implement or get save_tuning_sigma_ar_cpp
                         save_tuning_sigma_ar_cpp(ar_ParmNvaried, sigma_alltune_ParmNvaried, results_dir, rank_in_tune, n_ranks, N_parm);
                    }
                } // End if(tune_ranks[rank_in_tune] == 1)
            } // End for loop over ranks to tune

        } // End tuning step condition


        // Update stack counter for next iteration
        i_tmp_stack = i_next_stack;

    } // End of outer stack while loop


    // --- Return Results ---
    // Return the final state of the chains and related info
    return Rcpp::List::create(
        Rcpp::Named("sigma_RanksParm") = sigma_RanksParm,
        Rcpp::Named("transit_BetaParm") = transit_BetaParm,
        Rcpp::Named("logpost_all_ranks") = logpost_all_ranks,
        Rcpp::Named("i_accumul") = i_accumul,
        Rcpp::Named("i_accumul_accept") = i_accumul_accept
        // Add any other state variables you want to return
    );
}
