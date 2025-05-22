#' @useDynLib cppfiles, .registration = TRUE
#' @export save_the_seed_cpp
#' @export save_first_chain_cpp
#' @export save_init_parm_cpp
#' @export save_sigma_gauss_prop_cpp
#' @export save_ar_stack_cpp
#' @export save_the_batch_cpp
#' @export save_log_posterior_cpp
#' @export save_debug_gaussian_proposal_cpp
#' @export save_debug_stack_doswap_cpp
#' @export save_debug_stack_sequence_cpp
#' @export save_tuning_sigma_ar_cpp
#' @export save_the_batch_tune_cpp
#' @export r8_normal_01
#' @export r8_normal_ab
#' @export r8_logunif_ab
#' @export r8_unif_ab
#' @export i4_unif_ab
#' @export i4_unif_0a
#' @export normal_pdf
#' @export init_gaussian_proposal_cpp
#' @export bounce_inside_cpp
#' @export para_boundary_cpp
#' @export prior_para0_cpp
#' @export prior_para1_cpp
#' @export prior_para2_cpp
#' @export prior_para3_cpp
#' @export prior_para4_cpp
#' @export prior_para5_cpp
#' @export prior_para6_cpp
#' @export prior_para7_cpp
#' @export prior_para8_cpp
#' @export prior_para9_cpp
#' @export prior_para10_cpp
#' @export prior_para11_cpp
#' @export prior_para12_cpp
#' @export prior_para13_cpp
#' @export prior_para14_cpp
#' @export log_prior_cpp
#' @export save_debug_para_boundary_cpp
#' @export period_to_au_cpp
#' @export calc_osi_cpp
#' @export newton_solver_cpp
#' @export func_as_cpp
#' @export log_likelihood_cpp
#' @export logll_beta_cpp
#' @export run_mcmc_cpp
#' @keywords internal
"_PACKAGE"

# 你还可以在这里为导出的 C++ 函数创建 R 包装器（如果需要的话），
# 或者直接依赖 Rcpp 的自动包装。

# 如果你想确保所有 [[Rcpp::export]] 的函数都被导出到 R 用户，
# 可以在这个文件或者包含 Rcpp::export 的 .cpp 文件顶部（用 /// 注释）
# 添加 @export 标签，例如在 save.cpp 顶部：
# /// @export save_the_seed_cpp
# /// @export save_init_parm_cpp
# // [[Rcpp::export]]
# int save_the_seed_cpp(...) { ... }
#
# 或者，更简单的方式是保持你NAMESPACE里的 exportPattern，
# Rcpp 会为你导出的 C++ 函数创建同名的 R 包装函数，这些 R 函数
# 会因为 exportPattern 被导出。你需要确保 C++ 函数名以字母开头。
