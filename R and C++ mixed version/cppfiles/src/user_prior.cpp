#include <Rcpp.h>
#include <Rmath.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <random>

using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: init gaussian proposal
///////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
//
// set init gaussian proposal for the sampling
/// 
// [[Rcpp::export]]
NumericVector init_gaussian_proposal_cpp(const Rcpp::NumericVector& para_min, const Rcpp::NumericVector& para_max, double init_gp_ratio) {
  NumericVector ptr_sigma_prop(15);
  ptr_sigma_prop[0] = (para_max[0] - para_min[0]) * init_gp_ratio;
  ptr_sigma_prop[1] = (para_max[1] - para_min[1]) * init_gp_ratio;
  ptr_sigma_prop[2] = (para_max[2] - para_min[2]) * init_gp_ratio;
  ptr_sigma_prop[3] = (para_max[3] - para_min[3]) * init_gp_ratio;
  ptr_sigma_prop[4] = (para_max[4] - para_min[4]) * init_gp_ratio;
  ptr_sigma_prop[5] = (para_max[5] - para_min[5]) * init_gp_ratio;
  ptr_sigma_prop[6] = (para_max[6] - para_min[6]) * init_gp_ratio;
  ptr_sigma_prop[7] = (para_max[7] - para_min[7]) * init_gp_ratio;
  ptr_sigma_prop[8] = (para_max[8] - para_min[8]) * init_gp_ratio;
  ptr_sigma_prop[9] = (para_max[9] - para_min[9]) * init_gp_ratio;
  ptr_sigma_prop[10] = (para_max[10] - para_min[10]) * init_gp_ratio;
  ptr_sigma_prop[11] = (para_max[11] - para_min[11]) * init_gp_ratio;
  ptr_sigma_prop[12] = (para_max[12] - para_min[12]) * init_gp_ratio;
  ptr_sigma_prop[13] = (para_max[13] - para_min[13]) * init_gp_ratio;
  ptr_sigma_prop[14] = (para_max[14] - para_min[14]) * init_gp_ratio;
  return ptr_sigma_prop;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 6: check if a proposed point is within the range
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// TODO: add a new function for this ugly long boring function.
//
/// 
// [[Rcpp::export]]
double bounce_inside_cpp(double para, double para_min, double para_max) {
  if (para > para_max) {
    para = para_max - (para - para_max);
    para = std::max(para, para_min);
  }
  if (para < para_min) {
    para = para_min + (para_min - para);
    para = std::min(para, para_max);
  }
  return para;
}

/// 
// [[Rcpp::export]]
std::vector<double> para_boundary_cpp(std::vector<double> one_chain_new_cpp, // Pass by value
                                     const std::vector<double>& para_min_cpp,
                                     const std::vector<double>& para_max_cpp
                                     ) { 
  one_chain_new_cpp[0] = bounce_inside_cpp(one_chain_new_cpp[0], para_min_cpp[0], para_max_cpp[0]);
  one_chain_new_cpp[1] = bounce_inside_cpp(one_chain_new_cpp[1], para_min_cpp[1], para_max_cpp[1]);
  one_chain_new_cpp[2] = bounce_inside_cpp(one_chain_new_cpp[2], para_min_cpp[2], para_max_cpp[2]);
  one_chain_new_cpp[3] = bounce_inside_cpp(one_chain_new_cpp[3], para_min_cpp[3], para_max_cpp[3]);
  one_chain_new_cpp[4] = bounce_inside_cpp(one_chain_new_cpp[4], para_min_cpp[4], para_max_cpp[4]);
  one_chain_new_cpp[5] = bounce_inside_cpp(one_chain_new_cpp[5], para_min_cpp[5], para_max_cpp[5]);
  one_chain_new_cpp[6] = bounce_inside_cpp(one_chain_new_cpp[6], para_min_cpp[6], para_max_cpp[6]);
  one_chain_new_cpp[7] = bounce_inside_cpp(one_chain_new_cpp[7], para_min_cpp[7], para_max_cpp[7]);
  one_chain_new_cpp[8] = bounce_inside_cpp(one_chain_new_cpp[8], para_min_cpp[8], para_max_cpp[8]);
  one_chain_new_cpp[9] = bounce_inside_cpp(one_chain_new_cpp[9], para_min_cpp[9], para_max_cpp[9]);
  one_chain_new_cpp[10] = bounce_inside_cpp(one_chain_new_cpp[10], para_min_cpp[10], para_max_cpp[10]);
  one_chain_new_cpp[11] = bounce_inside_cpp(one_chain_new_cpp[11], para_min_cpp[11], para_max_cpp[11]);
  one_chain_new_cpp[12] = bounce_inside_cpp(one_chain_new_cpp[12], para_min_cpp[12], para_max_cpp[12]);
  one_chain_new_cpp[13] = bounce_inside_cpp(one_chain_new_cpp[13], para_min_cpp[13], para_max_cpp[13]);
  one_chain_new_cpp[14] = bounce_inside_cpp(one_chain_new_cpp[14], para_min_cpp[14], para_max_cpp[14]);
  return one_chain_new_cpp;
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////   Part 7: log prior
///////   NOTE: Write your own prior function of all parameters
///////
////////////////////////////////////////////
////////////////////////////////////////////
//
// NOTE: para_boundary ensures: ((a<=a_max) && (a>=a_min))
//
////////////////////////////////////////////
//
// in two-planet orbit retrieval:
// para0:  cos_i   p1
// para1:  ecc     p1
// para2:  anO     p1
// para3:  po      p1
// para4:  M0      p1
// para5:  mp      p1
// para6:  var_uke 
// para7:  period  p1
// para8:  cos_i   p2
// para9:  ecc     p2
// para10: anO     p2
// para11: po      p2
// para12: M0      p2
// para13: mp      p2
// para14: period  p2
//
/////////////////////////
/// 
// [[Rcpp::export]]
double prior_para0_cpp(double para0_min, double para0_max) {
  return 1.0 / (para0_max - para0_min);
}

/// 
// [[Rcpp::export]]
double prior_para1_cpp(double para1_min, double para1_max) {
  return 1.0 / (para1_max - para1_min);
}

/// 
// [[Rcpp::export]]
double prior_para2_cpp(double para2_min, double para2_max) {
  return 1.0 / (para2_max - para2_min);
}

/// 
// [[Rcpp::export]]
double prior_para3_cpp(double para3_min, double para3_max) {
  return 1.0 / (para3_max - para3_min);
}

/// 
// [[Rcpp::export]]
double prior_para4_cpp(double para4_min, double para4_max) {
  return 1.0 / (para4_max - para4_min);
}

/// 
// [[Rcpp::export]]
double prior_para5_cpp(double para5, double para5_min, double para5_max) {
  return 1.0 / (para5 * std::log(para5_max / para5_min));
}

/// 
// [[Rcpp::export]]
double prior_para6_cpp(double para6, double para6_min, double para6_max) {
  return 1.0 / (para6 * std::log(para6_max / para6_min));
}

/// 
// [[Rcpp::export]]
double prior_para7_cpp(double para7, double para7_min, double para7_max) {
  return 1.0 / (para7 * std::log(para7_max / para7_min));
}

/// 
// [[Rcpp::export]]
double prior_para8_cpp(double para8_min, double para8_max) {
  return 1.0 / (para8_max - para8_min);
}

/// 
// [[Rcpp::export]]
double prior_para9_cpp(double para9_min, double para9_max) {
  return 1.0 / (para9_max - para9_min);
}

/// 
// [[Rcpp::export]]
double prior_para10_cpp(double para10_min, double para10_max) {
  return 1.0 / (para10_max - para10_min);
}

/// 
// [[Rcpp::export]]
double prior_para11_cpp(double para11_min, double para11_max) {
  return 1.0 / (para11_max - para11_min);
}

/// 
// [[Rcpp::export]]
double prior_para12_cpp(double para12_min, double para12_max) {
  return 1.0 / (para12_max - para12_min);
}

/// 
// [[Rcpp::export]]
double prior_para13_cpp(double para13, double para13_min, double para13_max) {
  return 1.0 / (para13 * std::log(para13_max / para13_min));
}

/// 
// [[Rcpp::export]]
double prior_para14_cpp(double para14, double para14_min, double para14_max) {
  return 1.0 / (para14 * std::log(para14_max / para14_min));
}


//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////// Part 8: Combine all prior distributions 
///////
////////////////////////////////////////////
//
/// 
// [[Rcpp::export]]
double log_prior_cpp(const std::vector<double>& ptr_one_chain_cpp,
                    const std::vector<double>& para_min_cpp,
                    const std::vector<double>& para_max_cpp) {
  
  
  double log_prior = 0.0;

  log_prior += std::log(prior_para0_cpp(para_min_cpp[0], para_max_cpp[0]));
  log_prior += std::log(prior_para1_cpp(para_min_cpp[1], para_max_cpp[1]));
  log_prior += std::log(prior_para2_cpp(para_min_cpp[2], para_max_cpp[2]));
  log_prior += std::log(prior_para3_cpp(para_min_cpp[3], para_max_cpp[3]));
  log_prior += std::log(prior_para4_cpp(para_min_cpp[4], para_max_cpp[4]));
  log_prior += std::log(prior_para5_cpp(ptr_one_chain_cpp[5], para_min_cpp[5], para_max_cpp[5]));
  log_prior += std::log(prior_para6_cpp(ptr_one_chain_cpp[6], para_min_cpp[6], para_max_cpp[6]));
  log_prior += std::log(prior_para7_cpp(ptr_one_chain_cpp[7], para_min_cpp[7], para_max_cpp[7]));
  log_prior += std::log(prior_para8_cpp(para_min_cpp[8], para_max_cpp[8]));
  log_prior += std::log(prior_para9_cpp(para_min_cpp[9], para_max_cpp[9]));
  log_prior += std::log(prior_para10_cpp(para_min_cpp[10], para_max_cpp[10]));
  log_prior += std::log(prior_para11_cpp(para_min_cpp[11], para_max_cpp[11]));
  log_prior += std::log(prior_para12_cpp(para_min_cpp[12], para_max_cpp[12]));
  log_prior += std::log(prior_para13_cpp(ptr_one_chain_cpp[13], para_min_cpp[13], para_max_cpp[13]));
  log_prior += std::log(prior_para14_cpp(ptr_one_chain_cpp[14], para_min_cpp[14], para_max_cpp[14]));

  return log_prior;
  
}

////////////////////////////////////////////
////////////////////////////////////////////
/////  Save para debug
/////  NO need to change
////////////////////////////////////////////
////////////////////////////////////////////
//
//

/// 
// [[Rcpp::export]]
void save_debug_para_boundary_cpp(const std::string& p_s, double p, double p_min, double p_max, const std::string& result_dir) {
  std::string fname = result_dir + "/" + p_s + ".debug_para_boundary";
  std::ofstream con(fname, std::ios::app);
  if (con.is_open()) {
    con << std::scientific << p << " " << p_min << " " << p_max << "\n";
    con.close();
  } else {
    Rcpp::Rcerr << "Unable to open file: " << fname << "\n";
  }
}
