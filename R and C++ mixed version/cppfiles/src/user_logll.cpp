#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

/// 
// [[Rcpp::export]]
double period_to_au_cpp(double period_days, double ms_Msun) {
  const double PI = 3.14159265358;
  const double GG = 6.67259e-8;
  const double Msun = 1.9891e33;
  const double AU2cm = 1.4959787e13;

  double period_second = period_days * 24.0 * 3600.0;

  double a = period_second / (2.0 * PI);
  double b = GG * Msun * ms_Msun;
  double r = std::pow(a, 2.0) * b;
  double a_AU = std::pow(r, 1.0/3.0) / AU2cm;

  return a_AU;
}

/// 
// [[Rcpp::export]]
double calc_osi_cpp(double mp_Mearth, double ms_Msun, double a_AU, double d_pc) {
  double osi_alpha = 3.0 * mp_Mearth * std::pow(ms_Msun, -1.0) * a_AU * std::pow(d_pc, -1.0);
  return osi_alpha;
}

/// 
// [[Rcpp::export]]
double newton_solver_cpp(double EE_init, double e) {
  double tolerance = 1e-8;
  double eps = 1.0;
  double EE = EE_init;
  double M = EE_init;

  while (std::abs(eps) > tolerance) {
    double E1 = EE - (EE - e * std::sin(EE) - M) / (1.0 - e * std::cos(EE));
    eps = E1 - EE;
    EE = E1;
  }

  return EE;
}

/// 
// [[Rcpp::export]]
std::vector<std::vector<double>> func_as_cpp(const std::vector<double>& time_con, int Nline_time, double ecc, double osi, double cosi, double OmegaO, double M0, double omega, double per, std::vector<std::vector<double>>& Ntime_radec) {
  const double deg2rad = 0.0174532925;
  const double PI = 3.14159265358;

  omega *= deg2rad;
  M0 *= deg2rad;
  OmegaO *= deg2rad;

  double coso = std::cos(omega);
  double sino = std::sin(omega);
  double cosOg = std::cos(OmegaO);
  double sinOg = std::sin(OmegaO);

  double A = osi * (coso * cosOg - sino * sinOg * cosi);
  double B = osi * (coso * sinOg + sino * cosOg * cosi);
  double F = osi * (-sino * cosOg - coso * sinOg * cosi);
  double G = osi * (-sino * sinOg + coso * cosOg * cosi);

  for (int i = 0; i < Nline_time; ++i) {
    double EE = newton_solver_cpp(((2.0 * PI) / per * time_con[i] - M0), ecc);
    double X = std::cos(EE) - ecc;
    double Y = std::sqrt(1.0 - std::pow(ecc, 2.0)) * std::sin(EE);
    Ntime_radec[i][0] = B * X + G * Y; // C++ 索引从 0 开始
    Ntime_radec[i][1] = A * X + F * Y; // C++ 索引从 0 开始
  }

  return Ntime_radec;
}

/// 
// [[Rcpp::export]]
double log_likelihood_cpp(int Nline_time, const std::vector<std::vector<double>>& Ntime_radec, const std::vector<std::vector<double>>& data_radec, double var_uk) {
  const double PI2 = 6.28318530716;
  double log_llhd = 0.0;
  double sig_power = var_uk * var_uk;
  double AC_twice = 2.0 * std::log(std::pow(PI2, -0.5)) + 2.0 * std::log(std::pow(sig_power, -0.5));

  for (int i = 0; i < Nline_time; ++i) {
    double ra_once = -std::pow(data_radec[i][0] - Ntime_radec[i][0], 2.0) / (2.0 * sig_power); // C++ 索引从 0 开始
    double dec_once = -std::pow(data_radec[i][1] - Ntime_radec[i][1], 2.0) / (2.0 * sig_power); // C++ 索引从 0 开始
    log_llhd += ra_once + dec_once + AC_twice;
  }

  return log_llhd;
}

/// 
// [[Rcpp::export]]
double logll_beta_cpp(const std::vector<double>& ptr_one_chain_cpp, // Use C++ vector
                     int nline_data,
                     const std::vector<std::vector<double>>& data_cpp,    // Keep Rcpp - read only
                     int i_rank,
                     const std::vector<double>& beta_values_cpp) {
  
  
                     
  double ms_Msun = 1.0;
  double d_pc = 3.0;
  const int radec_int2 = 2;
  int ndim_data = 2;

  double cos_inc1 = ptr_one_chain_cpp[0]; // C++ 索引从 0 开始
  double ecc1 = ptr_one_chain_cpp[1];
  double an_Omg1 = ptr_one_chain_cpp[2];
  double p_omg1 = ptr_one_chain_cpp[3];
  double M0_1 = ptr_one_chain_cpp[4];
  double pl_m1 = ptr_one_chain_cpp[5];
  double var_uk = ptr_one_chain_cpp[6];
  double period1 = ptr_one_chain_cpp[7];
  double cos_inc2 = ptr_one_chain_cpp[8];
  double ecc2 = ptr_one_chain_cpp[9];
  double an_Omg2 = ptr_one_chain_cpp[10];
  double p_omg2 = ptr_one_chain_cpp[11];
  double M0_2 = ptr_one_chain_cpp[12];
  double pl_m2 = ptr_one_chain_cpp[13];
  double period2 = ptr_one_chain_cpp[14];

  double logll = 0.0;
  int Nline_time = nline_data / 2;
  std::vector<double> time_con(Nline_time);
  std::vector<std::vector<double>> data_radec(Nline_time, std::vector<double>(radec_int2));

  for (int i = 0; i < Nline_time; ++i) {
    time_con[i] = data_cpp[i*ndim_data][0]; // C++ 索引从 0 开始
    data_radec[i][0] = data_cpp[i*ndim_data][1];
    data_radec[i][1] = data_cpp[i*ndim_data+1][1];
  }

  double a_AU1 = period_to_au_cpp(period1, ms_Msun);
  double osi1 = calc_osi_cpp(pl_m1, ms_Msun, a_AU1, d_pc);

  double a_AU2 = period_to_au_cpp(period2, ms_Msun);
  double osi2 = calc_osi_cpp(pl_m2, ms_Msun, a_AU2, d_pc);

  std::vector<std::vector<double>> Ntime_radec1(Nline_time, std::vector<double>(radec_int2));
  std::vector<std::vector<double>> Ntime_radec2(Nline_time, std::vector<double>(radec_int2));

  Ntime_radec1 = func_as_cpp(time_con, Nline_time, ecc1, osi1, cos_inc1, an_Omg1, M0_1, p_omg1, period1, Ntime_radec1);
  Ntime_radec2 = func_as_cpp(time_con, Nline_time, ecc2, osi2, cos_inc2, an_Omg2, M0_2, p_omg2, period2, Ntime_radec2);

  std::vector<std::vector<double>> Ntime_radec(Nline_time, std::vector<double>(radec_int2));

  for (int i = 0; i < Nline_time; ++i) {
    Ntime_radec[i][0] = Ntime_radec1[i][0] + Ntime_radec2[i][0];
    Ntime_radec[i][1] = Ntime_radec1[i][1] + Ntime_radec2[i][1];
  }

  logll = log_likelihood_cpp(Nline_time, Ntime_radec, data_radec, var_uk);
  logll *= beta_values_cpp[i_rank]; // C++ 索引从 0 开始

  return logll;
  
}
