// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sim_epi_data
arma::vec sim_epi_data(double S0, double I0, double max_time, arma::vec beta_vec, double gamma_0, unsigned long user_seed);
RcppExport SEXP _BayesChange_sim_epi_data(SEXP S0SEXP, SEXP I0SEXP, SEXP max_timeSEXP, SEXP beta_vecSEXP, SEXP gamma_0SEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< double >::type I0(I0SEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_vec(beta_vecSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_epi_data(S0, I0, max_time, beta_vec, gamma_0, user_seed));
    return rcpp_result_gen;
END_RCPP
}
// psm
arma::mat psm(arma::mat M);
RcppExport SEXP _BayesChange_psm(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(psm(M));
    return rcpp_result_gen;
END_RCPP
}
// get_clust_VI
arma::vec get_clust_VI(arma::mat orders_mat);
RcppExport SEXP _BayesChange_get_clust_VI(SEXP orders_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type orders_mat(orders_matSEXP);
    rcpp_result_gen = Rcpp::wrap(get_clust_VI(orders_mat));
    return rcpp_result_gen;
END_RCPP
}
// detect_cp_uni
Rcpp::List detect_cp_uni(arma::vec data, int n_iterations, double q, double phi, double a, double b, double c, double par_theta_c, double par_theta_d, bool print_progress, unsigned long user_seed);
RcppExport SEXP _BayesChange_detect_cp_uni(SEXP dataSEXP, SEXP n_iterationsSEXP, SEXP qSEXP, SEXP phiSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP par_theta_cSEXP, SEXP par_theta_dSEXP, SEXP print_progressSEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type par_theta_c(par_theta_cSEXP);
    Rcpp::traits::input_parameter< double >::type par_theta_d(par_theta_dSEXP);
    Rcpp::traits::input_parameter< bool >::type print_progress(print_progressSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(detect_cp_uni(data, n_iterations, q, phi, a, b, c, par_theta_c, par_theta_d, print_progress, user_seed));
    return rcpp_result_gen;
END_RCPP
}
// detect_cp_multi
Rcpp::List detect_cp_multi(arma::mat data, int n_iterations, double q, double k_0, double nu_0, arma::mat phi_0, arma::vec m_0, double par_theta_c, double par_theta_d, double prior_var_gamma, bool print_progress, unsigned long user_seed);
RcppExport SEXP _BayesChange_detect_cp_multi(SEXP dataSEXP, SEXP n_iterationsSEXP, SEXP qSEXP, SEXP k_0SEXP, SEXP nu_0SEXP, SEXP phi_0SEXP, SEXP m_0SEXP, SEXP par_theta_cSEXP, SEXP par_theta_dSEXP, SEXP prior_var_gammaSEXP, SEXP print_progressSEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type k_0(k_0SEXP);
    Rcpp::traits::input_parameter< double >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_0(phi_0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< double >::type par_theta_c(par_theta_cSEXP);
    Rcpp::traits::input_parameter< double >::type par_theta_d(par_theta_dSEXP);
    Rcpp::traits::input_parameter< double >::type prior_var_gamma(prior_var_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type print_progress(print_progressSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(detect_cp_multi(data, n_iterations, q, k_0, nu_0, phi_0, m_0, par_theta_c, par_theta_d, prior_var_gamma, print_progress, user_seed));
    return rcpp_result_gen;
END_RCPP
}
// clust_cp_epi
Rcpp::List clust_cp_epi(arma::mat data, int n_iterations, int M, int B, int L, double gamma, double alpha, double q, double dt, double a0, double b0, double c0, double d0, double MH_var, double S0, double R0, double p, double coars, bool print_progress, unsigned long user_seed);
RcppExport SEXP _BayesChange_clust_cp_epi(SEXP dataSEXP, SEXP n_iterationsSEXP, SEXP MSEXP, SEXP BSEXP, SEXP LSEXP, SEXP gammaSEXP, SEXP alphaSEXP, SEXP qSEXP, SEXP dtSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP c0SEXP, SEXP d0SEXP, SEXP MH_varSEXP, SEXP S0SEXP, SEXP R0SEXP, SEXP pSEXP, SEXP coarsSEXP, SEXP print_progressSEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< double >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< double >::type MH_var(MH_varSEXP);
    Rcpp::traits::input_parameter< double >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type coars(coarsSEXP);
    Rcpp::traits::input_parameter< bool >::type print_progress(print_progressSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(clust_cp_epi(data, n_iterations, M, B, L, gamma, alpha, q, dt, a0, b0, c0, d0, MH_var, S0, R0, p, coars, print_progress, user_seed));
    return rcpp_result_gen;
END_RCPP
}
// clust_cp_uni
Rcpp::List clust_cp_uni(arma::mat data, int n_iterations, int B, int L, double gamma, double a, double b, double c, double q, double alpha_SM, double coars, bool print_progress, unsigned long user_seed);
RcppExport SEXP _BayesChange_clust_cp_uni(SEXP dataSEXP, SEXP n_iterationsSEXP, SEXP BSEXP, SEXP LSEXP, SEXP gammaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP qSEXP, SEXP alpha_SMSEXP, SEXP coarsSEXP, SEXP print_progressSEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_SM(alpha_SMSEXP);
    Rcpp::traits::input_parameter< double >::type coars(coarsSEXP);
    Rcpp::traits::input_parameter< bool >::type print_progress(print_progressSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(clust_cp_uni(data, n_iterations, B, L, gamma, a, b, c, q, alpha_SM, coars, print_progress, user_seed));
    return rcpp_result_gen;
END_RCPP
}
// clust_cp_multi
Rcpp::List clust_cp_multi(arma::cube data, int n_iterations, int B, int L, double gamma, double k_0, double nu_0, arma::mat phi_0, arma::vec m_0, double q, double alpha_SM, double coars, bool print_progress, unsigned long user_seed);
RcppExport SEXP _BayesChange_clust_cp_multi(SEXP dataSEXP, SEXP n_iterationsSEXP, SEXP BSEXP, SEXP LSEXP, SEXP gammaSEXP, SEXP k_0SEXP, SEXP nu_0SEXP, SEXP phi_0SEXP, SEXP m_0SEXP, SEXP qSEXP, SEXP alpha_SMSEXP, SEXP coarsSEXP, SEXP print_progressSEXP, SEXP user_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type k_0(k_0SEXP);
    Rcpp::traits::input_parameter< double >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_0(phi_0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_SM(alpha_SMSEXP);
    Rcpp::traits::input_parameter< double >::type coars(coarsSEXP);
    Rcpp::traits::input_parameter< bool >::type print_progress(print_progressSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type user_seed(user_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(clust_cp_multi(data, n_iterations, B, L, gamma, k_0, nu_0, phi_0, m_0, q, alpha_SM, coars, print_progress, user_seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesChange_sim_epi_data", (DL_FUNC) &_BayesChange_sim_epi_data, 6},
    {"_BayesChange_psm", (DL_FUNC) &_BayesChange_psm, 1},
    {"_BayesChange_get_clust_VI", (DL_FUNC) &_BayesChange_get_clust_VI, 1},
    {"_BayesChange_detect_cp_uni", (DL_FUNC) &_BayesChange_detect_cp_uni, 11},
    {"_BayesChange_detect_cp_multi", (DL_FUNC) &_BayesChange_detect_cp_multi, 12},
    {"_BayesChange_clust_cp_epi", (DL_FUNC) &_BayesChange_clust_cp_epi, 20},
    {"_BayesChange_clust_cp_uni", (DL_FUNC) &_BayesChange_clust_cp_uni, 13},
    {"_BayesChange_clust_cp_multi", (DL_FUNC) &_BayesChange_clust_cp_multi, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesChange(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
