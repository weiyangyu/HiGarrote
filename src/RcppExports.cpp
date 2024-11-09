// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_level
Rcpp::List get_level(Rcpp::DataFrame D, int p);
RcppExport SEXP _HiGarrote_get_level(SEXP DSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(get_level(D, p));
    return rcpp_result_gen;
END_RCPP
}
// contr_scale
Rcpp::NumericMatrix contr_scale(Rcpp::NumericMatrix x, int level_num);
RcppExport SEXP _HiGarrote_contr_scale(SEXP xSEXP, SEXP level_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type level_num(level_numSEXP);
    rcpp_result_gen = Rcpp::wrap(contr_scale(x, level_num));
    return rcpp_result_gen;
END_RCPP
}
// D_contr
Rcpp::DataFrame D_contr(Rcpp::DataFrame D, int p, Rcpp::IntegerVector mi, Rcpp::IntegerVector me_num, SEXP quali_id, SEXP quanti_eq_id);
RcppExport SEXP _HiGarrote_D_contr(SEXP DSEXP, SEXP pSEXP, SEXP miSEXP, SEXP me_numSEXP, SEXP quali_idSEXP, SEXP quanti_eq_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type mi(miSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type me_num(me_numSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quali_id(quali_idSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quanti_eq_id(quanti_eq_idSEXP);
    rcpp_result_gen = Rcpp::wrap(D_contr(D, p, mi, me_num, quali_id, quanti_eq_id));
    return rcpp_result_gen;
END_RCPP
}
// U_j_cpp
Rcpp::List U_j_cpp(Rcpp::List uni_level, int p, Rcpp::IntegerVector mi, SEXP quali_id, SEXP quanti_eq_id, SEXP quanti_ineq_id, SEXP quali_contr);
RcppExport SEXP _HiGarrote_U_j_cpp(SEXP uni_levelSEXP, SEXP pSEXP, SEXP miSEXP, SEXP quali_idSEXP, SEXP quanti_eq_idSEXP, SEXP quanti_ineq_idSEXP, SEXP quali_contrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type uni_level(uni_levelSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type mi(miSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quali_id(quali_idSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quanti_eq_id(quanti_eq_idSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quanti_ineq_id(quanti_ineq_idSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quali_contr(quali_contrSEXP);
    rcpp_result_gen = Rcpp::wrap(U_j_cpp(uni_level, p, mi, quali_id, quanti_eq_id, quanti_ineq_id, quali_contr));
    return rcpp_result_gen;
END_RCPP
}
// h_dist_cpp
Rcpp::List h_dist_cpp(Rcpp::NumericVector x, Rcpp::NumericMatrix m, bool two_level, bool quali);
RcppExport SEXP _HiGarrote_h_dist_cpp(SEXP xSEXP, SEXP mSEXP, SEXP two_levelSEXP, SEXP qualiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< bool >::type two_level(two_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type quali(qualiSEXP);
    rcpp_result_gen = Rcpp::wrap(h_dist_cpp(x, m, two_level, quali));
    return rcpp_result_gen;
END_RCPP
}
// h_j_cpp
Rcpp::List h_j_cpp(int p, Rcpp::List uni_level, Rcpp::List m_list, SEXP two_level_id, SEXP quali_id);
RcppExport SEXP _HiGarrote_h_j_cpp(SEXP pSEXP, SEXP uni_levelSEXP, SEXP m_listSEXP, SEXP two_level_idSEXP, SEXP quali_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type uni_level(uni_levelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type m_list(m_listSEXP);
    Rcpp::traits::input_parameter< SEXP >::type two_level_id(two_level_idSEXP);
    Rcpp::traits::input_parameter< SEXP >::type quali_id(quali_idSEXP);
    rcpp_result_gen = Rcpp::wrap(h_j_cpp(p, uni_level, m_list, two_level_id, quali_id));
    return rcpp_result_gen;
END_RCPP
}
// Psi_mat_cpp
arma::mat Psi_mat_cpp(const std::vector<arma::mat>& h_list_mat, const Rcpp::NumericVector& rho);
RcppExport SEXP _HiGarrote_Psi_mat_cpp(SEXP h_list_matSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type h_list_mat(h_list_matSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Psi_mat_cpp(h_list_mat, rho));
    return rcpp_result_gen;
END_RCPP
}
// initialize_NLLH_instance
void initialize_NLLH_instance(Rcpp::List h_list_mat, int n, int replicate, Rcpp::NumericVector y, double nugget, double epsilon, bool interpolate);
RcppExport SEXP _HiGarrote_initialize_NLLH_instance(SEXP h_list_matSEXP, SEXP nSEXP, SEXP replicateSEXP, SEXP ySEXP, SEXP nuggetSEXP, SEXP epsilonSEXP, SEXP interpolateSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type h_list_mat(h_list_matSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type replicate(replicateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type interpolate(interpolateSEXP);
    initialize_NLLH_instance(h_list_mat, n, replicate, y, nugget, epsilon, interpolate);
    return R_NilValue;
END_RCPP
}
// nllh_cpp_R
Rcpp::List nllh_cpp_R(Rcpp::NumericVector rho_lambda);
RcppExport SEXP _HiGarrote_nllh_cpp_R(SEXP rho_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho_lambda(rho_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(nllh_cpp_R(rho_lambda));
    return rcpp_result_gen;
END_RCPP
}
// nllh_GP_R
Rcpp::List nllh_GP_R(Rcpp::NumericVector rho);
RcppExport SEXP _HiGarrote_nllh_GP_R(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(nllh_GP_R(rho));
    return rcpp_result_gen;
END_RCPP
}
// rho_lambda_optim
Rcpp::List rho_lambda_optim(const Rcpp::NumericMatrix& ini_point, const Rcpp::List& h_list_mat, int n, int replicate, Rcpp::NumericVector y, double lambda_lb, double lambda_ub, double nugget, double epsilon, bool interpolate);
RcppExport SEXP _HiGarrote_rho_lambda_optim(SEXP ini_pointSEXP, SEXP h_list_matSEXP, SEXP nSEXP, SEXP replicateSEXP, SEXP ySEXP, SEXP lambda_lbSEXP, SEXP lambda_ubSEXP, SEXP nuggetSEXP, SEXP epsilonSEXP, SEXP interpolateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type ini_point(ini_pointSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_list_mat(h_list_matSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type replicate(replicateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lambda_lb(lambda_lbSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_ub(lambda_ubSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type interpolate(interpolateSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_lambda_optim(ini_point, h_list_mat, n, replicate, y, lambda_lb, lambda_ub, nugget, epsilon, interpolate));
    return rcpp_result_gen;
END_RCPP
}
// rho_optim_GP
Rcpp::List rho_optim_GP(const Rcpp::NumericMatrix& ini_point, const Rcpp::List& h_list_mat, int n, int replicate, Rcpp::NumericVector y, double nugget, double epsilon, bool interpolate);
RcppExport SEXP _HiGarrote_rho_optim_GP(SEXP ini_pointSEXP, SEXP h_list_matSEXP, SEXP nSEXP, SEXP replicateSEXP, SEXP ySEXP, SEXP nuggetSEXP, SEXP epsilonSEXP, SEXP interpolateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type ini_point(ini_pointSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_list_mat(h_list_matSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type replicate(replicateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type interpolate(interpolateSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_optim_GP(ini_point, h_list_mat, n, replicate, y, nugget, epsilon, interpolate));
    return rcpp_result_gen;
END_RCPP
}
// initialize_BETA_instance
void initialize_BETA_instance(Rcpp::List h_j_list, int p, Rcpp::List rho_list, Rcpp::IntegerVector mi);
RcppExport SEXP _HiGarrote_initialize_BETA_instance(SEXP h_j_listSEXP, SEXP pSEXP, SEXP rho_listSEXP, SEXP miSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type h_j_list(h_j_listSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rho_list(rho_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type mi(miSEXP);
    initialize_BETA_instance(h_j_list, p, rho_list, mi);
    return R_NilValue;
END_RCPP
}
// r_j_cpp_R
Rcpp::List r_j_cpp_R(Rcpp::List U_j_list, Rcpp::IntegerVector me_num);
RcppExport SEXP _HiGarrote_r_j_cpp_R(SEXP U_j_listSEXP, SEXP me_numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type U_j_list(U_j_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type me_num(me_numSEXP);
    rcpp_result_gen = Rcpp::wrap(r_j_cpp_R(U_j_list, me_num));
    return rcpp_result_gen;
END_RCPP
}
// beta_nng_cpp_R
Rcpp::NumericVector beta_nng_cpp_R(Rcpp::NumericMatrix U, Rcpp::NumericVector R, double lambda, int replicate, int n, Rcpp::NumericVector y, Rcpp::NumericMatrix Amat, double s2);
RcppExport SEXP _HiGarrote_beta_nng_cpp_R(SEXP USEXP, SEXP RSEXP, SEXP lambdaSEXP, SEXP replicateSEXP, SEXP nSEXP, SEXP ySEXP, SEXP AmatSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type U(USEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type replicate(replicateSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(beta_nng_cpp_R(U, R, lambda, replicate, n, y, Amat, s2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HiGarrote_get_level", (DL_FUNC) &_HiGarrote_get_level, 2},
    {"_HiGarrote_contr_scale", (DL_FUNC) &_HiGarrote_contr_scale, 2},
    {"_HiGarrote_D_contr", (DL_FUNC) &_HiGarrote_D_contr, 6},
    {"_HiGarrote_U_j_cpp", (DL_FUNC) &_HiGarrote_U_j_cpp, 7},
    {"_HiGarrote_h_dist_cpp", (DL_FUNC) &_HiGarrote_h_dist_cpp, 4},
    {"_HiGarrote_h_j_cpp", (DL_FUNC) &_HiGarrote_h_j_cpp, 5},
    {"_HiGarrote_Psi_mat_cpp", (DL_FUNC) &_HiGarrote_Psi_mat_cpp, 2},
    {"_HiGarrote_initialize_NLLH_instance", (DL_FUNC) &_HiGarrote_initialize_NLLH_instance, 7},
    {"_HiGarrote_nllh_cpp_R", (DL_FUNC) &_HiGarrote_nllh_cpp_R, 1},
    {"_HiGarrote_nllh_GP_R", (DL_FUNC) &_HiGarrote_nllh_GP_R, 1},
    {"_HiGarrote_rho_lambda_optim", (DL_FUNC) &_HiGarrote_rho_lambda_optim, 10},
    {"_HiGarrote_rho_optim_GP", (DL_FUNC) &_HiGarrote_rho_optim_GP, 8},
    {"_HiGarrote_initialize_BETA_instance", (DL_FUNC) &_HiGarrote_initialize_BETA_instance, 4},
    {"_HiGarrote_r_j_cpp_R", (DL_FUNC) &_HiGarrote_r_j_cpp_R, 2},
    {"_HiGarrote_beta_nng_cpp_R", (DL_FUNC) &_HiGarrote_beta_nng_cpp_R, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_HiGarrote(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
