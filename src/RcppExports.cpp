// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// predict_map_cpp
Rcpp::List predict_map_cpp(Eigen::MatrixXd data, Rcpp::List mcmc_sample, Eigen::MatrixXd dist_11, Eigen::MatrixXd dist_12, Eigen::MatrixXd dist_22, Rcpp::List params, int inner_reps, Rcpp::List args_progress, Rcpp::List args_functions, Rcpp::List args_misc);
RcppExport SEXP _genescaper_predict_map_cpp(SEXP dataSEXP, SEXP mcmc_sampleSEXP, SEXP dist_11SEXP, SEXP dist_12SEXP, SEXP dist_22SEXP, SEXP paramsSEXP, SEXP inner_repsSEXP, SEXP args_progressSEXP, SEXP args_functionsSEXP, SEXP args_miscSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mcmc_sample(mcmc_sampleSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type dist_11(dist_11SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type dist_12(dist_12SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type dist_22(dist_22SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type inner_reps(inner_repsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_progress(args_progressSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_misc(args_miscSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_map_cpp(data, mcmc_sample, dist_11, dist_12, dist_22, params, inner_reps, args_progress, args_functions, args_misc));
    return rcpp_result_gen;
END_RCPP
}
// wrangle_pairwise_Gst_cpp
Rcpp::List wrangle_pairwise_Gst_cpp(Rcpp::List freq_array);
RcppExport SEXP _genescaper_wrangle_pairwise_Gst_cpp(SEXP freq_arraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type freq_array(freq_arraySEXP);
    rcpp_result_gen = Rcpp::wrap(wrangle_pairwise_Gst_cpp(freq_array));
    return rcpp_result_gen;
END_RCPP
}
// sim_wrightfisher_cpp
Rcpp::List sim_wrightfisher_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
RcppExport SEXP _genescaper_sim_wrightfisher_cpp(SEXP argsSEXP, SEXP args_functionsSEXP, SEXP args_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_progress(args_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_wrightfisher_cpp(args, args_functions, args_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_genescaper_predict_map_cpp", (DL_FUNC) &_genescaper_predict_map_cpp, 10},
    {"_genescaper_wrangle_pairwise_Gst_cpp", (DL_FUNC) &_genescaper_wrangle_pairwise_Gst_cpp, 1},
    {"_genescaper_sim_wrightfisher_cpp", (DL_FUNC) &_genescaper_sim_wrightfisher_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_genescaper(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
