// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// after
NumericVector after(NumericVector x);
RcppExport SEXP _mixvlmc_after(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(after(x));
    return rcpp_result_gen;
END_RCPP
}
// before
NumericVector before(NumericVector x);
RcppExport SEXP _mixvlmc_before(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(before(x));
    return rcpp_result_gen;
END_RCPP
}
// forward_match_all_ctx_counts
List forward_match_all_ctx_counts(IntegerVector x, int nb_vals, int depth, Nullable<IntegerVector> nv_from);
RcppExport SEXP _mixvlmc_forward_match_all_ctx_counts(SEXP xSEXP, SEXP nb_valsSEXP, SEXP depthSEXP, SEXP nv_fromSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nb_vals(nb_valsSEXP);
    Rcpp::traits::input_parameter< int >::type depth(depthSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type nv_from(nv_fromSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_match_all_ctx_counts(x, nb_vals, depth, nv_from));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mixvlmc_after", (DL_FUNC) &_mixvlmc_after, 1},
    {"_mixvlmc_before", (DL_FUNC) &_mixvlmc_before, 1},
    {"_mixvlmc_forward_match_all_ctx_counts", (DL_FUNC) &_mixvlmc_forward_match_all_ctx_counts, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mixvlmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
