// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// init_sparse_matrix
SEXP init_sparse_matrix(const int nrow, const int ncol);
RcppExport SEXP _EigenSmallCountRounding_init_sparse_matrix(SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< const int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(init_sparse_matrix(nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// create_dummy_cpp
SEXP create_dummy_cpp(const Rcpp::CharacterMatrix groupings_inner, const Rcpp::CharacterMatrix groupings_publish);
RcppExport SEXP _EigenSmallCountRounding_create_dummy_cpp(SEXP groupings_innerSEXP, SEXP groupings_publishSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterMatrix >::type groupings_inner(groupings_innerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterMatrix >::type groupings_publish(groupings_publishSEXP);
    rcpp_result_gen = Rcpp::wrap(create_dummy_cpp(groupings_inner, groupings_publish));
    return rcpp_result_gen;
END_RCPP
}
// non_zeros_sparse_matrix
int non_zeros_sparse_matrix(SEXP mat_xptr);
RcppExport SEXP _EigenSmallCountRounding_non_zeros_sparse_matrix(SEXP mat_xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mat_xptr(mat_xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(non_zeros_sparse_matrix(mat_xptr));
    return rcpp_result_gen;
END_RCPP
}
// print_sparse_matrix
void print_sparse_matrix(Rcpp::XPtr<Eigen::SparseMatrix<double>> M);
RcppExport SEXP _EigenSmallCountRounding_print_sparse_matrix(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<Eigen::SparseMatrix<double>> >::type M(MSEXP);
    print_sparse_matrix(M);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EigenSmallCountRounding_init_sparse_matrix", (DL_FUNC) &_EigenSmallCountRounding_init_sparse_matrix, 2},
    {"_EigenSmallCountRounding_create_dummy_cpp", (DL_FUNC) &_EigenSmallCountRounding_create_dummy_cpp, 2},
    {"_EigenSmallCountRounding_non_zeros_sparse_matrix", (DL_FUNC) &_EigenSmallCountRounding_non_zeros_sparse_matrix, 1},
    {"_EigenSmallCountRounding_print_sparse_matrix", (DL_FUNC) &_EigenSmallCountRounding_print_sparse_matrix, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_EigenSmallCountRounding(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
