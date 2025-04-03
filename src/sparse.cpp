// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SparseMatrix; // maybe use int?
typedef Rcpp::XPtr<Eigen::SparseMatrix<double>> RSparseMatrix;

// [[Rcpp::export]]
SEXP init_sparse_matrix(const int nrow, const int ncol) {
  SparseMatrix *M = new SparseMatrix(nrow, ncol);

  RSparseMatrix out = Rcpp::XPtr<SparseMatrix>(M, true);

  return out;
}


// [[Rcpp::export]]
int non_zeros_sparse_matrix(SEXP mat_xptr) {
  Rcpp::XPtr<Eigen::SparseMatrix<double>> xptr(mat_xptr);
  Eigen::SparseMatrix<double> M = *xptr;

  return M.nonZeros();
}


// [[Rcpp::export]]
void print_sparse_matrix(Rcpp::XPtr<Eigen::SparseMatrix<double>> M) {
  Rcpp::Rcout << "Sparse matrix:" << std::endl;
  Rcpp::Rcout << *M << std::endl;
}
