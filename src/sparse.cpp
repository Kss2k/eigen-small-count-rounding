// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SparseMatrix; // maybe use int?
typedef Rcpp::XPtr<Eigen::SparseMatrix<double>> RSparseMatrix;

// [[Rcpp::export]]
Rcpp::XPtr<Eigen::SparseMatrix<double>> init_sparse_matrix(const int nrow, const int ncol) {
  SparseMatrix *M = new SparseMatrix(nrow, ncol);

  RSparseMatrix out = Rcpp::XPtr<SparseMatrix>(M);

  return out;
}


// [[Rcpp::export]]
void free_sparse_matrix(RSparseMatrix *M) {
  delete M;
}


// [[Rcpp::export]]
void print_sparse_matrix(Rcpp::XPtr<Eigen::SparseMatrix<double>> M) {
  Rcpp::Rcout << "Sparse matrix:" << std::endl;
  Rcpp::Rcout << *M << std::endl;
}
