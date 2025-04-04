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
SEXP create_dummy_cpp(const Rcpp::CharacterMatrix groupings_inner, 
                  const Rcpp::CharacterMatrix groupings_publish) {
  const int 
    n = groupings_inner.nrow(),
    m = groupings_publish.nrow(),
    k = groupings_publish.ncol();
 
  // create from triplets
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(n * m);

  // ugly as f...
  for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) {
    int match = 1;

    for (int l = 0; l < k; l++) {
      if (groupings_inner(i, l) != groupings_publish(j, l) & 
          groupings_publish(j, l) != "__total__") {
        match = 0;
        break;
      };
    }
    
    if (match) triplet_list.push_back(Eigen::Triplet<double>(i, j, 1));
  }
 
  // create the sparse matrix
  SparseMatrix *X = new SparseMatrix(n, m);
  X->setFromTriplets(triplet_list.begin(), triplet_list.end());

  return Rcpp::XPtr<SparseMatrix>(X, true);
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
