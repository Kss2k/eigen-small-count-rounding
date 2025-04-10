// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <unordered_map>
#include <string>
#include <random>
#include <vector>
#include <cmath>


// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SparseMatrixD; // maybe use int?
typedef Rcpp::XPtr<Eigen::SparseMatrix<double>> RSparseMatrix;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::VectorXi VectorXi;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::MatrixXi MatrixXi;

#define TOTAL_C -9

// [[Rcpp::export]]
SEXP init_sparse_matrix(const int nrow, const int ncol) {
  SparseMatrixD *M = new SparseMatrixD(nrow, ncol);

  RSparseMatrix out = Rcpp::XPtr<SparseMatrixD>(M, true);

  return out;
}


// SEXP create_dummy_cpp(const Rcpp::CharacterMatrix groupings_inner, 
//                   const Rcpp::CharacterMatrix groupings_publish) {
//   const int 
//     n = groupings_inner.nrow(),
//     m = groupings_publish.nrow(),
//     k = groupings_publish.ncol();
//  
//   // create from triplets
//   std::vector<Eigen::Triplet<double>> triplet_list;
//   triplet_list.reserve(n * m);
// 
//   // ugly as f...
//   for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) {
//     int match = 1;
// 
//     for (int l = 0; l < k; l++) {
//       if (groupings_inner(i, l) != groupings_publish(j, l) & 
//           groupings_publish(j, l) != "__total__") {
//         match = 0;
//         break;
//       };
//     }
//     
//     if (match) triplet_list.push_back(Eigen::Triplet<double>(i, j, 1));
//   }
//  
//   // create the sparse matrix
//   SparseMatrixD *X = new SparseMatrixD(n, m);
//   X->setFromTriplets(triplet_list.begin(), triplet_list.end());
// 
//   return Rcpp::XPtr<SparseMatrixD>(X, true);

// [[Rcpp::export]]
SEXP create_dummy_cpp(const Rcpp::IntegerMatrix& groupings_inner, 
                      const Rcpp::IntegerMatrix& groupings_publish) 
{
  const int n = groupings_inner.nrow();
  const int k = groupings_inner.ncol();  // must match groupings_publish.ncol()
  const int m = groupings_publish.nrow();

  if (groupings_publish.ncol() != k) {
    Rcpp::stop("groupings_inner and groupings_publish must have the same number of columns");
  }


  std::vector<Eigen::Triplet<double>> triplet_list;
  // If many rows match, you might want to carefully guess a capacity,
  // but n*m might be huge. We'll reserve something smaller if you want.
  triplet_list.reserve(std::min((size_t)n*m, (size_t)1000000)); // e.g., 1e6

  for (int i = 0; i < n; i++) {
    Rcpp::IntegerMatrix::ConstRow rowI = groupings_inner(i, Rcpp::_);
    
    for (int j = 0; j < m; j++) {
      Rcpp::IntegerMatrix::ConstRow rowJ = groupings_publish(j, Rcpp::_);

      bool match = true;
      for (int col = 0; col < k; col++) {
        int pub_val = rowJ[col];
        // If pub_val != -1 and also != rowI[col], we have mismatch
        if (pub_val != TOTAL_C && pub_val != rowI[col]) {
          match = false;
          break;
        }
      }
      if (match) {
        triplet_list.emplace_back(i, j, 1.0);
      }
    }
  }

  SparseMatrixD *X = new SparseMatrixD(n, m);
  X->setFromTriplets(triplet_list.begin(), triplet_list.end());

  // Return as external pointer
  return Rcpp::XPtr<SparseMatrixD>(X, true);
}


// [[Rcpp::export]]
int non_zeros_sparse_matrix(SEXP mat_xptr) {
  Rcpp::XPtr<SparseMatrixD> xptr(mat_xptr);
  SparseMatrixD M = *xptr;

  return M.nonZeros();
}


// [[Rcpp::export]]
void print_sparse_matrix(SEXP mat_xptr) {
  Rcpp::XPtr<SparseMatrixD> xptr(mat_xptr);
  SparseMatrixD M = *xptr;

  Rcpp::Rcout << "Sparse matrix:" << std::endl;
  Rcpp::Rcout << M << std::endl;
}


// [[Rcpp::export]]
void print_VectorXd(SEXP xptr) {
  Rcpp::XPtr<VectorXd> ptr(xptr);
  Rcpp::Rcout << *ptr << std::endl;
}


// [[Rcpp::export]]
void print_VectorXi(SEXP xptr) {
  Rcpp::XPtr<VectorXi> ptr(xptr);
  Rcpp::Rcout << *ptr << std::endl;
}


// [[Rcpp::export]]
SEXP calc_z(SEXP xptr_M, SEXP xptr_y) {
  Rcpp::XPtr<SparseMatrixD> ptr_M(xptr_M);
  Rcpp::XPtr<VectorXd> ptr_y(xptr_y);
  SparseMatrixD M = *ptr_M;
  VectorXd y = *ptr_y; 

  // z will have length m
  VectorXd z = M.transpose() * y;  // (m x n) * (n x 1) = (m x 1)
  
  return Rcpp::XPtr<VectorXd>(new VectorXd(z), true);
}


// [[Rcpp::export]]
SEXP copy_xptr_VectorXd(SEXP xptr) {
  Rcpp::XPtr<VectorXd> ptr(xptr);
  VectorXd *copy = new VectorXd(*ptr);

  return Rcpp::XPtr<VectorXd>(copy, true);
}


// [[Rcpp::export]]
SEXP calc_M(SEXP xptr_X, double b) {
  Rcpp::XPtr<SparseMatrixD> X(xptr_X);
  SparseMatrixD *M = new SparseMatrixD(b * (*X).transpose() * (*X));
  M->makeCompressed();

  return Rcpp::XPtr<SparseMatrixD>(M, true);
}


// [[Rcpp::export]]
Rcpp::NumericVector calc_c(SEXP mat_xptr, const Rcpp::NumericVector &z_i) {
  Rcpp::XPtr<SparseMatrixD> xptr(mat_xptr);
  SparseMatrixD M = *xptr;

  // Map y_i to an Eigen vector
  Eigen::Map<const VectorXd> Z(z_i.begin(), z_i.size());

  // z will have length m
  VectorXd c = M * Z;  // (m x n) * (n x 1) = (m x 1)

  return Rcpp::NumericVector(c.data(), c.data() + c.size());
}


// [[Rcpp::export]]
Rcpp::NumericVector as_numeric_vector(SEXP vec_xptr) {
  Rcpp::XPtr<VectorXd> xptr(vec_xptr);
  VectorXd vec = *xptr;
  return Rcpp::NumericVector(vec.data(), vec.data() + vec.size());
}


// [[Rcpp::export]]
SEXP as_xptr_vector(const Rcpp::NumericVector &vec) {
  VectorXd *v = new VectorXd(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    (*v)(i) = vec[i];
  }
  return Rcpp::XPtr<VectorXd>(v, true);
}


SparseMatrixD get_subsetter_matrix(const VectorXi &x) {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(x.sum());

  int j = 0;
  for (int i = 0; i < x.size(); i++) {
    if (x(i)) {
      triplets.emplace_back(i, j, 1.0);
      j++;
    }
  }

  // Create the sparse Matrix
  SparseMatrixD M(x.size(), j);
  M.setFromTriplets(triplets.begin(), triplets.end());

  return M;
}


// [[Rcpp::export]]
Rcpp::List reduce_X_y_cpp(SEXP Xptr,
                          SEXP y_i_ptr,
                          double b,
                          SEXP z_ptr)
{
  // 1) Convert Xptr -> X
  Rcpp::XPtr<SparseMatrixD> X(Xptr);
  const int n = X->rows();
  const int m = X->cols();

  // 2) Convert y_i_ptr -> an Eigen vector
  Rcpp::XPtr<VectorXd> Yptr(y_i_ptr);
  if (Yptr->size() != n) {
    Rcpp::stop("Size of y_i must match # of rows in X.");
  }

  // 3) Convert z_ptr -> an Eigen vector
  Rcpp::XPtr<VectorXd> Zptr(z_ptr);
  if (Zptr->size() != m) {
    Rcpp::stop("Size of z must match # of cols in X.");
  }

  VectorXd z_i = X->transpose() * (*Yptr);
  // Eigen::ArrayXd = z_i.array();
  if (!(z_i.array() < b).any()) {
    return Rcpp::List::create(
        Rcpp::Named("X_i") = R_NilValue,
        Rcpp::Named("y_i") = R_NilValue
        );
  }


  VectorXi z_less_b = (z_i.array() < b).cast<int>();
  SparseMatrixD X_sub_cols = *X * get_subsetter_matrix(z_less_b);

  VectorXi inner = VectorXi::Zero(n);
  for (int j = 0; j < X_sub_cols.outerSize(); j++) {
    for (SparseMatrixD::InnerIterator it(X_sub_cols, j); it; ++it) {
      // double val = it.value();
      int i = it.row();
      inner(i) = 1;
    }
  }

  SparseMatrixD inner_row_subsetter = get_subsetter_matrix(inner);
  SparseMatrixD X_sub_rows = inner_row_subsetter.transpose() * (*X);

  const int n_sub_rows = X_sub_rows.rows();
  VectorXi non_zero_var = VectorXi::Ones(m);
  for (int j = 0; j < X_sub_rows.outerSize(); j++) {
    int n_non_zero = 0;
    for (SparseMatrixD::InnerIterator it(X_sub_rows, j); it; ++it)
      n_non_zero++;
     
    if (n_non_zero == 0 || n_non_zero == n_sub_rows)
      non_zero_var(j) = 0;
  }


  SparseMatrixD non_v0_col_subsetter = get_subsetter_matrix(non_zero_var);

  SparseMatrixD X_i = inner_row_subsetter.transpose() * (*X) * non_v0_col_subsetter;
  VectorXd y_i = inner_row_subsetter.transpose() * (*Yptr);
  VectorXd z_e = non_v0_col_subsetter.transpose() * (*Zptr - z_i);


  return Rcpp::List::create(
    Rcpp::Named("X_i") = 
        Rcpp::XPtr<SparseMatrixD>(new SparseMatrixD(X_i), true),
    Rcpp::Named("y_i") = 
        Rcpp::XPtr<VectorXd>(new VectorXd(y_i), true),
    Rcpp::Named("mask_y") = 
        Rcpp::XPtr<VectorXi>(new VectorXi(inner), true),
    // Rcpp::Named("mask_X") = 
    //     Rcpp::XPtr<SparseMatrixD>(new SparseMatrixD(non_v0_col_subsetter), true),
    Rcpp::Named("z_e") = 
        Rcpp::XPtr<VectorXd>(new VectorXd(z_e), true)
  );
}


// [[Rcpp::export]]
void fill_vector_by_mask(SEXP xptr_x, SEXP xptr_y, SEXP xptr_mask) {
  Rcpp::XPtr<VectorXd> ptr_x(xptr_x);
  Rcpp::XPtr<VectorXd> ptr_y(xptr_y);
  Rcpp::XPtr<VectorXi> ptr_mask(xptr_mask);

  if (ptr_x->size() != ptr_mask->size())
    Rcpp::stop("Size of x and mask must match!");

  int j = 0;
  for (int i = 0; i < ptr_x->size(); i++)
    if ((*ptr_mask)(i)) (*ptr_x)(i) = (*ptr_y)(j++);
}


// [[Rcpp::export]]
int calc_n_b(SEXP xptr_y_i, SEXP xptr_y, SEXP xptr_y_rounded, int b) {
  Rcpp::XPtr<VectorXd> ptr_y_i(xptr_y_i);
  Rcpp::XPtr<VectorXd> ptr_y(xptr_y);
  Rcpp::XPtr<VectorXd> ptr_y_rounded(xptr_y_rounded);

  double s_y = ptr_y_i->sum();
  double s_e = (*ptr_y - *ptr_y_rounded).sum();
  int n_b = (int)std::round((s_y + s_e) / b);

  return n_b; 
}


/**
 * round_cells_cpp:
 *   - Xptr      : XPtr<SparseMatrix<double>>
 *   - y_ptr     : XPtr<VectorXd> (length n)
 *   - b         : double
 *   - n_b       : int
 *   - max_iter  : int
 *   - z_e_ptr   : XPtr<VectorXd> (length m)
 *   - seed      : int
 *
 * Always internally computes:
 *   z   = X^T * y + z_e
 *   y_i = 0-vector of length n
 *   M   = b * (X * X^T)  [dense n x n]
 *
 * Returns: XPtr<VectorXd> pointing to the final y_i (length n).
 */
// [[Rcpp::export]]
SEXP round_cells_cpp(SEXP Xptr,
                     SEXP y_ptr,
                     double b,
                     int n_b,
                     int max_iter,
                     SEXP z_e_ptr,
                     int seed)
{
  // 1) Get the sparse matrix X
  Rcpp::XPtr<SparseMatrixD> Xp(Xptr);
  SparseMatrixD &X = *Xp; 
  const int n = X.rows();
  const int m = X.cols();

  // 2) Get y
  Rcpp::XPtr<VectorXd> Yp(y_ptr);
  if (Yp->size() != n) {
    Rcpp::stop("Length of y must match # of rows in X.");
  }

  // 3) Get z_e
  Rcpp::XPtr<VectorXd> Zep(z_e_ptr);
  if (Zep->size() != m) {
    Rcpp::stop("Length of z_e must match # of cols in X.");
  }

  // 4) Compute z = X^T y + z_e
  VectorXd z = X.transpose() * (*Yp); // length m
  z += (*Zep); // now z is X^T y + z_e

  // 5) Create y_i (length n) = all 0
  VectorXd *Y_i = new VectorXd(n);
  Y_i->setZero();  // fill with 0

  // 6) Compute M = b * (X * X^T) as a dense matrix (n x n)
  MatrixXd M_mat = (X * X.transpose()).toDense();
  M_mat *= b;

  // 7) c_i = X * z  (length n)
  VectorXd c_i = X * z;

  // 8) RNG
  std::mt19937 rng(seed);

  // -------------------------------------------------------------------------
  // PART A: for (i in seq_len(n_b))
  // -------------------------------------------------------------------------
  int i = 0;
  for (; i < n_b; i++) {
    // find max c_i[r] among those with Y_i[r] == 0
    double max_val = -std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] == 0.0 && c_i[r] > max_val) {
        max_val = c_i[r];
      }
    }
    // gather all r that achieve this max_val
    std::vector<int> candidates;
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] == 0.0 && std::fabs(c_i[r] - max_val) < 1e-15) {
        candidates.push_back(r);
      }
    }
    if (candidates.empty()) {
      // no row with Y_i[r]==0 => break
      break;
    }

    // random tie-break
    int chosen;
    if (candidates.size() == 1) {
      chosen = candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
      chosen = candidates[dist(rng)];
    }

    // set Y_i[chosen] = b
    (*Y_i)[chosen] = b;

    // c_i = c_i - M[chosen, ]
    for (int r = 0; r < n; r++) {
      c_i[r] -= M_mat(chosen, r);
    }
  }

  // -------------------------------------------------------------------------
  // PART B: while (i < max_iter)
  // -------------------------------------------------------------------------
  double last_k_max = -std::numeric_limits<double>::infinity();
  // double last_k_min =  std::numeric_limits<double>::infinity();  // R code had it but never used it
  VectorXd last_y_i = *Y_i;  // copy current y_i

  while (i < max_iter) {
    i++;

    // k_min = min c_i[r] among r with Y_i[r] != 0
    double min_val =  std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] != 0.0 && c_i[r] < min_val) {
        min_val = c_i[r];
      }
    }
    // gather all r that achieve min_val
    std::vector<int> min_candidates;
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] != 0.0 && std::fabs(c_i[r] - min_val) < 1e-15) {
        min_candidates.push_back(r);
      }
    }
    if (min_candidates.empty()) {
      break; 
    }

    int chosen_min;
    if (min_candidates.size() == 1) {
      chosen_min = min_candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, min_candidates.size() - 1);
      chosen_min = min_candidates[dist(rng)];
    }

    // c_i += M[chosen_min, ]
    for (int r = 0; r < n; r++) {
      c_i[r] += M_mat(chosen_min, r);
    }

    // k_max = max c_i[r] among r with Y_i[r] == 0
    double max_val = -std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] == 0.0 && c_i[r] > max_val) {
        max_val = c_i[r];
      }
    }
    // gather all r that achieve max_val
    std::vector<int> max_candidates;
    for (int r = 0; r < n; r++) {
      if ((*Y_i)[r] == 0.0 && std::fabs(c_i[r] - max_val) < 1e-15) {
        max_candidates.push_back(r);
      }
    }
    if (max_candidates.empty()) {
      break;
    }

    int chosen_max;
    if (max_candidates.size() == 1) {
      chosen_max = max_candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, max_candidates.size() - 1);
      chosen_max = max_candidates[dist(rng)];
    }

    // if (k_min >= k_max) => break
    if (min_val >= max_val) {
      break;
    } else if (max_val <= last_k_max) {
      // y_i <- last_y_i
      *Y_i = last_y_i;
      break;
    }

    // c_i -= M[chosen_max, ]
    for (int r = 0; r < n; r++) {
      c_i[r] -= M_mat(chosen_max, r);
    }

    // y_i[chosen_min] = 0;  y_i[chosen_max] = b
    (*Y_i)[chosen_min] = 0.0;
    (*Y_i)[chosen_max] = b;

    last_k_max = max_val;
    last_y_i = *Y_i;
  }

  // Return the final y_i as an external pointer
  return Rcpp::XPtr<VectorXd>(Y_i, true);
}
