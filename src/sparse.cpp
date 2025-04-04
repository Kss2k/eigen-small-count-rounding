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
SEXP create_dummy_cpp(const Rcpp::CharacterMatrix& groupings_inner, 
                      const Rcpp::CharacterMatrix& groupings_publish) 
{
  const int n = groupings_inner.nrow();
  const int k = groupings_inner.ncol();  // must match groupings_publish.ncol()
  const int m = groupings_publish.nrow();

  if (groupings_publish.ncol() != k) {
    Rcpp::stop("groupings_inner and groupings_publish must have the same number of columns");
  }

  //-------------------------------------------------------------------------
  // 1) Collect unique strings from both matrices, assign them integer codes.
  //    We'll store "__total__" as a special code = -1.
  //-------------------------------------------------------------------------
  std::unordered_map<std::string,int> dict;
  dict.reserve(n * k + m * k); // just a rough reserve

  // We'll track the next code to assign:
  int nextCode = 0;

  auto getOrAssignCode = [&](const std::string &s) -> int {
    // Special rule for "__total__"
    if (s == "__total__") {
      return -1;  // wildcard code
    }
    // Otherwise, look up in dict
    auto it = dict.find(s);
    if (it == dict.end()) {
      // not found, assign new code
      int code = nextCode;
      dict[s] = code;
      nextCode++;
      return code;
    } else {
      return it->second;
    }
  };

  //-------------------------------------------------------------------------
  // 2) Build integer-coded matrices:
  //    inner_codes (n x k), publish_codes (m x k).
  //-------------------------------------------------------------------------
  // We'll store them in row-major order for convenience, but
  // an Eigen::MatrixXi with n rows, k cols is also possible.

  // Using a simple std::vector<int> for the 2D data:
  std::vector<int> inner_codes(n*k), publish_codes(m*k);

  // Fill inner_codes
  for (int i = 0; i < n; i++) {
    for (int col = 0; col < k; col++) {
      std::string s = Rcpp::as<std::string>(groupings_inner(i, col));
      inner_codes[i*k + col] = getOrAssignCode(s);
    }
  }

  // Fill publish_codes
  for (int j = 0; j < m; j++) {
    for (int col = 0; col < k; col++) {
      std::string s = Rcpp::as<std::string>(groupings_publish(j, col));
      publish_codes[j*k + col] = getOrAssignCode(s);
    }
  }

  //-------------------------------------------------------------------------
  // 3) Create triplets by comparing each row i in [0..n) with row j in [0..m).
  //    If for all col l in [0..k), either:
  //       publish_codes(j,l) == -1  (wildcard)
  //    or publish_codes(j,l) == inner_codes(i,l),
  //    => push triplet (i, j, 1).
  //-------------------------------------------------------------------------
  std::vector<Eigen::Triplet<double>> triplet_list;
  // If many rows match, you might want to carefully guess a capacity,
  // but n*m might be huge. We'll reserve something smaller if you want.
  triplet_list.reserve(std::min((size_t)n*m, (size_t)1000000)); // e.g., 1e6

  for (int i = 0; i < n; i++) {
    const int *rowI = &inner_codes[i*k];  // pointer to row i
    for (int j = 0; j < m; j++) {
      const int *rowJ = &publish_codes[j*k];  // pointer to row j

      bool match = true;
      for (int col = 0; col < k; col++) {
        int pub_val = rowJ[col];
        // If pub_val != -1 and also != rowI[col], we have mismatch
        if (pub_val != -1 && pub_val != rowI[col]) {
          match = false;
          break;
        }
      }
      if (match) {
        triplet_list.emplace_back(i, j, 1.0);
      }
    }
  }

  //-------------------------------------------------------------------------
  // 4) Build the (n x m) sparse matrix from the triplets and return as XPtr
  //-------------------------------------------------------------------------
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
Rcpp::NumericVector calc_z(SEXP mat_xptr, const Rcpp::NumericVector &y_i) {
  Rcpp::XPtr<SparseMatrixD> xptr(mat_xptr);
  SparseMatrixD M = *xptr;

  // Map y_i to an Eigen vector
  Eigen::Map<const Eigen::VectorXd> Y(y_i.begin(), y_i.size());

  // z will have length m
  Eigen::VectorXd z = M.transpose() * Y;  // (m x n) * (n x 1) = (m x 1)

  return Rcpp::NumericVector(z.data(), z.data() + z.size());
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
  Eigen::Map<const Eigen::VectorXd> Z(z_i.begin(), z_i.size());

  // z will have length m
  Eigen::VectorXd c = M * Z;  // (m x n) * (n x 1) = (m x 1)

  return Rcpp::NumericVector(c.data(), c.data() + c.size());
}


// [[Rcpp::export]]
Rcpp::NumericVector as_numeric_vector(SEXP vec_xptr) {
  Rcpp::XPtr<Eigen::VectorXd> xptr(vec_xptr);
  Eigen::VectorXd vec = *xptr;
  return Rcpp::NumericVector(vec.data(), vec.data() + vec.size());
}


// [[Rcpp::export]]
SEXP as_xptr_vector(const Rcpp::NumericVector &vec) {
  Eigen::VectorXd *v = new Eigen::VectorXd(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    (*v)(i) = vec[i];
  }
  return Rcpp::XPtr<Eigen::VectorXd>(v, true);
}


/**
 * reduce_X_y_cpp (Faster Version)
 *
 *  - Xptr      : External pointer to n x m SparseMatrix<double>
 *  - y_i_ptr   : External pointer to length-n VectorXd
 *  - b         : double threshold
 *  - z_ptr     : External pointer to length-m VectorXd
 *
 * Returns an R list with:
 *   - X_i   : XPtr<SparseMatrix<double>> (the reduced X)
 *   - y_i   : XPtr<VectorXd>            (the reduced y)
 *   - mask_y: R logical vector (length n)
 *   - mask_X: R logical vector (length m)
 *   - z_e   : XPtr<VectorXd>            (the subset of (z - z_i))
 */
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
  Rcpp::XPtr<Eigen::VectorXd> Yptr(y_i_ptr);
  if (Yptr->size() != n) {
    Rcpp::stop("Size of y_i must match # of rows in X.");
  }

  // 3) Convert z_ptr -> an Eigen vector
  Rcpp::XPtr<Eigen::VectorXd> Zptr(z_ptr);
  if (Zptr->size() != m) {
    Rcpp::stop("Size of z must match # of cols in X.");
  }

  // 4) Compute z_i = t(X) * y_i  (length = m)
  Eigen::VectorXd z_i = X->transpose() * (*Yptr);

  // 5) If !any(z_i < b), return (X=NULL, y=NULL)
  bool anyBelowB = false;
  for (int c = 0; c < m; ++c) {
    if (z_i[c] < b) {
      anyBelowB = true;
      break;
    }
  }
  if (!anyBelowB) {
    return Rcpp::List::create(
      Rcpp::Named("X_i") = R_NilValue,
      Rcpp::Named("y_i") = R_NilValue
    );
  }

  // -------------------------------------------------------------
  // Build colMaskB[c] = (z_i[c] < b).
  // Then rowMask[r] = true if row r has any positive entry
  // in columns c where colMaskB[c] = true.
  // -------------------------------------------------------------
  std::vector<bool> colMaskB(m, false);
  for (int c = 0; c < m; ++c) {
    colMaskB[c] = (z_i[c] < b);
  }

  std::vector<bool> rowMask(n, false);

  for (int c = 0; c < m; ++c) {
    if (!colMaskB[c]) continue;
    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      if (it.value() > 0.0) {
        rowMask[it.row()] = true;
      }
    }
  }

  // -------------------------------------------------------------
  // var0[c] = "all the same" among the rows where rowMask[r] is true.
  // We'll do a single pass over the column's nonzeros to discover
  // up to two distinct values (0 vs. nonzero).
  // -------------------------------------------------------------
  std::vector<bool> var0(m, false);

  // Precompute how many rows are "inner" => included
  int countIncluded = 0;
  for (int r = 0; r < n; ++r) {
    if (rowMask[r]) {
      countIncluded++;
    }
  }

  for (int c = 0; c < m; ++c) {
    // If 0 or 1 row is included, there's no variation
    if (countIncluded <= 1) {
      var0[c] = true;
      continue;
    }

    // We'll gather up to 2 distinct values among included rows
    // (0 is implicitly used for rows not in the column's nonzeros)
    std::vector<double> distinctVals;
    distinctVals.reserve(2);

    // Count how many included rows appear as nonzero
    int foundCount = 0;

    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      int r = it.row();
      if (!rowMask[r]) continue;  // skip rows not in "inner"
      double val = it.value();
      foundCount++;

      // Insert into distinctVals (up to 2):
      if (distinctVals.empty()) {
        distinctVals.push_back(val);
      } else if (distinctVals.size() == 1) {
        if (std::fabs(val - distinctVals[0]) > 1e-15) {
          // It's genuinely different
          distinctVals.push_back(val);
          // Now we have 2 distinct values => can break
          break;
        }
      } else {
        // distinctVals.size() == 2 => definitely variation
        break;
      }
    }

    // If we ended the loop with >2 distinct => variation
    // But we only allowed up to 2 in distinctVals, so check after
    if (distinctVals.size() > 1) {
      var0[c] = false;
      continue;
    }

    // Now see how many included rows we *didn't* see as nonzero
    int missing = countIncluded - foundCount;
    if (missing > 0) {
      // That means those rows have value 0
      if (distinctVals.empty()) {
        // Everything is 0
        // => all the same
        var0[c] = true;
      } else {
        // We have exactly 1 distinct value so far
        double stored = distinctVals[0];
        // If that value != 0 => we have 2 distinct => variation
        if (std::fabs(stored) > 1e-15) {
          var0[c] = false;
        } else {
          // stored == 0 => all 0 => still uniform
          var0[c] = true;
        }
      }
    } else {
      // missing=0 => all included rows were in the nonzeros
      // => we found exactly 0 or 1 distinct value
      //  - if 0 distinct => impossible, but let's handle gracefully
      //  - if 1 distinct => var0[c] = true
      if (distinctVals.empty()) {
        // means no included row was found => but missing=0 => contradictory
        // let's treat that as uniform => true
        var0[c] = true;
      } else {
        // exactly 1 distinct => uniform
        var0[c] = true;
      }
    }
  }

  // colMask[c] = !var0[c]
  std::vector<bool> colMask(m, false);
  for (int c = 0; c < m; ++c) {
    colMask[c] = !var0[c];
  }

  // -------------------------------------------------------------
  // Build the reduced matrix X_i = X[rowMask, colMask]
  // -------------------------------------------------------------
  int newN = 0;
  int newM = 0;
  std::vector<int> rowMap(n, -1), colMap(m, -1);

  {
    int idx = 0;
    for (int r = 0; r < n; ++r) {
      if (rowMask[r]) {
        rowMap[r] = idx++;
      }
    }
    newN = idx;
  }
  {
    int idx = 0;
    for (int c = 0; c < m; ++c) {
      if (colMask[c]) {
        colMap[c] = idx++;
      }
    }
    newM = idx;
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(X->nonZeros());

  for (int c = 0; c < m; ++c) {
    if (!colMask[c]) continue;
    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      int r = it.row();
      if (rowMask[r]) {
        triplets.emplace_back(rowMap[r], colMap[c], it.value());
      }
    }
  }

  SparseMatrixD *X_new = new SparseMatrixD(newN, newM);
  X_new->setFromTriplets(triplets.begin(), triplets.end());
  Rcpp::XPtr<SparseMatrixD> X_i_ptr(X_new, true);

  // -------------------------------------------------------------
  // Build new y_i (subset of original y_i)
  // -------------------------------------------------------------
  Eigen::VectorXd *Y_new = new Eigen::VectorXd(newN);
  {
    int idx = 0;
    for (int r = 0; r < n; ++r) {
      if (rowMask[r]) {
        (*Y_new)(idx++) = (*Yptr)(r);
      }
    }
  }
  Rcpp::XPtr<Eigen::VectorXd> Y_i_ptr(Y_new, true);

  // -------------------------------------------------------------
  // Build z_e = (z - z_i)[colMask]
  // -------------------------------------------------------------
  Eigen::VectorXd Z_diff = (*Zptr) - z_i; // length m
  Eigen::VectorXd *Z_e_new = new Eigen::VectorXd(newM);

  {
    int idx = 0;
    for (int c = 0; c < m; ++c) {
      if (colMask[c]) {
        (*Z_e_new)(idx++) = Z_diff[c];
      }
    }
  }
  Rcpp::XPtr<Eigen::VectorXd> Z_e_ptr(Z_e_new, true);

  // -------------------------------------------------------------
  // Build R logical vectors for rowMask, colMask
  // -------------------------------------------------------------
  Rcpp::LogicalVector mask_y(n), mask_X(m);
  for (int r = 0; r < n; ++r) {
    mask_y[r] = rowMask[r];
  }
  for (int c = 0; c < m; ++c) {
    mask_X[c] = colMask[c];
  }

  // -------------------------------------------------------------
  // Return as list
  // -------------------------------------------------------------
  return Rcpp::List::create(
    Rcpp::Named("X_i")    = X_i_ptr,
    Rcpp::Named("y_i")    = Y_i_ptr,
    Rcpp::Named("mask_y") = mask_y,
    Rcpp::Named("mask_X") = mask_X,
    Rcpp::Named("z_e")    = Z_e_ptr
  );
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
  Rcpp::XPtr<Eigen::VectorXd> Yp(y_ptr);
  if (Yp->size() != n) {
    Rcpp::stop("Length of y must match # of rows in X.");
  }

  // 3) Get z_e
  Rcpp::XPtr<Eigen::VectorXd> Zep(z_e_ptr);
  if (Zep->size() != m) {
    Rcpp::stop("Length of z_e must match # of cols in X.");
  }

  // 4) Compute z = X^T y + z_e
  Eigen::VectorXd z = X.transpose() * (*Yp); // length m
  z += (*Zep); // now z is X^T y + z_e

  // 5) Create y_i (length n) = all 0
  Eigen::VectorXd *Y_i = new Eigen::VectorXd(n);
  Y_i->setZero();  // fill with 0

  // 6) Compute M = b * (X * X^T) as a dense matrix (n x n)
  Eigen::MatrixXd M_mat = (X * X.transpose()).toDense();
  M_mat *= b;

  // 7) c_i = X * z  (length n)
  Eigen::VectorXd c_i = X * z;

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
  Eigen::VectorXd last_y_i = *Y_i;  // copy current y_i

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
  return Rcpp::XPtr<Eigen::VectorXd>(Y_i, true);
}
