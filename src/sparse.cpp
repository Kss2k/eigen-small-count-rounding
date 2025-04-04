// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <unordered_map>
#include <string>
#include <random>

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

/**
 * reduce_X_y_cpp:
 *   - Xptr: external pointer to an n x m sparse matrix (X).
 *   - y_i : length-n numeric vector.
 *   - b   : numeric threshold.
 *   - z   : length-m numeric vector.
 *
 * Returns a List with:
 *   - X_i   : new XPtr<SparseMatrix<double>> for reduced matrix
 *   - y_i   : reduced y_i
 *   - mask_y: logical row mask (the "inner" from R)
 *   - mask_X: logical column mask (the NOT of var0 from R)
 *   - z_e   : (z - z_i) subsetted by mask_X
 */
// [[Rcpp::export]]
Rcpp::List reduce_X_y_cpp(SEXP Xptr,
                          SEXP y_i,
                          double b,
                          SEXP &z)
{
  // 1) Convert Xptr -> X (a sparse matrix)
  Rcpp::XPtr<SparseMatrixD> X(Xptr);  
  const int n = X->rows(); // # of rows in X
  const int m = X->cols(); // # of cols in X

  // Safety checks (not strictly required but good practice):
  if (y_i.size() != n) {
    Rcpp::stop("Size of y_i must match # of rows in X.");
  }
  if (z.size() != m) {
    Rcpp::stop("Size of z must match # of cols in X.");
  }

  // 2) Compute z_i = t(X) * y_i
  //    X is n x m, so X^T is m x n. y_i is length n => result is length m
  // Map y_i to an Eigen vector
  // Eigen::Map<const Eigen::VectorXd> Y(y_i.begin(), n);
  Rcpp::XPtr<Eigen::VectorXd> y_i_ptr(y_i);
  Eigen::Vector Y = *y_i_ptr;

  // z_i will have length m
  Eigen::VectorXd z_i = X->transpose() * Y;  // (m x n) * (n x 1) = (m x 1)

  // 3) If !any(z_i < b), return list(X=NULL, y=NULL)
  bool anyBelowB = false;
  for (int col = 0; col < m; ++col) {
    if (z_i[col] < b) {
      anyBelowB = true;
      break;
    }
  }
  if (!anyBelowB) {
    // No column has z_i < b
    return Rcpp::List::create(
      Rcpp::Named("X_i") = R_NilValue,
      Rcpp::Named("y_i") = R_NilValue
      // If you want to return masks or z_e as NULL also, do so here
    );
  }

  // --------------------------------------------------------------------------
  // inner <- apply(X[, z_i < b], MARGIN=1, any(x>0))
  // 
  // We'll define colMaskB[c] = (z_i[c] < b).
  // Then rowMask[r] = true if row r has any positive entry in any col c
  // for which colMaskB[c] is true.
  // --------------------------------------------------------------------------
  std::vector<bool> colMaskB(m, false);
  for (int c = 0; c < m; ++c) {
    colMaskB[c] = (z_i[c] < b);
  }

  std::vector<bool> rowMask(n, false);

  // For each column c where colMaskB[c] is true, iterate over its non-zeros.
  // If a non-zero value > 0, mark rowMask[r] = true.
  for (int c = 0; c < m; ++c) {
    if (!colMaskB[c]) continue;
    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      if (it.value() > 0.0) {
        rowMask[it.row()] = true;
      }
    }
  }

  // --------------------------------------------------------------------------
  // var0 <- apply(X[inner, ], MARGIN=2, all(x == x[1]))
  //
  // For each column c, among the rows r where rowMask[r] is true ("inner"),
  // check if all the values are the same.
  // 
  // Because X is sparse, we have to remember: 
  //   - if a row r is not in the inner-iterator for column c, that value is 0.
  // --------------------------------------------------------------------------
  std::vector<bool> var0(m, false);

  for (int c = 0; c < m; ++c) {
    // We'll collect the values for all "inner" rows in this column.
    // If the row is not stored, that means the value is 0.
    // Then check if they are all the same.
    // 
    // If there are 0 or 1 'inner' rows, we consider that "all the same" => var0[c] = true.
    // (Because there's no variation with 0/1 points.)

    // Make a small map: row -> value
    // (We expect a column to have few nonzeros, so this is more efficient
    //  than scanning all rows.)
    std::unordered_map<int, double> colValues;
    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      colValues[it.row()] = it.value();
    }

    // Now iterate over all rows to find those where rowMask[r] = true
    double firstVal = 0.0;
    bool foundFirst = false;
    bool allSame = true;
    int countIncluded = 0;

    for (int r = 0; r < n; ++r) {
      if (!rowMask[r]) continue; // skip rows that are not "inner"
      countIncluded++;

      // value in X(r,c)
      double val = 0.0;
      auto search = colValues.find(r);
      if (search != colValues.end()) {
        val = search->second;
      }

      if (!foundFirst) {
        firstVal = val;
        foundFirst = true;
      } else {
        // Compare to firstVal
        if (val != firstVal) {
          allSame = false;
          break;
        }
      }
    }

    // If we have 0 or 1 included rows, or if we didn't break:
    //   => allSame remains true
    // => var0[c] = allSame
    // The R code sets var0=TRUE if "all x == x[1]"
    // so that means "no variation" => "allSame = true"
    if (countIncluded <= 1) {
      var0[c] = true;
    } else {
      var0[c] = allSame;
    }
  }

  // mask_X = !var0
  // We'll build that as well for clarity:
  std::vector<bool> colMask(m, true);
  for (int c = 0; c < m; ++c) {
    colMask[c] = !var0[c];
  }

  // --------------------------------------------------------------------------
  // Build the reduced matrix X_i = X[ rowMask, colMask ]
  // i.e. keep only rows where rowMask[r] is true and columns where colMask[c] is true.
  //
  // We'll do this by scanning all non-zero entries of X and collecting
  // triplets that match rowMask[r] & colMask[c].
  // --------------------------------------------------------------------------
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(X->nonZeros()); // upper bound

  // We'll also count how many rowMask & colMask are true (to size the new matrix).
  int newN = 0;
  int newM = 0;

  // Build an index mapping from old row -> new row index
  // and old col -> new col index
  std::vector<int> rowMap(n, -1);
  std::vector<int> colMap(m, -1);

  {
    int rCount = 0;
    for (int r = 0; r < n; ++r) {
      if (rowMask[r]) {
        rowMap[r] = rCount;
        rCount++;
      }
    }
    newN = rCount;
  }
  {
    int cCount = 0;
    for (int c = 0; c < m; ++c) {
      if (colMask[c]) {
        colMap[c] = cCount;
        cCount++;
      }
    }
    newM = cCount;
  }

  // Now iterate again
  for (int c = 0; c < m; ++c) {
    if (!colMask[c]) continue;
    for (SparseMatrixD::InnerIterator it(*X, c); it; ++it) {
      int r = it.row();
      if (!rowMask[r]) continue;
      // Keep it
      tripletList.emplace_back(rowMap[r], colMap[c], it.value());
    }
  }

  // Create a new sparse matrix from the triplets
  SparseMatrixD *X_new = new SparseMatrixD(newN, newM);
  X_new->setFromTriplets(tripletList.begin(), tripletList.end());

  // Wrap in an external pointer so R sees it as a "reference to a sparse matrix"
  Rcpp::XPtr<SparseMatrixD> X_i(X_new, true);

  // --------------------------------------------------------------------------
  // y_i <- y_i[ rowMask ]
  // We'll return that as a numeric vector
  // --------------------------------------------------------------------------
  Rcpp::NumericVector y_i_new(newN);
  {
    int idx = 0;
    for (int r = 0; r < n; ++r) {
      if (rowMask[r]) {
        y_i_new[idx] = y_i[r];
        idx++;
      }
    }
  }

  // --------------------------------------------------------------------------
  // z_e <- (z - z_i)[ colMask ]
  // i.e. a vector of length newM
  // --------------------------------------------------------------------------
  Rcpp::NumericVector z_e(newM);
  {
    int idx = 0;
    for (int c = 0; c < m; ++c) {
      if (colMask[c]) {
        z_e[idx] = (z[c] - z_i[c]);
        idx++;
      }
    }
  }

  // --------------------------------------------------------------------------
  // Also return the row/column masks:
  //    mask_y = inner (rowMask)
  //    mask_X = !var0 (colMask)
  //
  // We'll return them as logical vectors in R
  // --------------------------------------------------------------------------
  Rcpp::LogicalVector mask_y(n), mask_X(m);
  for (int r = 0; r < n; ++r) {
    mask_y[r] = rowMask[r];
  }
  for (int c = 0; c < m; ++c) {
    mask_X[c] = colMask[c];
  }

  // Finally, return everything as an R list
  return Rcpp::List::create(
    Rcpp::Named("X_i")    = X_i,
    Rcpp::Named("y_i")    = y_i_new,
    Rcpp::Named("mask_y") = mask_y,
    Rcpp::Named("mask_X") = mask_X,
    Rcpp::Named("z_e")    = z_e
  );
}


/*
 * round_cells_cpp
 * 
 * Translated from your R function:
 *
 *   round_cells <- function(X, y, z=NULL, y0=NULL, b, n_b, 
 *                           M=NULL, max.iter=1000, z_e, seed = 1234) { ... }
 *
 * Arguments:
 *  - Xptr     : External pointer to an (n x m) sparse matrix (X).
 *  - y        : NumericVector of length n.
 *  - z_       : (Optional) NumericVector of length m, if null => we compute z = t(X)*y + z_e
 *  - y0_      : (Optional) NumericVector of length n; if null => start with y_i=0
 *  - b        : double
 *  - n_b      : int
 *  - M_       : (Optional) NxN numeric matrix (if null => M = b * (X * X^T) ).
 *  - max_iter : maximum # of iterations in the second loop (default 1000)
 *  - z_e      : NumericVector of length m
 *  - seed     : int random seed
 *
 * Returns: A NumericVector (length n) = final y_i
 */

// [[Rcpp::export]]
Rcpp::NumericVector round_cells_cpp(
    SEXP Xptr, 
    const Rcpp::NumericVector &y, 
    double b, 
    int n_b,
    int max_iter ,
    const Rcpp::NumericVector &z_e,
    int seed,
    Rcpp::Nullable<Rcpp::NumericVector> z_=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> y0_=R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> M_=R_NilValue
)
{
  // 1) Get the sparse matrix X (n x m) from the external pointer
  Rcpp::XPtr<SparseMatrixD> Xp(Xptr);
  SparseMatrixD &X = *Xp; 
  const int n = X.rows();   // # rows
  const int m = X.cols();   // # cols

  // Basic checks
  if (y.size() != (size_t)n) {
    Rcpp::stop("Length of y must match the number of rows in X.");
  }
  if (z_e.size() != (size_t)m) {
    Rcpp::stop("Length of z_e must match the number of columns in X.");
  }

  // 2) Convert y to an Eigen vector
  Eigen::Map<const Eigen::VectorXd> Y_map(y.begin(), n);

  // 3) Prepare z (length m). If z is null => z = t(X)*y + z_e
  Eigen::VectorXd Z(m);
  if (z_.isNotNull()) {
    Rcpp::NumericVector z_in(z_);
    if ((int)z_in.size() != m) {
      Rcpp::stop("z has incorrect length (must be m).");
    }
    for (int col = 0; col < m; col++) {
      Z[col] = z_in[col];
    }
  } else {
    // Compute t(X)*y (which is length m), then add z_e
    Eigen::VectorXd XTy = X.transpose() * Y_map;  // (m x 1)
    for (int col = 0; col < m; col++) {
      Z[col] = XTy[col] + z_e[col];
    }
  }

  // 4) Define y_i. If y0 was given, use it; else all zeros.
  Rcpp::NumericVector y_i;
  if (y0_.isNotNull()) {
    Rcpp::NumericVector y0(y0_);
    if ((int)y0.size() != n) {
      Rcpp::stop("y0 must have length n.");
    }
    // Make a copy so we can modify it
    y_i = Rcpp::clone(y0);
  } else {
    // length n, fill with 0
    y_i = Rcpp::NumericVector(n, 0.0);
  }

  // 5) Define M (n x n). If given, just use it. Otherwise compute b * (X * X^T).
  Eigen::MatrixXd M_mat;  // NxN
  if (M_.isNotNull()) {
    Rcpp::NumericMatrix M_in(M_);
    if (M_in.nrow() != n || M_in.ncol() != n) {
      Rcpp::stop("M must be an n x n matrix.");
    }
    // Copy M_in into an Eigen dense matrix
    M_mat = Eigen::Map<const Eigen::MatrixXd>(M_in.begin(), n, n);
  } else {
    // M = b * (X * X^T) => n x n
    // X is n x m => X^T is m x n => product is n x n
    // We'll store in a dense matrix
    // This can be expensive for large n!
    Eigen::MatrixXd XXt = (X * X.transpose()).toDense();
    M_mat = b * XXt;
  }

  // 6) c_i = X %*% Z => length n
  Eigen::VectorXd c_i = X * Z;  // (n x 1)

  // 7) Set up RNG for sampling
  std::mt19937 rng(seed);

  // -------------------------------------------------------------------------
  // PART A: for (i in seq_len(n_b)) { ... }
  // -------------------------------------------------------------------------
  int i = 0;
  for (; i < n_b; i++) {
    // m <- max(c_i[y_i == 0])
    double max_val = -std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if (y_i[r] == 0.0 && c_i[r] > max_val) {
        max_val = c_i[r];
      }
    }

    // k <- which(c_i == m & y_i == 0)
    // gather the candidate indices
    std::vector<int> candidates;
    for (int r = 0; r < n; r++) {
      if (y_i[r] == 0.0 && std::fabs(c_i[r] - max_val) < 1e-15) {
        candidates.push_back(r);
      }
    }
    if (candidates.empty()) {
      // means no row with y_i[r]==0 => break out of loop
      break;
    }

    // If length(k)>1 => pick one at random
    int chosen;
    if (candidates.size() == 1) {
      chosen = candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
      chosen = candidates[dist(rng)];
    }

    // y_i[chosen] = b
    y_i[chosen] = b;

    // c_i = c_i - M[chosen, ]
    // row "chosen" of M => we subtract from each element of c_i
    for (int r = 0; r < n; r++) {
      c_i[r] -= M_mat(chosen, r);
    }
  }

  // -------------------------------------------------------------------------
  // PART B: while (i < max.iter) { i <- i+1; ... }
  // -------------------------------------------------------------------------
  double last_k_max = -std::numeric_limits<double>::infinity();
  double last_k_min =  std::numeric_limits<double>::infinity(); 
  // The R code sets last_k_min but never uses it. We still track it for completeness.
  
  // "last_y_i <- y_i" in R. We'll store a copy here:
  Rcpp::NumericVector last_y_i = Rcpp::clone(y_i);

  while (i < max_iter) {
    i++;

    // k_min <- min(c_i[y_i != 0])
    double min_val =  std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if (y_i[r] != 0.0 && c_i[r] < min_val) {
        min_val = c_i[r];
      }
    }

    // k_min_i <- which(c_i == min_val & y_i != 0)
    std::vector<int> min_candidates;
    for (int r = 0; r < n; r++) {
      if (y_i[r] != 0.0 && std::fabs(c_i[r] - min_val) < 1e-15) {
        min_candidates.push_back(r);
      }
    }
    if (min_candidates.empty()) {
      // No row with y_i != 0 => break
      break;
    }

    int chosen_min;
    if (min_candidates.size() == 1) {
      chosen_min = min_candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, min_candidates.size() - 1);
      chosen_min = min_candidates[dist(rng)];
    }

    // c_i <- c_i + M[chosen_min, ]
    for (int r = 0; r < n; r++) {
      c_i[r] += M_mat(chosen_min, r);
    }

    // k_max <- max(c_i[y_i == 0])
    double max_val = -std::numeric_limits<double>::infinity();
    for (int r = 0; r < n; r++) {
      if (y_i[r] == 0.0 && c_i[r] > max_val) {
        max_val = c_i[r];
      }
    }

    // k_max_i <- which(c_i == max_val & y_i == 0)
    std::vector<int> max_candidates;
    for (int r = 0; r < n; r++) {
      if (y_i[r] == 0.0 && std::fabs(c_i[r] - max_val) < 1e-15) {
        max_candidates.push_back(r);
      }
    }
    if (max_candidates.empty()) {
      // No row with y_i == 0 => break
      break;
    }

    int chosen_max;
    if (max_candidates.size() == 1) {
      chosen_max = max_candidates[0];
    } else {
      std::uniform_int_distribution<int> dist(0, max_candidates.size() - 1);
      chosen_max = max_candidates[dist(rng)];
    }

    // if (k_min >= k_max) break
    if (min_val >= max_val) {
      break;
    } 
    // else if (k_max <= last_k_max) { y_i <- last_y_i; break }
    else if (max_val <= last_k_max) {
      y_i = Rcpp::clone(last_y_i);
      break;
    }

    // c_i <- c_i - M[chosen_max, ]
    for (int r = 0; r < n; r++) {
      c_i[r] -= M_mat(chosen_max, r);
    }

    // y_i[chosen_min] <- 0
    // y_i[chosen_max] <- b
    y_i[chosen_min] = 0.0;
    y_i[chosen_max] = b;

    // last_k_max <- k_max
    last_k_max = max_val;
    // last_y_i <- y_i
    last_y_i = Rcpp::clone(y_i);
  }

  // Return final y_i
  return y_i;
}
