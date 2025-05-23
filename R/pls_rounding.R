TOTAL_C = -9
TOTAL_L = "__total__"
syms <- rlang::syms

create_dummy_old <- function(data, groupings) {
  n <- nrow(data)
  m <- nrow(groupings)
  k <- ncol(groupings)

  X <- matrix(0, nrow=n, ncol=m)

  dim_var <- colnames(groupings)
  data <- as.matrix(data[dim_var])

  for (i in seq_len(nrow(groupings))) {
    row <- matrix(groupings[i, ], nrow=n, ncol=k, byrow=TRUE)

    cond <- apply(data == row | row == TOTAL, MARGIN=1, FUN=all)
 
    X[cond, i] <- 1
  }

  X
}


create_dummy <- function(data, groupings) {

  dim_var <- colnames(groupings)
  groupings_inner <- as.matrix(data[dim_var])
  groupings_publish <- as.matrix(groupings)

  printf("Creating dummy matrix")
  create_dummy_cpp(
    groupings_inner=groupings_inner,
    groupings_publish=groupings_publish
  )
}


get_unique_and_total <- function(x) {
  c(as.character(unique(x)), TOTAL_L)
}


reduce_X_y <- function(X, y_i, b, z) {
  z_i <- t(X) %*% y_i

  if (!any(z_i < b)) 
    return(list(X=NULL, y=NULL))

  inner <- apply(X[, z_i < b, drop=FALSE], MARGIN=1, FUN=\(x) any(x > 0))
  var0  <- apply(X[inner, , drop=FALSE], MARGIN=2, FUN=\(x) all(x == x[1]))

  X_i <- X[inner, !var0] 
  y_i <- y_i[inner]
  z_e <- (z - z_i)[!var0]

  list(X_i=X_i, y_i=y_i, mask_y=inner, mask_X= !var0, z_e=z_e)
}


get_groupings <- function(data, dim_var, total=TRUE, exclude_no_total=FALSE) {
  # unique groupings in data, without totals

  printf("Getting groupings")
  groupings <- dplyr::group_by_at(data, dim_var) |>
    dplyr::summarize(.groups="drop")

  
  if (!total)
    return(as.matrix(groupings |> dplyr::arrange(!!!syms(dim_var))))

  for (v in dim_var) {
    groupings_total_v <- unique(groupings[dim_var != v])
    groupings_total_v[[v]] <- -9

    groupings <- rbind(groupings, groupings_total_v)
  }

  if (exclude_no_total)
    groupings <- groupings[apply(groupings == TOTAL_C, MARGIN=1, FUN=any), ]

  as.matrix(groupings |> dplyr::arrange(!!!syms(dim_var)))
}


printf <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}


get_levels <- function(factor) {
  c(levels(factor), TOTAL_L)
}


get_codes <- function(factor) {
  c(seq_along(levels(factor)), TOTAL_C)
}


pls_rounding <- function(data, dim_var, freq_var, round_base=3, total=TRUE, 
                         exclude_no_total=FALSE, seed = 1234, max_iter=1000) {
  set.seed(seed)

  if (length(freq_var) > 1)
    stop("Only one frequency variable is allowed.")

  if (any(freq_var < 0))
    stop("Frequency variable must be non-negative.")

  data <- data[data[[freq_var]] > 0, ] |> # zero-cells are not relevant for the algorithm
    dplyr::arrange(!!!syms(dim_var))
  factors <- list()

  for (v in dim_var) {
    factor <- as.factor(data[[v]])
    factors[[v]] <- list(levels=get_levels(factor), codes=get_codes(factor))
    data[[v]] <- as.integer(factor)
  }

  groupings_publish <- get_groupings(data, dim_var, total=total, exclude_no_total=exclude_no_total)
  groupings_inner <- as.matrix(data[dim_var])


  b <- round_base
  printf("Creating dummy matrix")
  X <- create_dummy_cpp(groupings_inner=groupings_inner, 
                        groupings_publish=groupings_publish)
  y <- as_xptr_vector(data[[freq_var]])
  z <- calc_z(X, y)

  y_rounded <- copy_xptr_VectorXd(y)
  y_i <- copy_xptr_VectorXd(y)

  printf("Rounding %s cells", length(y))
  i <- 0
  while (TRUE) {
    printf(" > iteration %s", (i <- i + 1))
    printf(" > reducing X")

    reduced <- reduce_X_y_cpp(X=X, y_i=y_rounded, b=b, z=z)

    printf(" > reduced X")
    if (is.null(reduced$X_i))
      break

    n_b <- calc_n_b(xptr_y_i=reduced$y_i, xptr_y=y, xptr_y_rounded=y_rounded, b=b)

    if (n_b == 0)
      break
    printf(" > rounding %s cells", n_b)

    y_i <- round_cells_cpp(X=reduced$X_i, y=reduced$y_i, b=b, z_e=reduced$z_e, n_b=n_b, 
                           max_iter=max_iter, seed=seed)
   
    # y_rounded[reduced$mask_y] <- as_numeric_vector(y_i)
    # print_VectorXd( y_i)
    # print_VectorXd( y_rounded)
    fill_vector_by_mask(xptr_x=y_rounded, xptr_y=y_i, xptr_mask=reduced$mask_y)
    # print_VectorXd( y_rounded)
  }


  original <- as_numeric_vector(calc_z(X, y))
  rounded <- as_numeric_vector(calc_z(X, y_rounded))

  as.data.frame(groupings_publish) |>
    dplyr::mutate(
      original = original,
      rounded = rounded,
      diff = rounded - original,
    )

}


round_cells <- function(X, y, z=NULL,y0=NULL, b, n_b, 
                        M=NULL, max.iter=1000, 
                        z_e, seed = 1234) {
  set.seed(seed)

  if (is.null(z))
    z <- as.vector(t(X) %*% y) + z_e
  if (is.null(y0))
    y_i <- rep(0, length(y))
  if (is.null(M))
    M <- b * (X %*% t(X))

  c_i <- X %*% z

  for (i in seq_len(n_b)) {
    m <- max(c_i[y_i == 0])

    k <- which(c_i == m & y_i == 0)

    if (length(k) > 1) k <- sample(k, 1)

    y_i[k] <- b

    c_i <- c_i - M[k, ]
  }


  last_k_max <- -Inf
  last_k_min <- Inf
  last_y_i <- NA

  while (i < max.iter) {
    i <- i + 1

    k_min <- min(c_i[y_i != 0]) 
    k_min_i <- which(c_i == k_min & y_i != 0)
    if (length(k_min_i) > 1) k_min_i <- sample(k_min_i, 1)

    c_i <- c_i + M[k_min_i, ]

    k_max <- max(c_i[y_i == 0])
    k_max_i <- which(c_i == k_max & y_i == 0)
    if (length(k_max_i) > 1) k_max_i <- sample(k_max_i, 1)

    if (k_min >= k_max)
      break
    else if (k_max <= last_k_max) {
      y_i <- last_y_i
      break
    }

    c_i <- c_i - M[k_max_i, ]

    y_i[k_min_i] <- 0
    y_i[k_max_i] <- b 

    last_k_max <- k_max
    last_y_i <- y_i
  }

  y_i
}
