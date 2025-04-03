create_dummy <- function(data, groupings) {
  n <- nrow(data)
  m <- nrow(groupings)
  k <- ncol(groupings)

  X <- matrix(0, nrow=n, ncol=m)

  dim_var <- colnames(groupings)
  data <- as.matrix(data[dim_var])

  for (i in seq_len(nrow(groupings))) {
    row <- matrix(groupings[i, ], nrow=n, ncol=k, byrow=TRUE)

    cond <- apply(data == row | row == "__total__", MARGIN=1, FUN=all)
 
    X[cond, i] <- 1
  }

  X
}


get_unique_and_total <- function(x) {
  c(as.character(unique(x)), "__total__")
}


reduce_X_y <- function(X, y, b) {
  z <- t(X) %*% y

  if (!any(z < b)) 
    return(list(X=NULL, y=NULL))

  inner <- apply(X[, z < b, drop=FALSE], MARGIN=1, FUN=\(x) any(x > 0))
  var0  <- apply(X[inner, , drop=FALSE], MARGIN=2, FUN=\(x) all(x == x[1]))

  X <- X[inner, !var0]
  y <- y[inner]

  list(X=X, y=y, mask_y = inner, mask_X = !var0)
}


get_groupings <- function(data, dim_var, total=TRUE, exclude_no_total=FALSE) {
  # unique groupings in data, without totals
  groupings <- dplyr::group_by_at(data, dim_var) |>
    dplyr::summarize(.groups="drop")

  if (!total)
    return(as.matrix(groupings))

  for (v in dim_var) {
    groupings_total_v <- unique(groupings[dim_var != v])
    groupings_total_v[[v]] <- "__total__"

    groupings <- rbind(groupings, groupings_total_v)
  }

  if (exclude_no_total)
    groupings <- groupings[apply(groupings == "__total__", MARGIN=1, FUN=any), ]

  as.matrix(groupings)
}


printf <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}


pls_rounding <- function(data, dim_var, freq_var, round_base=3, total=TRUE, 
                         exclude_no_total=FALSE, seed = 1234) {
  set.seed(seed)

  if (length(freq_var) > 1)
    stop("Only one frequency variable is allowed.")

  if (any(freq_var < 0))
    stop("Frequency variable must be non-negative.")

  data <- data[data[[freq_var]] > 0, ] # zero-cells are not relevant for the algorithm

  for (v in dim_var)
    data[[v]] <- as.character(data[[v]])

  groupings <- get_groupings(data, dim_var, total=total, exclude_no_total=exclude_no_total)

  b <- round_base
  X <- create_dummy(data, groupings)
  M <- b * (X %*% t(X))
  y <- data[[freq_var]]
  z <- t(X) %*% y
  y0 <- rep(0, length(y))

  c1 <- X %*% z
  y_rounded <- y

  X_i <- X
  y_i <- y

  last_mask <- rep(TRUE, length(y))

  while (TRUE) {
    reduced <- reduce_X_y(X=X, y=y_i, b=b)

    if (is.null(reduced$X))
      break

    y_rounded_sub <- round_cells(X=reduced$X, y=reduced$y, b=b)
    y_rounded[reduced$mask_y] <- y_rounded_sub
    last_mask[reduced$mask_y] <- FALSE

    y_i <- y_rounded
  }


  original <- t(X) %*% y
  rounded <- t(X) %*% y_rounded

  as.data.frame(groupings) |>
    dplyr::mutate(
      original = original,
      rounded = rounded,
      diff = original - rounded
    )

}


round_cells <- function(X, y, z=NULL,y0=NULL, b, M=NULL, max.iter=1000, 
                        seed = 1234) {
  set.seed(seed)

  if (is.null(z))
    z <- as.vector(t(X) %*% y)
  if (is.null(y0))
    y_i <- rep(0, length(y))
  if (is.null(M))
    M <- b * (X %*% t(X))

  c_i <- X %*% z

  nb <- round(sum(y) / b)

  for (i in seq_len(nb)) {
    printf("\rIteration %d", i)

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

    if (k_min >= k_max || (k_max == last_k_max))
      break
    else if (k_max < last_k_max) {
      y_i <- last_y_i
      break
    }

    c_i <- c_i - M[k_max_i, ]
    printf("\rIteration %d, kmax=%i", i, k_max)

    y_i[k_min_i] <- 0
    y_i[k_max_i] <- b 

    x <- as.vector(t(X) %*% y)
    x_hat <- as.vector(t(X) %*% y_i)
    printf("kmax = %d, kmin = %d, cov = %f", k_max, k_min, cov(x_hat, x))

    last_k_max <- k_max
    last_y_i <- y_i
  }

  y_i
}
