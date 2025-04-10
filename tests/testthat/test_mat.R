M <- \(...) matrix(..., byrow=TRUE)
D <- function(x) {
  d <- matrix(0, nrow=length(x), ncol=sum(x))

  j <- 1
  for (i in seq_along(x)) {
    if (x[i]) {
      d[i, j] <- 1
      j <- j + 1
    }
  }
  d
}

X <- M(c(1, 0, 0, 2, 3, 4,
         0, 5, 6, 0, 7, 0,
         8, 0, 0, 0, 0, 9,
         0, 10, 11, 0, 0, 12), nrow=4)
kx <- c(1, 1, 0, 0)
ky <- c(1, 0, 1, 1, 0, 1)
bx <- as.logical(kx)
by <- as.logical(ky)

X[bx, by]
X[bx, ]

t(t(X) %*% D(kx))
t(D(kx)) %*% X
t(D(kx)) %*% X %*% D(ky)
X %*% D(ky)

y <- c(1, 0, 0, 2, 3, 4)
x <- c(1, 0, 8, 0)
t(D(kx)) %*% x

