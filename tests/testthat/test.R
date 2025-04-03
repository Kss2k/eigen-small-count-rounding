library(tidyr)
library(dplyr)
devtools::load_all()

test_dat <- data.frame(
  row  = c("row1", "row2", "row3"), 
  col1 = c(6, 1, 0),
  col2 = c(0, 2, 1),
  col3 = c(1, 3, 1),
  col4 = c(3, 1, 0),
  col5 = c(4, 2, 2)
) |> pivot_longer(
  cols = starts_with("col"),
  names_to = "col",
  values_to = "n"
)

pls_rounding(test_dat, dim_var=c("col", "row"),
             freq_var="n", total=TRUE, exclude_no_total=TRUE, 
             round_base = 5)

X <- matrix(c(0, 1, 0, 1, 0, 0, 0,
              0, 0, 1, 1, 0, 0, 0,
              0, 0, 1, 0, 1, 0, 0,
              1, 0, 0, 0, 0, 1, 0,
              0, 1, 0, 0, 0, 1, 0,
              0, 0, 1, 0, 0, 0, 1), 
            nrow=6, ncol=7, byrow=TRUE)

y <- c(2, 1, 1, 3, 1, 2)

y0 <- rep(0, 6)
z1 <- t(X) %*% y
c1 <- X %*% z1
k <- 2
y0[k] <- 5
z2 <- t(X) %*% (y - y0)
c2 <- X %*% z2
k <- 4
y0[k] <- 5
