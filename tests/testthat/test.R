library(tidyr)
library(dplyr)
library(SmallCountRounding)

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


# generate a test datasets
set.seed(123)
indices <- 1:100
n <- 10000
col <- paste0("col", sample(indices, n, replace = TRUE))
row <- paste0("row", sample(indices, n, replace = TRUE))
n <- sample(0:10, n, replace = TRUE)
test_dat <- data.frame(
  row = row,
  col = col,
  n = n
) 

a <- pls_rounding(test_dat, dim_var=c("col", "row"),
             freq_var="n", total=TRUE, exclude_no_total=FALSE,  
             round_base = 5)
a2 <- lapply(a, FUN=\(x) {x[x=="__total__"] <- "Total"; x}) |>
  as.data.frame()
b <- PLSrounding(test_dat, 
            # dimVar=c("col", "row"),
            formula = ~ col * row,
            freqVar="n", 
            roundBase = 5)$publish
c <- left_join(x=a2, y=b, by=c("col", "row"))
cor(a$original, a$rounded)
cor(b$original, b$rounded)
