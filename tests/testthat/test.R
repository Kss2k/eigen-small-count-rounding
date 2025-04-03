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

a <- pls_rounding(test_dat, dim_var=c("col", "row"),
             freq_var="n", total=TRUE, exclude_no_total=TRUE, 
             round_base = 5)

b <- PLSrounding(test_dat, 
            # dimVar=c("col", "row"),
            formula = ~ col + row,
            freqVar="n", 
            roundBase = 5)$publish
