#' Performs a tidy two-sample wilcoxon rank-sum test
#' 
#' This is a wrapper around `stats::wilcoxon.test()`, but outputs results in 
#' a tidy format with calculation of different effect sizes.
#' 
#' @param x a numeric vector of data values for the first sample.
#' @param y a numeric vector of data values for the second sample.
#' @param ... further options to the `stats::wilcoxon.test()` function.
#' 
#' @return a data.frame with the following columns:
#' * n_x       number of observations in the first sample
#' * n_y       number of observations in the second sample
#' * x_less_y  number of pairwise comparisons where the values in x < y
#' * y_less_x  number of pairwise comparisons where the values in x > y
#' * x_equal_y number of pairwise comparisons where the values in x = y
#' * r         rank-biserial correlation
#' * vda       Vargas-Delayney A effect size measure (see `effsize::VD.A()` function)
#' 
wilcoxonRankSum <- function(x, y, ...){
  # Remove missing values
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  # Sample sizes
  n_x <- length(x)
  n_y <- length(y)
  
  # Make the test and tidy it into a table
  tidy_test <- wilcox.test(x, y, ...) %>% 
    broom::tidy()
  
  # Calculate effect sizes:
  # https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Effect_sizes
  
  # common language effect size (in both directions)
  x_less_y <- sum(as.vector(outer(x, y, "<")))/(n_x*n_y)
  y_less_x <- sum(as.vector(outer(x, y, ">")))/(n_x*n_y)
  x_equal_y <- sum(as.vector(outer(x, y, "==")))/(n_x*n_y)
  
  # Calculate rank-biserial correlation - this depends on the hypothesis
  # so using the output from the test
  r <- 1 - ((2*tidy_test$statistic)/(n_x*n_y))
  
  # Vargha and Delaney's A
  vda <- effsize::VD.A(x, y)$estimate
  
  # Output tidy data.frame
  tidy_test <- tidy_test %>% 
    mutate(n_x = n_x,
           n_y = n_y,
           x_less_y = x_less_y,
           y_less_x = y_less_x,
           x_equal_y = x_equal_y,
           r = r, 
           vda = vda)
}
