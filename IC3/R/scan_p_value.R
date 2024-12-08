#' Title: Calculate P-Value of scan statistic for Bernoulli random variables in 2D region 
#'
#' Description: This function calculates the p-value by performing a scan based on a matrix of communication statuses between cell pairs.
#'
#' @param A A numeric matrix where each row represents a cell pair. The first column represents the x-axis coordinates of the midpoint of this cell pair, the second column represents the y-axis coordinates, and the third column represents the communication status (1 for communication, 0 for no communication).
#' @param num_permutations Number of permutations to estimate the null distribution. Default is 1000.
#' @param min_size The ratio of the minimum value of the square window to the longest side of the slice
#' @param max_size The ratio of the maximum value of the square window to the longest side of the slice
#' @param sizeL The number of square side lengths to try to calculate (sampled evenly from minimum to maximum)
#' @param centerL The number of point center to try to calculate (sampled evenly from left to right or from bottle to top),if you choose centerL=10, we will try (10 times 10)100 possible centers.
#' @param max_size The ratio of the maximum value of the square window to the longest side of the slice
#' @return  5 parts. The first part is the statistic, the second part is the P value, the third part is the center of the optimal hotspot square, the fourth part is the side length of the optimal hotspot square. The fifth part is the null distribution result.
#' @examples
#' A <- matrix(c(1, 1, 1,
#'               2, 2, 0,
#'               3, 3, 1,
#'               4, 4, 0,
#'               5, 5, 0,
#'               6, 6, 0), 
#'             ncol = 3, byrow = TRUE)
#' scan_p_value(A, num_permutations = 1000)
#' @export
scan_p_value <- function(A, num_permutations = 1000,min_size = 0.02,max_size = 0.3,sizeL = 20,centerL = 10) {
  scan_statistic <- function(a, b, c, d) {  safe_log <- function(x) {
  if (x == 0) {
    return(0)
  } else {
    return(log(x))
  }
}
if (c == 0) {
  return(0)
}
if (d/c < b/a) {
  return(0)
} else {
  p1 <- d / c
  l1 <- d * safe_log(p1) + (c - d) * safe_log(1 - p1)
  p2 <- b / a
  l2 <- d * safe_log(p2) + (c - d) * safe_log(1 - p2)
  return(l1 - l2)
}
}

scan_statistic_from_matrix <- function(A, center, square_size) {
  x_center <- center[1]
  y_center <- center[2]
  x_min <- x_center - square_size / 2
  x_max <- x_center + square_size / 2
  y_min <- y_center - square_size / 2
  y_max <- y_center + square_size / 2
  a <- nrow(A)
  b <- sum(A[, 3] == 1)
  in_square <- A[, 1] >= x_min & A[, 1] <= x_max & A[, 2] >= y_min & A[, 2] <= y_max
  c <- sum(in_square)
  d <- sum(A[in_square, 3] == 1)
  return(scan_statistic(a, b, c, d))
}
scan_entire_region <- function(A, min_size ,max_size, sizeL, centerL) {
  x_min <- min(A[, 1])
  x_max <- max(A[, 1])
  y_min <- min(A[, 2])
  y_max <- max(A[, 2])  
  width <- x_max - x_min
  height <- y_max - y_min
  long_edge <- max(width, height)
  square_sizes <- seq(long_edge * min_size, long_edge * max_size, length.out =  sizeL)
  x_centers <- seq(x_min, x_max, length.out = centerL)
  y_centers <- seq(y_min, y_max, length.out = centerL)
  max_statistic <- -Inf
  best_center <- c(NA, NA)
  best_size <- NA
  for (square_size in square_sizes) {
    for (x_center in x_centers) {
      for (y_center in y_centers) {
        center <- c(x_center, y_center)
        current_stat <- scan_statistic_from_matrix(A, center, square_size)
        if (current_stat > max_statistic || 
            (current_stat == max_statistic && square_size < best_size)) {
          max_statistic <- current_stat
          best_center <- center
          best_size <- square_size
        }
      }
    }
  }
  
  return(list(statistic = max_statistic, center = best_center, size = best_size))
}
library(progress)
  observed_result <- scan_entire_region(A, min_size ,max_size, sizeL, centerL)
  observed_statistic <- observed_result$statistic
  best_center <- observed_result$center
  best_size <- observed_result$size
  
  num_points <- nrow(A)
  num_red_points <- sum(A[, 3] == 1)
  
  null_statistics <- numeric(num_permutations)
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed",
    total = num_permutations, clear = FALSE, width = 60
  )  
  for (i in 1:num_permutations) {
    pb$tick()  
    shuffled_labels <- sample(c(rep(1, num_red_points), rep(0, num_points - num_red_points)))
    A_shuffled <- A
    A_shuffled[, 3] <- shuffled_labels    
    null_statistics[i] <- scan_entire_region(A_shuffled, min_size ,max_size, sizeL, centerL)$statistic
  }
  
  p_value <- mean(null_statistics >= observed_statistic)
  
  return(list(observed_statistic = observed_statistic, 
              p_value = p_value, 
              best_center = best_center, 
              best_size = best_size, 
              null_statistics = null_statistics))
}
