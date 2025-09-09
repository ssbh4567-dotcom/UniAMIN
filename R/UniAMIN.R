#' Asymptotic Unimin test for testing equality means against umbrella-ordered alternatives in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_1 = mu_2 = ... = mu_k vs H_1:mu_1 <=.....<= mu_(h-1)<= mu_h >= mu_(h+1)>=....>= mu_k (at least one strict inequality), where mu_i represents the population means of the i-th treatment. The input consists of two variables: sample_data and significance_level. The output consists of the critical value, the UniAMIN test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @importFrom MASS mvrnorm
#' @import stats
#' @author Subha Halder

UniAMIN <- function(sample_data, significance_level, peak){
  # basic checks & cleanup
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_datasets <- length(sample_data)
  if (num_datasets < 3) stop("Need at least 3 groups (length(sample_data) >= 3).")
  if (!is.numeric(peak) || peak != as.integer(peak)) stop("'peak' must be an integer.")
  if (peak <= 1 || peak >= num_datasets) {
    stop(paste("Error: 'peak' must be greater than 1 and less than", num_datasets, "."))
  }

  set.seed(456)
  num_samples <- 100000
  n <- sapply(sample_data, length)
  proportions <- n / sum(n)
  var_data <- sapply(1:num_datasets, function(j) var(sample_data[[j]]))
  h <- as.integer(peak)
  s <- num_datasets - 1

  # compute b_sq and B
  b_sq <- numeric(s)
  for (p in 1:s) {
    b_sq[p] <- proportions[p] * var_data[p + 1] + proportions[p + 1] * var_data[p]
  }
  B <- diag(sqrt(b_sq))

  # prepare upper and lower diagonals of length s-1
  upper_diag <- numeric(max(0, s - 1))
  lower_diag <- numeric(max(0, s - 1))

  if (length(upper_diag) > 0) {
    # fill indices 1:(h-2) if any
    if ((h - 2) >= 1) {
      for (i in 1:(h - 2)) {
        upper_diag[i] <- - sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
        lower_diag[i] <- upper_diag[i]
      }
    }

    # center element at index (h-1)
    center_idx <- h - 1
    upper_diag[center_idx] <- sqrt(proportions[h - 1] * proportions[h + 1]) * var_data[h]
    lower_diag[center_idx] <- upper_diag[center_idx]

    # fill indices h:(s-1) if present
    if (h <= (s - 1)) {
      for (i in h:(s - 1)) {
        upper_diag[i] <- - sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
        lower_diag[i] <- upper_diag[i]
      }
    }
  }

  # build tridiagonal matrix using diag() — avoids recycling warnings
  diag_elements <- b_sq
  tridiag_matrix <- matrix(0, nrow = s, ncol = s)
  diag(tridiag_matrix) <- diag_elements
  if (length(upper_diag) > 0) diag(tridiag_matrix[-1, ]) <- upper_diag    # superdiagonal
  if (length(lower_diag) > 0) diag(tridiag_matrix[, -1]) <- lower_diag    # subdiagonal

  D <- solve(B) %*% tridiag_matrix %*% solve(B)

  # bootstrap for critical value (min)
  D_star_amin_values <- numeric(num_samples)
  for (k in 1:num_samples) {
    bootstrap_s <- MASS::mvrnorm(mu = rep(0, s), Sigma = D)
    D_star_amin_values[k] <- min(bootstrap_s)
  }
  quantile_value <- as.numeric(quantile(D_star_amin_values, probs = 1 - significance_level))

  # test statistic using your h-based split (Aunimin)
  Aunimin <- min(c(
    sapply(2:h, function(j) {
      (mean(sample_data[[j]]) - mean(sample_data[[j - 1]])) /
        sqrt(
          (var(sample_data[[j]]) / length(sample_data[[j]])) +
            (var(sample_data[[j - 1]]) / length(sample_data[[j - 1]]))
        )
    }),
    sapply((h + 1):num_datasets, function(j) {
      (-mean(sample_data[[j]]) + mean(sample_data[[j - 1]])) /
        sqrt(
          (var(sample_data[[j]]) / length(sample_data[[j]])) +
            (var(sample_data[[j - 1]]) / length(sample_data[[j - 1]]))
        )
    })
  ))

  if (Aunimin > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("UniAMIN Critical value:", quantile_value, 
               "; UniAMIN Test statistic:", Aunimin, "; Result:", result))
}
