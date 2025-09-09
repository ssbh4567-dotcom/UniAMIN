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
  # Check peak validity
  if (peak == 1 || peak == length(sample_data)) {
  stop(paste("Error: 'peak' must be greater than 1 and less than", length(sample_data), "."))
}
  set.seed(456)
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  proportions <- n / sum(n)
  var_data <- sapply(1:num_datasets, function(j) var(sample_data[[j]]))
  h <- peak
  s <- num_datasets - 1
  b_sq <- numeric(s)
  for (p in 1:s) {
    b_sq[p] <- proportions[p] * var_data[p+1] + proportions[p+1] * var_data[p]
  }
  B <- diag(sqrt(b_sq))
  if(num_datasets == 3){
    diag_elements <- b_sq
    upper_diag <- numeric(s - 1)
    lower_diag <- numeric(s - 1)
    upper_diag[h - 1] <- sqrt(proportions[h - 1] * proportions[h + 1]) * var_data[h]
    lower_diag[h - 1] <- upper_diag[h - 1]
    tridiag_matrix <- matrix(0, nrow = s, ncol = s)
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
  } else if (num_datasets == 4){
    diag_elements <- b_sq
    upper_diag <- numeric(s - 1)
    lower_diag <- numeric(s - 1)
    upper_diag[h - 1] <- sqrt(proportions[h - 1] * proportions[h + 1]) * var_data[h]
    lower_diag[h - 1] <- upper_diag[h - 1]
    for (i in h:(s - 1)) {
      upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
      lower_diag[i] <- upper_diag[i]
    }
    tridiag_matrix <- matrix(0, nrow = s, ncol = s)
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
  } else{
    diag_elements <- b_sq
    upper_diag <- numeric(s - 1)
    lower_diag <- numeric(s - 1)
    for (i in 1:(h - 2)) {
      upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
      lower_diag[i] <- upper_diag[i]
    }
    upper_diag[h - 1] <- sqrt(proportions[h - 1] * proportions[h + 1]) * var_data[h]
    lower_diag[h - 1] <- upper_diag[h - 1]
    for (i in h:(s - 1)) {
      upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
      lower_diag[i] <- upper_diag[i]
    }
    tridiag_matrix <- matrix(0, nrow = s, ncol = s)
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
  }
  D <- solve(B) %*% tridiag_matrix %*% solve(B)

  D_star_amin_values <- numeric(num_samples)
  for (k in 1:num_samples) {
    bootstrap_s <- MASS::mvrnorm(mu = rep(0, s), Sigma = D)
    D_star_amin_values[k] <- min(bootstrap_s)
  }
  quantile_value <- quantile(sort(D_star_amin_values), probs = 1-significance_level)
  Aunimin <- min(c(
    sapply(2:ceiling(num_datasets / 2), function(j) {
      (mean(sample_data[[j]]) - mean(sample_data[[j - 1]])) /
        sqrt(
          (var(sample_data[[j]]) / length(sample_data[[j]])) +
            (var(sample_data[[j - 1]]) / length(sample_data[[j - 1]]))
        )
    }),
    sapply((ceiling(num_datasets / 2) + 1):num_datasets, function(j) {
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
  return(paste("UniAMIN Critical value:", quantile_value, "; UniAMIN Test statistic:", Aunimin, "; Result:", result))
}




