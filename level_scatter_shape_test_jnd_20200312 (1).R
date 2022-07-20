# Post-hoc comparisons of two profiles

n_dims <- 3
I <- diag(n_dims)
J <- matrix(rep(1, n_dims ^ 2), nrow = n_dims, ncol = n_dims)
L <- I - 1 / n_dims * J

var_scatter <- function(theta, se) {
  sigma <- diag(as.vector(se) ^ 2)
  expected_sq_scatter <- sum(diag(L %*% sigma)) + t(theta) %*% L %*% theta
  var_sq_scatter <- 2 * sum(diag(L %*% sigma %*% L %*% sigma)) + 4 * t(theta) %*% L %*% sigma %*% L %*% theta
  return(1 / (4 * expected_sq_scatter) * var_sq_scatter)
}

compare_level <- function(theta1, theta2, se1, se2) {
  avg <- sapply(list(theta1, theta2), mean)
  var_avg <- sapply(list(se1, se2), function(se) sum(se ^ 2) / n_dims ^ 2)
  Z <- (avg[1] - avg[2]) / sqrt(var_avg[1] + var_avg[2])
  # JND_update_2020031: 2-sided p value instead
  # p_value <- pnorm(abs(Z), lower.tail = FALSE)
  p_value <- 2* (1- pnorm(abs(Z)))
  return(list(avg = avg, var_avg = var_avg, Z = Z, p_value = p_value))
}

compare_scatter <- function(theta1, theta2, se1, se2) {
  scatter <- sapply(list(theta1, theta2), function(theta) sqrt(t(theta) %*% L %*% theta))
  var_scatter1 <- var_scatter(theta1, se1)
  var_scatter2 <- var_scatter(theta2, se2)
  Z <- (scatter[1] - scatter[2]) / sqrt(var_scatter1 + var_scatter2)
  # Z <- (scatter[1] - scatter[2]) / sqrt(var_scatter(theta1, se1) + var_scatter(theta2, se2))
  # JND_update_2020031: 2-sided p value instead
  # p_value <- pnorm(abs(Z), lower.tail = FALSE)
  p_value <- 2* (1- pnorm(abs(Z)))
  return(list(scatter = scatter, var_scatter1 = var_scatter1, var_scatter2 = var_scatter2, Z = Z, p_value = p_value))
}

compare_shape <- function(theta1, theta2, se1, se2) {
  var_std_profile <- function(theta, se, k) {
    mu_sq_X <- (theta[k] - mean(theta)) ^ 2
    sigma_sq_X <- ((n_dims - 1) / n_dims * se[k]) ^ 2 + sum(se[-k] ^ 2) / n_dims ^ 2
    expected_sq_scatter <- sum(diag(L %*% diag(as.vector(se) ^ 2))) + t(theta) %*% L %*% theta
    mu_sq_S <- ((sqrt(expected_sq_scatter) - (t(theta) %*% L %*% theta) ^ (-3 / 2) / 8 * (t(theta) %*% L %*% theta - expected_sq_scatter) ^ 2) / sqrt(n_dims)) ^ 2
    sigma_sq_S <- var_scatter(theta, se)
    return(sigma_sq_X / mu_sq_S + mu_sq_X * sigma_sq_S / mu_sq_S ^ 2)
  }
  std_profile <- lapply(list(theta1, theta2), function(theta) (theta - mean(theta)) / as.vector(sqrt(t(theta) %*% L %*% theta / n_dims)))
  
  cov_mat <- matrix(0, nrow = n_dims, ncol = n_dims)
  for (k in 1:n_dims) {
    cov_mat[k, k] <- var_std_profile(theta1, se1, k) + var_std_profile(theta2, se2, k)
  }

  Z <- t(std_profile[[1]] - std_profile[[2]]) %*% solve(cov_mat) %*% (std_profile[[1]] - std_profile[[2]])
  # JND_update_2020031: 2-sided p value instead
  # p_value <- pnorm(abs(Z), lower.tail = FALSE)
  p_value <- 2* (1- pnorm(abs(Z)))
  return(list(std_profile1 = std_profile[[1]], std_profile2 = std_profile[[2]], Z = Z, p_value = p_value))
}
# JND_TODO_20200312: wait, shape should be pchisq? 
