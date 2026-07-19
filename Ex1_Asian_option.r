# Exercise 1: Asian Option with Control Variate
# Advanced Monte Carlo Methods - MFI 8302

library(ggplot2)
library(dplyr)

# Part 1: Geometric Asian Call Price (Analytical)
geometric_asian_price <- function(S0, K, r, sigma, T, m) {
  # Calculate parameters
  nu <- r - 0.5 * sigma^2
  mu_hat <- log(S0) + (nu * T / (2 * m)) * (m + 1)
  sigma2_hat <- (sigma^2 * T / (6 * m^2)) * (m + 1) * (2 * m + 1)
  sigma_hat <- sqrt(sigma2_hat)
  
  # Calculate d1 and d2
  d1 <- (mu_hat - log(K) + 0.5 * sigma2_hat) / sigma_hat
  d2 <- (mu_hat - log(K) - 0.5 * sigma2_hat) / sigma_hat
  
  # Price
  price <- exp(-r * T) * (exp(mu_hat + 0.5 * sigma2_hat) * pnorm(d1) - K * pnorm(d2))
  
  return(list(
    price = price,
    mu_hat = mu_hat,
    sigma_hat = sigma_hat,
    d1 = d1,
    d2 = d2
  ))
}

# Part 2: Arithmetic Asian Call Price (Plain Monte Carlo)
arithmetic_asian_mc <- function(S0, K, r, sigma, T, m, M) {
  dt <- T / m
  payoffs <- numeric(M)
  
  for (i in 1:M) {
    S <- S0
    sum_A <- 0
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      sum_A <- sum_A + S
    }
    
    avg <- sum_A / m
    payoffs[i] <- max(avg - K, 0)
  }
  
  price <- exp(-r * T) * mean(payoffs)
  se <- exp(-r * T) * sd(payoffs) / sqrt(M)
  
  return(list(price = price, se = se, payoffs = payoffs))
}

# Part 3: Arithmetic Asian with Geometric Control Variate
asian_control_variate <- function(S0, K, r, sigma, T, m, M) {
  dt <- T / m
  geo_price_known <- geometric_asian_price(S0, K, r, sigma, T, m)$price
  
  X <- numeric(M)  # Arithmetic payoffs
  Y <- numeric(M)  # Geometric payoffs
  
  for (i in 1:M) {
    S <- S0
    sum_A <- 0
    sum_G <- 0
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      sum_A <- sum_A + S
      sum_G <- sum_G + log(S)
    }
    
    avg_A <- sum_A / m
    avg_G <- exp(sum_G / m)
    
    X[i] <- max(avg_A - K, 0)
    Y[i] <- max(avg_G - K, 0)
  }
  
  # Compute optimal theta via regression
  theta_hat <- -cov(X, Y) / var(Y)
  
  # Controlled estimator
  Z <- X + theta_hat * (geo_price_known - Y)
  
  price_cv <- exp(-r * T) * mean(Z)
  se_cv <- exp(-r * T) * sd(Z) / sqrt(M)
  
  # Plain MC for comparison
  plain_price <- exp(-r * T) * mean(X)
  plain_se <- exp(-r * T) * sd(X) / sqrt(M)
  
  # Variance reduction factor
  vrf <- (plain_se / se_cv)^2
  
  return(list(
    price_cv = price_cv,
    se_cv = se_cv,
    theta_hat = theta_hat,
    plain_price = plain_price,
    plain_se = plain_se,
    vrf = vrf,
    correlation = cor(X, Y)
  ))
}

# Part 5: Repeat with different monitoring frequencies
analyze_m_frequencies <- function() {
  S0 <- 100; K <- 100; r <- 0.05; sigma <- 0.3; T <- 1
  M <- 1e5
  
  m_values <- c(12, 50, 250)
  results <- data.frame()
  
  for (m in m_values) {
    cat("\n=== m =", m, "===\n")
    
    # Geometric Asian price
    geo_result <- geometric_asian_price(S0, K, r, sigma, T, m)
    cat("Geometric Asian Price:", round(geo_result$price, 6), "\n")
    cat("  mu_hat:", round(geo_result$mu_hat, 6), "\n")
    cat("  sigma_hat:", round(geo_result$sigma_hat, 6), "\n")
    
    # Control variate results
    cv_result <- asian_control_variate(S0, K, r, sigma, T, m, M)
    cat("Arithmetic Asian Price (CV):", round(cv_result$price_cv, 6), "\n")
    cat("  SE (CV):", round(cv_result$se_cv, 6), "\n")
    cat("Arithmetic Asian Price (Plain):", round(cv_result$plain_price, 6), "\n")
    cat("  SE (Plain):", round(cv_result$plain_se, 6), "\n")
    cat("Optimal theta:", round(cv_result$theta_hat, 6), "\n")
    cat("Correlation:", round(cv_result$correlation, 6), "\n")
    cat("Variance Reduction Factor:", round(cv_result$vrf, 2), "\n")
    
    results <- rbind(results, data.frame(
      m = m,
      geo_price = geo_result$price,
      arith_cv = cv_result$price_cv,
      se_cv = cv_result$se_cv,
      arith_plain = cv_result$plain_price,
      se_plain = cv_result$plain_se,
      theta = cv_result$theta_hat,
      vrf = cv_result$vrf
    ))
  }
  
  return(results)
}

# Run the analysis
set.seed(42)
results <- analyze_m_frequencies()

# Plot convergence of VRF
ggplot(results, aes(x = m, y = vrf)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    title = "Variance Reduction Factor vs. Monitoring Frequency",
    x = "Number of Monitoring Dates (m)",
    y = "Variance Reduction Factor"
  ) +
  theme_minimal()

# Additional: Distribution comparison
distribution_analysis <- function() {
  S0 <- 100; K <- 100; r <- 0.05; sigma <- 0.3; T <- 1; m <- 12; M <- 1e5
  
  dt <- T / m
  X <- numeric(M)  # Arithmetic payoffs
  Y <- numeric(M)  # Geometric payoffs
  
  for (i in 1:M) {
    S <- S0
    sum_A <- 0
    sum_G <- 0
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      sum_A <- sum_A + S
      sum_G <- sum_G + log(S)
    }
    
    avg_A <- sum_A / m
    avg_G <- exp(sum_G / m)
    
    X[i] <- max(avg_A - K, 0)
    Y[i] <- max(avg_G - K, 0)
  }
  
  # Plot distributions
  df <- data.frame(
    Payoff = c(X, Y),
    Type = rep(c("Arithmetic", "Geometric"), each = M)
  )
  
  p <- ggplot(df, aes(x = Payoff, fill = Type)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Distribution of Arithmetic vs. Geometric Asian Payoffs",
      x = "Payoff",
      y = "Density"
    ) +
    theme_minimal()
  
  print(p)
  
  cat("\nCorrelation between payoffs:", cor(X, Y), "\n")
  cat("Mean Arithmetic Payoff:", mean(X), "\n")
  cat("Mean Geometric Payoff:", mean(Y), "\n")
}

distribution_analysis()