# Exercise 3: Computing Delta Three Ways
# Advanced Monte Carlo Methods - MFI 8302

library(ggplot2)
library(dplyr)
library(tidyr)

# Part 1: Analytical Black-Scholes Delta
bs_delta <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  delta <- pnorm(d1)
  return(delta)
}

bs_digital_delta <- function(S0, K, r, sigma, T) {
  d2 <- (log(S0/K) + (r - 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  # Digital call delta
  delta <- exp(-r * T) * dnorm(d2) / (S0 * sigma * sqrt(T))
  return(delta)
}

# European Call Payoff and its derivative
call_payoff <- function(S_T, K) {
  return(pmax(S_T - K, 0))
}

call_payoff_derivative <- function(S_T, K) {
  return(as.numeric(S_T > K))
}

digital_payoff <- function(S_T, K) {
  return(as.numeric(S_T > K))
}

# Part 1: Bump-and-Reprice Delta
bump_reprice_delta <- function(S0, K, r, sigma, T, M, h, payoff_func) {
  # Using Common Random Numbers (CRN)
  set.seed(42)  # For reproducibility
  
  # Generate common random numbers
  Z <- rnorm(M)
  
  deltas <- numeric(M)
  
  for (i in 1:M) {
    # Price at S0 + h
    S_T_plus <- (S0 + h) * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z[i])
    V_plus <- exp(-r * T) * payoff_func(S_T_plus, K)
    
    # Price at S0 - h
    S_T_minus <- (S0 - h) * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z[i])
    V_minus <- exp(-r * T) * payoff_func(S_T_minus, K)
    
    # Central difference
    deltas[i] <- (V_plus - V_minus) / (2 * h)
  }
  
  delta_hat <- mean(deltas)
  se <- sd(deltas) / sqrt(M)
  
  return(list(
    delta = delta_hat,
    se = se,
    deltas = deltas
  ))
}

# Part 2: Pathwise (IPA) Delta
pathwise_delta <- function(S0, K, r, sigma, T, M, payoff_func, payoff_deriv_func) {
  set.seed(42)
  
  deltas <- numeric(M)
  
  for (i in 1:M) {
    Z <- rnorm(1)
    S_T <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
    
    # Pathwise derivative
    dS_T_dS0 <- S_T / S0
    delta_path <- exp(-r * T) * payoff_deriv_func(S_T, K) * dS_T_dS0
    
    deltas[i] <- delta_path
  }
  
  delta_hat <- mean(deltas)
  se <- sd(deltas) / sqrt(M)
  
  return(list(
    delta = delta_hat,
    se = se,
    deltas = deltas
  ))
}

# Part 3: Likelihood Ratio Method Delta
lrm_delta <- function(S0, K, r, sigma, T, M, payoff_func) {
  set.seed(42)
  
  deltas <- numeric(M)
  nu <- r - 0.5 * sigma^2
  
  for (i in 1:M) {
    Z <- rnorm(1)
    S_T <- S0 * exp(nu * T + sigma * sqrt(T) * Z)
    
    # Score function w.r.t. S0
    score <- Z / (S0 * sigma * sqrt(T))
    
    # LRM delta
    delta_lrm <- exp(-r * T) * payoff_func(S_T, K) * score
    
    deltas[i] <- delta_lrm
  }
  
  delta_hat <- mean(deltas)
  se <- sd(deltas) / sqrt(M)
  
  return(list(
    delta = delta_hat,
    se = se,
    deltas = deltas
  ))
}

# Comprehensive Delta Comparison
delta_comparison <- function() {
  # Parameters
  S0 <- 100
  K <- 100
  r <- 0.05
  sigma <- 0.2
  T <- 0.5
  M <- 1e5
  h <- 1  # Bump size for central difference
  
  cat("=== European Call Delta ===\n")
  cat("Parameters: S0 =", S0, ", K =", K, ", r =", r, 
      ", sigma =", sigma, ", T =", T, "\n")
  cat("M =", M, ", h =", h, "\n\n")
  
  # Analytical Delta
  delta_bs <- bs_delta(S0, K, r, sigma, T)
  cat("1. Analytical BS Delta:", round(delta_bs, 6), "\n")
  
  # Bump-and-Reprice
  br_result <- bump_reprice_delta(S0, K, r, sigma, T, M, h, call_payoff)
  cat("2. Bump-and-Reprice Delta:", round(br_result$delta, 6), 
      "±", round(br_result$se, 6), "\n")
  cat("   Bias:", round(br_result$delta - delta_bs, 6), "\n")
  
  # Pathwise IPA
  ipa_result <- pathwise_delta(S0, K, r, sigma, T, M, call_payoff, call_payoff_derivative)
  cat("3. Pathwise IPA Delta:", round(ipa_result$delta, 6), 
      "±", round(ipa_result$se, 6), "\n")
  cat("   Bias:", round(ipa_result$delta - delta_bs, 6), "\n")
  
  # LRM
  lrm_result <- lrm_delta(S0, K, r, sigma, T, M, call_payoff)
  cat("4. LRM Delta:", round(lrm_result$delta, 6), 
      "±", round(lrm_result$se, 6), "\n")
  cat("   Bias:", round(lrm_result$delta - delta_bs, 6), "\n")
  
  # Compare SEs
  cat("\n5. Standard Error Comparison:\n")
  cat("   Bump-and-Reprice SE:", round(br_result$se, 6), "\n")
  cat("   Pathwise IPA SE:", round(ipa_result$se, 6), "\n")
  cat("   LRM SE:", round(lrm_result$se, 6), "\n")
  
  # Efficiency comparison
  cat("\n   Relative Efficiency (IPA vs BR):", round((br_result$se / ipa_result$se)^2, 2), "\n")
  cat("   Relative Efficiency (IPA vs LRM):", round((lrm_result$se / ipa_result$se)^2, 2), "\n")
  
  return(list(
    bs_delta = delta_bs,
    br = br_result,
    ipa = ipa_result,
    lrm = lrm_result
  ))
}

# Digital Call Delta Comparison
digital_delta_comparison <- function() {
  # Parameters
  S0 <- 100
  K <- 100
  r <- 0.05
  sigma <- 0.2
  T <- 0.5
  M <- 1e5
  h <- 0.5  # Smaller h for digital options
  
  cat("\n=== Digital Call Delta ===\n")
  cat("Parameters: S0 =", S0, ", K =", K, ", r =", r, 
      ", sigma =", sigma, ", T =", T, "\n")
  cat("M =", M, ", h =", h, "\n\n")
  
  # Analytical Delta
  delta_dig_bs <- bs_digital_delta(S0, K, r, sigma, T)
  cat("1. Analytical Digital Delta:", round(delta_dig_bs, 6), "\n")
  
  # Bump-and-Reprice (works for any payoff)
  br_result <- bump_reprice_delta(S0, K, r, sigma, T, M, h, digital_payoff)
  cat("2. Bump-and-Reprice Delta:", round(br_result$delta, 6), 
      "±", round(br_result$se, 6), "\n")
  cat("   Bias:", round(br_result$delta - delta_dig_bs, 6), "\n")
  
  # Pathwise IPA (FAILS for digital options - payoff not differentiable)
  cat("\n3. Pathwise IPA: NOT APPLICABLE (payoff not differentiable)\n")
  
  # LRM (works for any payoff)
  lrm_result <- lrm_delta(S0, K, r, sigma, T, M, digital_payoff)
  cat("4. LRM Delta:", round(lrm_result$delta, 6), 
      "±", round(lrm_result$se, 6), "\n")
  cat("   Bias:", round(lrm_result$delta - delta_dig_bs, 6), "\n")
  
  cat("\n5. Applicable Methods:\n")
  cat("   Bump-and-Reprice: Yes\n")
  cat("   Pathwise IPA: No (payoff discontinuous)\n")
  cat("   LRM: Yes\n")
  
  cat("\n6. SE Comparison:\n")
  cat("   Bump-and-Reprice SE:", round(br_result$se, 6), "\n")
  cat("   LRM SE:", round(lrm_result$se, 6), "\n")
  cat("   Relative Efficiency (BR vs LRM):", round((lrm_result$se / br_result$se)^2, 2), "\n")
  
  return(list(
    bs_delta = delta_dig_bs,
    br = br_result,
    lrm = lrm_result
  ))
}

# Convergence Analysis
convergence_analysis <- function() {
  S0 <- 100; K <- 100; r <- 0.05; sigma <- 0.2; T <- 0.5
  M_values <- c(1000, 5000, 10000, 50000, 100000)
  h <- 1
  
  results <- data.frame()
  
  for (M in M_values) {
    cat("\nM =", M, "\n")
    
    # Analytical
    delta_bs <- bs_delta(S0, K, r, sigma, T)
    
    # Bump-and-Reprice
    br <- bump_reprice_delta(S0, K, r, sigma, T, M, h, call_payoff)
    
    # IPA
    ipa <- pathwise_delta(S0, K, r, sigma, T, M, call_payoff, call_payoff_derivative)
    
    # LRM
    lrm <- lrm_delta(S0, K, r,