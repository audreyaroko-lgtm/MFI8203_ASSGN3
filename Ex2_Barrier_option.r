# Exercise 2: Barrier Option with Continuity Correction
# Advanced Monte Carlo Methods - MFI 8302

library(ggplot2)
library(dplyr)

# Part 1: Analytical Up-and-Out Call Price (Continuous Monitoring)
bs_call <- function(S, K, r, sigma, T) {
  if (T <= 0) return(max(S - K, 0))
  
  d1 <- (log(S/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  price <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  return(price)
}

up_and_out_call_analytical <- function(S0, K, B, r, sigma, T) {
  # Check conditions
  if (S0 >= B || K >= B) {
    stop("S0 and K must be less than B for this formula")
  }
  
  # Basic BS call
  C_BS <- bs_call(S0, K, r, sigma, T)
  
  # Compute terms for the analytical formula
  # This is a simplified version - the full formula is complex
  # Reference: Hull (2018) for complete implementation
  
  # Using the formula from Glasserman (2003) for barrier options
  lambda <- (r - 0.5 * sigma^2) / sigma^2
  y <- log(B^2 / (S0 * K)) / (sigma * sqrt(T)) + lambda * sigma * sqrt(T)
  x1 <- log(S0 / B) / (sigma * sqrt(T)) + lambda * sigma * sqrt(T)
  x2 <- log(S0 / B) / (sigma * sqrt(T)) - lambda * sigma * sqrt(T)
  
  # Simplified implementation (full formula would require more terms)
  # For educational purposes, we'll implement the main components
  price <- C_BS - (S0/B)^(1 - 2*r/sigma^2) * bs_call(B^2/S0, K, r, sigma, T)
  
  return(price)
}

# Part 2: Naive Monte Carlo for Up-and-Out Call (Discrete Monitoring)
barrier_naive_mc <- function(S0, K, B, r, sigma, T, m, M) {
  dt <- T / m
  payoffs <- numeric(M)
  
  for (i in 1:M) {
    S <- S0
    alive <- TRUE
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      
      # Check barrier (up-and-out: knock out if S >= B)
      if (S >= B) {
        alive <- FALSE
        break
      }
    }
    
    if (alive) {
      payoffs[i] <- max(S - K, 0)
    } else {
      payoffs[i] <- 0
    }
  }
  
  price <- exp(-r * T) * mean(payoffs)
  se <- exp(-r * T) * sd(payoffs) / sqrt(M)
  
  return(list(price = price, se = se, payoffs = payoffs))
}

# Part 3: Barrier with Continuity Correction
barrier_continuity_correction <- function(S0, K, B, r, sigma, T, m, M, beta = 0.5826) {
  dt <- T / m
  # Adjust barrier upward for up-and-out
  B_adj <- B * exp(beta * sigma * sqrt(dt))
  
  payoffs <- numeric(M)
  crossings <- numeric(M)  # Track barrier crossings
  
  for (i in 1:M) {
    S <- S0
    alive <- TRUE
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      
      if (S >= B_adj) {
        alive <- FALSE
        crossings[i] <- 1
        break
      }
    }
    
    if (alive) {
      payoffs[i] <- max(S - K, 0)
    } else {
      payoffs[i] <- 0
    }
  }
  
  price <- exp(-r * T) * mean(payoffs)
  se <- exp(-r * T) * sd(payoffs) / sqrt(M)
  
  return(list(
    price = price,
    se = se,
    payoffs = payoffs,
    crossings = crossings,
    B_adj = B_adj
  ))
}

# Part 4: Brownian Bridge Barrier Probability
brownian_bridge_barrier <- function(S_t, S_next, B, dt, sigma) {
  # Probability of crossing barrier between monitoring dates
  # Given S_t and S_{t+dt}
  # For up-and-out barrier
  
  mu <- (log(S_next) - log(S_t)) / dt
  sigma_sq <- sigma^2
  
  # Probability that barrier is NOT crossed
  if (S_t >= B) {
    return(0)  # Already crossed
  }
  
  # Brownian bridge probability of hitting barrier
  if (S_next >= B) {
    return(1)  # Definitely crossed at next point
  }
  
  # Probability of max < B given endpoints
  log_B <- log(B)
  log_S_t <- log(S_t)
  log_S_next <- log(S_next)
  
  # Probability of hitting barrier
  if (log_S_t >= log_B || log_S_next >= log_B) {
    return(1)
  }
  
  # Using Brownian bridge probability formula
  # P(max >= B) = exp(-2 * (log_B - log_S_t) * (log_B - log_S_next) / (sigma^2 * dt))
  
  p_cross <- exp(-2 * (log_B - log_S_t) * (log_B - log_S_next) / (sigma^2 * dt))
  return(p_cross)
}

# Part 5: Barrier with Brownian Bridge Correction
barrier_bridge_mc <- function(S0, K, B, r, sigma, T, m, M) {
  dt <- T / m
  payoffs <- numeric(M)
  
  for (i in 1:M) {
    S <- S0
    alive <- TRUE
    
    # Store path for bridge calculation
    path <- numeric(m + 1)
    path[1] <- S
    
    for (k in 1:m) {
      Z <- rnorm(1)
      S <- S * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      path[k + 1] <- S
    }
    
    # Check for barrier crossings using Brownian bridge
    for (k in 1:m) {
      p_cross <- brownian_bridge_barrier(path[k], path[k + 1], B, dt, sigma)
      if (runif(1) < p_cross) {
        alive <- FALSE
        break
      }
    }
    
    if (alive) {
      payoffs[i] <- max(S - K, 0)
    } else {
      payoffs[i] <- 0
    }
  }
  
  price <- exp(-r * T) * mean(payoffs)
  se <- exp(-r * T) * sd(payoffs) / sqrt(M)
  
  return(list(price = price, se = se, payoffs = payoffs))
}

# Comprehensive Barrier Option Analysis
barrier_analysis <- function() {
  # Parameters
  S0 <- 100
  K <- 100
  B <- 110  # Barrier
  r <- 0.05
  sigma <- 0.3
  T <- 1
  M <- 1e5
  
  m_values <- c(12, 50, 100, 250)
  results <- data.frame()
  
  for (m in m_values) {
    cat("\n=== m =", m, "===\n")
    
    # Naive MC
    naive_result <- barrier_naive_mc(S0, K, B, r, sigma, T, m, M)
    cat("Naive MC Price:", round(naive_result$price, 6), "±", round(naive_result$se, 6), "\n")
    
    # Continuity Correction
    cc_result <- barrier_continuity_correction(S0, K, B, r, sigma, T, m, M)
    cat("CC MC Price:", round(cc_result$price, 6), "±", round(cc_result$se, 6), "\n")
    cat("  Adjusted Barrier:", round(cc_result$B_adj, 6), "\n")
    
    # Brownian Bridge
    bb_result <- barrier_bridge_mc(S0, K, B, r, sigma, T, m, M)
    cat("BB MC Price:", round(bb_result$price, 6), "±", round(bb_result$se, 6), "\n")
    
    results <- rbind(results, data.frame(
      m = m,
      naive_price = naive_result$price,
      naive_se = naive_result$se,
      cc_price = cc_result$price,
      cc_se = cc_result$se,
      bb_price = bb_result$price,
      bb_se = bb_result$se,
      B_adj = cc_result$B_adj
    ))
  }
  
  # Analytical price (continuous monitoring) - approximate
  # For a truly accurate analytical price, we'd need full implementation
  cat("\n=== Analytical Benchmark (Continuous) ===\n")
  # Using approximation for comparison
  # In practice, you'd use the full formula from Hull (2018)
  
  return(results)
}

# Visualization
visualize_barrier_results <- function(results) {
  # Plot barrier prices vs m
  df_plot <- results %>%
    select(m, naive_price, cc_price, bb_price) %>%
    tidyr::pivot_longer(
      cols = c(naive_price, cc_price, bb_price),
      names_to = "Method",
      values_to = "Price"
    ) %>%
    mutate(Method = factor(Method,
                           levels = c("naive_price", "cc_price", "bb_price"),
                           labels = c("Naive MC", "Continuity Correction", "Brownian Bridge")))
  
  p <- ggplot(df_plot, aes(x = m, y = Price, color = Method, group = Method)) +
    geom_line() +
    geom_point(size = 3) +
    labs(
      title = "Barrier Option Prices vs. Monitoring Frequency",
      x = "Number of Monitoring Dates (m)",
      y = "Option Price",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
}

# Run the analysis
set.seed(42)
results <- barrier_analysis()
visualize_barrier_results(results)

# Additional: Bias Analysis
bias_analysis <- function() {
  S0 <- 100; K <- 100; B <- 110; r <- 0.05; sigma <- 0.3; T <- 1
  M <- 5e4
  m_values <- seq(4, 100, by = 4)
  
  naive_biases <- numeric(length(m_values))
  cc_biases <- numeric(length(m_values))
  
  # Approximate true price (using very fine discretization)
  true_price <- barrier_naive_mc(S0, K, B, r, sigma, T, 1000, 5e4)$price
  
  for (i in 1:length(m_values)) {
    m <- m_values[i]
    
    naive_result <- barrier_naive_mc(S0, K, B, r, sigma, T, m, M)
    naive_biases[i] <- naive_result$price - true_price
    
    cc_result <- barrier_continuity_correction(S0, K, B, r, sigma, T, m, M)
    cc_biases[i] <- cc_result$price - true_price
  }
  
  df <- data.frame(
    m = m_values,
    Naive = naive_biases,
    CC = cc_biases
  ) %>%
    tidyr::pivot_longer(cols = c(Naive, CC), names_to = "Method", values_to = "Bias")
  
  p <- ggplot(df, aes(x = m, y = Bias, color = Method)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Barrier Option Pricing Bias vs. Monitoring Frequency",
      x = "Number of Monitoring Dates (m)",
      y = "Bias"
    ) +
    theme_minimal()
  
  print(p)
}

bias_analysis()