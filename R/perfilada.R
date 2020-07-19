library(magrittr)

perfiled <- function(f,
                     kicks,
                     names_par,
                     par_int,
                     data,
                     ...){
  
  names_par_int <- names_par[par_int]
  names_par_pert <- names_par[!par_int]
  
  v <- double(length = length(kicks))
  
  emv <- NULL
  pert <- NULL
  
  body(f) %<>% as.list %>%
    append(quote(if (is.list(var))
      eval(parse(
        text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")
      ))), length(body(f)) - 1L) %>%
    as.call %>%
    as.expression
  
  myoptim <-
    function(...)
      tryCatch(
        expr = optim(...),
        error = function(e)
          NA
      )
  
  p_log_lik <- function(a, x) {
    
    fn <- function(par, x) {
      -sum(log(f(par, x, var = list(names_par_int, a))))
    }
    
    result_pert <- myoptim(par = kicks, fn = fn, x = data, ...)
    
    if(is.na(result_pert)) return(result_pert)
    
    pert <<- result_pert$par[!par_int]
    
    v[par_int] <- a
    
    -sum(log(f(par = v, x, var = list(names_par_pert, pert))))
  }
  
  result <- myoptim(par = kicks[par_int], fn = p_log_lik, x = data, ...)
  
  if(is.na(result)) return(result)
  
  emv[par_int] <- result$par
  emv[!par_int] <- pert
  result$par <- emv
  
  result
  
}

# Distribuição Lindley Ponderada ------------------------------------------

pdf_lp <- function(par, x, var = NULL) {
  mu <- par[1L]
  beta <- par[2L]

  mu^(beta + 1) / ((mu + beta) * gamma(beta)) * x^(beta - 1) * (1 + x) * exp(-mu * x)
}

rlp <- function(n, mu, beta) {
  p <- rbinom(n, size = 1, prob = mu / (mu + beta))
  p * rgamma(n, shape = beta, rate = mu) +
    (1 - p) * rgamma(n, shape = beta + 1, rate = mu)
}


set.seed(0L)

data <- rlp(n = 10000, mu = 5.5, beta = 4.5)

log_lik <- function(par, x) {
  -sum(log(pdf_lp(par, x)))
}

# Resultado por meio da varossimilhança perfilada:
resultado_perf <- 
  perfiled(
  f = pdf_lp,
  kicks = c(1, 1),
  names_par = c("mu", "beta"),
  par_int = c(F, T),
  data = data,
  method = "BFGS"
)

beta <- resultado_perf$par[2] # Estimativa de interesse
mu_perf_numerico <- resultado_perf$par[1]



# Por meio da teoria da dissertação, calculando a EMV perfilada de mu, em que beta é o parâmetro de interesse.
# A função abaixo implementa a estimativa de mu em função do parâmetro de interesse beta.
emvp_mu <- function(data, beta) {
  n <- length(data)
  t0 <- sum(data)
  
  (beta * (n - t0) + sqrt(beta^2 * (n + t0)^2 + 4 * n * beta * t0)) / (2 * t0)
}
mu_perf_teorico <- emvp_mu(data = data, beta = beta)

# Comparando os resultados:

mu_perf_numerico
mu_perf_teorico



# Distribuição Exp-Weibull ------------------------------------------------
pdf_ew <- function(par, x, var = NULL) {
  alpha <- par[1]
  sigma <- par[2]
  theta <- par[3]
  
  if (is.list(var))
    eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
  
  alpha * theta / sigma * (1 - exp(-(x / sigma) ^ alpha)) ^ (theta - 1) *
    exp(-(x / sigma) ^ alpha) * (x / sigma) ^ (alpha - 1)
}

rew <- function(n, alpha, sigma, theta) {
  u <- runif(n, 0, 1)
  sigma * (-log(1 - u ^ (1 / theta))) ^ (1 / alpha)
}


set.seed(0)
data <- rew(n = 5000L, alpha = 0.5, sigma = 1.5, theta = 1.5) 