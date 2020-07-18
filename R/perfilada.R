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

# Unrestricted log-likelihood -----------------------------------

set.seed(0)
data <- rew(n = 5000L, alpha = 0.5, sigma = 1.5, theta = 1.5) 

func <- function(f,
                 kicks,
                 names_par,
                 par_int,
                 data){
  
  names_par_int <- names_par[par_int]
  names_par_pert <- names_par[!par_int]
  
  p_log_lik <- function(a, x) {
    
    fn <- function(par, x) {
      -sum(log(f(par, x, var = list(names_par_int, a))))
    }
    
    pert <- optim(par = kicks, fn = fn, x = data, method = "Nelder-Mead")$par[!par_int]
    
    -sum(log(f(par = a, x, var = list(names_par_pert, pert))))
    
  }
  
  optim(par = kicks[par_int], fn = p_log_lik, x = data, method = "Nelder-Mead")
  
}

func(
  f = pdf_ew,
  kicks = c(1, 1, 1),
  names_par = c("alpha", "sigma", "theta"),
  par_int = c(T, F, T),
  data = data
)



fn1 <- function(par, x) {
  -sum(log(pdf_ew(par, x)))
}

optim(par = c(1,1,1), fn = fn1, x = data)


# Distribuição Lindley Ponderada ------------------------------------------

pdf_lp <- function(par, x, var = NULL) {
  mu <- par[1L]
  beta <- par[2L]
  
  if (is.list(var))
    eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
  
  mu^(beta + 1) / ((mu + beta) * gamma(beta)) * x^(beta - 1) * (1 + x) * exp(-mu * x)
}

rlp <- function(n, mu, beta) {
  p <- mu / (mu + beta)
  
  p * rgamma(n = n, shape = beta, scale = mu) + (1 - p) * rgamma(n = n, shape = beta + 1, scale = mu)
}

data <- rlp(n = 1000, mu = 1.5, beta = 5.5)

log_lik <- function(par, x) {
  -sum(log(pdf_lp(par, x)))
}


optim(par = c(1, 1), fn = log_lik, x = data, method = "Nelder-Mead")

func(
  f = pdf_lp,
  kicks = c(1, 1),
  names_par = c("mu", "beta"),
  par_int = c(T, F),
  data = data
)