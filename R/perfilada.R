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


data <- rew(n = 5000L, alpha = 3.5, sigma = 1.5, theta = 1.5) 

func <- function(kicks,
                 names_par,
                 par_int,
                 data){
  
  names_par_int <- names_par[par_int]
  names_par_pert <- names_par[!par_int]
  
  p_log_lik <- function(a, x) {
    
    f <- function(par, x) {
      -sum(log(pdf_ew(par, x, var = list(names_par_int, a))))
    }
    
    pert <- optim(par = kicks, fn = f, x = data, method = "Nelder-Mead")$par[-(1L:sum(par_int))]
    
    -sum(log(pdf_ew(par = a, x, var = list(names_par_pert, pert))))
    
  }
  
  optim(par = kicks[par_int], fn = p_log_lik, x = data, method = "Nelder-Mead")
  
}

func(
  kicks = c(1, 1, 1),
  names_par = c("alpha", "sigma", "theta"),
  par_int = c(T, T, F),
  data = data
)