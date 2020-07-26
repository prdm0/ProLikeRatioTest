rlindleyp <- function(n = 1L, mu, beta, ...) {
  
  # Densidade da Lindley-Ponderada:
  pdf_lp <- function(y)
    mu^(beta + 1) / ((mu + beta) * gamma(beta)) * y^(beta - 1) * (1 + y) * exp(-mu * y)

  # Precisamos escolher c adequadaramente para que se rejeite a menor quantidade
  # possivel. Isto é, precisamos escolher o menor valor de c, uma vez que
  # P(aceitar) = 1/c.
  
  f <- function(y) {
    -1 * pdf_lp(y) / (1/100)
  }
  
  # Valor de c ótimo:  
  c <-  abs(optim(par = 1L, fn = f, method = "BFGS")$value)
  
  one_step <- function(i) {
    repeat{
      y <- runif(n = 1, 0, 100)
      u <- runif(n = 1L)
      cond <- pdf_lp(y) / (c * 1/100)
      
      if(u < cond) 
        break
    }
    y
  }
  
  sapply(X = 1L:n, FUN = one_step)
}

rlp <- function(n, mu, beta) {
  p <- rbinom(n, size = 1, prob = mu / (mu + beta))
  p * rgamma(n, shape = beta, rate = mu) +
    (1 - p) * rgamma(n, shape = beta + 1, rate = mu)
}


# Vendo o grafico ---------------------------------------------------------
set.seed(0L)
mu <- 0.2
beta <- 0.3

dados <- rlindleyp(n = 500, mu = mu, beta = beta)
#dados <- rlp(n = 500, mu = mu, beta = beta)
# Densidade da Lindley-Ponderada:
pdf_lp <- function(y, mu, beta)
  mu^(beta + 1) / ((mu + beta) * gamma(beta)) * y^(beta - 1) * (1 + y) * exp(-mu * y)

x <- seq(0, max(ceiling(dados)), length.out = 500L)
hist(dados, xlab = "Dados", probability = TRUE, ylab = "Densidade", main = "")
lines(x = x, y = pdf_lp(y = x, mu = mu, beta = beta))