#' @title Teste da razão de verossimilhança generalizada
#' @author Pedro Rafael D. Marinho
#' @description Calcula a estatística da razão de verossimilhança generalizada.
#' @details A função recebe como argumento uma função densidade de probabilidade ou uma função de probabilidade
#' que é passada como argumento à \code{f}. O objeto passado como argumento de \code{f} deverá ser implementado conforme os
#' exemplos. Note que se a função for implementada segundo os exemplos, não haverá a necessidade de implementar a função sob a
#' hipótese nula, visto que são as mesmas funções. Para especificar os parâmetros que serão fixados, i.e, para especificar
#' a distribuição sob a hipótese nula, utiliza-se o argumento \code{par0} que receberá uma lista formada por dois vetores.
#' O primeiro vetor da lista de verá ser um vetor de strings com os nomes das variáveis que deseja-se fixar e o segundo vetor
#' deverá conter os valores que serão atribuidos à cada uma das variáveis passada ao primeiro vetor.
#' @param f Função densidade de probabilidade considerada no teste da razão de verossimilhança. Essa função deverá ser
#' implementada conforme os exemplos.
#' @param data Um conjunto de dados que será considerado para a realização do teste.
#' @param par0 Uma lista contendo como primeiro elemento um vetor com o nomes das variáveis que serão fixadas como hipótese nula
#' (variáveis aos quais desejamos testar) e um vetor com os valores fixados para cada uma das respectivas variáveis.
#' @param par1 Lista com dois elementos, sendo o primeiro um vetor com os nomes das variáveis que receberão valores fixos sob a
#' hipótese alternativa e o segundo elemendo é um outro vetor com os valores impostos às variáveis.
#' @param kicks Chutes iniciais que desejamos considerar no método de otimização.
#' @param ... Lista de argumentos adicionais que serão passados à função \code{optim()} otimizada para otimização. Por exemplo,
#' será possível escolher o método de otimização a ser utilizado.
#' @importFrom stats optim
#' @importFrom stats qchisq
#' @return Valor da estatística do teste da razão de verossimilhança generalizada.
#' @examples
#'  pdf_w <- function(par, x, var = NULL){
#'      alpha <- par[1]
#'      beta <- par[2]
#'      if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
#'      dweibull(x, shape = alpha, scale = beta)
#'   }
#'
#'   rw <- function(n = 1L, alpha, beta){
#'    rweibull(n = n, shape = alpha, scale = beta)
#'   }
#'
#'   data <- rw(n = 100L, alpha = 1, beta = 1)
#'
#'   lrt(f = pdf_w, data = data, kicks = c(1, 1), par0 = list("beta", 1))
#' @export
lrt <- function(f, data, kicks, par0 = NULL, ...) {
  if (is.null(par0))
    stop("Informar uma lista informando o parâmetro e o valor sob a hipótese nula.")

  # Log-Likelihood under the null hypothesis. -----------------------------------
  log_lik_h0 <- function(par, x) {
    -sum(log(f(par, x, var = par0)))
  }

  # Unrestricted log-likelihood. -----------------------------------
  log_lik <- function(par, x) {
    -sum(log(f(par, x)))
  }

  myoptim <-
    function(...)
      tryCatch(
        expr = optim(...),
        error = function(e)
          NA
      )

  par_h0 <- myoptim(par = kicks, fn = log_lik_h0, x = data, ...)

  if (!is.list(par_h0) || par_h0$convergence != 0L)
    return(NA)

  par_h <- myoptim(par = kicks, fn = log_lik, x = data, ...)

  if (!is.list(par_h) || par_h$convergence != 0L)
    return(NA)

  lambda <-
    2 * (log_lik_h0(par = par_h0$par, x = data) - log_lik(par_h$par, x = data))

  lambda[lambda < 0] <- 0

  # Estatística de razão de verossimilhança:
  lambda
}

#' @title Simulações de Monte-Carlo para o teste da razão de verossimilhança generalizado.
#' @author Pedro Rafael D. Marinho
#' @description Realiza uma única iteração de um procedimento de Monte-Carlo para o teste da razão de verossimilhança generalizado. Dado um nível
#' de significância, será retornado 1 (um) se a estatística de teste está acima do quantil da distribuição qui-quadrado e 0 (zero),
#' caso contrário.
#' @param N Número de réplicas de Monte-Carlo a ser considerada.
#' @param n Tamanho da amostra a ser considerada.
#' @param sig Nível de significância adotado.
#' @param f Função densidade de probabilidade considerada no teste. Essa função deverá ser implementada conforme o exemplo abaixo.
#' @param q Função responsável pela geração de observações de uma variável aleatório com função densidade passada para \code{f}.
#' @param kicks Vetor com os chutes iniciais utilizados para a otimização.
#' @param par0 Lista com dois elementos, sendo o primeiro um vetor com os nomes das variáveis que receberão valores fixos sob a
#' hipótese nula e o segundo elemendo é um outro vetor com os valores impostos às variáveis.
#' @param ncores Número de núcleos a ser considerado. Por padrão, \code{ncores = 1L}.
#' @param bilateral Se \code{TRUE}, retorna os quantis para um teste bilateral. O padrão considera
#' \code{bilateral  = FALSE}.
#' @param p Valor utilizado para controlar o parâmetro da Qui-quadrado inf.
#' @param step Tamanho do passo da integral numeérica responsável pela obtenção dos quantis da Qui-Quadrado inf.
#' O padrão considera \code{step = 1e-3}.
#' @param ... Lista de argumetos que serão passados para a função passada à \code{q}.
#' @importFrom tibble as_tibble
#' @importFrom stats dchisq pchisq
#' @return Retornará 0 (zero) se a estatística calculado não estiver acima do quantil da distribuição qui-quadrado e 1 (um),
#' caso contrário.
#' @examples
#' pdf_ew <- function(par, x, var = NULL){
#' alpha <- par[1]
#' sigma <- par[2]
#' theta <- par[3]
#'
#' if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
#'
#' alpha * theta / sigma * (1 - exp(-(x / sigma) ^ alpha)) ^ (theta - 1) *
#'   exp(-(x / sigma) ^ alpha) * (x / sigma) ^ (alpha - 1)
#' }
#'
#' rew <- function(n, alpha, sigma, theta){
#'   u <- runif(n, 0, 1)
#'   sigma * (-log(1 - u ^ (1 / theta))) ^ (1 / alpha)
#' }
#'
#' set.seed(1L, kind = "L'Ecuyer-CMRG")
#'
#' tictoc::tic()
#' result <- mc(N = 100L,
#'              n = 50L,
#'              sig = 0.05,
#'              f = pdf_ew,
#'              q = rew,
#'              kicks = c(1, 1, 1),
#'              par0 = list("theta", 1),
#'              ncores = 1L,
#'              p = 0.5,
#'              bilateral = FALSE,
#'              step = 0.001, alpha = 1, sigma = 1, theta = 1)
#' tictoc::toc()
#' @export
# Simulação de Monte-Carlo ------------------------------------------------
mc <- function(N = 1L,
               n = 50L,
               sig = 0.05,
               f,
               q,
               kicks,
               par0,
               ncores = 1L,
               p,
               bilateral = FALSE,
               step = 1e-3,
               ...) {
  c <- p * log(n) # Regra para escolhar de c.

  # Função densidade de probabilidade Qui-Quadrado inf.
  fdp_chisq_inf <- function(par, x) {
    k <- par[1]
    c <- par[2]

    indet <- 1 - pchisq(q = x, df = k)

    # Corrigindo problema de indeterminação que poderá ocorrer.
    indet[indet == 0] <- .Machine$double.eps

    dchisq(x = x, df = k) * (1 - (1 - pchisq(q = x, df = k)) ^ c + c * pchisq(q = x, df = k) *
                               (indet) ^ (c - 1))
  }

  quantile_chisq <- function(sig, bilateral = FALSE) {
    ifelse(
      bilateral == FALSE,
      result <- qchisq(p = 1 - sig, df =  length(par0[[1]])),
      result <-
        list(
          q1 = qchisq(p = sig / 2, df = length(par0[[1]])),
          q2 = qchisq(p = 1 - sig / 2, df = length(par0[[1]]))
        )
    )
    result
  }

  # Quantil da distribuição qui-quadrado.
  q_chisq <- quantile_chisq(sig, bilateral)

  # Quantil obtido da distribuição qui-quadrado inf.
  q_inf <-
    est_q(
      fn = fdp_chisq_inf,
      alpha = sig,
      bilateral = bilateral,
      step = step,
      c = c,
      k = length(par0[[1]])
    )

  mc_one_step <- function(i) {
    # Selecionando uma amostra que não gere erro nos chutes iniciais ----------
    repeat {
      amostra <- q(n, ...)
      result <- lrt(
        f = f,
        data = amostra,
        kicks = kicks,
        par0 = par0
      )
      if (!is.na(result))
        break
    }

    if (bilateral == TRUE) {
      sucess <- ifelse(result > q_chisq$q2 || result < q_chisq$q1, 1L, 0L)
      sucess_inf <-
        ifelse(result > q_inf$q2 || result < q_inf$q1, 1L, 0L)
      return(list(
        result = result,
        sucess = sucess,
        sucess_inf = sucess_inf
      ))
    } else {
      sucess <- ifelse(result > q_chisq, 1L, 0L)
      sucess_inf <- ifelse(result > q_inf, 1L, 0L)
      return(list(
        result = result,
        sucess = sucess,
        sucess_inf = sucess_inf
      ))
    }
  } # End mc_one_step().

  result_vector <-
    unlist(pbmcapply::pbmclapply(
      X = 1L:N,
      FUN = mc_one_step,
      mc.cores = ncores
    ))

  result <-
    as_tibble(matrix(result_vector, byrow = TRUE, ncol = 3L))

  names(result) <- c("lambda", "sucess_chisq", "sucess_inf")

  list(
    result = result,
    prop_chisq = mean(result$sucess_chisq),
    prop_inf = mean(result$sucess_inf)
  )
}

#' @title Encontra o quantil da distribuição Qui-Quadrado inf.
#' @author Pedro Rafael D. Marinho
#' @description  Encontra o quantil da distribuição Qui-Quadrado inf.
#' @details A função recebe como argumento uma função densidade de probabilidade ou uma função de probabilidade
#' que é passada como argumento à \code{f}. O objeto passado como argumento de \code{f} deverá ser implementado conforme os
#' exemplos. Note que se a função for implementada segundo os exemplos, não haverá a necessidade de implementar a função sob a
#' hipótese nula, visto que são as mesmas funções. Para especificar os parâmetros que serão fixados, i.e, para especificar
#' a distribuição sob a hipótese nula, utiliza-se o argumento \code{par0} que receberá uma lista formada por dois vetores.
#' O primeiro vetor da lista de verá ser um vetor de strings com os nomes das variáveis que deseja-se fixar e o segundo vetor
#' deverá conter os valores que serão atribuidos à cada uma das variáveis passada ao primeiro vetor.
#' @param fn Recebe a distribuição Qui-Quadrado inf.
#' @param alpha Nível de significância adotado.
#' @param bilateral Se \code{TRUE}, retorna os quantis para um teste bilateral. O padrão considera
#' \code{bilateral  = FALSE}.
#' @param step Tamanho do passo da integral numeérica responsável pela obtenção dos quantis da Qui-Quadrado inf.
#' O padrão considera \code{step = 1e-4}.
#' @param c Parâmetro da distribuição Qui-Quadrado inf.
#' @param k Parâmetro da distribuição Qui-Quadrado inf.
#' @importFrom stats integrate
#' @examples
#' fdp_chisq_inf <- function(par, x) {
#'   k <- par[1]
#'   c <- par[2]
#'   dchisq(x = x, df = k) * (1 - (1 - pchisq(q = x, df = k)) ^ c + c * pchisq(q = x, df = k) *
#'                         (1 - pchisq(q = x, df = k)) ^ (c - 1))
#' }
#' est_q(fn = fdp_chisq_inf, alpha = 0.05, bilateral = FALSE, c = 1, k = 1)
#' @export
est_q <- function(fn,
                  alpha = 0.05,
                  bilateral = FALSE,
                  step = 1e-4,
                  c,
                  k) {
  seq_q <- seq(from = step, to = 1e3L, by = step)

  test_q1 <- function(q) {
    integrate(
      f = fn,
      lower = 0,
      upper = q,
      par = c(k, c)
    )$value
  }

  test_q2 <- function(q) {
    integrate(
      f = fn,
      lower = q,
      upper = Inf,
      par = c(k, c)
    )$value
  }

  if (bilateral == TRUE) {
    for (i in  seq_q) {
      q1 <- test_q1(i)
      if (q1 >= alpha / 2)
        break
    }

    for (j in  seq_q) {
      q2 <- test_q2(j)
      if (q2 <= alpha / 2)
        break
    }

    return(list(q1 = i, q2 = j))

  } else {
    for (j in  seq_q) {
      q2 <- test_q2(j)
      if (q2 <= alpha)
        break
    }
    return(j)
  }
}

bootstraping <- function(B = 250L, f, kicks, par0, sample_true) {
  # A função resample() obtem uma amostra com reposição de x:
  # resample <- function(x) {
  #   n <- length(x)
  #   # Selecionando observações uniformemente distribuídas em x.
  #   # Poderia ser utilizado a função sample().
  #   x[floor(n * runif(n = n, min = 0, max = 1) + 1L)]
  # }

  # A função boot() calcula uma statística em uma única amostra
  # bootstrap:

  vec_t_boot <- double(length = B)

  for (b in 1L:B) {
    repeat {
      t_boot <- lrt(
        f,
        data = sample(
          x = sample_true,
          size = length(sample_true),
          replace = TRUE
        ),
        kicks = kicks,
        par0 = par0
      )
      if (!is.na(t_boot))
        break

    }

    vec_t_boot[b] <- t_boot
  }

  vec_t_boot

  # purrr::map_dbl(.x = 1L:B,
  #                .f = boot,
  #                sample_true = sample_true)
}

#' @title Função para o cálculo do poder do teste.
#' @param N Número de réplicas de Monte-Carlo.
#' @param B Número de réplicas de bootstrap.
#' @param n Tamanho da amostra.
#' @param sig Nível de significância.
#' @param q Função quantílica para a geração da amostra.
#' @param kicks Chutes iniciais.
#' @param par0 Lista de parâmetros para serem fixados na hipótese nula.
#' @param ncores Quantidade de cores utilizados na paralelização.
#' @param p Valor do parâmetro \strong{p}.
#' @param step Passos para das integrais calculadas para obtenção do quantil da distribuiçao Chisq-Inf.
#' @param ... Parâmetros para geração da amostra verdadeira.
#' @examples
#'# pdf_ew <- function(par, x, var = NULL) {
#'#   alpha <- par[1]
#'#   sigma <- par[2]
#'#   theta <- par[3]
#'#
#'#   if (is.list(var))
#'#     eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
#'#
#'#   alpha * theta / sigma * (1 - exp(-(x / sigma) ^ alpha)) ^ (theta - 1) *
#'#     exp(-(x / sigma) ^ alpha) * (x / sigma) ^ (alpha - 1)
#'# }
#'#
#'# rew <- function(n, alpha, sigma, theta){
#'#   u <- runif(n, 0, 1)
#'#   sigma * (-log(1 - u ^ (1 / theta))) ^ (1 / alpha)
#'# }
#'#
#'# set.seed(1L, kind = "L'Ecuyer-CMRG")
#'#
#'# tictoc::tic()
#'# power_test(N = 1e4L,
#'#            B = 250L,
#'#            n = 150L,
#'#            f = pdf_ew,
#'#            sig = 0.05,
#'#            q = rew,
#'#            kicks = c(1, 1, 1),
#'#            par0 = list(c("alpha", "theta"), c(1.5, 1.7)),
#'#            ncores = 8L,
#'#            p = 0.5,
#'#            step = 1e-3,
#'#            alpha = 1.5,
#'#            sigma = 1.5,
#'#            theta = 1.7
#'#
#'# )
#'# tictoc::toc()
#' @export
power_test <-
  function(N = 1e3L,
           B = 250L,
           n,
           f,
           sig,
           q,
           kicks,
           par0,
           ncores,
           p = 0.5,
           step = 1e-3,
           ...) {
    fdp_chisq_inf <- function(par, x) {
      k <- par[1]
      c <- par[2]
      dchisq(x = x, df = k) * (1 - (1 - pchisq(q = x, df = k)) ^ c + c * pchisq(q = x, df = k) *
                                 (1 - pchisq(q = x, df = k)) ^ (c - 1))
    }

    lrt_boot <- function(x) {
      lrt(
        f = f,
        data = x,
        kicks = kicks,
        par0 = par0
      )
    }

    # Um passo de Monte-Carlo
    mc_one_step <- function(i) {
      repeat {
        sample_true <- q(n = n, ...)
        t <- lrt(f,
                 data = sample_true,
                 kicks = kicks,
                 par0 = par0)

        if (!is.na(t))
          break

      }

      q0 <- function(q, n, par0, ...) {
        ndot <- ...length()

        var <-
          eval(parse(text = paste0(
            "names(", as.character(deparse(substitute(list(
              ...
            )))), ")"
          )))

        for (i in 1:ndot) {
          eval(parse(text = paste0(var[i], " <- ", ...elt(i))))
        }

        eval(parse(text = paste(par0[[1]], "<-", par0[[2]])))
        return(eval(body(q)))
      }

      sample_0 <- q0(q, n, par0, ...)

      q_boot <-
        quantile(bootstraping(
          B = B,
          f = f,
          kicks = kicks,
          par0 = par0,
          sample_true = sample_0
        ),
        prob = 1 - sig)

      sucess_theoretical <- ifelse(t > q_boot, 1L, 0L)
      sucess_chisc <-
        ifelse(t > qchisq(p = 1 - sig, df = length(par0[[1]])), 1L, 0L)

      sucess_chisc_inf <-
        ifelse(
          t > est_q(
            fdp_chisq_inf,
            alpha = sig,
            bilateral = FALSE,
            step = step,
            c = p * log(n),
            k = length(par0[[1]])
          ),
          1L,
          0L
        )

      c(sucess_theoretical, sucess_chisc, sucess_chisc_inf)

    }

    result_list <-
      unlist(pbmcapply::pbmclapply(
        X = 1L:N,
        FUN = mc_one_step,
        mc.cores = ncores
      ))

    matrix_result <- matrix(result_list, ncol = 3L, byrow = TRUE)
    colnames(matrix_result) <-
      c("theoretical", "chisq", "chisq_inf")
    apply(X = matrix_result, MARGIN = 2L, FUN = mean)

}
