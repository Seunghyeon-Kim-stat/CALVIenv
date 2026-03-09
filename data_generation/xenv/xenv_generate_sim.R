
library(LaplacesDemon)

predictor_generate_par <- function(r, p, m) {

  mu_X_tru <- runif(p, -10, 10)
  mu_Y_tru <- runif(r, -10, 10)
  eta_tru <- matrix(runif(m * r, min = 5, max = 10), nrow = m, ncol = r)
  A <- A_tru <- matrix(runif(m * (p - m), min = -1, max = 1), nrow = p-m, ncol = m)
  CA <- rbind(diag(1, m), A)
  DA <- rbind(-t(A), diag(1, p-m))

  CAtCA <- t(CA) %*% CA
  Gamma <- CA %*% sqrtmatinv(CAtCA)

  DAtDA <- t(DA) %*% DA
  Gamma0 <- DA %*% sqrtmatinv(DAtDA)

  Omega_tru <- 5*rinvwishart(nu = m + 2, S = diag(50, m))
  Omega0_tru <- rinvwishart(nu = p - m + 2, S = diag(0.1, p-m))

  Sigma_Y_X_tru <- rinvwishart(nu = r + 1, S = diag(5, r))

  Sigma_X1 <- Gamma %*% Omega_tru %*% t(Gamma)
  Sigma_X2 <- Gamma0 %*% Omega0_tru %*% t(Gamma0)

  Sigma_X_tru <- Sigma_X1 + Sigma_X2

  # reparameterization

  eta_tilde_tru <- sqrtmatinv(CAtCA) %*% eta_tru
  Omega_tilde_tru <- sqrtmat(CAtCA) %*% Omega_tru %*% sqrtmat(CAtCA)
  Omega0_tilde_tru <- sqrtmat(DAtDA) %*% Omega0_tru %*% sqrtmat(DAtDA)

  # calculate beta

  beta_tru <- t(CA %*% eta_tilde_tru)

  list(
    mu_X_tru = mu_X_tru,
    mu_Y_tru = mu_Y_tru,
    beta_tru = beta_tru,
    eta_tilde_tru = eta_tilde_tru,
    Omega_tilde_tru = Omega_tilde_tru,
    Omega0_tilde_tru = Omega0_tilde_tru,
    A_tru = A_tru,
    Gamma_tru = Gamma,
    Gamma0_tru = Gamma0,
    Sigma_X_tru = Sigma_X_tru,
    Sigma_Y_X_tru = Sigma_Y_X_tru
  )
}

predictor_generate_data <- function(mu_X_tru,
                                    mu_Y_tru,
                                    beta_tru,
                                    Sigma_X_tru,
                                    Sigma_Y_X_tru,
                                    n, r, p, ...){

  X_muX <- matrix(rnorm(n*p), n, p) %*% sqrtmat(Sigma_X_tru)
  X <-  X_muX + tcrossprod(rep(1, n), mu_X_tru)

  Y <- matrix(rnorm(n*r), n, r) %*% sqrtmat(Sigma_Y_X_tru) +
    X_muX %*% t(beta_tru) +
    tcrossprod(rep(1, n), mu_Y_tru)

  list(X = X, Y = Y)
}

# n = 100
# r = 3
# p = 15
# m_true = 5
# 
# set.seed(100)
# all_pars <- predictor_generate_par(r, p, m_true)
# all_pars$A_tru
# 
# generate_multiple_datasets <- function(n_reps, all_pars) {
#   replicate(n_reps, {
#     do.call(predictor_generate_data, c(all_pars, list(n = n, r = r, p = p)))
#   }, simplify = FALSE)
# }
# 
# datasets <- generate_multiple_datasets(100, all_pars)
# 
# saveRDS(datasets, file = paste0("simulation/datasets_n",n,"_m",m_true,".rds"))

# datasets <- readRDS(file = paste0("D:/Dropbox/Project/LaplaceVI/Simulation/predictor datasets/datasets_n",50,"_m",m_true,".rds"))
