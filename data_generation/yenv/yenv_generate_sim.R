
set_norm_keep_mean0 <- function(A, target_norm) {
  if (target_norm < 0) stop("target_norm must be >= 0")
  A0 <- A - mean(A)                 # 평균 0으로 센터링
  n0 <- sqrt(sum(A0^2))             # 현재 Frobenius 노름
  if (n0 == 0) {
    stop("Centered matrix has zero norm; cannot rescale.")
  }
  s  <- target_norm / n0            # 스케일 계수
  X  <- s * A0                       # 평균 0, 노름 target_norm
  dim(X) <- dim(A)
  X
}

response_generate_par <- function(r, p, u, Anorm = NULL) {

  mu_tru <- runif(r, 0, 10)
  eta_tru <- matrix(runif(u * p, min = 0, max = 10), nrow = u, ncol = p)
  if (is.null(Anorm)){
    A <- A_tru <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)
  }else{
    A <- A_tru <- set_norm_keep_mean0(
      matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u),
      Anorm
    )
  }

  CA <- rbind(diag(1, u), A)
  DA <- rbind(-t(A), diag(1, r-u))

  CAtCA <- t(CA) %*% CA
  Gamma <- CA %*% sqrtmatinv(CAtCA)

  DAtDA <- t(DA) %*% DA
  Gamma0 <- DA %*% sqrtmatinv(DAtDA)

  Omega_tru <- diag(sort(runif(u, 0, 1), decreasing = TRUE),
                    ncol = u, nrow = u)
  Omega0_tru <- diag(sort(runif(r-u, 5, 10), decreasing = TRUE),
                     ncol = r-u, nrow = r-u)

  Sigma1 <- Gamma %*% Omega_tru %*% t(Gamma)
  Sigma2 <- Gamma0 %*% Omega0_tru %*% t(Gamma0)

  Sigma_tru <- Sigma1 + Sigma2

  # reparameterization

  eta_tilde_tru <- sqrtmat(CAtCA) %*% eta_tru
  Omega_tilde_tru <- sqrtmat(CAtCA) %*% Omega_tru %*% sqrtmat(CAtCA)
  Omega0_tilde_tru <- sqrtmat(DAtDA) %*% Omega0_tru %*% sqrtmat(DAtDA)

  # calculate beta

  beta_tru <- CA %*% solve(CAtCA) %*% eta_tilde_tru

  list(
    mu_tru = mu_tru,
    beta_tru = beta_tru,
    eta_tilde_tru = eta_tilde_tru,
    Omega_tilde_tru = Omega_tilde_tru,
    Omega0_tilde_tru = Omega0_tilde_tru,
    A_tru = A_tru,
    Gamma_tru = Gamma,
    Gamma0_tru = Gamma0,
    Sigma_tru = Sigma_tru
  )
}

response_generate_data <- function(mu_tru, beta_tru,
                          Sigma_tru, mux=3, sigmax=1,
                          n, r, p, ...){

  X <- matrix(rnorm(n*p, mean = mux, sd = sigmax),
              nrow = n, ncol = p)
  eps <- matrix(rnorm(n*r), n, r) %*% sqrtmat(Sigma_tru)
  Y <- tcrossprod(rep(1, n), mu_tru) + X %*% t(beta_tru) + eps

  list(X = X, Y = Y)
}

# set.seed(100)
# n <- 1000
# r <- 20
# p <- 7
# u_true <- 5
# 
# all_pars <- response_generate_par(r, p, u_true)
# all_pars$beta_tru
# 
# generate_multiple_datasets <- function(n, n_reps, all_pars) {
#   replicate(n_reps, {
#     do.call(response_generate_data, c(all_pars, list(n = n, r = r, p = p)))
#   }, simplify = FALSE)
# }
# 
# datasets <- generate_multiple_datasets(n, 100, all_pars)
# saveRDS(datasets, file = paste0("/home/jupyter-tmdgus4970/2025_CALVI/simulation/datasets_n",n,"_u",u_true,".rds"))
