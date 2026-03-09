
#install.packages("Renvlp")
library(Renvlp)

response_get_prior <- function(r, p, u, nuvalue, nu0value, Avalue, Mvalue){

  #mu0 <- matrix(0, nrow = r, ncol = 1)
  #Sigma0 <- 1e6 * diag(1, nrow = r)

  nu <- u + nuvalue
  psi <- 1e-06

  nu0 <- (r - u) + nu0value
  psi0 <- 1e-06

  B0 <- matrix(0, nrow = r, ncol = p)

  M <- 1 / Mvalue * diag(1, nrow = p)

  A0 <- matrix(0, nrow = r-u, ncol = u)
  U_A <- Avalue * diag(1, nrow = r-u)
  V_A <- Avalue * diag(1, nrow = u)

  return(list(nu = nu,
              psi = psi,
              nu0 = nu0,
              psi0 = psi0,
              B0 = B0,
              M = M,
              A0 = A0,
              U_A = U_A,
              V_A = V_A))
}

response_get_init <- function(n, r, p, u){

  hatA <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)
  HA <- diag(1e-6, u*(r-u), u*(r-u))

  CA <- rbind(diag(1, u), hatA)
  DA <- rbind(-t(hatA), diag(1, r-u))

  CAtCA <- t(CA) %*% CA
  Gamma <- CA %*% sqrtmatinv(CAtCA)

  DAtDA <- t(DA) %*% DA
  Gamma0 <- DA %*% sqrtmatinv(DAtDA)

  nu_q <- n+p+u
  Psi_q <- diag(sort((nu_q-u-1)*runif(u, 0, 1), decreasing = TRUE),
                ncol = u, nrow = u)
  invPsi_q <- solve(Psi_q)

  nu0_q <- n+r-u
  Psi0_q <- diag(sort((nu0_q-(r-u)-1)*runif(r-u, 5, 10), decreasing = TRUE),
                 ncol = r-u, nrow = r-u)
  invPsi0_q <- solve(Psi0_q)

  eta_q <- sqrtmat(CAtCA)%*%matrix(runif(u * p, min = 0, max = 10), nrow = u, ncol = p)
  U_q <- 1e-6 * diag(1, nrow = u)
  V_q <- 1e-6 * diag(1, nrow = p)

  mu_q <- runif(r, 0, 10)
  Sigma_q <- Gamma %*% rinvwishart(nu_q, Psi_q) %*% t(Gamma) +
    Gamma0 %*% rinvwishart(nu0_q, Psi0_q) %*% t(Gamma0)

  CA_start <- rbind(diag(u), hatA)
  DA_start <- rbind(-t(hatA), diag(r-u))

  return(list(hatA = hatA, CA = CA_start, DA = DA_start, HA = HA,
              nu_q = nu_q, invPsi_q = invPsi_q,
              nu0_q = nu0_q, invPsi0_q = invPsi0_q,
              eta_q = eta_q, U_q = U_q, V_q = V_q,
              mu_q = mu_q, Sigma_q = Sigma_q))
}

response_get_init_MLE <- function(X, Y, n, r, p, u){

  Y <- data.matrix(Y)
  tYY <- crossprod(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  YctYc <- crossprod(Yc)
  Y_bar <- attr(Yc,"scaled:center")

  X <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X_bar <- attr(Xc,"scaled:center")

  if (u == 0){

    m <- env(X, Y, u)

    gamma <- m$Gamma

    # A MLE
    A_start <- 0
    CA_start <- 0
    DA_start <- diag(1, r)
    HA <- 0

    DAtDA_start <- t(DA_start) %*% DA_start

    gamma_start <- m$Gamma
    gamma0_start <- m$Gamma0

    # Omega MLE

    tilde_Omega <- 0
    nu_q <- n+p+u
    invPsi_q <- 0

    # eta MLE
    eta_start <- 0

    U_q <- 0
    V_q <- 0

    # Omega0 MLE

    Sigma2 <- m$Gamma0 %*% m$Omega0 %*% t(m$Gamma0)
    Omega0 <- t(gamma0_start) %*% Sigma2 %*% gamma0_start
    tilde_Omega0 <- sqrtmat(DAtDA_start) %*% Omega0 %*% sqrtmat(DAtDA_start)
    nu0_q <- n+r-u
    Psi0_q <- (nu0_q-(r-u)-1) * tilde_Omega0
    invPsi0_q <- solve(Psi0_q)

    # mu MLE
    mu_start <- m$mu

    # Sigma MLE
    Sigma_start <- diag(1e-6, nrow = r)

    Sigma <- m$Sigma
  }
  else if (u < r){

    m <- env(X, Y, u)

    gamma <- m$Gamma
    G1 <- as.matrix(gamma[1:u, ])
    G2 <- gamma[-(1:u), ]

    A_start <- G2 %*% solve(G1)

    HA <- diag(1e-6, u*(r-u), u*(r-u))

    Z <- matrix(0, nrow=(r-u), ncol=u)
    K <- rbind(diag(1, u), Z)
    L <- rbind(t(Z), diag(1, r-u))

    CA_start <- K + L %*% A_start
    DA_start <- L - K %*% t(A_start)

    CAtCA_start <- t(CA_start) %*% CA_start
    DAtDA_start <- t(DA_start) %*% DA_start

    gamma_start <- CA_start %*% sqrtmatinv(CAtCA_start)
    gamma_t <- t(gamma_start)
    gamma0_start <- DA_start %*% sqrtmatinv(DAtDA_start)
    gamma0_t <- t(gamma0_start)

    # Omega MLE
    Sigma1 <- m$Gamma %*% m$Omega %*% t(m$Gamma)
    Omega <- t(gamma_start) %*% Sigma1 %*% gamma_start
    tilde_Omega <- sqrtmat(CAtCA_start) %*% Omega %*% sqrtmat(CAtCA_start)
    nu_q <- n+p+u
    Psi_q <- (nu_q-u-1) * tilde_Omega
    invPsi_q <- solve(Psi_q)

    # eta MLE
    eta_start <- sqrtmat(CAtCA_start) %*% gamma_t %*% m$beta

    U_q <- Psi_q / nu_q
    V_q <- solve(XctXc)

    # Omega0 MLE
    Sigma2 <- m$Gamma0 %*% m$Omega0 %*% t(m$Gamma0)
    Omega0 <- t(gamma0_start) %*% Sigma2 %*% gamma0_start
    tilde_Omega0 <- sqrtmat(DAtDA_start) %*% Omega0 %*% sqrtmat(DAtDA_start)
    nu0_q <- n+r-u
    Psi0_q <- (nu0_q-(r-u)-1) * tilde_Omega0
    invPsi0_q <- solve(Psi0_q)

    # mu MLE
    mu_start <- m$mu + CA_start%*%solve(CAtCA_start)%*%eta_start%*%X_bar

    # Sigma MLE
    Sigma_start <- diag(1e-6, nrow = r)

    Sigma <- m$Sigma

  }
  else if (u == r){

    m <- Renvlp::env(X, Y, u)

    gamma <- m$Gamma

    # A MLE
    A_start <- 0
    CA_start <- diag(1, r)
    DA_start <- 0
    HA <- 0

    CAtCA_start <- t(CA_start) %*% CA_start

    gamma_start <- gamma
    gamma0_start <- 0

    # Omega MLE

    Sigma1 <- m$Gamma %*% m$Omega %*% t(m$Gamma)
    Omega <- t(gamma_start) %*% Sigma1 %*% gamma_start
    tilde_Omega <- sqrtmat(CAtCA_start) %*% Omega %*% sqrtmat(CAtCA_start)
    nu_q <- n+p+u
    Psi_q <- (nu_q-u-1) * tilde_Omega
    invPsi_q <- solve(Psi_q)

    # eta MLE
    eta_start <- m$beta

    U_q <- Psi_q / nu_q
    V_q <- solve(XctXc)

    # Omega0 MLE

    tilde_Omega0 <- 0
    nu0_q <- n+r-u
    invPsi0_q <- 0

    # mu MLE
    mu_start <- m$mu + eta_start%*%X_bar

    # Sigma MLE
    Sigma_start <- diag(1e-6, nrow = r)

    Sigma <- m$Sigma

  }

  return(list(hatA = A_start, CA = CA_start, DA = DA_start,HA = HA,
              nu_q = nu_q, invPsi_q = invPsi_q,
              nu0_q = nu0_q, invPsi0_q = invPsi0_q,
              eta_q = eta_start, U_q = U_q, V_q = V_q,
              mu_q = mu_start, Sigma_q = Sigma_start,
              Sigma = Sigma))
}

add_gaussian_noise <- function(mat, sd = 10) {
  # check if symmetric pd using cholesky

  check_pd <- tryCatch(chol(mat), error = function(e) e)

  if (is.null(mat)) {
    NULL
  } else if (is(check_pd, "error")) {
    mat + rnorm(length(mat), 0, sd)
  } else {
    crossprod(check_pd + rnorm(length(check_pd), 0, sd))
  }
}

response_get_init_MLE_with_noise <- function(X, Y, n, r, p, u, sd = 1) {
  init <- response_get_init_MLE(X, Y, n, r, p, u)

  if (u > 0 && u < r) {
    noisy_hatA <- add_gaussian_noise(init$hatA, sd = sd)

    Z <- matrix(0, nrow = (r - u), ncol = u)
    K <- rbind(diag(1, u), Z)
    L <- rbind(t(Z), diag(1, r - u))

    CA_new <- K + L %*% noisy_hatA
    DA_new <- L - K %*% t(noisy_hatA)

    init$hatA <- noisy_hatA
    init$CA   <- CA_new
    init$DA   <- DA_new
  }

  init$eta_q   <- add_gaussian_noise(init$eta_q, sd = sd)
  init$mu_q    <- add_gaussian_noise(init$mu_q, sd = sd)

  return(init)
}


response_get_init_MAP <- function(X, Y, n, r, p, u,
                                  nu, psi, nu0, psi0, B0, M, A0, U_A, V_A,
                                  BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, ...){

  # --- data prep ---
  Y  <- data.matrix(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  YctYc <- crossprod(Yc)
  Y_bar <- attr(Yc,"scaled:center")

  X  <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X_bar <- attr(Xc,"scaled:center")

  # posterior mean of B under ridge prior (B0, M)
  hatB0 <- t(tXcY_MtB0) %*% inv_XctXc_M

  # \tilde{G} = Y'Y + B0 M B0' - \hat{B0}(X'X+M)\hat{B0}'
  tildeG <- YctYc + BMB - hatB0 %*% XctXc_M %*% t(hatB0)

  # handy zero-mats for degenerate dims
  zmat <- function(m,n) matrix(0, m, n)

  if (u == 0) {
    ## ---------- full reduct: u = 0 ----------
    # A: (r-u) x u = r x 0
    A_start  <- zmat(r, 0)
    CA_start <- zmat(r, 0)        # r x 0
    DA_start <- diag(1, r)        # r x r
    HA       <- zmat(0, 0)        # u*(r-u) = 0

    # Omega (u x u) — empty
    nu_q    <- n + p + nu
    invPsi_q <- zmat(0, 0)

    # Omega0 (r-u x r-u) = r x r
    Omega0   <- t(DA_start) %*% (YctYc + psi0 * diag(1, r)) %*% DA_start / (n + nu0 + (r - u) + 1)
    nu0_q    <- n + nu0
    Psi0_q   <- (nu0_q - (r - u) - 1) * Omega0
    invPsi0_q <- solve(Psi0_q)

    # eta: (u x p) = 0 x p
    eta_start <- zmat(0, p)
    U_q <- zmat(0, 0)
    V_q <- solve(XctXc_M)

    # mu, Sigma
    mu_start    <- Y_bar
    Sigma_start <- solve(DA_start %*% solve(Omega0) %*% t(DA_start))  # = Omega0

    return(list(hatA = A_start, CA = CA_start, DA = DA_start, HA = HA,
                nu_q = nu_q, invPsi_q = invPsi_q,
                nu0_q = nu0_q, invPsi0_q = invPsi0_q,
                eta_q = eta_start, U_q = U_q, V_q = V_q,
                mu_q = mu_start, Sigma_q = Sigma_start))
  }

  if (u == r) {
    ## ---------- full envelope: u = r ----------
    # A: (r-u) x u = 0 x r
    A_start  <- zmat(0, r)
    CA_start <- diag(1, r)        # r x r
    DA_start <- zmat(r, 0)        # r x 0
    HA       <- zmat(0, 0)

    # Omega (r x r)
    Omega   <- t(CA_start) %*% (tildeG + psi * diag(1, r)) %*% CA_start / (n + p + nu + u + 1)
    nu_q    <- n + p + nu
    Psi_q   <- (nu_q - u - 1) * Omega
    invPsi_q <- solve(Psi_q)

    # Omega0 — empty
    nu0_q     <- n + nu0
    invPsi0_q <- zmat(0, 0)

    # eta: r x p  (CA=I ⇒ eta = hatB0)
    eta_start <- hatB0
    U_q <- Omega
    V_q <- solve(XctXc_M)

    # mu, Sigma
    mu_start    <- Y_bar
    Sigma_start <- solve(CA_start %*% solve(Omega) %*% t(CA_start))  # = Omega

    return(list(hatA = A_start, CA = CA_start, DA = DA_start, HA = HA,
                nu_q = nu_q, invPsi_q = invPsi_q,
                nu0_q = nu0_q, invPsi0_q = invPsi0_q,
                eta_q = eta_start, U_q = U_q, V_q = V_q,
                mu_q = mu_start, Sigma_q = Sigma_start))
  }

  ## ---------- interior case: 0 < u < r ----------
  m <- MAPCA(X=X, Y=Y, u=u, A0=A0, K=U_A, L=V_A,
             Psi=psi * diag(1, u), Psi0=psi0 * diag(1, r - u),
             e=B0, M=M, nu=nu, nu0=nu0,
             maxiter = 100, method = "BFGS", ...)

  CA <- m$CA
  gamma <- CA %*% sqrtmatinv(crossprod(CA))
  A <- find_A_from_gamma(gamma)

  CA_start <- rbind(diag(1, u), A)
  DA_start <- rbind(-t(A), diag(1, r - u))
  HA <- diag(1e-6, u * (r - u))

  # Omega (u x u)
  Omega   <- t(CA_start) %*% (tildeG + psi * diag(1, r)) %*% CA_start / (n + p + nu + u + 1)
  nu_q    <- n + p + nu
  Psi_q   <- (nu_q - u - 1) * Omega
  invPsi_q <- solve(Psi_q)

  # Omega0 (r-u x r-u)
  Omega0   <- t(DA_start) %*% (YctYc + psi0 * diag(1, r)) %*% DA_start / (n + nu0 + (r - u) + 1)
  nu0_q    <- n + nu0
  Psi0_q   <- (nu0_q - (r - u) - 1) * Omega0
  invPsi0_q <- solve(Psi0_q)

  # eta: u x p
  eta_start <- t(CA_start) %*% hatB0

  U_q <- Omega
  V_q <- solve(XctXc_M)

  mu_start    <- Y_bar
  Sigma_start <- solve(CA_start %*% solve(Omega) %*% t(CA_start) +
                         DA_start %*% solve(Omega0) %*% t(DA_start))

  beta <- CA_start %*% solve(crossprod(CA_start)) %*% eta_start

  return(list(hatA = A, CA = CA_start, DA = DA_start, HA = HA,
              nu_q = nu_q, invPsi_q = invPsi_q,
              nu0_q = nu0_q, invPsi0_q = invPsi0_q,
              eta_q = eta_start, U_q = U_q, V_q = V_q,
              mu_q = mu_start, Sigma_q = Sigma_start,
              beta = beta))
}
