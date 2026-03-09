
elmwise_mean_in_list <- function(lst) {
  n_lst <- length(lst)
  Reduce("+", lst) / n_lst
}

elmwise_sd_in_list <- function(lst) {
  lst_orig <- lst
  lst <- lst[!sapply(lst, is.null)]
  n_lst <- length(lst)
  xbar <- elmwise_mean_in_list(lst)
  x2bar <- elmwise_mean_in_list(
    lapply(lst, "^", 2)
  )
  # browser()
  sqrt(pmax(matrix(0, nrow(xbar), ncol(xbar)), (x2bar - xbar^2)) * n_lst / (n_lst - 1))
}



autotune_param <- function(draw_param_list,
                           tune_param,
                           accpt_list,
                           tune_nterm,
                           tune.incr,
                           tune.accpt.prop.lower,
                           tune.accpt.prop.upper,
                           adjust_scale,
                           ...) {
  # browser()

  n_full <- draw_param_list %>%
    .[!sapply(., is.null)] %>%
    length()

  ave_accpt_last_nterm <- accpt_list[max(1, n_full - tune_nterm + 1):n_full] %>%
    elmwise_mean_in_list()

  idx_increase <- which(ave_accpt_last_nterm > tune.accpt.prop.upper)
  idx_decrease <- which(ave_accpt_last_nterm < tune.accpt.prop.lower)

  if (length(idx_increase) > 0) {
    tune_param[idx_increase] <- tune_param[idx_increase] * (1 + tune.incr)
  }

  if (length(idx_decrease) > 0) {
    tune_param[idx_decrease] <- tune_param[idx_decrease] * (1 - tune.incr)
  }


  if (adjust_scale) {
    # last nterms iterations of draw_param
    # each column of draw_param_j_last_nterm is
    # from draw_param single interation of draw_param
    par_sd <- elmwise_sd_in_list(
      draw_param_list[(n_full - tune_nterm + 1):n_full]
    )
    par_sd_mean <- mean(par_sd)
    tune_param_mean <- mean(tune_param)
    if (par_sd_mean > 0) {
      tune_param <-
        0.5 * par_sd / par_sd_mean * tune_param_mean +
        0.5 * tune_param
    }
  }

  tune_param
}


autotune_SCAM <- function(draw_param_list,
                          epsilon = 0.05,
                          prev_var = NULL,
                          prev_mean = NULL,
                          ...) {
  if (is.null(prev_mean) | is.null(prev_var)) {
    curr_mean <- elmwise_mean_in_list(draw_param_list)
    curr_sd <- elmwise_sd_in_list(draw_param_list)
    curr_var <- curr_sd^2
  } else {
    # browser()
    n <- length(draw_param_list)
    curr_mean <- (n - 1) / n * prev_mean + 1 / n * draw_param_list[[n]]
    curr_var <- (n - 2) / (n - 1) * prev_var + prev_mean^2 +
      1 / (n - 1) * draw_param_list[[n]]^2 - n / (n - 1) * curr_mean^2
    curr_sd <- sqrt(curr_var)
  }
  tau <- curr_sd * 2.38
  attr(tau, "curr_mean") <- curr_mean
  attr(tau, "curr_var") <- curr_var

  tau
}


rwmh_colwise_A <- function(A_start,
                           lpd_start,
                           lpd_func,
                           tau,
                           random_scan = TRUE,
                           elementwise = FALSE,
                           random_update_fraction = 1,
                           n_update_max = 100,
                           n_rwmh_iter = 1,
                           ...) {
  dims <- dim(A_start)
  u <- dims[2]
  r <- sum(dims)
  A <- A_start

  count <- 0



  for (iter in 1:n_rwmh_iter) {
    accpt.A <- matrix(0, dims[1], dims[2])
    lpd_A <- lpd_func(A = A, ...)

    if (!elementwise) {
      if (random_scan) {
        one_to_u <- sample(1:u)
      } else {
        one_to_u <- 1:u
      }

      for (j in one_to_u) {
        count <- count + 1
        if (count > n_update_max) break
        # if (autotune_size == "single") tau_curr <- tau
        tau_curr <- tau[, j]
        A_j_star <- A[, j] + rnorm(r - u, 0, tau_curr)
        A_star <- A
        A_star[, j] <- A_j_star
        lpd_A_star <- lpd_func(A = A_star, ...)

        if (log(runif(1)) < (c(lpd_A_star) - c(lpd_A))) {
          # accept
          A <- A_star
          lpd_A <- lpd_A_star
          accpt.A[, j] <- 1
        }
      }

      # browser()
    } else if (elementwise) {
      n_full <- length(A)
      n_update <- (n_full * random_update_fraction) %>%
        ceiling() %>%
        pmin(n_full)


      update_set <- if (n_update < n_full) {
        sample(1:n_full, n_update)
      } else if (random_scan) {
        sample(1:n_full)
      } else {
        1:n_full
      }


      # browser()

      if (n_update < n_full) {
        accpt.A[-(update_set)] <- 0
      }

      for (j in update_set) {
        count <- count + 1
        if (count > n_update_max) break
        # browser()
        # if (autotune_size == "single") tau_curr <- tau
        tau_curr <- tau[j]
        A_j_star <- A[j] + rnorm(1, 0, tau_curr)
        A_star <- A
        A_star[j] <- A_j_star
        lpd_A_star <- lpd_func(A = A_star, ...)

        if (log(runif(1)) < (c(lpd_A_star) - c(lpd_A))) {
          # accept
          A <- A_star
          lpd_A <- lpd_A_star
          accpt.A[j] <- 1
        }
      }
    }
  }

  # browser()

  out <- list(
    A = A,
    lpd = lpd_A,
    accpt = accpt.A
  )

  out
}


xenv_get_prior <- function(r, p, m){

  # muX0 <- 0
  # SigmaX0 <- 1e+06
  #
  # muY0 <- 0
  # SigmaY0 <- 1e+06

  nuY <- r
  psiY <- 1e-06

  nuX1 <- m
  psiX1 <- 1e-06

  nuX0 <- p-m
  psiX0 <- 1e-06

  A0 <- matrix(0, nrow = (p-m), ncol = m)
  U_A <- 1e6 * diag(1, nrow = (p-m))
  V_A <- 1e6 * diag(1, nrow = m)

  psi_eta <- 1e10
  B0 <- matrix(0, nrow = p, ncol = r)
  M <- 1 / 1e6 * diag(1, nrow = r)

  return(list(nuY = nuY,
              psiY = psiY,
              nuX1 = nuX1,
              psiX1 = psiX1,
              nuX0 = nuX0,
              psiX0 = psiX0,
              A0 = A0,
              U_A = U_A,
              V_A = V_A,
              psi_eta = psi_eta,
              B0 = B0,
              M = M))
}

# get_init <- function(n, r, p, u){
#
#   hatA <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)
#   HA <- diag(1e-6, u*(r-u), u*(r-u))
#
#   CA <- rbind(diag(1, u), hatA)
#   DA <- rbind(-t(hatA), diag(1, r-u))
#
#   CAtCA <- t(CA) %*% CA
#   Gamma <- CA %*% sqrtmatinv(CAtCA)
#
#   DAtDA <- t(DA) %*% DA
#   Gamma0 <- DA %*% sqrtmatinv(DAtDA)
#
#   nu_q <- n+p+u
#   Psi_q <- diag(sort((nu_q-u-1)*runif(u, 0, 1), decreasing = TRUE),
#                 ncol = u, nrow = u)
#   invPsi_q <- solve(Psi_q)
#
#   nu0_q <- n+r-u
#   Psi0_q <- diag(sort((nu0_q-(r-u)-1)*runif(r-u, 5, 10), decreasing = TRUE),
#                  ncol = r-u, nrow = r-u)
#   invPsi0_q <- solve(Psi0_q)
#
#   eta_q <- sqrtmat(CAtCA)%*%matrix(runif(u * p, min = 0, max = 10), nrow = u, ncol = p)
#   U_q <- 1e-6 * diag(1, nrow = u)
#   V_q <- 1e-6 * diag(1, nrow = p)
#
#   mu_q <- runif(r, 0, 10)
#   Sigma_q <- Gamma %*% rinvwishart(nu_q, Psi_q) %*% t(Gamma) +
#     Gamma0 %*% rinvwishart(nu0_q, Psi0_q) %*% t(Gamma0)
#
#   return(list(hatA = hatA, HA = HA,
#               nu_q = nu_q, invPsi_q = invPsi_q,
#               nu0_q = nu0_q, invPsi0_q = invPsi0_q,
#               eta_q = eta_q, U_q = U_q, V_q = V_q,
#               mu_q = mu_q, Sigma_q = Sigma_q))
# }

xenv_get_init_MLE <- function(X, Y, n, r, p, m){

  # library(Renvlp)

  Y <- data.matrix(Y)
  tYY <- crossprod(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  YctYc <- crossprod(Yc)
  Y_bar <- attr(Yc,"scaled:center")

  X <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X_bar <- attr(Xc,"scaled:center")

  if (m ==0){

    xm <- Renvlp::xenv(X, Y, m)

    A_start <- matrix(0, (p-m), m)
    HA <- diag(0, (p-m)*m, (p-m)*m)
    CA_start <- rbind(diag(1, m), A_start)
    DA_start <- rbind(-t(A_start), diag(1, p-m))

    # SigmaYcX
    SigmaY <- xm$SigmaYcX
    nuY_q <- n + r
    PsiY_q <- (nuY_q-r-1) * SigmaY
    invPsiY_q <- solve(PsiY_q)

    # eta MLE
    eta_start <- matrix(0, m, r)
    U_q <- 1e-2 * diag(1, nrow = m)
    V_q <- 1e-2 * diag(1, nrow = r)
    Sigmaeta_q <- kronecker(V_q, U_q)

    # Omega MLE
    # SigmaX1 <- matrix(0, p, p)
    # tilde_OmegaX1 <- t(CA_start) %*% SigmaX1 %*% CA_start
    nuX1_q <- n + m - r
    # PsiX1_q <- (nuX1_q-m-1) * tilde_OmegaX1
    invPsiX1_q <- matrix(0, m, m)

    # Omega0 MLE
    SigmaX2 <- xm$Gamma0 %*% xm$Omega0 %*% t(xm$Gamma0)
    tilde_OmegaX0 <- t(DA_start) %*% SigmaX2 %*% DA_start
    nuX0_q <- n + p - m
    PsiX0_q <- (nuX0_q-(p-m)-1) * tilde_OmegaX0
    invPsiX0_q <- solve(PsiX0_q)

    # muX MLE
    muX_start <- X_bar
    SigmaX_start <- 1e-4 * diag(1, nrow = p)

    # muY MLE
    muY_start <- xm$mu + t(eta_start) %*% t(CA_start) %*% muX_start
    SigmaY_start <- 1e-4 * diag(1, nrow = r)

    # Sigma MLE
    Sigma_start <- xm$Sigma

  }
  else if(m < p){

    xm <- Renvlp::xenv(X, Y, m)
    gamma <- xm$Gamma
    G1 <- as.matrix(gamma[1:m, ])
    G2 <- gamma[-(1:m), ]

    # A MLE
    A_start <- G2 %*% solve(G1)
    HA <- diag(1e-6, (p-m)*m, (p-m)*m)

    CA_start <- rbind(diag(1, m), A_start)
    DA_start <- rbind(-t(A_start), diag(1, p-m))
    CAtCA_start <- t(CA_start) %*% CA_start
    DAtDA_start <- t(DA_start) %*% DA_start

    CAtCA_start_sqrtinv <- sqrtmatinv(CAtCA_start)

    gamma_start <- CA_start %*% CAtCA_start_sqrtinv
    gamma_t <- t(gamma_start)

    # SigmaYcX
    SigmaY <- xm$SigmaYcX
    nuY_q <- n + r
    PsiY_q <- (nuY_q-r-1) * SigmaY
    invPsiY_q <- solve(PsiY_q)

    # eta MLE
    eta_start <- CAtCA_start_sqrtinv %*% gamma_t %*% xm$beta
    U_q <- CAtCA_start_sqrtinv %*% (1e-2 * diag(1, nrow = m)) %*% CAtCA_start_sqrtinv
    V_q <- 1e-2 * diag(1, nrow = r)
    Sigmaeta_q <- kronecker(V_q, U_q)

    # Omega MLE
    SigmaX1 <- xm$Gamma %*% xm$Omega %*% t(xm$Gamma)
    tilde_OmegaX1 <- t(CA_start) %*% SigmaX1 %*% CA_start
    nuX1_q <- n + m - r
    PsiX1_q <- (nuX1_q-m-1) * tilde_OmegaX1
    invPsiX1_q <- solve(PsiX1_q)

    # Omega0 MLE
    SigmaX2 <- xm$Gamma0 %*% xm$Omega0 %*% t(xm$Gamma0)
    tilde_OmegaX0 <- t(DA_start) %*% SigmaX2 %*% DA_start
    nuX0_q <- n + p - m
    PsiX0_q <- (nuX0_q-(p-m)-1) * tilde_OmegaX0
    invPsiX0_q <- solve(PsiX0_q)

    # muX MLE
    muX_start <- X_bar
    SigmaX_start <- 1e-4 * diag(1, nrow = p)

    # muY MLE
    muY_start <- xm$mu + t(eta_start) %*% t(CA_start) %*% muX_start
    SigmaY_start <- 1e-4 * diag(1, nrow = r)

    # Sigma MLE
    Sigma_start <- xm$Sigma

  }
  else if(m == p){

    xm <- Renvlp::xenv(X, Y, m)
    gamma <- xm$Gamma
    G1 <- as.matrix(gamma[1:m, ])
    G2 <- gamma[-(1:m), ]

    # A MLE
    A_start <- G2 %*% solve(G1)
    HA <- diag(0, (p-m)*m, (p-m)*m)

    CA_start <- rbind(diag(1, m), A_start)
    DA_start <- rbind(-t(A_start), diag(1, p-m))
    CAtCA_start <- t(CA_start) %*% CA_start
    DAtDA_start <- t(DA_start) %*% DA_start

    CAtCA_start_sqrtinv <- sqrtmatinv(CAtCA_start)

    gamma_start <- CA_start %*% CAtCA_start_sqrtinv
    gamma_t <- t(gamma_start)

    # SigmaYcX
    SigmaY <- xm$SigmaYcX
    nuY_q <- n + r
    PsiY_q <- (nuY_q-r-1) * SigmaY
    invPsiY_q <- solve(PsiY_q)

    # eta MLE
    eta_start <- CAtCA_start_sqrtinv %*% gamma_t %*% xm$beta
    U_q <- CAtCA_start_sqrtinv %*% (1e-2 * diag(1, nrow = m)) %*% CAtCA_start_sqrtinv
    V_q <- 1e-2 * diag(1, nrow = r)
    Sigmaeta_q <- kronecker(V_q, U_q)

    # Omega MLE
    SigmaX1 <- xm$Gamma %*% xm$Omega %*% t(xm$Gamma)
    tilde_OmegaX1 <- t(CA_start) %*% SigmaX1 %*% CA_start
    nuX1_q <- n + m - r
    PsiX1_q <- (nuX1_q-m-1) * tilde_OmegaX1
    invPsiX1_q <- solve(PsiX1_q)

    # Omega0 MLE
    # SigmaX2 <- xm$Gamma0 %*% xm$Omega0 %*% t(xm$Gamma0)
    # tilde_OmegaX0 <- t(DA_start) %*% SigmaX2 %*% DA_start
    nuX0_q <- n + p - m
    # PsiX0_q <- (nuX0_q-(p-m)-1) * tilde_OmegaX0
    invPsiX0_q <- matrix(0, (p-m), (p-m))

    # muX MLE
    muX_start <- X_bar
    SigmaX_start <- 1e-4 * diag(1, nrow = p)

    # muY MLE
    muY_start <- xm$mu + t(eta_start) %*% t(CA_start) %*% muX_start
    SigmaY_start <- 1e-4 * diag(1, nrow = r)

    # Sigma MLE
    Sigma_start <- xm$Sigma

  }


  return(list(hatA = A_start,  CA = CA_start, DA = DA_start, HA = HA,
              nuY_q = nuY_q, invPsiY_q = invPsiY_q,
              nuX1_q = nuX1_q, invPsiX1_q = invPsiX1_q,
              nuX0_q = nuX0_q, invPsiX0_q = invPsiX0_q,
              eta_q = eta_start, Sigmaeta_q = Sigmaeta_q,
              muX_q = muX_start, SigmaX_q = SigmaX_start,
              muY_q = muY_start, SigmaY_q = SigmaY_start))
}



get_gamma_predenv <- function(X, Y, m) {
  ols <- lm(X ~ Y - 1)
  beta_hat <- t(ols$coefficients)
  Sx <- cov(X)
  # Sy <- cov(Y)
  # Sxy <- cov(X, Y)
  # SxCy <- Sx - Sxy %*% solve_chol(Sy) %*% t(Sxy)

  E1 <- eigen(Sx)$vectors
  E1_betahat <- t(E1) %*% beta_hat
  E1_betahat_rowss <- rowSums(E1_betahat^2)
  order_E1_betahat <- order(E1_betahat_rowss, decreasing = TRUE)

  #
  #   E2 <- eigen(SxCy)$vectors
  #   E2_betahat <- t(E2) %*% beta_hat
  #   E2_betahat_rowss <- rowSums(E2_betahat^2)
  #   order_E2_betahat <- order(E2_betahat_rowss, decreasing = TRUE)

  # if (max(E1_betahat_rowss[order_E1_betahat[1:m]]) >
  #     max(E2_betahat_rowss[order_E2_betahat[1:m]])) {
  #   # replace max by sum
  #   Gamma <- E1[, order_E1_betahat[1:m], drop=FALSE]
  # } else {
  #   Gamma <- E2[, order_E2_betahat[1:m], drop=FALSE]
  # }
  Gamma <- E1[, order_E1_betahat[1:m], drop = FALSE]

  Gamma
}

# form:
# 1 = xenv_tol, if not broken. if broken, form changed to 2
# 2 = get_gamma_predenv

start_point_predenv <- function(X, Y, m, form = 1) {
  dim_Y <- dim(Y)
  r <- dim_Y[2]
  n <- dim_Y[1]
  p <- ncol(X)
  muX <- Xbar <- colSums(X) / n
  muY <- Ybar <- colSums(Y) / n
  Zbar <- c(Xbar, Ybar)


  Y.c <- Y - tcrossprod(rep(1, n), Ybar)

  X.c <- X - tcrossprod(rep(1, n), Xbar)

  start_penv <- tryCatch(
    xenv_tol(X = X, Y = Y, u = m),
    error = function(e) e
  )


  # browser()
  # this gamma is temporary, will be replaced after calculating A
  if (is(start_penv, "error") | form == 2) {
    gamma <- get_gamma_predenv(X = X.c, Y = Y.c, m = m)
  } else if (form == 1) {
    gamma <- start_penv$Gamma
  }

  # beta <- start_penv$be

  SigX <- cov(X)
  SigXY <- cov(X, Y)
  SigY <- cov(Y)


  if (m >= 1 & m < p) {
    A <- find_A_from_gamma(gamma)
    gamma_gamma0 <- find_gammas_from_A(A)
    gamma <- gamma_gamma0$gamma
    gamma0 <- gamma_gamma0$gamma0
    # # gamma <- start_penv$Gamma
    # # this gamma is temporary, will be replaced after calculating A
    # G1 <- as.matrix(gamma[1:m, , drop = FALSE])
    # # check if G1 is invertible - else reorganize the predictors
    # if (abs(det(G1)) < 1e-7) {
    #   gamma.t <- t(gamma)
    #   X.order <- qr(gamma.t, tol = 1e-7)$pivot
    #   X <- X[, X.order]
    #   gamma <- gamma[X.order, ]
    # }
    #
    # G2 <- gamma[-(1:m), ]
    # A <- A.start <- G2 %*% solve(G1)
    # CA <- CA.start <- rbind(diag(1, m), A.start)
    # DA <- DA.start <- rbind(-t(A.start), diag(1, p - m))
    # CAtCA <- t(CA) %*% CA
    # DAtDA <- t(DA) %*% DA
    # gamma <- gamma.start <- CA %*% sqrtmatinv(CAtCA)
    gamma.t <- t(gamma)
    # gamma0 <- gamma0.start <- DA %*% sqrtmatinv(DAtDA)
    gamma0.t <- t(gamma0)


    Omega <- Omega.start <- gamma.t %*% SigX %*% gamma
    Omega.inv <- solve_chol(Omega)
    Omega0 <- Omega0.start <- gamma0.t %*% SigX %*% gamma0
    Omega0.inv <- solve_chol(Omega0)
    SigmaX <-
      gamma %*% Omega %*% gamma.t + gamma0 %*% Omega0 %*% gamma0.t

    eta <- eta.start <- Omega.inv %*% gamma.t %*% SigXY
    eta.t <- t(eta)
    # beta_start <- gamma %*% eta

  } else if (m == p) {
    # browser()
    tmp <- lm.fit(x = cbind(1, X), y = Y)
    beta <- as.matrix(coef(tmp))[-1, , drop = FALSE]
    eta <- unname(beta)
    gamma <- diag(1, nrow = p)
    gamma0 <- diag(1, nrow = 0)
    Omega <- cov(X)
    Omega0 <- diag(1, nrow = 0)
    # SigmaYcX <- start_penv$SigmaYcX
    A <- matrix(0, m, p-m)
    SigmaX <- SigX
  } else {
    beta <- matrix(0, p, r)
    eta <- matrix(0, m, r)
    gamma <- diag(1, nrow = 0)
    gamma0 <- diag(1, nrow = p)
    Omega0 <- cov(X)
    Omega <- diag(1, nrow = 0)
    A <- matrix(0, m, p-m)
    SigmaX <- SigX
  }

  eta.t <- t(eta)

  # SigmaYcX <- SigmaYcX.start <- start_penv$SigmaYcX
  SigmaYcX <- SigmaYcX.start <- SigY - eta.t %*% Omega %*% eta
  SigmaYcX.inv <- solve_chol(SigmaYcX)


  list(
    muX = muX, muY = muY, gamma = gamma, gamma0 = gamma0,
    eta = eta, Omega = Omega, Omega0 = Omega0, SigmaX = SigmaX,
    SigmaYcX = SigmaYcX, A = A
  )
}

# 기존 버젼
#
# calc_gamma_lpd_pred <- function(gamma, gamma0,
#                                 XctXc,
#                                 Xc, Yc,
#                                 # Omega.inv_plus_eta.Syx.inv.eta.t,
#                                 Omega.inv,
#                                 SigmaYcX.inv,
#                                 eta,
#                                 Omega0.inv
#                                 # XctYc.SigmaYcX.inv.eta.t
# ) {
#   -0.5 * (
#     sum(
#       Omega.inv * (crossprod(gamma, XctXc) %*% gamma)
#     ) +
#       sum(
#         Omega0.inv * (crossprod(gamma0, XctXc) %*% gamma0)
#       ) +
#
#       # sum(
#       #   SigmaYcX.inv * crossprod(Yc - (Xc %*% gamma %*% eta))
#       # )
#
#       sum(
#         SigmaYcX.inv * (
#           # crossprod(Xc %*% gamma %*% eta) -
#           #   2 * crossprod(Xc %*% gamma %*% eta, Yc)
#           crossprod(Yc - (Xc %*% gamma %*% eta))
#         )
#       )
#     # (-2) * sum(gamma * XctYc.SigmaYcX.inv.eta.t)
#   )
# }

# prior 바꾼 버젼

calc_gamma_lpd_pred <- function(invCAtCA,
                                M.half,
                                gamma, gamma0,
                                XctXc,
                                Xc, Yc,
                                # Omega.inv_plus_eta.Syx.inv.eta.t,
                                Omega.inv,
                                SigmaYcX.inv,
                                eta,
                                Omega0.inv
                                # XctYc.SigmaYcX.inv.eta.t
) {
  -0.5 * (
    sum(
      Omega.inv * (crossprod(gamma, XctXc) %*% gamma)
    ) +
      sum(
        Omega0.inv * (crossprod(gamma0, XctXc) %*% gamma0)
      ) +

      # sum(
      #   SigmaYcX.inv * crossprod(Yc - (Xc %*% gamma %*% eta))
      # )

      sum(
        SigmaYcX.inv * (
          # crossprod(Xc %*% gamma %*% eta) -
          #   2 * crossprod(Xc %*% gamma %*% eta, Yc)
          crossprod(Yc - (Xc %*% gamma %*% eta))
        )
      ) +

      sum(
        Omega.inv * tcrossprod(invCAtCA %*% eta %*% M.half)
      )
    # (-2) * sum(gamma * XctYc.SigmaYcX.inv.eta.t)
  )
}

calc_lpiror_A_pred <- function(A, A0, K.half.inv, L.half.inv) {
  -0.5 * sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
}


lpd_A_pred_semi_marginal <- function(A,
                                     XctXc,
                                     Xc,
                                     Yc,
                                     Omega.inv,
                                     SigmaYcX.inv,
                                     eta,
                                     Omega0.inv,
                                     A0,
                                     K.half.inv,
                                     L.half.inv) {
  gamma_gamma0 <- find_gammas_from_A(A)
  m <- dim(A)[2]
  r <- dim(Yc)[2]

  invCAtCA <-
    tryCatch({
      solve_chol(crossprod( rbind(diag(m), A) ))
    }, error = function(e) {
      message("Error computing invCAtCA_end")
      browser()
      stop(e)
    })

  M <- 1e-06 * diag(r)
  M.half <- diag(sqrt(diag(M)), nrow = r)

  lpd_gamma <- calc_gamma_lpd_pred(
    invCAtCA = invCAtCA,
    M.half = M.half,
    gamma = gamma_gamma0$gamma,
    gamma0 = gamma_gamma0$gamma0,
    XctXc = XctXc,
    Xc = Xc,
    Yc = Yc,
    Omega.inv = Omega.inv,
    SigmaYcX.inv = SigmaYcX.inv,
    eta = eta,
    Omega0.inv = Omega0.inv
  )

  lprior_A <- calc_lpiror_A_pred(
    A,
    A0 = A0,
    K.half.inv = K.half.inv,
    L.half.inv = L.half.inv
  )

  lpd_A <- lpd_gamma + lprior_A

  attr(lpd_A, "gamma_gamma0") <- gamma_gamma0

  lpd_A
}

###########################
# 이 함수 통채로 바꿔야함 #
###########################

sample_A_pred_metrop_manifold_proposal <- function(m, n, p, r,
                                                   Omega.inv_plus_eta.Syx.inv.eta.t,
                                                   Omega.inv,
                                                   Omega0.inv,
                                                   SigmaYcX.inv,
                                                   eta,
                                                   XctYc.SigmaYcX.inv.eta.t,
                                                   current_A,
                                                   current_gamma_gamma0,
                                                   current_det_jacob_A_over_gamma_gamma0,
                                                   XctXc,
                                                   Xc,
                                                   Yc,
                                                   YctYc = crossprod(Yc),
                                                   col_pairs_updt_set,
                                                   A0,
                                                   L.half.inv,
                                                   K.half.inv,
                                                   jacobian_only_gamma = FALSE,
                                                   adjust_jacobian = FALSE,
                                                   nbreaks_trapezoid_A_prop = 5,
                                                   browser = FALSE,
                                                   random_scan = TRUE,
                                                   adjust_gb22 = TRUE,
                                                   iter = 1,
                                                   n_update_max = 100,
                                                   ...) {
  # browser()

  log_dens_const <- (-0.5) * sum(SigmaYcX.inv * YctYc)

  lprior_A_start <- calc_lpiror_A_pred(
    A = current_A,
    A0 = A0,
    K.half.inv = K.half.inv,
    L.half.inv = L.half.inv
  )

  gamma_start <- gamma <- current_gamma_gamma0$gamma
  gamma0_start <- gamma0 <- current_gamma_gamma0$gamma0
  A_start <- find_A_from_gamma(gamma_start)

  invCAtCA <-
    tryCatch({
      solve_chol(crossprod( rbind(diag(m), A_start) ))
    }, error = function(e) {
      message("Error computing invCAtCA_end")
      browser()
      stop(e)
    })

  M <- 1e-06 * diag(r)
  M.half <- diag(sqrt(diag(M)), nrow = r)

  # convert from A to Omat (track the jacobian)
  Omegastar.inv_eig <- eigen(Omega.inv_plus_eta.Syx.inv.eta.t, symmetric = TRUE)
  Omega0.inv_eig <- eigen(Omega0.inv, symmetric = TRUE)
  log_jacob_A_over_gamma_gamma0_start <- current_det_jacob_A_over_gamma_gamma0

  # needed for the GB22 distribution
  lambda <- c(Omegastar.inv_eig$values, Omega0.inv_eig$values)

  lpd_gamma_gamma0_start <- calc_gamma_lpd_pred(
    invCAtCA = invCAtCA,
    M.half = M.half,
    gamma = gamma,
    gamma0 = gamma0,
    XctXc = XctXc,
    # Omega.inv_plus_eta.Syx.inv.eta.t = Omega.inv_plus_eta.Syx.inv.eta.t,
    Xc = Xc, Yc = Yc,
    # Omega.inv_plus_eta.Syx.inv.eta.t,
    Omega.inv = Omega.inv,
    SigmaYcX.inv = SigmaYcX.inv,
    eta = eta,
    Omega0.inv = Omega0.inv
  )

  lpd_A_start <- lpd_gamma_gamma0_start +
    # Jacobian A -> (gamma, gamma0) -> A
    log_jacob_A_over_gamma_gamma0_start +
    lprior_A_start

  Omat <- cbind(
    gamma %*% Omegastar.inv_eig$vectors,
    gamma0 %*% Omega0.inv_eig$vectors
  )

  Omat_old <- Omat


  transition_adj <- 0
  # all_transition <- c()

  tau_i_all <- matrix(0, p, p)
  tau_i_all[, 1:m] <- XctYc.SigmaYcX.inv.eta.t %*% Omegastar.inv_eig$vectors

  counter <- 1
  # Omat_new <- Omat

  update_set <- seq_len(length(col_pairs_updt_set)) %>%
    {
      if (random_scan) sample(., replace = FALSE) else .
    }

  count <- 0
  transition_adj_all <- c()

  for (idx in update_set) {
    count <- count + 1
    if (count > n_update_max) break

    i <- col_pairs_updt_set[[idx]][1]
    j <- col_pairs_updt_set[[idx]][2]
    H_i <- XctXc * lambda[i] * (0.5)
    H_j <- XctXc * lambda[j] * (0.5)
    tau_i <- tau_i_all[, i, drop = FALSE]
    tau_j <- tau_i_all[, j, drop = FALSE]

    # if (count == 2) browser()


    calc_A_B_C_D_params <- function(Oij) {
      list(
        A = crossprod(Oij, H_i) %*% Oij,
        B = crossprod(Oij, H_j) %*% Oij,
        C = crossprod(Oij, tau_i),
        D = crossprod(Oij, tau_j)
      )
    }


    Oij_old <- Omat[, c(i, j)]
    params_given_old <- calc_A_B_C_D_params(Oij_old)

    log_dens_diff <- -Inf

    Z_GB22_gen <- tryCatch(
      rgmatrixbingham22_4param_an(
        A = params_given_old$A,
        B = params_given_old$B,
        C = params_given_old$C,
        D = params_given_old$D,
        log_den = TRUE
      ),
      error = function(e) e
    )

    Z_GB22 <- Z_GB22_gen$Z
    if (is(Z_GB22, "error")) browser()
    Omat[, c(i, j)] <- Oij_new <- Oij_old %*% Z_GB22

    params_given_new <- calc_A_B_C_D_params(Oij_new)


    # if (count == 3) browser()

    log_den_proposal <- mapply(
      function(this_params, this_Z) {
        dgmatrixbingham22_4param_an(
          Z = this_Z, # (Z_GB22), #t(Z_GB22),
          A = this_params$A, # 0.5 * crossprod(Oij_new, H_i) %*% Oij_new,
          B = this_params$B, # 0.5 * crossprod(Oij_new, H_j) %*% Oij_new,
          C = this_params$C,
          D = this_params$D,
          log_den = TRUE,
          # nbreak = n_breaks
        )
      },
      list(
        new_given_old = params_given_old,
        old_given_new = params_given_new
      ),
      list(
        new_given_old = Z_GB22,
        old_given_new = t(Z_GB22)
      ),
      SIMPLIFY = FALSE
    )

    log_dens_diff <- log_den_proposal$old_given_new - log_den_proposal$new_given_old

    # browser()

    if (is.na(transition_adj) | is.na(log_dens_diff) | is.infinite(log_dens_diff)) browser()

    if (is.na(log_dens_diff)) {
      Omat[, c(i, j)] <- Oij_old
      next
    } else {
      transition_adj <- transition_adj + log_dens_diff
    } # counter <- counter + 1
    # all_transition <- c(all_transition, log_dens_diff)

    transition_adj_all <- c(transition_adj_all, log_dens_diff)
  }



  # if (iter > 50) browser()

  gamma_end <- tryCatch(
    tcrossprod(
      Omat[, 1:m, drop = FALSE],
      Omegastar.inv_eig$vectors
    ),
    error = function(e) e
  )

  # if (is(gamma_end, "error")) browser()

  gamma0_end <- tcrossprod(
    Omat[, -(1:m), drop = FALSE],
    Omega0.inv_eig$vectors
  )


  # if (nrow(current_A) >= ncol(current_A)) {
  # r-u >= A
  A_end <- find_A_from_gamma(gamma = gamma_end)
  # } else {
  # A_end <- find_A_from_gamma0(gamma0 = gamma0_end)
  # }

  lprior_A_end <- calc_lpiror_A_pred(
    A = A_end,
    A0 = A0,
    K.half.inv = K.half.inv,
    L.half.inv = L.half.inv
  )

  # tmp_tmp <- lapply(
  #   list(old_A = current_A,
  #        new_A = A_end),
  #   function(this_A) {
  #     lapply(
  #       c(only_gamma = TRUE, both_gamma_gamma0 = FALSE),
  #       function(xx)
  #         find_gammas_from_A(
  #           this_A,
  #           jacobian = TRUE,
  #           log = TRUE,
  #           jacobian_only_gamma = xx
  #         )
  #     )
  #   }
  # )

  # tmp_tmp$old_A %>% sapply("[[", "det_Jacobian_Omat_over_A")
  # tmp_tmp$new_A %>% sapply("[[", "det_Jacobian_Omat_over_A")

  # gamma gamma0 in the "A" formulation
  gamma_gamma0_A_end <- find_gammas_from_A(
    A_end,
    jacobian = adjust_jacobian,
    log = TRUE,
    jacobian_only_gamma = jacobian_only_gamma
  )

  # gamma_gamma0_A_end$det_Jacobian_Omat_over_A
  # jacobian_gamma_A_numeric(
  #   A_end,
  #   jacobian_only_gamma = !jacobian_only_gamma#,
  #   # method = "simple"
  # )

  if (adjust_jacobian) {
    log_jacob_A_over_gamma_gamma0_end <- -gamma_gamma0_A_end$det_Jacobian_Omat_over_A
  } else {
    log_jacob_A_over_gamma_gamma0_end <- 0
  }
  # lpds of (gamma, gamma0) in the "A" formulation
  # and the general steifel formulation
  gamma_gamma0_end_list <- list(
    list(
      gamma = gamma_end,
      gamma0 = gamma0_end
    ),
    gamma_gamma0_A_end
  )

  invCAtCA_end <-
    tryCatch({
      solve_chol(crossprod( rbind(diag(m), A_end) ))
    }, error = function(e) {
      message("Error computing invCAtCA_end")
      browser()
      stop(e)
    })

  tmp_lpd_gamma_gamma0_vec <- sapply(
    gamma_gamma0_end_list,
    function(this_gamma_gamma0) {
      calc_gamma_lpd_pred(
        invCAtCA = invCAtCA_end,
        M.half = M.half,
        gamma = this_gamma_gamma0$gamma,
        gamma0 = this_gamma_gamma0$gamma0,
        XctXc = XctXc,
        # Omega.inv_plus_eta.Syx.inv.eta.t = Omega.inv_plus_eta.Syx.inv.eta.t,
        Xc = Xc, Yc = Yc,
        # Omega.inv_plus_eta.Syx.inv.eta.t,
        Omega.inv = Omega.inv,
        SigmaYcX.inv = SigmaYcX.inv,
        eta = eta,
        Omega0.inv = Omega0.inv
      )
    }
  )
  # take the maximum
  max_indic <- which.max(tmp_lpd_gamma_gamma0_vec)[1]
  lpd_gamma_gamma0_end <- tmp_lpd_gamma_gamma0_vec[max_indic]
  gamma_gamma0_end_return <- gamma_gamma0_end_list[[max_indic]]

  # browser()

  lpd_A_end <- lpd_gamma_gamma0_end +
    # Jacobian (gamma, gamma0) -> A
    log_jacob_A_over_gamma_gamma0_end +
    lprior_A_end


  # transition_adj <- 0

  metrop_log_ratio <- (lpd_A_end - lpd_A_start + transition_adj)

  # if (browser) browser()

  # do metropolis accept reject
  check_accpt <- tryCatch(
    log(runif(1)) < metrop_log_ratio,
    error = function(e) e
  )

  if (is(check_accpt, "error")) browser()

  computation <- tryCatch(
    {
      if (check_accpt) {
        # accept
        A <- A_end
        lpd_A <- lpd_A_end
        accpt <- 1
        gamma_gamma0 <- gamma_gamma0_A_end
        det_jacob <- log_jacob_A_over_gamma_gamma0_end
      } else {
        # reject
        A <- current_A
        lpd_A <- lpd_A_start
        accpt <- 0
        gamma_gamma0 <- current_gamma_gamma0
        det_jacob <- log_jacob_A_over_gamma_gamma0_start
      }
    },
    error = function(e) e
  )

  if (is(computation, "error")) browser()

  attr(lpd_A, "gamma_gamma0") <- gamma_gamma0_end_return

  # browser()

  list(
    A = A,
    lpd = lpd_A,
    accpt = accpt,
    det_jacob_A_over_gamma_gamma0 = det_jacob # ,
    # gamma_gamma0 = gamma_gamma0_end_return
  )
}

#' #' @rdname Benvlp_MC_pred
#' #' @inheritParams Benvlp_MC_pred
#' #' @param m predictor envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' #' @export
#' Benvlp_MC_pred_gibbs <- Benvlp_MC_pred

# Benvlp_MC_pred_gibbs_devel <- function(Y, X,
#                                        m,
#                                        n.iter = 1000,
#                                        tau = 0.01,
#                                        autotune = TRUE,
#                                        tune.accpt.prop.lower = 0.2,
#                                        tune.accpt.prop.upper = 0.3,
#                                        tune.incr = 0.05,
#                                        burnin.prop = 0.5,
#                                        tune.burnin.prop = 0.75,
#                                        tune_nterm = 25,
#                                        tune_adj_scale_nterm = 50,
#                                        checkstep = FALSE,
#                                        starting = "mle",
#                                        show_progress = TRUE,
#                                        n.chains = 1,
#                                        chains_parallel = FALSE,
#                                        A_proposal = "rwmh",
#                                        jacobian_only_gamma = FALSE,
#                                        posterior_gamma_space = TRUE,
#                                        nbreaks_trapezoid_A_prop = 5,
#                                        n_iter_A_rwmh = 1,
#                                        rwmh_elementwise = TRUE,
#                                        rwmh_elementwise_fraction = 1,
#                                        n_A_update_max = Inf,
#                                        ...) {
#   X <- data.matrix(X)
#   Y <- data.matrix(Y)
#
#   stopifnot(A_proposal %in% c("stiefel-experimental", "rwmh", "rwmh_marginal"))
#
#
#   dim_X <- dim(X)
#   dim_Y <- dim(Y)
#
#   n <- dim_X[1]
#   p <- dim_X[2]
#   r <- dim_Y[2]
#
#
#   Y <- Y.orig <- data.matrix(Y)
#   Yc <- scale(Y, center = TRUE, scale = FALSE)
#   Y.bar <- attr(Yc, "scaled:center")
#   YctYc <- crossprod(Yc)
#
#
#   X <- X.orig <- data.matrix(X)
#   Xc <- scale(X, center = TRUE, scale = FALSE)
#   X.bar <- attr(Xc, "scaled:center")
#   XctXc <- crossprod(Xc)
#
#   XctYc <- crossprod(Xc, Yc)
#
#   Z <- cbind(X, Y)
#
#   Z.bar <- c(X.bar, Y.bar)
#
#   n.burnin <- ceiling(n.iter * burnin.prop)
#   n.iter.final <- n.iter - n.burnin
#   n.tune <- ceiling(n.burnin * tune.burnin.prop)
#   burnin_iter <- seq_len(n.burnin)
#   # log_2pi <- log(2*pi)
#
#   if (m < 0 || m > p || m != as.integer(m)) {
#     stop(paste0("\'m\' must be an integer between 0 and ", p, " (p)"))
#   }
#
#   X.order <- 1:p # changed if G1 is singular
#
#
#   A_exists <- m < p & m > 0
#
#   if (length(tau) == 1) {
#     tau <- matrix(tau, p - m, m)
#   }
#
#
#
#   # prior parameters
#   PsiY <- diag(1, r) * 1e-6
#   nuY <- r
#   PsiX <- diag(1, m) * 1e-6
#   nuX <- m
#   PsiX0 <- diag(p - m) * 1e-6
#   nuX0 <- p - m
#
#   # 기존 prior
#
#   e <- matrix(0, m, r)
#   is_e_zero <- TRUE
#   M <- 1e-6 * diag(1, nrow = m)
#   M.half <- diag(sqrt(diag(M)), nrow = m)
#   M.inv <- diag(1 / diag(M), nrow = m)
#   M.e <- M %*% e
#
#   # 변경된 prior
#
#   # e <- matrix(0, p, r)
#   # is_e_zero <- TRUE
#   # M <- 1e-6 * diag(1, nrow = r)
#   # M.half <- diag(sqrt(diag(M)), nrow = r)
#   # M.inv <- diag(1 / diag(M), nrow = r)
#   # M.e <- M %*% t(e)
#   # psi_eta <- 1e10
#
#   if (A_exists) {
#     A0 <- matrix(0, p - m, m)
#     K <- 1e6 * diag(1, p - m)
#     K.inv <- 1e-6 * diag(1, p - m)
#     K.half.inv <- 1e-3 * diag(1, p - m)
#     L <- 1e6 * diag(1, m)
#     L.inv <- 1e-6 * diag(1, m)
#     L.half.inv <- 1e-3 * diag(1, m)
#   }
#
#   start_penv <- start_point_predenv(X, Y, m, 1)
#
#   # basically attach(start_penv) but does not
#   # affect the environment
#   {
#     for (x in names(start_penv)) {
#       assign(x, start_penv[[x]])
#     }
#   }
#   SigmaYcX.inv <-
#     tryCatch({
#       solve_chol(SigmaYcX)
#     }, error = function(e) {
#       message("Error computing SigmaYcX.inv")
#       browser()
#       stop(e)
#     })
#   Omega.inv <-
#     tryCatch({
#       solve_chol(Omega)
#     }, error = function(e) {
#       message("Error computing Omega.inv")
#       browser()
#       stop(e)
#     })
#   Omega0.inv <-
#     tryCatch({
#       solve_chol(Omega0)
#     }, error = function(e) {
#       message("Error computing Omega0.inv")
#       browser()
#       stop(e)
#     })
#
#   gamma_gamma0 <- find_gammas_from_A(A)
#
#   log_jacob_A_over_gamma_gamma0_end <-
#     det_jacob_A_over_gamma_gamma0 <- 0
#
#   if (m > 0) {
#     eta.t <- t(eta)
#     gamma.t <- t(gamma)
#   }
#
#
#   log_sqrt_2pi <- log(sqrt(2 * pi))
#
#   # browser()
#
#   runMC <- function(chain_no = 1) {
#     runif(chain_no)
#     one_n <- rep(1, n)
#     # e.check <- crossprod(Xc, Yc) + e
#
#     # needs to be recalucated at every iteration if A_exists
#     # if (m > 0) {
#     #   # M.check <- M + crossprod(Xc %*% gamma)
#     #   # M.check_inv <- solve_chol(M.check)
#     #   # M.check_half <- sqrtmat(M.check)
#     #   M.tilde <- M + crossprod(X_muX %*% gamma)
#     #   e.tilde <- e + crossprod(X_muX, Y_muY)
#     #   M.tilde_inv <- solve_chol(M.tilde)
#     #   H.tilde <- YctYc +
#     #     crossprod(sqrtmatinv(M) %*% crossprod(gamma, e)) -
#     #     crossprod(sqrtmatinv(M.check) %*% crossprod(gamma, e.check))
#     #   eta.mean <- M.check_inv %*% crossprod(gamma, e.check)
#     # } else {
#     #   H.tilde <- YctYc
#     # }
#
#
#     if (m == 0) {
#       eta <- 0
#     }
#
#
#     A_exists <- m > 0 & m < p
#
#     muZ.all <- muX.all <- muY.all <- SigmaYcX.all <-
#       SigmaX.all <- SigmaX.inv.all <-
#       Omega.all <- Omega0.all <- eta.all <- beta.all <-
#       gamma.all <- gamma0.all <- A.all <- accpt.A.all <-
#       llik.contri.all <- vector("list", n.iter)
#
#     # llik.all2 <- llik.all3 <-
#     lpd.A.all <- accpt.A.ave <- llik.all <- lpd.full.all <- rep(0, n.iter)
#
#
#     DeltaZ <- matrix(0, p + r, p + r) # just initializes, will be replaced at every iteration
#
#
#     col_pairs_all <- combn(1:p, 2)
#     # if (u <= r-u) {
#     #   col_pairs_all <- combn(1:u, 2)
#     # } else {
#     #   col_pairs_all <- combn((u+1):r, 2)
#     # }
#     # col_pairs_all <- t(cbind(1:r, c(2:r, 1)))
#
#     col_pairs_list <- lapply(
#       1:ncol(col_pairs_all),
#       function(jj) {
#         col_pairs_all[, jj]
#       }
#     )
#
#     if (m <= p - m) {
#       col_pairs_list <- col_pairs_list[
#         sapply(col_pairs_list, function(x) any(x %in% 1:m))
#       ]
#     } else {
#       col_pairs_list <- col_pairs_list[
#         sapply(col_pairs_list, function(x) any(x %in% (m + 1):p))
#       ]
#     }
#
#     pb <- txtProgressBar(min = 1, max = n.iter, style = 3)
#
#     for (iter in 1:n.iter) {
#       # browser()
#
#       # generate A, if 0 < m < p
#       if (A_exists & (iter == 1 | iter %% n_iter_A_rwmh == 0)) {
#         if (A_proposal == "stiefel-experimental") {
#           SigmaYcX.inv.eta.t <- tcrossprod(SigmaYcX.inv, eta)
#           Omega.inv_plus_eta.Syx.inv.eta.t.t <- Omega.inv +
#             eta %*% SigmaYcX.inv.eta.t
#           XctYc.SigmaYcX.inv.eta.t <- XctYc %*% SigmaYcX.inv.eta.t
#           mh_step <- sample_A_pred_metrop_manifold_proposal(
#             m = m,
#             n = n,
#             p = p,
#             r = r,
#             Omega.inv_plus_eta.Syx.inv.eta.t = Omega.inv_plus_eta.Syx.inv.eta.t.t,
#             Omega.inv = Omega.inv,
#             Omega0.inv = Omega0.inv,
#             SigmaYcX.inv = SigmaYcX.inv,
#             eta = eta,
#             XctYc.SigmaYcX.inv.eta.t = XctYc.SigmaYcX.inv.eta.t,
#             current_A = A,
#             current_gamma_gamma0 = gamma_gamma0,
#             current_det_jacob_A_over_gamma_gamma0 = det_jacob_A_over_gamma_gamma0,
#             XctXc = XctXc,
#             Xc = Xc,
#             Yc = Yc,
#             col_pairs_updt_set = col_pairs_list,
#             A0 = A0,
#             L.half.inv = L.half.inv,
#             K.half.inv = K.half.inv,
#             browser = (iter >= 1000),
#             jacobian_only_gamma = jacobian_only_gamma,
#             adjust_jacobian = !posterior_gamma_space,
#             nbreaks_trapezoid_A_prop = nbreaks_trapezoid_A_prop,
#             iter = iter,
#             n_update_max = n_A_update_max
#           )
#
#           det_jacob_A_over_gamma_gamma0 <- mh_step$det_jacob_A_over_gamma_gamma0
#         } else if (A_proposal == "rwmh") {
#           lpd_A_curr <- lpd_A_pred_semi_marginal(
#             A,
#             XctXc = XctXc,
#             Xc = Xc,
#             Yc = Yc,
#             Omega.inv = Omega.inv,
#             SigmaYcX.inv = SigmaYcX.inv,
#             eta = eta,
#             Omega0.inv = Omega0.inv,
#             A0 = A0,
#             K.half.inv = K.half.inv,
#             L.half.inv = L.half.inv
#           )
#
#           mh_step <- rwmh_colwise_A(
#             A_start = A,
#             lpd_start = lpd_A_curr,
#             tau = tau,
#             elementwise = rwmh_elementwise,
#             n_update_max = n_A_update_max,
#             n_rwmh_iter = n_iter_A_rwmh,
#             random_update_fraction = rwmh_elementwise_fraction,
#             lpd_func = lpd_A_pred_semi_marginal,
#             XctXc = XctXc,
#             Xc = Xc,
#             Yc = Yc,
#             Omega.inv = Omega.inv,
#             SigmaYcX.inv = SigmaYcX.inv,
#             eta = eta,
#             Omega0.inv = Omega0.inv,
#             A0 = A0,
#             K.half.inv = K.half.inv,
#             L.half.inv = L.half.inv
#           )
#         } else if (A_proposal == "rwmh_marginal") {
#
#         }
#
#
#         # autotune in burnin to make accpt prob between
#         # tune.accpt.prop.lower and tune.accpt.prop.upper
#         if (autotune &
#             iter >= tune_nterm &
#             iter %% 5 == 0 &
#             iter <= n.tune &
#             grepl("rwmh", A_proposal)) {
#           # nterm <- tune_nterm
#           # accpt.mean <- mean(accpt.A.ave[(iter-nterm+1):iter])
#           # if (accpt.mean > tune.accpt.prop.upper) {
#           #   tau <- tau * (1 + tune.incr)
#           # } else if (accpt.mean < tune.accpt.prop.lower) {
#           #   tau <- tau * (1 - tune.incr)
#           # }
#           tau <- autotune_param(
#             draw_param_list = A.all[1:iter],
#             tune_param = tau,
#             accpt_list = accpt.A.all[1:iter],
#             tune_nterm = tune_nterm,
#             tune.incr = tune.incr,
#             tune.accpt.prop.lower = tune.accpt.prop.lower,
#             tune.accpt.prop.upper = tune.accpt.prop.upper,
#             adjust_scale = ((iter %% tune_adj_scale_nterm) == 0)
#           )
#         }
#
#
#         A.all[[iter]] <- A <- mh_step$A
#         accpt.A.all[[iter]] <- mh_step$accpt # tcrossprod(rep(1, p-m), mh_step$accpt)
#         accpt.A.ave[iter] <- mean(mh_step$accpt)
#         lpd_A <- mh_step$lpd
#         lpd.A.all[iter] <- c(lpd_A)
#         gamma_gamma0 <- attr(lpd_A, "gamma_gamma0")
#         gamma.all[[iter]] <- gamma <- gamma_gamma0$gamma
#         gamma0.all[[iter]] <- gamma0 <- gamma_gamma0$gamma0
#       }
#
#       # browser()
#       ################################
#       # 이 부분 다시 정리해서 바꾸기 #
#       ################################
#       invCAtCA <- solve_chol(crossprod(rbind(diag(m), A)))
#       # generate eta
#       if (m > 0) {
#
#         # 기존 prior
#
#         M.tilde <- M + crossprod(Xc %*% gamma)
#         e.tilde <- crossprod(gamma, XctYc) + M.e
#         M.tilde_inv <- solve_chol(M.tilde)
#         eta.mean <- M.tilde_inv %*% e.tilde
#         eta.all[[iter]] <- eta <-
#           rMatrixNormal(eta.mean, M.tilde_inv, SigmaYcX)
#
#         # 변경된 prior
#
#         # eta.Sig <- solve_chol(
#         #   kronecker(M, 1/psi_eta * invCAtCA %*% Omega.inv %*% invCAtCA) + kronecker(SigmaYcX.inv, t(gamma) %*% XctXc %*% gamma)
#         # )
#         # eta.mean <- eta.Sig %*% matrix((SigmaYcX.inv %*% t(XctYc) %*% gamma), m*r, 1) # + M.e %*% gamma %*% Omega.inv %*% invCAtCA
#         # eta.all[[iter]] <- eta <-
#         #   matrix(rmvnorm(1, eta.mean, eta.Sig), m, r)
#
#       } else {
#         eta.all[[iter]] <- eta <- NULL
#       }
#
#       # browser()
#
#       ################################
#       # 이 부분 다시 정리해서 바꾸기 #
#       ################################
#
#       # generate SigmaYcX
#       if (m > 0) {
#
#         # 기존 prior
#         PsiY.tilde <- PsiY +
#           crossprod(Yc - Xc %*% gamma %*% eta) +
#           crossprod(M.half %*% (eta - e))
#
#         # 변경된 prior
#         # PsiY.tilde <- PsiY +
#         #   crossprod(Yc - Xc %*% gamma %*% eta)
#
#       } else {
#         PsiY.tilde <- PsiY + YctYc
#       }
#       SigmaYcX.all[[iter]] <- SigmaYcX <-
#         rinvwish(r, Phi = PsiY.tilde, nu = (nuY + n))
#       SigmaYcX.inv <- solve_chol(SigmaYcX)
#
#       # browser()
#
#       ################################
#       # 이 부분 다시 정리해서 바꾸기 #
#       ################################
#
#       # generate Omega if m > 0
#       if (m > 0) {
#
#         # 기존 prior
#         PsiX.tilde <- PsiX + crossprod(Xc %*% gamma)
#         Omega.all[[iter]] <-
#           Omega <- rinvwish(m, Phi = PsiX.tilde, nu = nuX + n - 1)
#         Omega.inv <- solve_chol(Omega)
#
#         # 변경된 prior
#         # PsiX.tilde <- PsiX + crossprod(Xc %*% gamma) + psi_eta * (invCAtCA %*% eta - t(gamma) %*% e) %*% M %*% t((invCAtCA %*% eta - t(gamma) %*% e))
#         # Omega.all[[iter]] <-
#         #   Omega <- rinvwish(m, Phi = PsiX.tilde, nu = (nuX + n - r))
#         # Omega.inv <- solve_chol(Omega)
#
#       } else {
#         Omega.all[[iter]] <- Omega <- Omega.inv <- NULL
#       }
#
#       # browser()
#
#       # generate Omega0 if m < p
#       if (m < p) {
#         PsiX0.tilde <- PsiX0 + crossprod(Xc %*% gamma0)
#         Omega0.all[[iter]] <-
#           Omega0 <- rinvwish((p - m), Phi = PsiX0.tilde, nu = (nuX0 + n - 1))
#         Omega0.inv <- solve_chol(Omega0)
#       } else {
#         Omega0.all[[iter]] <- Omega0 <- Omega0.inv <- NULL
#       }
#
#       # calculate Sigma and Sigma.inv
#       if (m > 0) {
#         SigmaX.inv.1 <- gamma %*% tcrossprod(Omega.inv, gamma)
#         SigmaX.inv.0 <- gamma0 %*% tcrossprod(Omega0.inv, gamma0)
#         SigmaX.inv.all[[iter]] <- SigmaX.inv <-
#           SigmaX.inv.1 + SigmaX.inv.0
#
#         SigmaX.all[[iter]] <- SigmaX <-
#           gamma %*% tcrossprod(Omega, gamma) +
#           gamma0 %*% tcrossprod(Omega0, gamma0)
#       } else {
#         SigmaX.inv.0 <- Omega0.inv
#         SigmaX.inv.1 <- NULL
#         SigmaX.inv.all[[iter]] <- SigmaX.inv <-
#           SigmaX.inv.0
#
#         SigmaX.all[[iter]] <- SigmaX <- Omega0
#       }
#
#       # calculate beta
#       if (m > 0) {
#         beta.all[[iter]] <- beta <- gamma %*% eta
#       } else if (m == 0) {
#         beta.all[[iter]] <- beta <- matrix(0, nrow = p, ncol = r)
#       }
#
#       # browser()
#
#       # generate muX and muY
#       # generate muZ - then extract muX and muY
#       # SigmaX <- gamma %*% tcrossprod(Omega, gamma) +
#       #   gamma0 %*% tcrossprod(Omega0, gamma0)
#
#       DeltaZ[1:p, 1:p] <- SigmaX
#       if (m > 0) {
#         Omega.eta <- Omega %*% eta
#         DeltaZ[1:p, (p + 1):(p + r)] <- gamma %*% Omega.eta
#         DeltaZ[(p + 1):(p + r), (p + 1):(p + r)] <- SigmaYcX + crossprod(eta, Omega.eta)
#       } else {
#         # DeltaZ[1:p, (p+1):(p+r)] <- 0
#         # ^not need, initialized at 0
#         DeltaZ[(p + 1):(p + r), (p + 1):(p + r)] <- SigmaYcX
#       }
#
#       DeltaZ[(p + 1):(p + r), 1:p] <- t(DeltaZ[1:p, (p + 1):(p + r)])
#       DeltaZ.n <- DeltaZ / n
#
#       # browser()
#
#       if (!all(eigen(DeltaZ.n)$values > 0)){
#         DeltaZ.n <- make_positive_definite(DeltaZ.n)
#       }
#
#       muZ.all[[iter]] <- muZ <- rmvnorm(dim = p + r, mu = Z.bar, sigma = DeltaZ.n)
#       muX.all[[iter]] <- muX <- muZ[1:p]
#       muY.all[[iter]] <- muY <- muZ[(p + 1):(r + p)]
#
#       X_muX <- X - tcrossprod(one_n, muX)
#       Y_muY <- Y - tcrossprod(one_n, muY)
#
#       # browser()
#
#       # # if A_exists,  recalculate M.check, M.check_inv, and H.tilde --
#       # if (A_exists) {
#       #   M.check <- M + crossprod(Xc %*% gamma)
#       #   M.check_inv <- solve_chol(M.check)
#       #   M.check_half <- sqrtmat(M.check)
#       #   H.tilde <- YctYc +
#       #     crossprod(M.half %*% crossprod(gamma, e)) -
#       #     crossprod(sqrtmatinv(M.check) %*% crossprod(gamma, e.check))
#       # }
#
#       log_det_SigmaYcX <- log(det(SigmaYcX))
#       if (m > 0) log_det_Omega <- log(det(Omega))
#       if (m < p) log_det_Omega0 <- log(det(Omega0))
#       resi <- Y_muY - X_muX %*% beta
#       llik.contri.1 <- rowSums((resi %*% SigmaYcX.inv) * resi)
#       llik.contri.2 <- rowSums((X_muX %*% SigmaX.inv) * X_muX)
#
#
#       # llikelihood contribution
#       if (A_exists) {
#         llik.contri.all[[iter]] <- llik.contri <-
#           -0.5 * n * (p + r) * log_sqrt_2pi +
#           -0.5 * (log_det_SigmaYcX + log_det_Omega + log_det_Omega0 +
#                     llik.contri.1 + llik.contri.2)
#       }
#       # else if (m == p) {
#       #   llik.contri.all[[iter]] <- llik.contri <-
#       #     - 0.5 * n * (p + r) * log_sqrt_2pi +
#       #     - 0.5 * (log_det_SigmaYcX + log_det_Omega +
#       #                llik.contri.1 + llik.contri.2)
#       # }
#       else {
#         llik.contri.all[[iter]] <- llik.contri <-
#           -0.5 * n * (p + r) * log_sqrt_2pi +
#           -0.5 * (log_det_SigmaYcX + log_det_Omega0 +
#                     llik.contri.1 + llik.contri.2)
#       }
#
#       # total log likelihood
#       llik.all[[iter]] <- llik <- sum(llik.contri)
#
#
#       # # log posterior density
#       # if (A_exists) {
#       #   lpd.full.all[[iter]] <- lpd_full <-
#       #     llik - -0.5 * (
#       #       (nuY + r + p + 1) * log_det_SigmaYcX +
#       #         sum(SigmaYcX.inv * PsiY) +
#       #         sum(SigmaYcX.inv *
#       #           crossprod(M.half %*% (eta))) +
#       #         (nuX + m + 1) * log_det_Omega +
#       #         sum(Omega.inv * PsiX) +
#       #         (nuX0 + p - m + 1) * log_det_Omega0 +
#       #         sum(Omega0.inv * PsiX0) +
#       #         sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
#       #     )
#       # } else if (m == p) {
#       #   lpd.full.all[[iter]] <- lpd_full <-
#       #     llik - -0.5 * (
#       #       (nuY + r + p + 1) * log_det_SigmaYcX +
#       #         sum(SigmaYcX.inv * PsiY) +
#       #         sum(SigmaYcX.inv *
#       #           crossprod(M.half %*% (eta))) +
#       #         (nuX + m + 1) * log_det_Omega +
#       #         sum(Omega.inv * PsiX)
#       #     )
#       # } else {
#       #   lpd.full.all[[iter]] <- lpd_full <-
#       #     llik - -0.5 * (
#       #       (nuY + r + p + 1) * log_det_SigmaYcX +
#       #         sum(SigmaYcX.inv * PsiY) +
#       #         (nuX0 + p - m + 1) * log_det_Omega0 +
#       #         sum(Omega0.inv * PsiX0)
#       #     )
#       # }
#
#
#       # tune A, if exists & proposal is rwmh
#       # if (A_exists & autotune &
#       #     iter >= tune_nterm & iter %% 5 == 0 & iter <= n.tune &
#       #     A_proposal != "stiefel") {
#       #
#       #
#       #   # accpt_A_col <- Reduce(
#       #   #   '+',
#       #   #   accpt.A.all[(iter-tune_nterm+1):iter]
#       #   # ) / tune_nterm
#       #   #
#       #   # for(j in 1:m) {
#       #   #   ave_accpt <- accpt_A_col[j]
#       #   #
#       #   #   if (ave_accpt > tune.accpt.prop.upper) {
#       #   #     tau[, j] <- tau[, j] * (1 + tune.incr)
#       #   #   } else if (ave_accpt < tune.accpt.prop.lower) {
#       #   #     tau[, j] <- tau[, j] * (1 - tune.incr)
#       #   #   }
#       #   #
#       #   #   if (iter %% tune_nterm == 0) {
#       #   #     # last nterms iterations of A[, j]
#       #   #     # each column of A_j_last_nterm is
#       #   #     # from a single interation of A[, j]
#       #   #     A_j_last_nterm <- do.call(
#       #   #       cbind,
#       #   #       lapply(A.all[(iter-tune_nterm+1):iter],
#       #   #              function(x) x[, j]))
#       #   #     par_sd <- apply(A_j_last_nterm, 1, sd)
#       #   #     par_sd_mean <- sum(par_sd) / (p - m)
#       #   #     tau_j_mean <- sum(tau[, j]) /(p - m)
#       #   #
#       #   #     if(par_sd_mean > 0)
#       #   #       tau[, j] <- 0.5*par_sd/par_sd_mean*tau_j_mean +
#       #   #       0.5*tau[, j]
#       #   #   }
#       #   # }
#       # }
#
#
#       # if (iter == 100) browser()
#
#       if (show_progress) {
#         utils::setTxtProgressBar(pb, iter)
#       }
#     }
#
#     MC <- list( # "n.iter" = n.iter, "X" = X, "Y" = Y, "X.order" = X.order, "m" = m,
#       # "burnin" = burnin,
#       muX = muX.all[-burnin_iter],
#       muY = muY.all[-burnin_iter],
#       eta = eta.all[-burnin_iter],
#       SigmaYcX = SigmaYcX.all[-burnin_iter],
#       Omega = Omega.all[-burnin_iter],
#       Omega0 = Omega0.all[-burnin_iter],
#       A = A.all[-burnin_iter],
#       beta = beta.all[-burnin_iter],
#       lpd.A = lpd.A.all[-burnin_iter],
#       gamma = gamma.all[-burnin_iter],
#       gamma0 = gamma0.all[-burnin_iter],
#       accpt.A.ave = accpt.A.ave[-burnin_iter],
#       # lpd.full = lpd.full.all[-burnin_iter],
#       llik = llik.all[-burnin_iter],
#       # "llik2" = llik.all2[-burnin_iter],
#       # "llik3" = llik.all3[-burnin_iter],
#       llik.contri = llik.contri.all[-burnin_iter],
#       accpt.A.all = accpt.A.all[-burnin_iter],
#       tau = tau,
#       prior_param = list()
#     )
#
#     MC
#   }
#
#   lapply_ <- function(...) {
#     if (n.chains == 1 | !chains_parallel) {
#       lapply(...)
#     } else {
#       future.apply::future_apply(...,
#                                  future.seed = TRUE
#       )
#     }
#   }
#
#   all_MCs <- lapply_(
#     1:n.chains,
#     function(j) {
#       runMC(chain_no = j)
#     }
#   )
#
#
#   # reshape the final list
#   # each parameter, e.g. Omega is an
#   # n.chain long list of n.iter.final long lists.
#
#   vars <- setdiff(names(all_MCs[[1]]), "prior_param")
#
#   out_vars <- lapply(
#     vars,
#     function(x) {
#       out <- lapply(all_MCs, "[[", x)
#       names(out) <- paste0("Chain_", 1:n.chains)
#       out
#     }
#   )
#   names(out_vars) <- vars
#
#   out_res <- c(
#     list(
#       n.iter = n.iter,
#       n.iter.final = n.iter.final,
#       X = X.orig,
#       Y = Y, X.order = X.order,
#       m = m, burnin = n.burnin,
#       tune = n.tune,
#       n.chains = n.chains,
#       n_par = r + r * (r + 1) / 2 + p + p * (p + 1) / 2 + r * m,
#       prior_param = all_MCs[[1]]$prior_param
#     ),
#     out_vars
#   )
#
#   class(out_res) <- c("Benvlp", "Benvlp_pred")
#
#   out_res
# }
#

Benvlp_MC_pred_gibbs_devel <- function(Y, X,
                                       m,
                                       n.iter = 1000,
                                       tau = 0.01,
                                       autotune = TRUE,
                                       tune.accpt.prop.lower = 0.2,
                                       tune.accpt.prop.upper = 0.3,
                                       tune.incr = 0.05,
                                       burnin.prop = 0.5,
                                       tune.burnin.prop = 0.75,
                                       tune_nterm = 25,
                                       tune_adj_scale_nterm = 50,
                                       checkstep = FALSE,
                                       starting = "mle",
                                       show_progress = TRUE,
                                       n.chains = 1,
                                       chains_parallel = FALSE,
                                       A_proposal = "rwmh",
                                       jacobian_only_gamma = FALSE,
                                       posterior_gamma_space = TRUE,
                                       nbreaks_trapezoid_A_prop = 5,
                                       n_iter_A_rwmh = 1,
                                       rwmh_elementwise = TRUE,
                                       rwmh_elementwise_fraction = 1,
                                       n_A_update_max = Inf,
                                       ...) {

  tt1 <- Sys.time()

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  stopifnot(A_proposal %in% c("stiefel-experimental", "rwmh", "rwmh_marginal"))


  dim_X <- dim(X)
  dim_Y <- dim(Y)

  n <- dim_X[1]
  p <- dim_X[2]
  r <- dim_Y[2]


  Y <- Y.orig <- data.matrix(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  Y.bar <- attr(Yc, "scaled:center")
  YctYc <- crossprod(Yc)


  X <- X.orig <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  X.bar <- attr(Xc, "scaled:center")
  XctXc <- crossprod(Xc)

  XctYc <- crossprod(Xc, Yc)

  Z <- cbind(X, Y)

  Z.bar <- c(X.bar, Y.bar)

  n.burnin <- ceiling(n.iter * burnin.prop)
  n.iter.final <- n.iter - n.burnin
  n.tune <- ceiling(n.burnin * tune.burnin.prop)
  burnin_iter <- seq_len(n.burnin)
  # log_2pi <- log(2*pi)

  if (m < 0 || m > p || m != as.integer(m)) {
    stop(paste0("\'m\' must be an integer between 0 and ", p, " (p)"))
  }

  X.order <- 1:p # changed if G1 is singular


  A_exists <- m < p & m > 0

  if (length(tau) == 1) {
    tau <- matrix(tau, p - m, m)
  }



  # prior parameters
  PsiY <- diag(1, r) * 1e-6
  nuY <- r
  PsiX <- diag(1, m) * 1e-6
  nuX <- m
  PsiX0 <- diag(p - m) * 1e-6
  nuX0 <- p - m

  # 기존 prior

  # e <- matrix(0, m, r)
  # is_e_zero <- TRUE
  # M <- 1e-6 * diag(1, nrow = m)
  # M.half <- diag(sqrt(diag(M)), nrow = m)
  # M.inv <- diag(1 / diag(M), nrow = m)
  # M.e <- M %*% e

  # 변경된 prior

  e <- matrix(0, p, r)
  is_e_zero <- TRUE
  M <- 1e-6 * diag(1, nrow = r)
  M.half <- diag(sqrt(diag(M)), nrow = r)
  M.inv <- diag(1 / diag(M), nrow = r)
  M.e <- M %*% t(e)
  psi_eta <- 1e10

  if (A_exists) {
    A0 <- matrix(0, p - m, m)
    K <- 1e6 * diag(1, p - m)
    K.inv <- 1e-6 * diag(1, p - m)
    K.half.inv <- 1e-3 * diag(1, p - m)
    L <- 1e6 * diag(1, m)
    L.inv <- 1e-6 * diag(1, m)
    L.half.inv <- 1e-3 * diag(1, m)
  }

  start_penv <- start_point_predenv(X, Y, m, 1)

  # basically attach(start_penv) but does not
  # affect the environment
  {
    for (x in names(start_penv)) {
      assign(x, start_penv[[x]])
    }
  }
  SigmaYcX.inv <-
    tryCatch({
      solve_chol(SigmaYcX)
    }, error = function(e) {
      message("Error computing SigmaYcX.inv")
      browser()
      stop(e)
    })

  Omega.inv <-
    tryCatch({
      if (dim(Omega)[1] == 0){
        NULL
      }else{
        solve_chol(Omega)
      }
    }, error = function(e) {
      message("Error computing Omega.inv")
      browser()
      stop(e)
    })
  Omega0.inv <-
    tryCatch({
      if (dim(Omega0)[1] == 0){
        NULL
      }else{
        solve_chol(Omega0)
      }
    }, error = function(e) {
      message("Error computing Omega0.inv")
      browser()
      stop(e)
    })

  gamma_gamma0 <- find_gammas_from_A(A)

  log_jacob_A_over_gamma_gamma0_end <-
    det_jacob_A_over_gamma_gamma0 <- 0

  if (m > 0) {
    eta.t <- t(eta)
    gamma.t <- t(gamma)
  }


  log_sqrt_2pi <- log(sqrt(2 * pi))

  # browser()

  runMC <- function(chain_no = 1) {
    runif(chain_no)
    one_n <- rep(1, n)
    # e.check <- crossprod(Xc, Yc) + e

    # needs to be recalucated at every iteration if A_exists
    # if (m > 0) {
    #   # M.check <- M + crossprod(Xc %*% gamma)
    #   # M.check_inv <- solve_chol(M.check)
    #   # M.check_half <- sqrtmat(M.check)
    #   M.tilde <- M + crossprod(X_muX %*% gamma)
    #   e.tilde <- e + crossprod(X_muX, Y_muY)
    #   M.tilde_inv <- solve_chol(M.tilde)
    #   H.tilde <- YctYc +
    #     crossprod(sqrtmatinv(M) %*% crossprod(gamma, e)) -
    #     crossprod(sqrtmatinv(M.check) %*% crossprod(gamma, e.check))
    #   eta.mean <- M.check_inv %*% crossprod(gamma, e.check)
    # } else {
    #   H.tilde <- YctYc
    # }


    if (m == 0) {
      eta <- 0
    }


    A_exists <- m > 0 & m < p

    muZ.all <- muX.all <- muY.all <- SigmaYcX.all <-
      SigmaX.all <- SigmaX.inv.all <-
      Omega.all <- Omega0.all <- eta.all <- beta.all <-
      gamma.all <- gamma0.all <- A.all <- accpt.A.all <-
      llik.contri.all <- vector("list", n.iter)

    # llik.all2 <- llik.all3 <-
    lpd.A.all <- accpt.A.ave <- llik.all <- lpd.full.all <- rep(0, n.iter)


    DeltaZ <- matrix(0, p + r, p + r) # just initializes, will be replaced at every iteration


    col_pairs_all <- combn(1:p, 2)
    # if (u <= r-u) {
    #   col_pairs_all <- combn(1:u, 2)
    # } else {
    #   col_pairs_all <- combn((u+1):r, 2)
    # }
    # col_pairs_all <- t(cbind(1:r, c(2:r, 1)))

    col_pairs_list <- lapply(
      1:ncol(col_pairs_all),
      function(jj) {
        col_pairs_all[, jj]
      }
    )

    if (m <= p - m) {
      col_pairs_list <- col_pairs_list[
        sapply(col_pairs_list, function(x) any(x %in% 1:m))
      ]
    } else {
      col_pairs_list <- col_pairs_list[
        sapply(col_pairs_list, function(x) any(x %in% (m + 1):p))
      ]
    }

    pb <- txtProgressBar(min = 1, max = n.iter, style = 3)

    for (iter in 1:n.iter) {
      # browser()

      # generate A, if 0 < m < p
      if (A_exists & (iter == 1 | iter %% n_iter_A_rwmh == 0)) {
        if (A_proposal == "stiefel-experimental") {

          SigmaYcX.inv.eta.t <- tcrossprod(SigmaYcX.inv, eta)
          Omega.inv_plus_eta.Syx.inv.eta.t.t <- Omega.inv +
            eta %*% SigmaYcX.inv.eta.t
          XctYc.SigmaYcX.inv.eta.t <- XctYc %*% SigmaYcX.inv.eta.t
          mh_step <- sample_A_pred_metrop_manifold_proposal(
            m = m,
            n = n,
            p = p,
            r = r,
            Omega.inv_plus_eta.Syx.inv.eta.t = Omega.inv_plus_eta.Syx.inv.eta.t.t,
            Omega.inv = Omega.inv,
            Omega0.inv = Omega0.inv,
            SigmaYcX.inv = SigmaYcX.inv,
            eta = eta,
            XctYc.SigmaYcX.inv.eta.t = XctYc.SigmaYcX.inv.eta.t,
            current_A = A,
            current_gamma_gamma0 = gamma_gamma0,
            current_det_jacob_A_over_gamma_gamma0 = det_jacob_A_over_gamma_gamma0,
            XctXc = XctXc,
            Xc = Xc,
            Yc = Yc,
            col_pairs_updt_set = col_pairs_list,
            A0 = A0,
            L.half.inv = L.half.inv,
            K.half.inv = K.half.inv,
            browser = (iter >= 1000),
            jacobian_only_gamma = jacobian_only_gamma,
            adjust_jacobian = !posterior_gamma_space,
            nbreaks_trapezoid_A_prop = nbreaks_trapezoid_A_prop,
            iter = iter,
            n_update_max = n_A_update_max
          )

          det_jacob_A_over_gamma_gamma0 <- mh_step$det_jacob_A_over_gamma_gamma0
        } else if (A_proposal == "rwmh") {
          lpd_A_curr <- lpd_A_pred_semi_marginal(
            A,
            XctXc = XctXc,
            Xc = Xc,
            Yc = Yc,
            Omega.inv = Omega.inv,
            SigmaYcX.inv = SigmaYcX.inv,
            eta = eta,
            Omega0.inv = Omega0.inv,
            A0 = A0,
            K.half.inv = K.half.inv,
            L.half.inv = L.half.inv
          )

          mh_step <- rwmh_colwise_A(
            A_start = A,
            lpd_start = lpd_A_curr,
            tau = tau,
            elementwise = rwmh_elementwise,
            n_update_max = n_A_update_max,
            n_rwmh_iter = n_iter_A_rwmh,
            random_update_fraction = rwmh_elementwise_fraction,
            lpd_func = lpd_A_pred_semi_marginal,
            XctXc = XctXc,
            Xc = Xc,
            Yc = Yc,
            Omega.inv = Omega.inv,
            SigmaYcX.inv = SigmaYcX.inv,
            eta = eta,
            Omega0.inv = Omega0.inv,
            A0 = A0,
            K.half.inv = K.half.inv,
            L.half.inv = L.half.inv
          )
        } else if (A_proposal == "rwmh_marginal") {

        }


        # autotune in burnin to make accpt prob between
        # tune.accpt.prop.lower and tune.accpt.prop.upper
        if (autotune &
            iter >= tune_nterm &
            iter %% 5 == 0 &
            iter <= n.tune &
            grepl("rwmh", A_proposal)) {
          # nterm <- tune_nterm
          # accpt.mean <- mean(accpt.A.ave[(iter-nterm+1):iter])
          # if (accpt.mean > tune.accpt.prop.upper) {
          #   tau <- tau * (1 + tune.incr)
          # } else if (accpt.mean < tune.accpt.prop.lower) {
          #   tau <- tau * (1 - tune.incr)
          # }
          tau <- autotune_param(
            draw_param_list = A.all[1:iter],
            tune_param = tau,
            accpt_list = accpt.A.all[1:iter],
            tune_nterm = tune_nterm,
            tune.incr = tune.incr,
            tune.accpt.prop.lower = tune.accpt.prop.lower,
            tune.accpt.prop.upper = tune.accpt.prop.upper,
            adjust_scale = ((iter %% tune_adj_scale_nterm) == 0)
          )
        }


        A.all[[iter]] <- A <- mh_step$A
        accpt.A.all[[iter]] <- mh_step$accpt # tcrossprod(rep(1, p-m), mh_step$accpt)
        accpt.A.ave[iter] <- mean(mh_step$accpt)
        lpd_A <- mh_step$lpd
        lpd.A.all[iter] <- c(lpd_A)
        gamma_gamma0 <- attr(lpd_A, "gamma_gamma0")
        gamma.all[[iter]] <- gamma <- gamma_gamma0$gamma
        gamma0.all[[iter]] <- gamma0 <- gamma_gamma0$gamma0
      }

      # browser()
      ################################
      # 이 부분 다시 정리해서 바꾸기 #
      ################################
      # browser()
      invCAtCA <- solve_chol(crossprod(rbind(diag(m), A)))
      # generate eta
      if (m > 0) {

        # 기존 prior

        # M.tilde <- M + crossprod(Xc %*% gamma)
        # e.tilde <- crossprod(gamma, XctYc) + M.e
        # M.tilde_inv <- solve_chol(M.tilde)
        # eta.mean <- M.tilde_inv %*% e.tilde
        # eta.all[[iter]] <- eta <-
        #   rMatrixNormal(eta.mean, M.tilde_inv, SigmaYcX)

        # 변경된 prior

        eta.Sig <- solve_chol(
          kronecker(M, 1/psi_eta * invCAtCA %*% Omega.inv %*% invCAtCA) + kronecker(SigmaYcX.inv, t(gamma) %*% XctXc %*% gamma)
        )
        eta.mean <- eta.Sig %*% matrix(t(SigmaYcX.inv %*% t(XctYc) %*% gamma), m*r, 1) # + M.e %*% gamma %*% Omega.inv %*% invCAtCA
        eta.all[[iter]] <- eta <-
          matrix(rmvnorm(1, eta.mean, eta.Sig), m, r)

      } else {
        eta.all[[iter]] <- eta <- NULL
      }

      # browser()

      ################################
      # 이 부분 다시 정리해서 바꾸기 #
      ################################

      # generate SigmaYcX
      if (m > 0) {

        # 기존 prior
        # PsiY.tilde <- PsiY +
        #   crossprod(Yc - Xc %*% gamma %*% eta) +
        #   crossprod(M.half %*% (eta - e))

        # 변경된 prior
        PsiY.tilde <- PsiY +
          crossprod(Yc - Xc %*% gamma %*% eta)

      } else {
        PsiY.tilde <- PsiY + YctYc
      }
      SigmaYcX.all[[iter]] <- SigmaYcX <-
        rinvwish(r, Phi = PsiY.tilde, nu = (nuY + n))
      SigmaYcX.inv <- solve_chol(SigmaYcX)

      # browser()

      ################################
      # 이 부분 다시 정리해서 바꾸기 #
      ################################

      # generate Omega if m > 0
      if (m > 0) {

        # 기존 prior
        # PsiX.tilde <- PsiX + crossprod(Xc %*% gamma)
        # Omega.all[[iter]] <-
        #   Omega <- rinvwish(m, Phi = PsiX.tilde, nu = nuX + n - 1)
        # Omega.inv <- solve_chol(Omega)

        # 변경된 prior
        PsiX.tilde <- PsiX + crossprod(Xc %*% gamma) + psi_eta * (invCAtCA %*% eta - t(gamma) %*% e) %*% M %*% t((invCAtCA %*% eta - t(gamma) %*% e))
        Omega.all[[iter]] <-
          Omega <- rinvwish(m, Phi = PsiX.tilde, nu = (nuX + n - r))
        Omega.inv <- solve_chol(Omega)

      } else {
        Omega.all[[iter]] <- Omega <- Omega.inv <- NULL
      }

      # browser()

      # generate Omega0 if m < p
      if (m < p) {
        PsiX0.tilde <- PsiX0 + crossprod(Xc %*% gamma0)
        Omega0.all[[iter]] <-
          Omega0 <- rinvwish((p - m), Phi = PsiX0.tilde, nu = (nuX0 + n - 1))
        Omega0.inv <- solve_chol(Omega0)
      } else {
        Omega0.all[[iter]] <- Omega0 <- Omega0.inv <- NULL
      }

      # calculate Sigma and Sigma.inv
      if (m > 0) {
        SigmaX.inv.1 <- gamma %*% tcrossprod(Omega.inv, gamma)
        SigmaX.inv.0 <- gamma0 %*% tcrossprod(Omega0.inv, gamma0)
        SigmaX.inv.all[[iter]] <- SigmaX.inv <-
          SigmaX.inv.1 + SigmaX.inv.0

        SigmaX.all[[iter]] <- SigmaX <-
          gamma %*% tcrossprod(Omega, gamma) +
          gamma0 %*% tcrossprod(Omega0, gamma0)
      } else {
        SigmaX.inv.0 <- Omega0.inv
        SigmaX.inv.1 <- NULL
        SigmaX.inv.all[[iter]] <- SigmaX.inv <-
          SigmaX.inv.0

        SigmaX.all[[iter]] <- SigmaX <- Omega0
      }

      # calculate beta
      if (m > 0) {
        beta.all[[iter]] <- beta <- gamma %*% eta
      } else if (m == 0) {
        beta.all[[iter]] <- beta <- matrix(0, nrow = p, ncol = r)
      }

      # browser()

      # generate muX and muY
      # generate muZ - then extract muX and muY
      # SigmaX <- gamma %*% tcrossprod(Omega, gamma) +
      #   gamma0 %*% tcrossprod(Omega0, gamma0)

      DeltaZ[1:p, 1:p] <- SigmaX
      if (m > 0) {
        Omega.eta <- Omega %*% eta
        DeltaZ[1:p, (p + 1):(p + r)] <- gamma %*% Omega.eta
        DeltaZ[(p + 1):(p + r), (p + 1):(p + r)] <- SigmaYcX + crossprod(eta, Omega.eta)
      } else {
        # DeltaZ[1:p, (p+1):(p+r)] <- 0
        # ^not need, initialized at 0
        DeltaZ[(p + 1):(p + r), (p + 1):(p + r)] <- SigmaYcX
      }

      DeltaZ[(p + 1):(p + r), 1:p] <- t(DeltaZ[1:p, (p + 1):(p + r)])
      DeltaZ.n <- DeltaZ / n

      #browser()

      if (!all(eigen(DeltaZ.n)$values > 0)){
        DeltaZ.n <- make_positive_definite(DeltaZ.n)
      }

      muZ.all[[iter]] <- muZ <- rmvnorm(dim = p + r, mu = Z.bar, sigma = DeltaZ.n)
      muX.all[[iter]] <- muX <- muZ[1:p]
      muY.all[[iter]] <- muY <- muZ[(p + 1):(r + p)]

      X_muX <- X - tcrossprod(one_n, muX)
      Y_muY <- Y - tcrossprod(one_n, muY)

      # browser()

      # # if A_exists,  recalculate M.check, M.check_inv, and H.tilde --
      # if (A_exists) {
      #   M.check <- M + crossprod(Xc %*% gamma)
      #   M.check_inv <- solve_chol(M.check)
      #   M.check_half <- sqrtmat(M.check)
      #   H.tilde <- YctYc +
      #     crossprod(M.half %*% crossprod(gamma, e)) -
      #     crossprod(sqrtmatinv(M.check) %*% crossprod(gamma, e.check))
      # }

      log_det_SigmaYcX <- log(det(SigmaYcX))
      if (m > 0) log_det_Omega <- log(det(Omega))
      if (m < p) log_det_Omega0 <- log(det(Omega0))
      resi <- Y_muY - X_muX %*% beta
      llik.contri.1 <- rowSums((resi %*% SigmaYcX.inv) * resi)
      llik.contri.2 <- rowSums((X_muX %*% SigmaX.inv) * X_muX)


      # llikelihood contribution
      if (A_exists) {
        llik.contri.all[[iter]] <- llik.contri <-
          -0.5 * (p + r) * log_sqrt_2pi +
          -0.5 * (log_det_SigmaYcX + log_det_Omega + log_det_Omega0 +
                    llik.contri.1 + llik.contri.2)
      }
      # else if (m == p) {
      #   llik.contri.all[[iter]] <- llik.contri <-
      #     - 0.5 * n * (p + r) * log_sqrt_2pi +
      #     - 0.5 * (log_det_SigmaYcX + log_det_Omega +
      #                llik.contri.1 + llik.contri.2)
      # }
      else {
        llik.contri.all[[iter]] <- llik.contri <-
          -0.5 * n * (p + r) * log_sqrt_2pi +
          -0.5 * (log_det_SigmaYcX + log_det_Omega0 +
                    llik.contri.1 + llik.contri.2)
      }

      # total log likelihood
      llik.all[[iter]] <- llik <- sum(llik.contri)


      # # log posterior density
      # if (A_exists) {
      #   lpd.full.all[[iter]] <- lpd_full <-
      #     llik - -0.5 * (
      #       (nuY + r + p + 1) * log_det_SigmaYcX +
      #         sum(SigmaYcX.inv * PsiY) +
      #         sum(SigmaYcX.inv *
      #           crossprod(M.half %*% (eta))) +
      #         (nuX + m + 1) * log_det_Omega +
      #         sum(Omega.inv * PsiX) +
      #         (nuX0 + p - m + 1) * log_det_Omega0 +
      #         sum(Omega0.inv * PsiX0) +
      #         sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
      #     )
      # } else if (m == p) {
      #   lpd.full.all[[iter]] <- lpd_full <-
      #     llik - -0.5 * (
      #       (nuY + r + p + 1) * log_det_SigmaYcX +
      #         sum(SigmaYcX.inv * PsiY) +
      #         sum(SigmaYcX.inv *
      #           crossprod(M.half %*% (eta))) +
      #         (nuX + m + 1) * log_det_Omega +
      #         sum(Omega.inv * PsiX)
      #     )
      # } else {
      #   lpd.full.all[[iter]] <- lpd_full <-
      #     llik - -0.5 * (
      #       (nuY + r + p + 1) * log_det_SigmaYcX +
      #         sum(SigmaYcX.inv * PsiY) +
      #         (nuX0 + p - m + 1) * log_det_Omega0 +
      #         sum(Omega0.inv * PsiX0)
      #     )
      # }


      # tune A, if exists & proposal is rwmh
      # if (A_exists & autotune &
      #     iter >= tune_nterm & iter %% 5 == 0 & iter <= n.tune &
      #     A_proposal != "stiefel") {
      #
      #
      #   # accpt_A_col <- Reduce(
      #   #   '+',
      #   #   accpt.A.all[(iter-tune_nterm+1):iter]
      #   # ) / tune_nterm
      #   #
      #   # for(j in 1:m) {
      #   #   ave_accpt <- accpt_A_col[j]
      #   #
      #   #   if (ave_accpt > tune.accpt.prop.upper) {
      #   #     tau[, j] <- tau[, j] * (1 + tune.incr)
      #   #   } else if (ave_accpt < tune.accpt.prop.lower) {
      #   #     tau[, j] <- tau[, j] * (1 - tune.incr)
      #   #   }
      #   #
      #   #   if (iter %% tune_nterm == 0) {
      #   #     # last nterms iterations of A[, j]
      #   #     # each column of A_j_last_nterm is
      #   #     # from a single interation of A[, j]
      #   #     A_j_last_nterm <- do.call(
      #   #       cbind,
      #   #       lapply(A.all[(iter-tune_nterm+1):iter],
      #   #              function(x) x[, j]))
      #   #     par_sd <- apply(A_j_last_nterm, 1, sd)
      #   #     par_sd_mean <- sum(par_sd) / (p - m)
      #   #     tau_j_mean <- sum(tau[, j]) /(p - m)
      #   #
      #   #     if(par_sd_mean > 0)
      #   #       tau[, j] <- 0.5*par_sd/par_sd_mean*tau_j_mean +
      #   #       0.5*tau[, j]
      #   #   }
      #   # }
      # }


      # if (iter == 100) browser()

      if (show_progress) {
        utils::setTxtProgressBar(pb, iter)
      }
    }

    MC <- list( # "n.iter" = n.iter, "X" = X, "Y" = Y, "X.order" = X.order, "m" = m,
      # "burnin" = burnin,
      muX = muX.all[-burnin_iter],
      muY = muY.all[-burnin_iter],
      eta = eta.all[-burnin_iter],
      SigmaYcX = SigmaYcX.all[-burnin_iter],
      Omega = Omega.all[-burnin_iter],
      Omega0 = Omega0.all[-burnin_iter],
      A = A.all[-burnin_iter],
      beta = beta.all[-burnin_iter],
      lpd.A = lpd.A.all[-burnin_iter],
      gamma = gamma.all[-burnin_iter],
      gamma0 = gamma0.all[-burnin_iter],
      accpt.A.ave = accpt.A.ave[-burnin_iter],
      # lpd.full = lpd.full.all[-burnin_iter],
      llik = llik.all[-burnin_iter],
      # "llik2" = llik.all2[-burnin_iter],
      # "llik3" = llik.all3[-burnin_iter],
      llik.contri = llik.contri.all[-burnin_iter],
      accpt.A.all = accpt.A.all[-burnin_iter],
      tau = tau,
      prior_param = list()
    )

    MC
  }

  lapply_ <- function(...) {
    if (n.chains == 1 | !chains_parallel) {
      lapply(...)
    } else {
      future.apply::future_apply(...,
                                 future.seed = TRUE
      )
    }
  }

  all_MCs <- lapply_(
    1:n.chains,
    function(j) {
      runMC(chain_no = j)
    }
  )


  # reshape the final list
  # each parameter, e.g. Omega is an
  # n.chain long list of n.iter.final long lists.

  vars <- setdiff(names(all_MCs[[1]]), "prior_param")

  out_vars <- lapply(
    vars,
    function(x) {
      out <- lapply(all_MCs, "[[", x)
      names(out) <- paste0("Chain_", 1:n.chains)
      out
    }
  )
  names(out_vars) <- vars

  tt2 <- Sys.time()

  total_time <- as.numeric(difftime(tt2, tt1, units = "secs"))

  out_res <- c(
    list(
      n.iter = n.iter,
      n.iter.final = n.iter.final,
      X = X.orig,
      Y = Y, X.order = X.order,
      m = m, burnin = n.burnin,
      tune = n.tune,
      n.chains = n.chains,
      n_par = r + r * (r + 1) / 2 + p + p * (p + 1) / 2 + r * m,
      prior_param = all_MCs[[1]]$prior_param,
      total_time = total_time
    ),
    out_vars
  )

  class(out_res) <- c("Benvlp", "Benvlp_pred")

  out_res
}

# samplecode <- Benvlp_MC_pred_gibbs_devel(
#   X = X,
#   Y = Y,
#   m = 2,
#   n.iter = 10000,
#   burnin.prop = 1/3,
#   n.chains = 1,
#   strating = "map",
#   A_proposal = "rwmh", #stiefel-experimental
#   jacobian_only_gamma = FALSE,
#   compute_llik = TRUE
# )
