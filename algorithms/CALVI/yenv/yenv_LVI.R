response_LVI <- function(
    X, Y, u,
    maxiter       = 200,
    n_iter        = 160,
    show_progress = FALSE,
    tol           = 1e-6,
    make_positive = TRUE,
    init          = c("map","mle","random"),
    compute_llik  = FALSE,
    # prior / 기타
    nuvalue = 0, nu0value = 0, Avalue = 1e6, Mvalue = 1e6,
    ...
){
  ## --------- 안전장치 & 기본 준비 ---------
  init <- match.arg(init)
  X <- data.matrix(X); Y <- data.matrix(Y)
  stopifnot(nrow(X) == nrow(Y))
  n <- nrow(X); p <- ncol(X); r <- ncol(Y)
  
  ## positive-definite 보정 헬퍼
  .is_pos_def <- function(M) {
    if (!is.matrix(M)) return(FALSE)
    if (!isSymmetric(M)) return(FALSE)
    ev <- tryCatch(eigen(M, only.values = TRUE)$values, error=function(e) return(NA))
    if (any(is.na(ev))) return(FALSE)
    all(ev > 0)
  }
  .make_pos_def <- function(M, eps = 1e-8) {
    if (!is.matrix(M)) stop("Matrix required")
    M <- (M + t(M))/2
    ev <- eigen(M)
    ev$values[ev$values < eps] <- eps
    M2 <- ev$vectors %*% (diag(ev$values, nrow(M))) %*% t(ev$vectors)
    (M2 + t(M2))/2
  }
  
  ## 진행바: 텍스트 바만 사용 (headless에서 안전)
  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  ## 공통 전처리
  tYY   <- crossprod(Y)
  Yc    <- scale(Y, center = TRUE, scale = FALSE)
  Y_bar <- attr(Yc, "scaled:center")
  
  Xc    <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X_bar <- attr(Xc, "scaled:center")
  
  Z <- matrix(0, nrow = (r - u), ncol = u)
  K <- rbind(diag(1, u), Z)
  L <- rbind(t(Z), diag(1, r - u))
  if ((r - u) == 1L || u == 1L) {
    ## K_{r-u,u} = I
    KK <- diag((r - u) * u)
  } else {
    KK <- commutation.matrix(r - u, u)
  }
  prior_pars <- response_get_prior(r, p, u, nuvalue, nu0value, Avalue, Mvalue)
  B0 <- prior_pars$B0
  M  <- prior_pars$M
  
  BMB          <- B0 %*% M %*% t(B0)
  tXcY_MtB0    <- t(Xc) %*% Y + M %*% t(B0)
  XctXc_M      <- t(XctXc) + M
  inv_XctXc_M  <- solve(XctXc_M)
  
  onevector <- matrix(1, n, 1)
  
  ## 초기값
  if (init == "mle") {
    variational_pars <- init_pars <- response_get_init_MLE(X, Y, n, r, p, u)
  } else if (init == "random") {
    variational_pars <- init_pars <- response_get_init_MLE_with_noise(X, Y, n, r, p, u, sd = 10)
  } else { # "map"
    variational_pars <- init_pars <- do.call(
      response_get_init_MAP,
      c(list(X = X, Y = Y, u = u, n = n, r = r, p = p,
             BMB = BMB, tXcY_MtB0 = tXcY_MtB0, XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M),
        prior_pars)
    )
  }
  
  Lvaluelist <- list()
  best_Lvalue <- -Inf
  i <- 1L
  
  start_time <- Sys.time()
  
  while (i < n_iter) {
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    
    G <- crossprod(Y - onevector %*% t(variational_pars[["mu_q"]])) + n * variational_pars[["Sigma_q"]]
    
    if (u == 0) {
      ## A, eta는 의미 없으므로 0/단위 구성
      variational_pars[["hatA"]] <- matrix(0, (r - u), u)
      variational_pars[["HA"]]   <- matrix(0, (r - u) * u, (r - u) * u)
      variational_pars[["CA"]]   <- matrix(0, r, u)
      variational_pars[["DA"]]   <- diag(1, r)
      
      variational_pars[["nu_q"]]    <- n + p
      variational_pars[["invPsi_q"]] <- matrix(0, u, u)
      variational_pars[["eta_q"]]   <- matrix(0, u, p)
      variational_pars[["U_q"]]     <- matrix(0, u, u)
      variational_pars[["V_q"]]     <- matrix(0, p, p)
      
      Lvalue0 <- do.call(response_convergence,
                         c(prior_pars, variational_pars,
                           list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      
      tilde_Omega0_q <- do.call(response_UpdateOmega0,
                                c(prior_pars, variational_pars,
                                  list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                       G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                       XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["nu0_q"]] <- tilde_Omega0_q[[1]]
      tmpinvPsi0_q <- tilde_Omega0_q[[2]]
      if (make_positive && !.is_pos_def(tmpinvPsi0_q)) tmpinvPsi0_q <- .make_pos_def(tmpinvPsi0_q)
      variational_pars[["invPsi0_q"]] <- tmpinvPsi0_q
      
      tilde_mu_q <- do.call(response_Updatemu,
                            c(prior_pars, variational_pars,
                              list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                   G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                   XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["mu_q"]]    <- tilde_mu_q[[1]]
      variational_pars[["Sigma_q"]] <- tilde_mu_q[[2]]
      
      G <- crossprod(Y - onevector %*% t(variational_pars[["mu_q"]])) + n * variational_pars[["Sigma_q"]]
      
      Lvalue <- do.call(response_convergence,
                        c(prior_pars, variational_pars,
                          list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                               G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                               XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      
      Lvaluelist <- c(Lvaluelist, Lvalue)
      if (abs(Lvalue0 - Lvalue) < tol * abs(Lvalue0)) break
      Lvalue0 <- Lvalue; i <- i + 1L
      
    } else if (u < r) {
      Lvalue0 <- do.call(response_convergence,
                         c(prior_pars, variational_pars,
                           list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      
      hatA <- do.call(response_UpdateA.optim,
                      c(list(maxiter = maxiter, method = "BFGS"),
                        prior_pars, variational_pars,
                        list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                             G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                             XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L, KK = KK)))
      variational_pars[["hatA"]] <- hatA[[1]]
      tmpHA <- hatA[[3]]
      if (make_positive && !.is_pos_def(tmpHA)) tmpHA <- .make_pos_def(tmpHA)
      variational_pars[["HA"]] <- tmpHA
      variational_pars[["CA"]] <- hatA[[4]]
      variational_pars[["DA"]] <- hatA[[5]]
      
      tilde_eta_q <- do.call(response_Updateeta,
                             c(prior_pars, variational_pars,
                               list(X = X, Y = Y, Xc = Xc, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                    G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                    XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["eta_q"]] <- tilde_eta_q[[1]]
      variational_pars[["U_q"]]   <- tilde_eta_q[[2]]
      variational_pars[["V_q"]]   <- tilde_eta_q[[3]]
      
      tilde_Omega_q <- do.call(response_UpdateOmega,
                               c(prior_pars, variational_pars,
                                 list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                      G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                      XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["nu_q"]] <- tilde_Omega_q[[1]]
      tmpinvPsi_q <- tilde_Omega_q[[2]]
      if (make_positive && !.is_pos_def(tmpinvPsi_q)) tmpinvPsi_q <- .make_pos_def(tmpinvPsi_q)
      variational_pars[["invPsi_q"]] <- tmpinvPsi_q
      
      tilde_Omega0_q <- do.call(response_UpdateOmega0,
                                c(prior_pars, variational_pars,
                                  list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                       G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                       XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["nu0_q"]] <- tilde_Omega0_q[[1]]
      tmpinvPsi0_q <- tilde_Omega0_q[[2]]
      if (make_positive && !.is_pos_def(tmpinvPsi0_q)) tmpinvPsi0_q <- .make_pos_def(tmpinvPsi0_q)
      variational_pars[["invPsi0_q"]] <- tmpinvPsi0_q
      
      tilde_mu_q <- do.call(response_Updatemu,
                            c(prior_pars, variational_pars,
                              list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                   G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                   XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["mu_q"]]    <- tilde_mu_q[[1]]
      variational_pars[["Sigma_q"]] <- tilde_mu_q[[2]]
      
      G <- crossprod(Y - onevector %*% t(variational_pars[["mu_q"]])) + n * variational_pars[["Sigma_q"]]
      
      Lvalue <- do.call(response_convergence,
                        c(prior_pars, variational_pars,
                          list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                               G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                               XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      Lvaluelist <- c(Lvaluelist, Lvalue)
      if (abs(Lvalue0 - Lvalue) < tol * abs(Lvalue0)) break
      Lvalue0 <- Lvalue; i <- i + 1L
      
    } else { # u == r
      variational_pars[["hatA"]] <- matrix(0, (r - u), u)
      variational_pars[["HA"]]   <- matrix(0, (r - u) * u, (r - u) * u)
      variational_pars[["CA"]]   <- diag(1, r)
      variational_pars[["DA"]]   <- matrix(0, r, (r - u))
      
      variational_pars[["nu0_q"]]    <- n
      variational_pars[["invPsi0_q"]] <- matrix(0, (r - u), (r - u))
      
      Lvalue0 <- do.call(response_convergence,
                         c(prior_pars, variational_pars,
                           list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      
      tilde_eta_q <- do.call(response_Updateeta,
                             c(prior_pars, variational_pars,
                               list(X = X, Y = Y, Xc = Xc, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                    G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                    XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["eta_q"]] <- tilde_eta_q[[1]]
      variational_pars[["U_q"]]   <- tilde_eta_q[[2]]
      variational_pars[["V_q"]]   <- tilde_eta_q[[3]]
      
      tilde_Omega_q <- do.call(response_UpdateOmega,
                               c(prior_pars, variational_pars,
                                 list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                      G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                      XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["nu_q"]] <- tilde_Omega_q[[1]]
      tmpinvPsi_q <- tilde_Omega_q[[2]]
      if (make_positive && !.is_pos_def(tmpinvPsi_q)) tmpinvPsi_q <- .make_pos_def(tmpinvPsi_q)
      variational_pars[["invPsi_q"]] <- tmpinvPsi_q
      
      tilde_mu_q <- do.call(response_Updatemu,
                            c(prior_pars, variational_pars,
                              list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                                   G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                                   XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      variational_pars[["mu_q"]]    <- tilde_mu_q[[1]]
      variational_pars[["Sigma_q"]] <- tilde_mu_q[[2]]
      
      G <- crossprod(Y - onevector %*% t(variational_pars[["mu_q"]])) + n * variational_pars[["Sigma_q"]]
      
      Lvalue <- do.call(response_convergence,
                        c(prior_pars, variational_pars,
                          list(X = X, Y = Y, X_bar = X_bar, Y_bar = Y_bar, u = u, n = n, r = r, p = p,
                               G = G, BMB = BMB, tXcY_MtB0 = tXcY_MtB0,
                               XctXc_M = XctXc_M, inv_XctXc_M = inv_XctXc_M, K = K, L = L)))
      Lvaluelist <- c(Lvaluelist, Lvalue)
      if (abs(Lvalue0 - Lvalue) < tol * abs(Lvalue0)) break
      Lvalue0 <- Lvalue; i <- i + 1L
    }
  }
  
  end_time <- Sys.time()
  last_Lvalue <- Lvalue
  
  # beta 구성
  if (u == 0) {
    variational_pars[["beta"]] <- matrix(0, r, p)
  } else {
    variational_pars[["beta"]] <- variational_pars[["CA"]] %*%
      solve(crossprod(variational_pars[["CA"]])) %*% variational_pars[["eta_q"]]
  }
  
  # (선택) 대략적 log-likelihood 계산
  if (isTRUE(compute_llik)) {
    llik_vec <- yenv_llik(
      X, Y, N = 1,
      mu_q       = variational_pars$mu_q,
      Sigma_q    = variational_pars$Sigma_q,
      eta_q      = variational_pars$eta_q,
      U_q        = variational_pars$U_q,
      V_q        = variational_pars$V_q,
      nu_q       = variational_pars$nu_q,
      invPsi_q   = variational_pars$invPsi_q,
      nu0_q      = variational_pars$nu0_q,
      invPsi0_q  = variational_pars$invPsi0_q,
      hatA       = variational_pars$hatA,
      HA         = variational_pars$HA
    )
    llik_mc <- llik_vec
    LVI_BIC <- -2 * llik_mc + log(n) * (r + p*u + r*(r + 1)/2)
  } else {
    llik_mc <- NA
    LVI_BIC <- NA
  }
  
  # out <- calc_curvature_bounds(
  #   A = variational_pars$hatA,
  #   n = n,
  #   G = G,
  #   BMB = BMB, psi = prior_pars$psi, psi0 = prior_pars$psi0,
  #   nu = prior_pars$nu, nu0 = prior_pars$nu0,
  #   nu_q = variational_pars$nu_q, invPsi_q = variational_pars$invPsi_q,
  #   nu0_q = variational_pars$nu0_q, invPsi0_q = variational_pars$invPsi0_q
  # )
  
  total_run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    variational_pars,
    init_pars,
    total_run_time,
    last_Lvalue,
    Lvaluelist,
    llik_mc,
    LVI_BIC
    # out$lambda_min_1n_Jn,
    # out$lambda_max_1n_d2Phi
  )
}
