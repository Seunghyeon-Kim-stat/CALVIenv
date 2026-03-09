
library(progress)
library(matrixcalc)
library(dplyr)

response_UpdateA.optim <- function(maxiter = 100, method = "BFGS",
                                   nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, HA0=NULL, # prior
                                   hatA, CA, DA, HA,
                                   nu_q, invPsi_q,
                                   nu0_q, invPsi0_q,
                                   eta_q, U_q, V_q,
                                   mu_q, Sigma_q, # variational parameters
                                   X, Y, X_bar, Y_bar, u, n, p, r, # data
                                   G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, K, L, KK,... # common terms
                                   ){

  # common part

  C1 <- nu_q * invPsi_q
  C2 <- G + BMB + psi*diag(r)
  C3 <- nu0_q * invPsi0_q
  C4 <- G + psi0*diag(r)
  C5 <- nu_q * invPsi_q %*% eta_q %*% tXcY_MtB0

  o_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
  o_term2 <- -1/2 * tr( C3 %*% t(DA) %*% C4 %*% DA )
  o_term3 <- tr( C5 %*% CA )

  invU_A <- solve(U_A)
  invV_A <- solve(V_A)

  o_term4 <- -1/2 * tr(invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))

  o_logdet <- ((2*n + nu + nu0) / 2) * logdetfunction(diag(u) + crossprod(hatA))

  obj1 <- o_term1 + o_term2 + o_term3 + o_term4 + o_logdet
  obj0 <- obj1 <- obj1 / n

  i <- 1

  while (i < maxiter) {

    ftol <- 1e-3

    for (m in (u + 1):r) {
      j <- m

      fA <- function(x){

        hatA[j-u,] <- x
        CA <- rbind(diag(u), hatA)
        DA <- rbind(-t(hatA), diag(r-u))

        f_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
        f_term2 <- -1/2 * tr( C3 %*% t(DA) %*% C4 %*% DA )
        f_term3 <- tr( C5 %*% CA )
        f_term4 <- -1/2 * tr(invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))

        f_logdet <- ((2*n + nu + nu0) / 2) * logdetfunction(diag(u) + crossprod(hatA))

        fA_value <- f_term1 + f_term2 + f_term3 + f_term4 + f_logdet

        return(-fA_value)
      }

      gA <- function(x){

        hatA[j-u,] <- x
        CA <- rbind(diag(u), hatA)
        DA <- rbind(-t(hatA), diag(r-u))

        d_term1 <- - t(L) %*% C2 %*% CA %*% C1
        d_term2 <- C3 %*% t(DA) %*% C4 %*% K
        d_term3 <- t(C5 %*% L)
        d_term4 <- - invU_A %*% (hatA - A0) %*% invV_A

        d_logdet <- (2*n + nu + nu0) * solve((diag(r-u) + hatA %*% t(hatA)), hatA)

        gA_value <- (d_term1 + d_term2 + d_term3 + d_term4 + d_logdet)[j-u,]

        return(-gA_value)
      }

      res <- stats::optim(
        hatA[j-u, ],
        fn = fA,
        gr = gA,
        method = "BFGS"
      )
      hatA[j-u, ] <- res$par

    }

    CA <- rbind(diag(u), hatA)
    DA <- rbind(-t(hatA), diag(r-u))

    o_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
    o_term2 <- -1/2 * tr( C3 %*% t(DA) %*% C4 %*% DA )
    o_term3 <- tr( C5 %*% CA )
    o_term4 <- -1/2 * tr(invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))

    o_logdet <- ((2*n + nu + nu0) / 2) * logdetfunction(diag(u) + crossprod(hatA))

    obj2 <- o_term1 + o_term2 + o_term3 + o_term4 + o_logdet
    obj2 <- obj2 / n

    if (abs(obj0 - obj2) < ftol * abs(obj0)) {
      break
    } else {
      obj0 <- obj2
      i <- i + 1
    }

  }
  # Calculate Hessian
  
  invJ <- solve(diag(r-u) + tcrossprod(hatA))

  h_term1 <- - kronecker(t(C1), t(L) %*% C2 %*% L)
  h_term2 <- - kronecker(t(K) %*% t(C4) %*% K, C3)
  h_term4 <- - kronecker(t(invV_A), invU_A)

  #h_logdet <- kronecker((diag(u) - t(hatA) %*% invJ %*% hatA), invJ) - tcrossprod(vec(invJ %*% hatA))
  
  h_logdet <- if ((r-u) == 0L) {
    
    ## r == u : dimension 0 x 0
    matrix(0, 0, 0)
    
  } else {
    
    IA <- diag(u) - t(hatA) %*% invJ %*% hatA      # u x u
    B  <- invJ %*% hatA                            # (r-u) x u
    
    term1 <- kronecker(IA, invJ)                   # (u*ru) x (u*ru)
    term2 <- kronecker(t(B), B)                    # (u*ru) x (ru*u)
    
    if ((r-u) == 1L || u == 1L) {
      ## commutation.matrix(ru, u) = I
      term1 - term2
    } else {
      term1 - term2 %*% KK
    }
  }
  
  HA <- h_term1 + h_term2 + h_term4 + (2*n + nu + nu0)*h_logdet
  HA <- -solve(HA)
  HA <- (HA + t(HA)) / 2

  CA <- rbind(diag(u), hatA)
  DA <- rbind(-t(hatA), diag(r-u))

  return(list(hatA, FA=obj2, HA, CA, DA))
}

response_UpdateOmega <- function(nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, # prior
                                 hatA, CA, DA, HA,
                                 nu_q, invPsi_q,
                                 nu0_q, invPsi0_q,
                                 eta_q, U_q, V_q,
                                 mu_q, Sigma_q, # variational distribution
                                 X, Y, X_bar, Y_bar, u, n, p, r, # data
                                 G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, K, L, ... # common terms
                                 ){

  LPL <- t(L) %*% (G+BMB+psi*diag(r)) %*% L

  nu_q <- n+p+nu

  tilde_Psi_q <- divide_matrix(HA, LPL) + t(CA) %*% (G+BMB+psi*diag(r)) %*% CA -
    2 * eta_q %*% tXcY_MtB0 %*% CA + tr(XctXc_M %*% V_q) * U_q + eta_q %*% XctXc_M %*% t(eta_q)

  invPsi_q <- solve(tilde_Psi_q)
  invPsi_q <- (t(invPsi_q) + invPsi_q) / 2

  return(list(nu_q, invPsi_q))

}

response_UpdateOmega0 <- function(nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, # prior
                                  hatA, CA, DA, HA,
                                  nu_q, invPsi_q,
                                  nu0_q, invPsi0_q,
                                  eta_q, U_q, V_q,
                                  mu_q, Sigma_q, # variational distribution
                                  X, Y, X_bar, Y_bar, u, n, p, r, # data
                                  G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, K, L, ... # common terms
                                  ){

  KPK <- t(K)%*%(G+psi0*diag(r))%*%K

  nu0_q <- n+nu0

  tilde_Psi0_q <- divide_matrix(HA, KPK, commutation = TRUE) + t(DA)%*%(G+psi0*diag(r))%*%DA

  invPsi0_q <- solve(tilde_Psi0_q)
  invPsi0_q <- (t(invPsi0_q) + invPsi0_q) / 2

  return(list(nu0_q, invPsi0_q))

}

response_Updateeta <- function(nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, # prior
                               hatA, CA, DA, HA,
                               nu_q, invPsi_q,
                               nu0_q, invPsi0_q,
                               eta_q, U_q, V_q,
                               mu_q, Sigma_q, # variational distribution
                               X, Y, X_bar, Y_bar, u, n, p, r, # data
                               G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, K, L, ... # common terms
                               ){

  U_q <- solve(invPsi_q)/(nu_q)
  V_q <- inv_XctXc_M

  tilde_eta_q <- t(tXcY_MtB0 %*% CA) %*% V_q

  return(list(tilde_eta_q, U_q, V_q))

}

response_Updatemu <- function(nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, # prior
                              hatA, CA, DA, HA,
                              nu_q, invPsi_q,
                              nu0_q, invPsi0_q,
                              eta_q, U_q, V_q,
                              mu_q, Sigma_q, # variational distribution
                              X, Y, X_bar, Y_bar, u, n, p, r, # data
                              G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M, K, L, ... # common terms
                              ){

  if (u ==0){

    # S1 <- n * nu_q * (CA%*%invPsi_q%*%t(CA))
    S2 <- n * nu0_q * (DA%*%invPsi0_q%*%t(DA))

    inv_tilde_Sigma_q <- S2 # + solve(Sigma0)
    tilde_Sigma_q <- solve(inv_tilde_Sigma_q)

    tilde_mu_q <- Y_bar

  }
  else if(u < r){

    LULTtrVinvPsi_q <- L %*% divide_matrix(HA, invPsi_q, commutation=TRUE) %*% t(L)
    KVKTtrUinvPsi0_q <- K %*% divide_matrix(HA, invPsi0_q) %*% t(K)

    S1 <- n * nu_q * (CA%*%invPsi_q%*%t(CA) + LULTtrVinvPsi_q)
    S2 <- n * nu0_q * (DA%*%invPsi0_q%*%t(DA) + KVKTtrUinvPsi0_q)

    inv_tilde_Sigma_q <- S1 + S2 # + solve(Sigma0)
    tilde_Sigma_q <- solve(inv_tilde_Sigma_q)

    tilde_mu_q <- Y_bar

  }
  else if (u == r) {

    S1 <- n * nu_q * (CA%*%invPsi_q%*%t(CA))
    # S2 <- n * nu0_q * (DA%*%invPsi0_q%*%t(DA))

    inv_tilde_Sigma_q <- S1 # + S2 + solve(Sigma0)
    tilde_Sigma_q <- solve(inv_tilde_Sigma_q)

    tilde_mu_q <- Y_bar

  }

  return(list(tilde_mu_q, tilde_Sigma_q))

}

response_convergence <- function(nu, psi, nu0, psi0, B0, M, A0, U_A=NULL, V_A=NULL, HA0=NULL, # prior
                                 hatA, CA, DA, HA,
                                 nu_q, invPsi_q,
                                 nu0_q, invPsi0_q,
                                 eta_q, U_q, V_q,
                                 mu_q, Sigma_q, # variational distribution
                                 X, Y, X_bar, Y_bar, u, n, p, r, # data
                                 G, BMB, tXcY_MtB0, XctXc_M, inv_XctXc_M,  K, L,... # common terms
                                 ){

  Aprior <- matrix_normal_logconst(U_A, V_A)

  const <-
    Aprior +
    invWishart_logconst(nu, psi, u) +
    invWishart_logconst(nu0, psi0, r-u) +
    - 0.5 * u * p * log(2 * pi) +
    0.5 * u * logdetfunction(M)

  C1 <- nu_q * invPsi_q
  C2 <- G + BMB + psi*diag(r)
  C3 <- nu0_q * invPsi0_q
  C4 <- G + psi0*diag(r)
  C5 <- nu_q * invPsi_q %*% eta_q %*% tXcY_MtB0

  f_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
  f_term2 <- -1/2 * tr( C3 %*% t(DA) %*% C4 %*% DA )
  f_term3 <- tr( C5 %*% CA )

  if (u < r & u > 0){
    invU_A <- solve(U_A)
    invV_A <- solve(V_A)
    f_term4 <- -1/2 * tr(invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))
  }else{
    f_term4 <- 0
  }

  f_logdet <- ((2*n + nu + nu0) / 2) * logdetfunction(diag(r-u) + tcrossprod(hatA))

  fA <-  - 0.5 * (nu_q + u + 1) * expectation_loginvWishart(nu_q, invPsi_q) -
    0.5 * (nu0_q + r - u + 1) * expectation_loginvWishart(nu0_q, invPsi0_q) -
    0.5 * nu_q * tr((tr(t(XctXc_M) %*% V_q) * U_q + eta_q%*%XctXc_M%*%t(eta_q))%*%invPsi_q) +
    f_term1 + f_term2 + f_term3 + f_term4 + f_logdet

  Entropy_mu_q <- normal_entropy(mu_q, Sigma_q)

  Entropy_eta_q <- matrix_normal_entropy(M = eta_q, U = U_q, V = V_q)

  Entropy_Omega_q <- invWishart_entropy(nu_q, invPsi_q)

  Entropy_Omega0_q <- invWishart_entropy(nu0_q, invPsi0_q)

  Entropy_A_q <- normal_entropy(hatA, HA)

  Lvalue <- const + fA - 1/2 * (r-u) * u +
    Entropy_mu_q + Entropy_eta_q + Entropy_Omega_q + Entropy_Omega0_q + Entropy_A_q

  return(Lvalue)

}

yenv_llik <- function(X, Y, N,
                      mu_q, Sigma_q,
                      eta_q, U_q, V_q,
                      nu_q, invPsi_q,
                      nu0_q, invPsi0_q,
                      hatA, HA) {

  X <- data.matrix(X)
  Y <- data.matrix(Y)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X)
  r <- ncol(Y)
  u <- ncol(hatA)

  log_weights <- numeric(N)

  for (iter in 1:N) {

    if (u == 0){
      # reduced model
      # non eta, A, Omega, X

      # Draw from variational posterior

      # mu <- matrix(mvtnorm::rmvnorm(n = 1, mean = mu_q, sigma = Sigma_q), r, 1)
      # Omega0 <- CholWishart::rInvWishart(n = 1, df = nu0_q, Sigma = solve(invPsi0_q))[,,1]

      mu <- mu_q
      Omega0 <- solve(invPsi0_q)/(nu0_q - r + u - 1)

      DA <- rbind(diag(1, r - u))

      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(mu)

      tr_term2 <- sum(diag((Y_1u %*% DA) %*% solve(Omega0) %*% t(Y_1u %*% DA)))

      log_weights[iter] <- -0.5 * n * r * log(2 * pi) -
        0.5 * n * logdetfunction(Omega0) -
        0.5 * tr_term2

    }else if(u == r){
      # full model
      # non A, Omega0

      # Draw from variational posterior
      U_q <- make_positive_definite(U_q)
      V_q <- make_positive_definite(V_q)

      # mu <- matrix(mvtnorm::rmvnorm(n = 1, mean = mu_q, sigma = Sigma_q), r, 1)
      # eta <- matrixNormal::rmatnorm(s = 1, M = eta_q, U = U_q, V = V_q)
      # Omega <- CholWishart::rInvWishart(n = 1, df = nu_q, Sigma = solve(invPsi_q))[,,1]

      mu <- mu_q
      eta <- eta_q
      Omega <- solve(invPsi_q)/(nu_q - u - 1)

      CA <- rbind(diag(1, u))

      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(mu)
      Xcteta <- Xc %*% t(eta)

      tr_term1 <- sum(diag((Y_1u %*% CA - Xcteta) %*% solve(Omega) %*% t(Y_1u %*% CA - Xcteta)))

      log_weights[iter] <- -0.5 * n * r * log(2 * pi) -
        0.5 * n * logdetfunction(Omega) -
        0.5 * tr_term1

    }else{
      # Draw from variational posterior
      U_q <- make_positive_definite(U_q)
      V_q <- make_positive_definite(V_q)

      # mu <- matrix(mvtnorm::rmvnorm(n = 1, mean = mu_q, sigma = Sigma_q), r, 1)
      # eta <- matrixNormal::rmatnorm(s = 1, M = eta_q, U = U_q, V = V_q)
      # Omega <- CholWishart::rInvWishart(n = 1, df = nu_q, Sigma = solve(invPsi_q))[,,1]
      # Omega0 <- CholWishart::rInvWishart(n = 1, df = nu0_q, Sigma = solve(invPsi0_q))[,,1]
      # A <- matrix(mvtnorm::rmvnorm(n = 1, mean = matrix(hatA, (r-u)*u, 1), sigma = HA), r - u, u)

      mu <- mu_q
      eta <- eta_q
      Omega <- solve(invPsi_q)/(nu_q - u - 1)
      Omega0 <- solve(invPsi0_q)/(nu0_q - r + u - 1)
      A <- hatA

      CA <- rbind(diag(1, u), A)
      DA <- rbind(-t(A), diag(1, r - u))

      log_det_Iu_AA <- logdetfunction(diag(u) + crossprod(A))
      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(mu)
      Xcteta <- Xc %*% t(eta)

      tr_term1 <- sum(diag((Y_1u %*% CA - Xcteta) %*% solve(Omega) %*% t(Y_1u %*% CA - Xcteta)))
      tr_term2 <- sum(diag((Y_1u %*% DA) %*% solve(Omega0) %*% t(Y_1u %*% DA)))

      log_weights[iter] <- -0.5 * n * r * log(2 * pi) -
        0.5 * n * logdetfunction(Omega) -
        0.5 * n * logdetfunction(Omega0) +
        n * log_det_Iu_AA -
        0.5 * tr_term1 -
        0.5 * tr_term2
    }

  }

  return(log_weights)
}

