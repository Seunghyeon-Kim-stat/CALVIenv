
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

