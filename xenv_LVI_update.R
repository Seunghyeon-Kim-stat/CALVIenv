
library(progress)
library(matrixcalc)
library(dplyr)

xenv_A_LVI <- function(maxiter = 100, method = "BFGS",
                       nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                       hatA, CA, DA, HA,
                       nuY_q, invPsiY_q,
                       nuX1_q, invPsiX1_q,
                       nuX0_q, invPsiX0_q,
                       eta_q, Sigmaeta_q,
                       muX_q, SigmaX_q,
                       muY_q, SigmaY_q, # variational distribution
                       X, Y, X_bar, Y_bar, m, n, p, r, # data
                       GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
                       ){

  # common part
  # Gp

  if (m == 0){

    hatA <- matrix(0, (p-m), m)
    CA <- rbind(diag(m), hatA)
    DA <- rbind(-t(hatA), diag(p-m))
    HA <- matrix(0, (p-m)*m, (p-m)*m)

  }

  else if (m < p){

    invU_A <- solve(U_A)
    invV_A <- solve(V_A)

    ftol <- 1e-3

    C1 <- nuY_q * eta_q %*% invPsiY_q %*% t(eta_q) + divide_matrix(Sigmaeta_q, nuY_q * invPsiY_q, commutation = TRUE)
    C2 <- GpX
    C3 <- nuX1_q * invPsiX1_q
    C4 <- GpX + 1/psi_eta * BMB + psiX1 * diag(p)
    C5 <- nuX0_q * invPsiX0_q
    C6 <- GpX + psiX0 * diag(p)
    C7 <- nuY_q * eta_q %*% invPsiY_q %*% tYmuXmu + 1/psi_eta * nuX1_q * invPsiX1_q %*% eta_q %*% M %*% t(B0)

    o_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
    o_term2 <- -1/2 * tr( C3 %*% t(CA) %*% C4 %*% CA )
    o_term3 <- -1/2 * tr( C5 %*% t(DA) %*% C6 %*% DA )
    o_term4 <- tr( C7 %*% CA)
    o_term5 <- -1/2 * tr( invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))
    logdet <- ((2*n + nuX1 + nuX0) / 2) * logdetfunction(diag(p-m) + tcrossprod(hatA))

    obj1 <- o_term1 + o_term2 + o_term3 + o_term4 + o_term5 + logdet
    obj1 <- obj1 / n

    i <- 1

    while (i < maxiter) {
      for (k in (m + 1):p) {
        j <- k

        fA <- function(x){

          hatA[j-m,] <- x
          CA <- rbind(diag(m), hatA)
          DA <- rbind(-t(hatA), diag(p-m))

          f_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
          f_term2 <- -1/2 * tr( C3 %*% t(CA) %*% C4 %*% CA )
          f_term3 <- -1/2 * tr( C5 %*% t(DA) %*% C6 %*% DA )
          f_term4 <- tr( C7 %*% CA)
          f_term5 <- -1/2 * tr( invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))
          f_logdet <- ((2*n + nuX1 + nuX0) / 2) * logdetfunction(diag(p-m) + tcrossprod(hatA))

          fA_value <- f_term1 + f_term2 + f_term3 + f_term4 + f_term5 + f_logdet

          return(-fA_value)
        }

        gA <- function(x){

          hatA[j-m,] <- x
          CA <- rbind(diag(m), hatA)
          DA <- rbind(-t(hatA), diag(p-m))

          d_term1 <- - t(L) %*% C2 %*% CA %*% C1
          d_term2 <- - t(L) %*% C4 %*% CA %*% C3
          d_term3 <- - C5 %*% t(DA) %*% C6 %*% K
          d_term4 <- t(L) %*% t(C7)
          d_term5 <- - invU_A %*% (hatA - A0) %*% invV_A

          d_logdet <- (2*n + nuX1 + nuX0) * solve((diag(p-m) + hatA %*% t(hatA)), hatA)

          # solve(A) %*% B -> A X = B
          # X = solve(A, B)

          gA_value <- (d_term1 + d_term2 + d_term3 + d_term4 + d_term5 + d_logdet)[j-m,]

          return(-gA_value)
        }

        res <- stats::optim(
          hatA[j-m, ],
          fn = fA,
          gr = gA,
          method = "BFGS"
        )
        hatA[j-m, ] <- res$par
      }

      CA <- rbind(diag(m), hatA)
      DA <- rbind(-t(hatA), diag(p-m))

      o_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
      o_term2 <- -1/2 * tr( C3 %*% t(CA) %*% C4 %*% CA )
      o_term3 <- -1/2 * tr( C5 %*% t(DA) %*% C6 %*% DA )
      o_term4 <- tr( C7 %*% CA)
      o_term5 <- -1/2 * tr( invV_A %*% t(hatA - A0) %*% invU_A %*% (hatA - A0))
      logdet <- ((2*n + nuX1 + nuX0) / 2) * logdetfunction(diag(m) + crossprod(hatA))

      obj2 <- o_term1 + o_term2 + o_term3 + o_term4 + o_term5 + logdet
      obj2 <- obj2 / n

      if (abs(obj1 - obj2) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj2
        i <- i + 1
      }

    }

    # Calculate Hessian

    invJ <- solve(diag(p-m) + tcrossprod(hatA))

    h_term1 <- - kronecker(C1, t(L) %*% C2 %*% L)
    h_term2 <- - kronecker(C3, t(L) %*% C4 %*% L)
    h_term3 <- - kronecker(t(K) %*% C6 %*% K, C5)
    h_term5 <- - kronecker(t(invV_A), invU_A)

    h_logdet <- kronecker((diag(m) - t(hatA) %*% invJ %*% hatA), invJ) - tcrossprod(vec(invJ %*% hatA))

    HA <- h_term1 + h_term2 + h_term3 + h_term5 + (2*n + nuX1 + nuX0)*h_logdet
    HA <- -solve(HA)
    HA <- (HA + t(HA)) / 2
    HA <- make_positive_definite(HA)

    CA <- rbind(diag(m), hatA)
    DA <- rbind(-t(hatA), diag(p-m))
  }

  else if (m == p){

    hatA <- matrix(0, (p-m), m)
    CA <- rbind(diag(m), hatA)
    DA <- rbind(-t(hatA), diag(p-m))
    HA <- matrix(0, (p-m)*m, (p-m)*m)

  }

  return(list(hatA, HA, CA, DA))
}

xenv_SigmaYcX_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                              hatA, CA, DA, HA,
                              nuY_q, invPsiY_q,
                              nuX1_q, invPsiX1_q,
                              nuX0_q, invPsiX0_q,
                              eta_q, Sigmaeta_q,
                              muX_q, SigmaX_q,
                              muY_q, SigmaY_q, # variational distribution
                              X, Y, X_bar, Y_bar, m, n, p, r, # data
                              GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
                              ){

  onevector <- matrix(1, n, 1)
  tCAtXmuXmuCA <- t(CA) %*% GpX %*% CA

  if(m ==0){

    S1 <- tCAtXmuXmuCA + divide_matrix(HA, t(L) %*% GpX %*% L)

    nuY_q <- n + nuY

    PsiY_q <- GpY + psiY * diag(r) - 2 * tYmuXmu %*% CA %*% eta_q + t(eta_q) %*% S1 %*% eta_q + divide_matrix(Sigmaeta_q, S1)
    invPsiY_q <- solve(PsiY_q)
    invPsiY_q <- (t(invPsiY_q) + invPsiY_q) / 2

  }
  else if(m < p){

    S1 <- tCAtXmuXmuCA + divide_matrix(HA, t(L) %*% GpX %*% L)

    nuY_q <- n + nuY

    PsiY_q <- GpY + psiY * diag(r) - 2 * tYmuXmu %*% CA %*% eta_q + t(eta_q) %*% S1 %*% eta_q + divide_matrix(Sigmaeta_q, S1)
    invPsiY_q <- solve(PsiY_q)
    invPsiY_q <- (t(invPsiY_q) + invPsiY_q) / 2

  }
  else if(m == p){

    S1 <- tCAtXmuXmuCA + divide_matrix(HA, t(L) %*% GpX %*% L)

    nuY_q <- n + nuY

    PsiY_q <- GpY + psiY * diag(r) - 2 * tYmuXmu %*% CA %*% eta_q + t(eta_q) %*% S1 %*% eta_q + divide_matrix(Sigmaeta_q, S1)
    invPsiY_q <- solve(PsiY_q)
    invPsiY_q <- (t(invPsiY_q) + invPsiY_q) / 2

  }


  return(list(nuY_q, invPsiY_q))
}

xenv_muX_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                         hatA, CA, DA, HA,
                         nuY_q, invPsiY_q,
                         nuX1_q, invPsiX1_q,
                         nuX0_q, invPsiX0_q,
                         eta_q, Sigmaeta_q,
                         muX_q, SigmaX_q,
                         muY_q, SigmaY_q, # variational distribution
                         X, Y, X_bar, Y_bar, m, n, p, r, # data
                         GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  muX_q <- X_bar

  S2 <- nuX1_q * invPsiX1_q
  S3 <- nuX0_q * invPsiX0_q
  S4 <- nuY_q * invPsiY_q
  S5 <- eta_q %*% S4 %*% t(eta_q) + divide_matrix(Sigmaeta_q, S4, commutation = TRUE)

  if (m ==0){

    invSigmaX_q <-
      n * (
        CA %*% (S2+S5) %*% t(CA) + #L %*% divide_matrix(HA, S2+S5, commutation = TRUE) %*% t(L) +
          DA %*% S3 %*% t(DA) #+ K %*% divide_matrix(HA, S3) %*% t(K)
      )

  }
  else if(m < p){

    invSigmaX_q <-
      n * (
        CA %*% (S2+S5) %*% t(CA) + L %*% divide_matrix(HA, S2+S5, commutation = TRUE) %*% t(L) +
          DA %*% S3 %*% t(DA) + K %*% divide_matrix(HA, S3) %*% t(K)
      )

  }
  else if(m == p){

    invSigmaX_q <-
      n * (
        CA %*% (S2+S5) %*% t(CA) + #L %*% divide_matrix(HA, S2+S5, commutation = TRUE) %*% t(L) +
          DA %*% S3 %*% t(DA)# + K %*% divide_matrix(HA, S3) %*% t(K)
      )

  }

  SigmaX_q <- solve(invSigmaX_q)
  SigmaX_q <- (t(SigmaX_q) + SigmaX_q) / 2

  return(list(muX_q, SigmaX_q))
}

xenv_muY_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                         hatA, CA, DA, HA,
                         nuY_q, invPsiY_q,
                         nuX1_q, invPsiX1_q,
                         nuX0_q, invPsiX0_q,
                         eta_q, Sigmaeta_q,
                         muX_q, SigmaX_q,
                         muY_q, SigmaY_q, # variational distribution
                         X, Y, X_bar, Y_bar, m, n, p, r, # data
                         GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  muY_q <- Y_bar

  SigmaY_q <- solve(invPsiY_q) / (n * nuY_q)
  SigmaY_q <- (t(SigmaY_q) + SigmaY_q) / 2

  return(list(muY_q, SigmaY_q))
}


xenv_eta_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                         hatA, CA, DA, HA,
                         nuY_q, invPsiY_q,
                         nuX1_q, invPsiX1_q,
                         nuX0_q, invPsiX0_q,
                         eta_q, Sigmaeta_q,
                         muX_q, SigmaX_q,
                         muY_q, SigmaY_q, # variational distribution
                         X, Y, X_bar, Y_bar, m, n, p, r, # data
                         GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  if (m == 0){

    eta_q <- matrix(0, m, r)
    Sigmaeta_q <- matrix(0, m*r, m*r)

  }
  else if (m < p){

    tCAtXmuXmuCA <- t(CA) %*% GpX %*% CA

    S1 <- tCAtXmuXmuCA + divide_matrix(HA, t(L) %*% GpX %*% L)

    invSigmaeta_q <- kronecker(nuY_q * invPsiY_q, S1) + kronecker(M, nuX1_q/psi_eta * invPsiX1_q)
    invSigmaeta_q <- (t(invSigmaeta_q) + invSigmaeta_q) / 2

    if(!is.positive.definite(invSigmaeta_q)){
      invSigmaeta_q <- make_positive_definite(invSigmaeta_q)
    }

    Sigmaeta_q <- chol2inv(chol(invSigmaeta_q))
    Sigmaeta_q <- (t(Sigmaeta_q) + Sigmaeta_q) / 2

    Q <- nuY_q * invPsiY_q %*% tYmuXmu %*% CA + nuX1_q/psi_eta * M %*% t(B0) %*% CA %*% invPsiX1_q
    eta_q <- Sigmaeta_q %*% vec(t(Q))
    eta_q <- matrix(eta_q, m, r)

  }
  else if (m == p){

    tCAtXmuXmuCA <- t(CA) %*% GpX %*% CA

    S1 <- tCAtXmuXmuCA # + divide_matrix(HA, t(L) %*% GpX %*% L)

    invSigmaeta_q <- kronecker(nuY_q * invPsiY_q, S1) + kronecker(M, nuX1_q/psi_eta * invPsiX1_q)
    invSigmaeta_q <- (t(invSigmaeta_q) + invSigmaeta_q) / 2

    if(!is.positive.definite(invSigmaeta_q)){
      invSigmaeta_q <- make_positive_definite(invSigmaeta_q)
    }

    Sigmaeta_q <- chol2inv(chol(invSigmaeta_q))
    Sigmaeta_q <- (t(Sigmaeta_q) + Sigmaeta_q) / 2

    Q <- nuY_q * invPsiY_q %*% tYmuXmu %*% CA + nuX1_q/psi_eta * M %*% t(B0) %*% CA %*% invPsiX1_q
    eta_q <- Sigmaeta_q %*% vec(t(Q))
    eta_q <- matrix(eta_q, m, r)

  }

  return(list(eta_q, Sigmaeta_q))
}

xenv_Omega_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                           hatA, CA, DA, HA,
                           nuY_q, invPsiY_q,
                           nuX1_q, invPsiX1_q,
                           nuX0_q, invPsiX0_q,
                           eta_q, Sigmaeta_q,
                           muX_q, SigmaX_q,
                           muY_q, SigmaY_q, # variational distribution
                           X, Y, X_bar, Y_bar, m, n, p, r, # data
                           GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  nuX1_q <- n + nuX1 - r

  if (m == 0){

    invPsiX1_q <- matrix(0, m, m)

  }
  else if (m < p){

    S5 <- (GpX + 1/psi_eta * BMB + psiX1 * diag(p))

    PsiX1_q <-
      t(CA) %*% S5 %*% CA + divide_matrix(HA, t(L) %*% S5 %*% L) +
      eta_q %*% (1/psi_eta * M) %*% t(eta_q) + divide_matrix(Sigmaeta_q, 1/psi_eta * M, commutation = TRUE)

    invPsiX1_q <- solve(PsiX1_q)
    invPsiX1_q <- (t(invPsiX1_q) + invPsiX1_q) / 2

  }
  else if (m == p){

    S5 <- (GpX + 1/psi_eta * BMB + psiX1 * diag(p))

    PsiX1_q <-
      t(CA) %*% S5 %*% CA + divide_matrix(HA, t(L) %*% S5 %*% L) +
      eta_q %*% (1/psi_eta * M) %*% t(eta_q) + divide_matrix(Sigmaeta_q, 1/psi_eta * M, commutation = TRUE)

    invPsiX1_q <- solve(PsiX1_q)
    invPsiX1_q <- (t(invPsiX1_q) + invPsiX1_q) / 2

  }

  return(list(nuX1_q, invPsiX1_q))
}

xenv_Omega0_LVI <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                            hatA, CA, DA, HA,
                            nuY_q, invPsiY_q,
                            nuX1_q, invPsiX1_q,
                            nuX0_q, invPsiX0_q,
                            eta_q, Sigmaeta_q,
                            muX_q, SigmaX_q,
                            muY_q, SigmaY_q, # variational distribution
                            X, Y, X_bar, Y_bar, m, n, p, r, # data
                            GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  nuX0_q <- n + nuX0

  if (m == 0){

    S5 <- (GpX + psiX0 * diag(p))

    PsiX0_q <-
      t(DA) %*% S5 %*% DA + divide_matrix(HA, t(K) %*% S5 %*% K, commutation = TRUE)

    invPsiX0_q <- solve(PsiX0_q)
    invPsiX0_q <- (t(invPsiX0_q) + invPsiX0_q) / 2

  }
  else if (m < p){

    S5 <- (GpX + psiX0 * diag(p))

    PsiX0_q <-
      t(DA) %*% S5 %*% DA + divide_matrix(HA, t(K) %*% S5 %*% K, commutation = TRUE)

    invPsiX0_q <- solve(PsiX0_q)
    invPsiX0_q <- (t(invPsiX0_q) + invPsiX0_q) / 2

  }
  else if(m == p){

    invPsiX0_q <- matrix(0, (p-m), (p-m))

  }

  return(list(nuX0_q, invPsiX0_q))
}

xenv_convergence <- function(nuY, psiY, nuX1, psiX1, nuX0, psiX0, A0, U_A, V_A, psi_eta, B0, M,# prior
                             hatA, CA, DA, HA,
                             nuY_q, invPsiY_q,
                             nuX1_q, invPsiX1_q,
                             nuX0_q, invPsiX0_q,
                             eta_q, Sigmaeta_q,
                             muX_q, SigmaX_q,
                             muY_q, SigmaY_q, # variational distribution
                             X, Y, X_bar, Y_bar, m, n, p, r, # data
                             GpX, GpY, BMB, tYmuXmu, K, L,... # common terms
){

  onevector <- matrix(1, n, 1)

  const <-
    matrix_normal_logconst(U_A, V_A) +
    invWishart_logconst(nuY, psiY, r) +
    invWishart_logconst(nuX1, psiX1, m) +
    invWishart_logconst(nuX0, psiX0, (p-m)) -
    0.5 * (m*r + n*r + n*p) * log(2 * pi) +
    0.5 * m * logdetfunction(M) -
    0.5 * m * r * log(psi_eta)

  extra_part <-
    - 0.5 * (nuY_q + r + 1) * expectation_loginvWishart(nuY_q, invPsiY_q) -
    0.5 * (nuX1_q + m + 1) * expectation_loginvWishart(nuX1_q, invPsiX1_q) -
    0.5 * (nuX0_q + p - m + 1) * expectation_loginvWishart(nuX0_q, invPsiX0_q) -
    0.5 * tr( nuY_q * (GpY + psiY * diag(r)) %*% invPsiY_q ) -
    0.5 * (nuX1_q / psi_eta) * tr( ( t(eta_q) %*% invPsiX1_q %*% eta_q + divide_matrix(Sigmaeta_q, invPsiX1_q) ) %*% M)

  C1 <- nuY_q * eta_q %*% invPsiY_q %*% t(eta_q) + divide_matrix(Sigmaeta_q, nuY_q * invPsiY_q, commutation = TRUE)
  C2 <- GpX
  C3 <- nuX1_q * invPsiX1_q
  C4 <- GpX + 1/psi_eta * BMB + psiX1 * diag(p)
  C5 <- nuX0_q * invPsiX0_q
  C6 <- GpX + psiX0 * diag(p)
  C7 <- nuY_q * eta_q %*% invPsiY_q %*% tYmuXmu + 1/psi_eta * nuX1_q * invPsiX1_q %*% eta_q %*% M %*% t(B0)

  f_term1 <- -1/2 * tr( C1 %*% t(CA) %*% C2 %*% CA )
  f_term2 <- -1/2 * tr( C3 %*% t(CA) %*% C4 %*% CA )
  f_term3 <- -1/2 * tr( C5 %*% t(DA) %*% C6 %*% DA )
  f_term4 <- tr( C7 %*% CA)

  if (m == 0){
    f_term5 <- 0
  }
  else if(m < p){
    f_term5 <- -1/2 * tr( solve(V_A) %*% t(hatA - A0) %*% solve(U_A) %*% (hatA - A0))
  }
  else if(m == p){
    f_term5 <- 0
  }

  f_logdet <- ((2*n + nuX1 + nuX0) / 2) * logdetfunction(diag(p-m) + tcrossprod(hatA))

  fA <- f_term1 + f_term2 + f_term3 + f_term4 + f_term5 + f_logdet

  Entropy_muX_q <- normal_entropy(muX_q, SigmaX_q)

  Entropy_muY_q <- normal_entropy(muY_q, SigmaY_q)

  Entropy_eta_q <- normal_entropy(eta_q, Sigmaeta_q)

  Entropy_SigmaYcX <- invWishart_entropy(nuY_q, invPsiY_q)

  Entropy_Omega_q <- invWishart_entropy(nuX1_q, invPsiX1_q)

  Entropy_Omega0_q <- invWishart_entropy(nuX0_q, invPsiX0_q)

  Entropy_A_q <- normal_entropy(hatA, HA)

  Lvalue <- const + extra_part + fA - 1/2 * (p-m) * m +
    Entropy_muX_q + Entropy_muY_q + Entropy_eta_q + Entropy_SigmaYcX + Entropy_Omega_q + Entropy_Omega0_q + Entropy_A_q

  return(Lvalue)

}

xenv_llik <- function(X, Y, N,
                      nuY_q, invPsiY_q,
                      nuX1_q, invPsiX1_q,
                      nuX0_q, invPsiX0_q,
                      eta_q, Sigmaeta_q,
                      muX_q, SigmaX_q,
                      muY_q, SigmaY_q,
                      hatA, HA # variational distribution
                      ) {

  X <- data.matrix(X)
  Y <- data.matrix(Y)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  m <- ncol(hatA)

  log_weights <- numeric(N)

  for (iter in 1:N) {

    # Draw from variational posterior
    # muX <- matrix(mvtnorm::rmvnorm(n = 1, mean = muX_q, sigma = SigmaX_q), p, 1)
    # muY <- matrix(mvtnorm::rmvnorm(n = 1, mean = muY_q, sigma = SigmaY_q), r, 1)
    # eta <- matrix(mvtnorm::rmvnorm(n = 1, mean = eta_q, sigma = Sigmaeta_q), m, r)
    # OmegaX1 <- CholWishart::rInvWishart(n = 1, df = nuX1_q, Sigma = solve(invPsiX1_q))[,,1]
    # OmegaX0 <- CholWishart::rInvWishart(n = 1, df = nuX0_q, Sigma = solve(invPsiX0_q))[,,1]
    # SigmaY <- CholWishart::rInvWishart(n = 1, df = nuY_q, Sigma = solve(invPsiY_q))[,,1]
    # A <- matrix(mvtnorm::rmvnorm(n = 1, mean = matrix(hatA, (p-m)*m, 1), sigma = HA), p - m, m)

    if (m == 0){
      # reduced model
      # non eta, A, OmegaX1, X

      muX <- muX_q
      muY <- muY_q
      OmegaX0 <- solve(invPsiX0_q)/(nuX0_q - p + m - 1)
      SigmaY <- solve(invPsiY_q)/(nuY_q - r - 1)
      DA <- rbind(diag(1, p - m))

      X_1u <- X - matrix(rep(1,n), ncol = 1) %*% t(muX)
      XutXu <- crossprod(X_1u)
      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(muY)

      tr_term1 <- sum(diag((Y_1u) %*% solve(SigmaY) %*% t(Y_1u)))
      tr_term3 <- sum(diag(t(DA) %*% XutXu %*% DA %*% solve(OmegaX0)))

      log_weights[iter] <- -0.5 * (n * r + n * p) * log(2 * pi) -
        0.5 * n * logdetfunction(SigmaY) -
        0.5 * n * logdetfunction(OmegaX0) -
        0.5 * tr_term3

    }else if(m == p){
      # full model
      # non A, Omega0

      muX <- muX_q
      muY <- muY_q
      eta <- matrix(eta_q, m, r)
      OmegaX1 <- solve(invPsiX1_q)/(nuX1_q - m - 1)
      SigmaY <- solve(invPsiY_q)/(nuY_q - r - 1)
      CA <- rbind(diag(1, m))

      X_1u <- X - matrix(rep(1,n), ncol = 1) %*% t(muX)
      XutXu <- crossprod(X_1u)
      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(muY)

      tr_term1 <- sum(diag((Y_1u - X_1u %*% CA %*% eta) %*% solve(SigmaY) %*% t(Y_1u - X_1u %*% CA %*% eta)))
      tr_term2 <- sum(diag(t(CA) %*% XutXu %*% CA %*% solve(OmegaX1)))

      log_weights[iter] <- -0.5 * (n * r + n * p) * log(2 * pi) -
        0.5 * n * logdetfunction(SigmaY) -
        0.5 * n * logdetfunction(OmegaX1) -
        0.5 * tr_term1 -
        0.5 * tr_term2
    }else{

      muX <- muX_q
      muY <- muY_q
      eta <- matrix(eta_q, m, r)
      OmegaX1 <- solve(invPsiX1_q)/(nuX1_q - m - 1)
      OmegaX0 <- solve(invPsiX0_q)/(nuX0_q - p + m - 1)
      SigmaY <- solve(invPsiY_q)/(nuY_q - r - 1)
      A <- hatA

      CA <- rbind(diag(1, m), A)
      DA <- rbind(-t(A), diag(1, p - m))

      log_det_Iu_AA <- logdetfunction(diag(m) + crossprod(A))
      X_1u <- X - matrix(rep(1,n), ncol = 1) %*% t(muX)
      XutXu <- crossprod(X_1u)
      Y_1u <- Y - matrix(rep(1, n), ncol = 1) %*% t(muY)

      tr_term1 <- sum(diag((Y_1u - X_1u %*% CA %*% eta) %*% solve(SigmaY) %*% t(Y_1u - X_1u %*% CA %*% eta)))
      tr_term2 <- sum(diag(t(CA) %*% XutXu %*% CA %*% solve(OmegaX1)))
      tr_term3 <- sum(diag(t(DA) %*% XutXu %*% DA %*% solve(OmegaX0)))

      log_weights[iter] <- -0.5 * (n * r + n * p) * log(2 * pi) -
        0.5 * n * logdetfunction(SigmaY) -
        0.5 * n * logdetfunction(OmegaX1) -
        0.5 * n * logdetfunction(OmegaX0) +
        n * log_det_Iu_AA -
        0.5 * tr_term1 -
        0.5 * tr_term2 -
        0.5 * tr_term3
    }

  }

  return(log_weights)
}

