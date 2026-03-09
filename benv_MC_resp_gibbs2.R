
library(dplyr)
library(rstan)
library(magrittr)
library(tcltk)
library(Matrix)
library(BAMBI)

GE <- fuMatrixGE <- function(A) {
  # Gaussian elimination, p must be less than or equal to n
  a <- dim(A)
  n <- a[1]
  p <- a[2]
  idx <- rep(0, p)
  res.idx <- 1:n

  i <- 1
  while (i <= p) {
    tmp <- max(abs(A[res.idx, i]))
    Stmp <- setdiff(which(abs(A[, i]) == tmp), idx)
    idx[i] <- Stmp[1]
    res.idx <- setdiff(res.idx, idx[i])
    for (j in 1:(n - i)) {
      A[res.idx[j], ] <- A[res.idx[j], ] - A[res.idx[j], i] / A[idx[i], i] * A[idx[i], ]
    }
    i <- i + 1
  }
  c(idx, res.idx)
}

env <- function(X, Y, u, asy = TRUE, init = NULL) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)

  if (a[1] != nrow(X)) stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) stop("u must be an integer between 0 and r.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")

  sigY <- stats::cov(Y) * (n - 1) / n
  sigYX <- stats::cov(Y, X) * (n - 1) / n
  sigX <- stats::cov(X) * (n - 1) / n
  invsigX <- chol2inv(chol(sigX))
  betaOLS <- sigYX %*% invsigX

  U <- tcrossprod(betaOLS, sigYX)
  M <- sigY - U

  if (qr(M)$rank < r){

    sigY <- stats::cov(Y) * (n - 1) / n + 1e-6 * diag(r)
    sigYX <- stats::cov(Y, X) * (n - 1) / n
    sigX <- stats::cov(X) * (n - 1) / n + 1e-6 * diag(p)
    invsigX <- chol2inv(chol(sigX))
    betaOLS <- sigYX %*% invsigX

    U <- tcrossprod(betaOLS, sigYX)
    M <- sigY - U

  }

  tmp <- envMU(M, U, u)

  if (!is.null(init)) {
    if (nrow(init) != r || ncol(init) != u) stop("The initial value should have r rows and u columns.")
    tmp0 <- qr.Q(qr(init), complete = TRUE)
    tmp$Gammahat <- as.matrix(tmp0[, 1:u])
    tmp$Gamma0hat <- as.matrix(tmp0[, (u + 1):r])
  }

  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  objfun <- tmp$objfun
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL

  if (u == 0) {
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    muhat <- colMeans(Y)
    betahat <- matrix(0, r, p)
    Sigmahat <- sigY
    loglik <- -n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) ratio <- matrix(1, r, p)
  } else if (u == r) {
    etahat <- betaOLS
    Omegahat <- M
    Omega0hat <- NULL
    muhat <- colMeans(Y) - betaOLS %*% colMeans(X)
    betahat <- betaOLS
    Sigmahat <- M
    loglik <- -n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- matrix(1, r, p)
    }
  } else {
    etahat <- crossprod(Gammahat, betaOLS)
    betahat <- Gammahat %*% etahat
    muhat <- colMeans(Y) - betahat %*% colMeans(X)
    Omegahat <- crossprod(Gammahat, M) %*% Gammahat
    Omega0hat <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
    Sigma1 <- Gammahat %*% tcrossprod(Omegahat, Gammahat)
    Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
    loglik <- -n * r / 2 * (log(2 * pi) + 1) - n / 2 * objfun
    if (asy == T) {
      covMatrix <- kronecker(invsigX, M)
      asyFm <- matrix(sqrt(diag(covMatrix)), nrow = r)
      invOmega0hat <- chol2inv(chol(Omega0hat))
      invOmegahat <- chol2inv(chol(Omegahat))
      temp <- kronecker(etahat %*% tcrossprod(sigX, etahat), invOmega0hat) + kronecker(invOmegahat, Omega0hat) + kronecker(Omegahat, invOmega0hat) - 2 * kronecker(diag(u), diag(r - u))
      temp2 <- kronecker(t(etahat), Gamma0hat)
      covMatrix <- kronecker(invsigX, Sigma1) + temp2 %*% chol2inv(chol(temp)) %*% t(temp2)
      asySE <- matrix(sqrt(diag(covMatrix)), nrow = r)
      ratio <- asyFm / asySE
    }
  }

  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, mu = muhat, beta = betahat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, loglik = loglik, n = n, covMatrix = covMatrix, asySE = asySE, ratio = ratio, CA = tmp$CA, idx = tmp$idx))
}

envMU <- function(M, U, u) {
  dimM <- dim(M)
  dimU <- dim(U)
  r <- dimM[1]

  if (dimM[1] != dimM[2] & dimU[1] != dimU[2]) stop("M and U should be square matrices.")
  if (dimM[1] != dimU[1]) stop("M and U should have the same dimension.")
  if (qr(M)$rank < r) stop("M should be positive definite.")
  if (u > r & u < 0) stop("u should be between 0 and r.")


  if (u == 0) {
    Gammahat <- NULL
    Gamma0hat <- diag(r)
    MU <- M + U
    tmp.MU <- eigen(MU)
    objfun <- sum(log(tmp.MU$values))
    GEidx <- 1:r
    Ginit <- NULL
  } else if (u == r) {
    Gammahat <- diag(r)
    Gamma0hat <- NULL
    tmp.M <- eigen(M)
    objfun <- sum(log(tmp.M$values))
    GEidx <- 1:r
    Ginit <- diag(r)
  } else if (u == 1) {
    maxiter <- 100
    ftol <- 1e-3


    MU <- M + U
    tmp.MU <- eigen(MU)
    invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, "*") %*% t(tmp.MU$vectors)
    invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), "*") %*% t(tmp.MU$vectors)

    midmatrix <- U
    startv <- function(a) t(a) %*% midmatrix %*% a
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])

    # 		if (qr(MU)$rank == r) {
    eig1 <- eigen(t(init) %*% M %*% init)
    eig2 <- eigen(t(init) %*% invMU %*% init)
    obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))

    midmatrix <- invMU2 %*% tcrossprod(U, invMU2)
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1]])
    e1 <- eigen(t(init.MU) %*% M %*% init.MU)
    e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
    obj2 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj2 < obj1) {
      init <- init.MU
      obj1 <- obj2
    }

    # 			if (qr(M)$rank == r) {
    tmp.M <- eigen(M)
    midmatrix <- U
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])
    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj3 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj3 < obj1) {
      init <- init.M
      obj1 <- obj3
    }

    invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), "*") %*% t(tmp.M$vectors)
    midmatrix <- invM2 %*% tcrossprod(U, invM2)
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1]])

    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj4 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj4 < obj1) {
      init <- init.M
      obj1 <- obj4
    }
    # 			}
    # 		}

    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1], ])

    i <- 1
    while (i < maxiter) {
      fobj <- function(x) {
        T1 <- crossprod(x, x)
        T2 <- crossprod(x, M) %*% x
        T3 <- crossprod(x, invMU) %*% x
        -2 * log(T1) + log(T2) + log(T3)
      }

      gobj <- function(x) {
        W1 <- crossprod(x, x)
        T1 <- x / as.vector(W1)
        W2 <- crossprod(x, M) %*% x
        T2 <- M %*% x / as.vector(W2)
        W3 <- crossprod(x, invMU) %*% x
        T3 <- invMU %*% x / as.vector(W3)
        -2 * T1 + T2 + T3
      }

      res <- stats::optim(Ginit, fobj, gobj, method = "BFGS")
      g <- as.matrix(res$par)
      a <- qr(g)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u + 1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
  } else if (u == r - 1) {
    maxiter <- 100
    ftol <- 1e-3

    MU <- M + U
    tmp.MU <- eigen(MU)
    invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, "*") %*% t(tmp.MU$vectors)
    invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), "*") %*% t(tmp.MU$vectors)

    midmatrix <- U
    startv <- function(a) t(a) %*% midmatrix %*% a
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])

    # 	  if (qr(MU)$rank == r) {
    eig1 <- eigen(t(init) %*% M %*% init)
    eig2 <- eigen(t(init) %*% invMU %*% init)
    obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))

    midmatrix <- invMU2 %*% tcrossprod(U, invMU2)
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
    e1 <- eigen(t(init.MU) %*% M %*% init.MU)
    e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
    obj2 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj2 < obj1) {
      init <- init.MU
      obj1 <- obj2
    }

    # 	    if (qr(M)$rank == r) {
    tmp.M <- eigen(M)
    midmatrix <- U
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj3 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj3 < obj1) {
      init <- init.M
      obj1 <- obj3
    }

    invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), "*") %*% t(tmp.M$vectors)
    midmatrix <- invM2 %*% tcrossprod(U, invM2)
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])

    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj4 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj4 < obj1) {
      init <- init.M
      obj1 <- obj4
    }
    # 	    }
    # 	  }

    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])

    j <- GEidx[r]

    g <- as.matrix(Ginit[j, ])
    t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
    t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]

    GUGt2 <- g + t2
    GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]

    GVGt2 <- g + t3
    GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j]

    invC1 <- chol2inv(chol(GUG))
    invC2 <- chol2inv(chol(GVG))

    fobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2
      T3 <- invC2 %*% tmp3
      -2 * log(1 + sum(x^2)) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
    }

    gobj <- function(x) {
      tmp2 <- x + t2
      tmp3 <- x + t3
      T2 <- invC1 %*% tmp2
      T3 <- invC2 %*% tmp3
      -4 * x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))
    }

    i <- 1
    while (i < maxiter) {
      res <- stats::optim(Ginit[j, ], fobj, gobj, method = "BFGS")
      Ginit[j, ] <- res$par

      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }

    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u + 1):r, drop = FALSE]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
  } else {
    maxiter <- 100
    ftol <- 1e-3

    MU <- M + U
    tmp.MU <- eigen(MU)
    invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, "*") %*% t(tmp.MU$vectors)
    invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), "*") %*% t(tmp.MU$vectors)

    midmatrix <- U
    startv <- function(a) t(a) %*% midmatrix %*% a
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])

    # 		if (qr(MU)$rank == r) {
    eig1 <- eigen(t(init) %*% M %*% init)
    eig2 <- eigen(t(init) %*% invMU %*% init)
    obj1 <- sum(log(eig1$values)) + sum(log(eig2$values))

    midmatrix <- invMU2 %*% tcrossprod(U, invMU2)
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
    e1 <- eigen(t(init.MU) %*% M %*% init.MU)
    e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
    obj2 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj2 < obj1) {
      init <- init.MU
      obj1 <- obj2
    }

    # 			if (qr(M)$rank == r) {
    tmp.M <- eigen(M)
    midmatrix <- U
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj3 <- sum(log(e1$values)) + sum(log(e2$values))

    if (obj3 < obj1) {
      init <- init.M
      obj1 <- obj3
    }

    invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), "*") %*% t(tmp.M$vectors)
    midmatrix <- invM2 %*% tcrossprod(U, invM2)
    tmp2.M <- apply(tmp.M$vectors, 2, startv)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])

    e1 <- eigen(t(init.M) %*% M %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj4 <- sum(log(e1$values)) + sum(log(e2$values))
    if (obj4 < obj1) {
      init <- init.M
      obj1 <- obj4
    }
    # 			}
    # 		}

    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])


    GUG <- crossprod(Ginit, (M %*% Ginit))
    GVG <- crossprod(Ginit, (invMU %*% Ginit))


    t4 <- crossprod(Ginit[GEidx[(u + 1):r], ], Ginit[GEidx[(u + 1):r], ]) + diag(u)
    i <- 1
    while (i < maxiter) {
      for (j in GEidx[(u + 1):r]) {
        g <- as.matrix(Ginit[j, ])
        t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
        t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]

        GUGt2 <- g + t2
        GUG <- GUG - tcrossprod(GUGt2, GUGt2) * M[j, j]

        GVGt2 <- g + t3
        GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invMU[j, j]

        t4 <- t4 - tcrossprod(g, g)
        invC1 <- chol2inv(chol(GUG))
        invC2 <- chol2inv(chol(GVG))
        invt4 <- chol2inv(chol(t4))

        fobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2
          T3 <- invC2 %*% tmp3
          -2 * log(1 + x %*% T1) + log(1 + M[j, j] * crossprod(tmp2, T2)) + log(1 + invMU[j, j] * crossprod(tmp3, T3))
        }

        gobj <- function(x) {
          tmp2 <- x + t2
          tmp3 <- x + t3
          T1 <- invt4 %*% x
          T2 <- invC1 %*% tmp2
          T3 <- invC2 %*% tmp3
          -4 * T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / M[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3))
        }

        res <- stats::optim(Ginit[j, ], fobj, gobj, method = "BFGS")
        Ginit[j, ] <- res$par
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        GUGt2 <- g + t2
        GUG <- GUG + tcrossprod(GUGt2, GUGt2) * M[j, j]

        GVGt2 <- g + t3
        GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invMU[j, j]
      }
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      e1 <- eigen(t(Gammahat) %*% M %*% Gammahat)
      e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)
      obj5 <- sum(log(e1$values)) + sum(log(e2$values))
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }

    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u + 1):r]
    objfun <- obj5 + sum(log(tmp.MU$values))
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
  }
  return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat, objfun = objfun, idx = GEidx, CA = Ginit))
}

MAPCA <- function(X, Y, u, A0, K=NULL, L=NULL, HA0=NULL, Psi, Psi0, e, M, nu, nu0,
                  maxiter = 100, method = "Nelder-Mead", ...) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  e <- as.matrix(e)
  Xc <- scale(X, scale = FALSE)
  Yc <- scale(Y, scale = FALSE)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)

  invK <- chol2inv(chol(K))
  invL <- chol2inv(chol(L))

  sumY <- crossprod(Yc)
  tmp <- crossprod(Xc) + M
  tmp2 <- crossprod(Yc, Xc) + e %*% M
  echeck <- tmp2 %*% chol2inv(chol(tmp))
  Gtilde <- sumY + e %*% tcrossprod(M, e) - tcrossprod(tmp2, echeck)
  fit <- env(X, Y, u)
  A <- fit$CA[fit$idx[(u + 1):r], ]
  GEidx <- fit$idx
  Ginit <- fit$CA

  if (u == 0) {
    A <- NULL
    Ginit <- NULL
    G0init <- diag(r)
    objfun <- NULL
  } else if (u == r) {
    A <- NULL
    Ginit <- diag(r)
    G0init <- NULL
    objfun <- NULL
  } else if (u == r - 1) {
    maxiter = 100
    ftol <- 1e-3
    A <- t(as.matrix(Ginit[GEidx[(u + 1):r], ]))
    G0init <- matrix(0, r, r - u)
    G0init[GEidx[1:u], ] <- -t(A)
    G0init[GEidx[(u + 1):r], ] <- diag(r - u)
    GG12 <- sqrtmat(crossprod(Ginit))
    G0G012 <- sqrtmat(crossprod(G0init))
    S1 <- (n + nu - 1) / 2 * log(det(crossprod(Ginit, Gtilde) %*% Ginit + crossprod(GG12, Psi) %*% GG12))
    S2 <- -(n + nu - 1) / 2 * log(det(crossprod(Ginit)))
    S3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(G0init)))
    S4 <- (n + nu0 - 1) / 2 * log(det(crossprod(G0init, sumY) %*% G0init + crossprod(G0G012, Psi0) %*% G0G012))
    S5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))
    obj1 <- S1 + S2 + S3 + S4 + S5

    i <- 1
    while (i < maxiter) {
      m <- r
      j <- GEidx[m]

      fobj <- function(x) {
        CA <- Ginit
        CA[j, ] <- x
        A <- CA[GEidx[(u + 1):r], ]
        DA <- matrix(0, r, r - u)
        DA[GEidx[1:u], ] <- -t(A)
        DA[GEidx[(u + 1):r], ] <- diag(r - u)
        CACA12 <- sqrtmat(crossprod(CA))
        DADA12 <- sqrtmat(crossprod(DA))
        T1 <- (n + nu - 1) / 2 * log(det(crossprod(CA, Gtilde) %*% CA + crossprod(CACA12, Psi) %*% CACA12))
        T2 <- -(n + nu - 1) / 2 * log(det(crossprod(CA)))
        T3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(DA)))
        T4 <- (n + nu0 - 1) / 2 * log(det(crossprod(DA, sumY) %*% DA + crossprod(DADA12, Psi0) %*% DADA12))
        T5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))
        T1 + T2 + T3 + T4 + T5
      }

      gobj <- function(x) {
        CA <- Ginit
        CA[j, ] <- x
        A <- CA[GEidx[(u + 1):r], ]
        DA <- matrix(0, r, r - u)
        DA[GEidx[1:u], ] <- -t(A)
        DA[GEidx[(u + 1):r], ] <- diag(r - u)
        CACAinv <- chol2inv(chol(crossprod(CA)))
        DADAinv <- chol2inv(chol(crossprod(DA)))
        CACA12 <- sqrtmat(crossprod(CA))
        DADA12 <- sqrtmat(crossprod(DA))
        tmp <- kronecker(diag(u), t(CA))
        tmp2 <- kronecker(diag(r - u), t(DA))
        tmp3 <- kronecker(diag(u), crossprod(CA, Gtilde))
        tmp4 <- kronecker(diag(r - u), crossprod(DA, sumY))
        idxuu <- c(matrix(c(1:(u^2)), nrow = u, byrow = T))
        idxruru <- c(matrix(c(1:((r - u)^2)), nrow = r - u, byrow = T))
        idxca <- j + r * (0:(u - 1))
        idxa <- (m - u) + (r - u) * (0:(u - 1))
        idxda <- (r * (m - u - 1) + 1):(r * (m - u - 1) + u)
        T1 <- -(n + nu - 1) / 2 * t(as.matrix(c(CACAinv))) %*% (tmp + tmp[idxuu, ])
        T1 <- T1[idxca]
        T2 <- -(n + nu - 1) / 2 * t(as.matrix(c(DADAinv))) %*% (tmp2 + tmp2[idxruru, ])
        T2 <- -T2[idxda]
        T3 <- (n + nu - 1) / 2 * t(as.matrix(c(chol2inv(chol(crossprod(CA, Gtilde) %*% CA + crossprod(CACA12, Psi) %*% CACA12))))) %*% (tmp3 + tmp3[idxuu, ] + (kronecker(diag(u), CACA12 %*% Psi) + kronecker(CACA12 %*% Psi, diag(u))) %*% chol2inv(chol(kronecker(diag(u), CACA12) + kronecker(CACA12, diag(u)))) %*% (tmp + tmp[idxuu, ]))
        T3 <- T3[idxca]
        T4 <- (n + nu0 - 1) / 2 * t(as.matrix(c(chol2inv(chol(crossprod(DA, sumY) %*% DA + crossprod(DADA12, Psi0) %*% DADA12))))) %*% (tmp4 + tmp4[idxruru, ] + (kronecker(diag(r - u), DADA12 %*% Psi0) + kronecker(DADA12 %*% Psi0, diag(r - u))) %*% chol2inv(chol(kronecker(diag(r - u), DADA12) + kronecker(DADA12, diag(r - u)))) %*% (tmp2 + tmp2[idxruru, ]))
        T4 <- T4[idxda]
        T5 <- t(c(A - A0)) %*% kronecker(invL, invK)
        T5 <- T5[idxa]
        T1 + T2 + T3 + T4 + T5
      }

      res <- stats::optim(
        Ginit[j, ],
        fn = fobj,
        gr = if (method == "BFGS") gobj else NULL,
        method = method
      )
      Ginit[j, ] <- res$par

      A <- Ginit[GEidx[(u + 1):r], ]
      G0init <- matrix(0, r, r - u)
      G0init[GEidx[1:u], ] <- -t(A)
      G0init[GEidx[(u + 1):r], ] <- diag(r - u)
      GG12 <- sqrtmat(crossprod(Ginit))
      G0G012 <- sqrtmat(crossprod(G0init))
      S1 <- (n + nu - 1) / 2 * log(det(crossprod(Ginit, Gtilde) %*% Ginit + crossprod(GG12, Psi) %*% GG12))
      S2 <- -(n + nu - 1) / 2 * log(det(crossprod(Ginit)))
      S3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(G0init)))
      S4 <- (n + nu0 - 1) / 2 * log(det(crossprod(G0init, sumY) %*% G0init + crossprod(G0G012, Psi0) %*% G0G012))
      S5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))
      obj2 <- S1 + S2 + S3 + S4 + S5

      if (abs(obj1 - obj2) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj2
        i <- i + 1
      }
      objfun <- obj1
    }
  } else {
    maxiter = 100
    ftol <- 1e-3
    A <- Ginit[GEidx[(u + 1):r], ]
    G0init <- matrix(0, r, r - u)
    G0init[GEidx[1:u], ] <- -t(A)
    G0init[GEidx[(u + 1):r], ] <- diag(r - u)
    GG12 <- sqrtmat(crossprod(Ginit))
    G0G012 <- sqrtmat(crossprod(G0init))
    S1 <- (n + nu - 1) / 2 * log(det(crossprod(Ginit, Gtilde) %*% Ginit + crossprod(GG12, Psi) %*% GG12))
    S2 <- -(n + nu - 1) / 2 * log(det(crossprod(Ginit)))
    S3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(G0init)))
    S4 <- (n + nu0 - 1) / 2 * log(det(crossprod(G0init, sumY) %*% G0init + crossprod(G0G012, Psi0) %*% G0G012))
    S5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))
    obj1 <- S1 + S2 + S3 + S4 + S5

    i <- 1
    while (i < maxiter) {

      for (m in (u + 1):r) {
        j <- GEidx[m]

        fobj <- function(x) {
          CA <- Ginit
          CA[j, ] <- x
          A <- CA[GEidx[(u + 1):r], ]
          DA <- matrix(0, r, r - u)
          DA[GEidx[1:u], ] <- -t(A)
          DA[GEidx[(u + 1):r], ] <- diag(r - u)
          CACA12 <- sqrtmat(crossprod(CA))
          DADA12 <- sqrtmat(crossprod(DA))
          T1 <- (n + nu - 1) / 2 * log(det(crossprod(CA, Gtilde) %*% CA + crossprod(CACA12, Psi) %*% CACA12))
          T2 <- -(n + nu - 1) / 2 * log(det(crossprod(CA)))
          T3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(DA)))
          T4 <- (n + nu0 - 1) / 2 * log(det(crossprod(DA, sumY) %*% DA + crossprod(DADA12, Psi0) %*% DADA12))
          T5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))
          T1 + T2 + T3 + T4 + T5
        }

        gobj <- function(x) {
          CA <- Ginit
          CA[j, ] <- x
          A <- CA[GEidx[(u + 1):r], ]
          DA <- matrix(0, r, r - u)
          DA[GEidx[1:u], ] <- -t(A)
          DA[GEidx[(u + 1):r], ] <- diag(r - u)
          CACAinv <- chol2inv(chol(crossprod(CA)))
          DADAinv <- chol2inv(chol(crossprod(DA)))
          CACA12 <- sqrtmat(crossprod(CA))
          DADA12 <- sqrtmat(crossprod(DA))
          tmp <- kronecker(diag(u), t(CA))
          tmp2 <- kronecker(diag(r - u), t(DA))
          tmp3 <- kronecker(diag(u), crossprod(CA, Gtilde))
          tmp4 <- kronecker(diag(r - u), crossprod(DA, sumY))
          idxuu <- c(matrix(c(1:(u^2)), nrow = u, byrow = T))
          idxruru <- c(matrix(c(1:((r - u)^2)), nrow = r - u, byrow = T))
          idxca <- j + r * (0:(u - 1))
          idxa <- (m - u) + (r - u) * (0:(u - 1))
          idxda <- (r * (m - u - 1) + 1):(r * (m - u - 1) + u)
          T1 <- -(n + nu - 1) / 2 * t(as.matrix(c(CACAinv))) %*% (tmp + tmp[idxuu, ])
          T1 <- T1[idxca]
          T2 <- -(n + nu - 1) / 2 * t(as.matrix(c(DADAinv))) %*% (tmp2 + tmp2[idxruru, ])
          T2 <- -T2[idxda]
          T3 <- (n + nu - 1) / 2 * t(as.matrix(c(chol2inv(chol(crossprod(CA, Gtilde) %*% CA + crossprod(CACA12, Psi) %*% CACA12))))) %*% (tmp3 + tmp3[idxuu, ] + (kronecker(diag(u), CACA12 %*% Psi) + kronecker(CACA12 %*% Psi, diag(u))) %*% chol2inv(chol(kronecker(diag(u), CACA12) + kronecker(CACA12, diag(u)))) %*% (tmp + tmp[idxuu, ]))
          T3 <- T3[idxca]
          T4 <- (n + nu0 - 1) / 2 * t(as.matrix(c(chol2inv(chol(crossprod(DA, sumY) %*% DA + crossprod(DADA12, Psi0) %*% DADA12))))) %*% (tmp4 + tmp4[idxruru, ] + (kronecker(diag(r - u), DADA12 %*% Psi0) + kronecker(DADA12 %*% Psi0, diag(r - u))) %*% chol2inv(chol(kronecker(diag(r - u), DADA12) + kronecker(DADA12, diag(r - u)))) %*% (tmp2 + tmp2[idxruru, ]))
          T4 <- T4[idxda]
          T5 <- t(c(A - A0)) %*% kronecker(invL, invK)
          T5 <- T5[idxa]
          T1 + T2 + T3 + T4 + T5
        }

        res <- stats::optim(
          Ginit[j, ],
          fn = fobj,
          gr = if (method == "BFGS") gobj else NULL,
          method = method
        )
        Ginit[j, ] <- res$par
      }
      A <- Ginit[GEidx[(u + 1):r], ]
      G0init <- matrix(0, r, r - u)
      G0init[GEidx[1:u], ] <- -t(A)
      G0init[GEidx[(u + 1):r], ] <- diag(r - u)
      GG12 <- sqrtmat(crossprod(Ginit))
      G0G012 <- sqrtmat(crossprod(G0init))
      S1 <- (n + nu - 1) / 2 * log(det(crossprod(Ginit, Gtilde) %*% Ginit + crossprod(GG12, Psi) %*% GG12))
      S2 <- -(n + nu - 1) / 2 * log(det(crossprod(Ginit)))
      S3 <- -(n + nu0 - 1) / 2 * log(det(crossprod(G0init)))
      S4 <- (n + nu0 - 1) / 2 * log(det(crossprod(G0init, sumY) %*% G0init + crossprod(G0G012, Psi0) %*% G0G012))

      S5 <- 0.5 * sum(diag(invK %*% (A - A0) %*% invL %*% t(A - A0)))

      obj2 <- S1 + S2 + S3 + S4 + S5

      if (abs(obj1 - obj2) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj2
        i <- i + 1
      }
      objfun <- obj1
    }
  }

  return(list(A = A, CA = Ginit, DA = G0init, G.tilde = Gtilde, Xc = Xc, Yc = Yc, e.check = echeck))
}

Benvlp_resp_map <- function(X, Y, u,
                            max.iter = 100,
                            method = "BFGS",
                            ...) {

  X <- data.matrix(X); Y <- data.matrix(Y)
  n <- nrow(X); p <- ncol(X); r <- ncol(Y)

  dots <- list(...); pp <- dots$prior_params

  if (is.null(pp)) pp <- list()

  # ---- default priors ----
  if (is.null(pp$e))    pp$e    <- matrix(0, r, p)
  if (is.null(pp$M))    pp$M    <- (1/1e6) * diag(rchisq(p, 1), nrow = p)
  if (is.null(pp$Psi))  pp$Psi  <- diag(1, u) / 1e6
  if (is.null(pp$nu))   pp$nu   <- u
  if (is.null(pp$Psi0)) pp$Psi0 <- diag(1, r - u) / 1e6
  if (is.null(pp$nu0))  pp$nu0  <- (r - u)
  if (is.null(pp$A0))   pp$A0   <- matrix(0, r - u, u)
  if (is.null(pp$K))    pp$K    <- 1e6 * diag(1, r - u)
  if (is.null(pp$L))    pp$L    <- 1e6 * diag(1, u)

  # sanity
  if (u < 0 || u > r) stop("u must be in [0, r].")

  if (u > 0 && u < r) {
    tmp <- do.call(
      MAPCA,
      c(list(X = X, Y = Y, u = u, maxiter = max.iter, method = method), pp)
    )

    CA      <- tmp$CA
    gamma   <- CA %*% sqrtmatinv(crossprod(CA))
    A       <- find_A_from_gamma(gamma)
    e.check <- tmp$e.check
    Yc      <- tmp$Yc
    G.tilde <- tmp$G.tilde

    gg      <- find_gammas_from_A(A)
    gamma   <- gg$gamma
    gamma0  <- gg$gamma0

    Omega  <- (crossprod(gamma,  G.tilde) %*% gamma  + pp$Psi ) / (n + p + pp$nu  + u     + 1)
    Omega0 <- (crossprod(gamma0, crossprod(Yc)) %*% gamma0 + pp$Psi0) / (n + p + pp$nu0 + (r-u) + 1)
    eta    <- crossprod(gamma, e.check)
    beta   <- gamma %*% eta
    Sigma  <- gamma %*% tcrossprod(Omega, gamma) + gamma0 %*% tcrossprod(Omega0, gamma0)

  } else if (u == 0) {
    # no envelope
    beta   <- matrix(0, r, p)
    gamma  <- NULL
    A      <- NULL
    eta    <- NULL
    gamma0 <- diag(r)

    Yc     <- scale(Y, center = TRUE, scale = FALSE)
    Omega0 <- (crossprod(Yc) + pp$Psi0) / (n + p + pp$nu0 + r + 1)
    Omega  <- NULL
    Sigma  <- Omega0

  } else { # u == r
    Xc   <- scale(X, center = TRUE, scale = FALSE)
    Yc   <- scale(Y, center = TRUE, scale = FALSE)
    sumY <- crossprod(Yc)
    tmp  <- crossprod(Xc) + pp$M
    tmp2 <- crossprod(Yc, Xc) + pp$e %*% pp$M

    e.check <- tmp2 %*% solve_chol(tmp)
    G.tilde <- sumY + pp$e %*% tcrossprod(pp$M, pp$e) - tcrossprod(tmp2, e.check)

    A      <- diag(1, r)      # dummy (not used)
    gamma  <- diag(1, r)
    gamma0 <- NULL
    eta    <- e.check
    beta   <- e.check
    Omega  <- (crossprod(gamma, G.tilde) %*% gamma + pp$Psi) / (n + p + pp$nu + r + 1)
    Omega0 <- NULL
    Sigma  <- Omega
  }

  mu <- colMeans(Y) - c(beta %*% colMeans(X))

  list(
    mu = mu, eta = eta, beta = beta, A = A,
    gamma = gamma, gamma0 = gamma0,
    Omega = Omega, Omega0 = Omega0,
    Sigma = Sigma, Sigma.inv = solve_chol(Sigma),
    prior_params = pp
  )
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
  if (length(dims) == 0) dims <- c(0L, 0L)
  ru <- dims[1]; u <- dims[2]

  # u=0 또는 r-u=0 => A가 빈 행렬이면 그대로 반환
  if (ru == 0L || u == 0L) {
    lp <- if (!missing(lpd_start) && !is.null(lpd_start)) {
      lpd_start
    } else {
      # lpd를 계산할 수 있으면 계산, 아니면 NA
      tryCatch(lpd_func(A = A_start, ...), error = function(e) NA_real_)
    }
    return(list(A = A_start, lpd = lp, accpt = matrix(0, nrow = ru, ncol = u)))
  }

  A <- A_start
  count <- 0L

  for (iter in seq_len(n_rwmh_iter)) {
    accpt.A <- matrix(0, ru, u)
    lpd_A   <- if (!missing(lpd_start) && !is.null(lpd_start) && iter == 1L) {
      lpd_start
    } else {
      lpd_func(A = A, ...)
    }

    if (!elementwise) {
      one_to_u <- if (random_scan) sample.int(u) else seq_len(u)
      for (j in one_to_u) {
        count <- count + 1L
        if (count > n_update_max) break

        tau_curr <- if (is.matrix(tau)) tau[, j] else if (length(tau) == ru) tau else rep(tau, ru)
        A_j_star <- A[, j] + rnorm(ru, 0, tau_curr)
        A_star   <- A; A_star[, j] <- A_j_star

        lpd_A_star <- lpd_func(A = A_star, ...)
        if (log(runif(1)) < (c(lpd_A_star) - c(lpd_A))) {
          A <- A_star; lpd_A <- lpd_A_star; accpt.A[, j] <- 1
        }
      }

    } else {
      n_full   <- length(A)
      n_update <- ceiling(n_full * random_update_fraction)
      n_update <- min(n_update, n_full)

      update_set <- if (n_update < n_full) {
        sample.int(n_full, n_update)
      } else if (random_scan) {
        sample.int(n_full)
      } else {
        seq_len(n_full)
      }

      if (n_update < n_full) accpt.A[-update_set] <- 0
      for (j in update_set) {
        count <- count + 1L
        if (count > n_update_max) break

        tau_curr <- if (length(tau) == n_full) tau[j] else if (length(tau) == 1L) tau else stop("tau length mismatch.")
        A_star   <- A; A_star[j] <- A[j] + rnorm(1, 0, tau_curr)

        lpd_A_star <- lpd_func(A = A_star, ...)
        if (log(runif(1)) < (c(lpd_A_star) - c(lpd_A))) {
          A <- A_star; lpd_A <- lpd_A_star; accpt.A[j] <- 1
        }
      }
    }
  }

  list(A = A, lpd = lpd_A, accpt = accpt.A)
}



log_I0 <- function(kappa) {
  if (kappa < 100) {
    kappa + log(besselI(kappa, 0, expon.scaled = TRUE))
  } else {
    kappa + Bessel::besselIasym(kappa, 0, expon.scaled = TRUE, log = TRUE)
  }
}

rgmatrixbingham22 <- function(A, B) {
  a <- -(A[1, 1] + B[2, 2])
  b <- -(A[2, 2] + B[1, 1])
  c <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])

  r <- sqrt((a - b)^2 + c^2)
  alpha <- atan2(y = c, x = a - b)

  theta <- tryCatch(
    BAMBI::minuspi_to_pi(
      BAMBI::rvm(n = 1, mu = alpha, kappa = r / 2)
    ) / 2,
    error = function(e) e
  )

  if (is(theta, "error")) browser()

  s <- ifelse(runif(1) < 0.5, 1, -1)

  Z <- matrix(
    c(
      cos(theta), s * sin(theta),
      sin(theta), -s * cos(theta)
    ),
    ncol = 2,
    byrow = TRUE
  )

  Z
}

dgmatrixbingham22 <- function(Z, A, B, C = 0, D = 0, log = TRUE, adjust = TRUE) {
  a <- -(A[1, 1] + B[2, 2])
  b <- -(A[2, 2] + B[1, 1])
  c <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])

  r <- sqrt((a - b)^2 + c^2)
  # alpha <- atan2(y = c, x = a-b)
  kappa <- r / 2


  out1 <- Z[1, 1]^2 * A[1, 1] +
    Z[2, 1]^2 * A[2, 2] +
    Z[1, 1] * Z[2, 1] * (A[1, 2] + A[2, 1])

  out2 <- Z[1, 2]^2 * B[1, 1] +
    Z[2, 2]^2 * B[2, 2] +
    Z[1, 2] * Z[2, 2] * (B[1, 2] + B[2, 1])


  log_bessel_kappa <- log_I0(kappa)

  out <- -(out1 + out2) -
    log_bessel_kappa -
    log(2 * pi) -
    (a + b) / 2

  if (!adjust) {
    out <- out + (a + b) / 2
  }

  if (!log) {
    out <- exp(out)
  }

  out
}

log_sum_exp <- function(x) {
  xm <- max(x)
  log(sum(exp(x - xm))) + xm
}

# vectorized in theta, not in s
lden_num_matrixbingham22_4param <- function(theta = 0.5, s = 1, A, B, C, D) {
  a2 <- -(A[1, 1] + B[2, 2])
  b2 <- -(A[2, 2] + B[1, 1])
  c2 <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])
  r2 <- sqrt((a2 - b2)^2 + c2^2)
  # alpha <- atan2(y = c, x = a-b)
  kappa2 <- r2 / 2
  alpha2 <- atan2(y = c2, x = a2 - b2)

  a1 <- C[1] - D[1] * s
  b1 <- C[2] + D[2] * s
  kappa1 <- sqrt(a1^2 + b1^2)
  alpha1 <- atan2(y = b1, x = a1)

  (kappa2 * cos(2 * theta - alpha2)) + (kappa1 * cos(theta - alpha1))
}

lden_num_matrixbingham22_4param_part2 <- function(theta = 0.5, s = 1, A, B, C, D) {
  a2 <- -(A[1, 1] + B[2, 2])
  b2 <- -(A[2, 2] + B[1, 1])
  c2 <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])
  r2 <- sqrt((a2 - b2)^2 + c2^2)
  # alpha <- atan2(y = c, x = a-b)
  kappa2 <- r2 / 2
  alpha2 <- atan2(y = c2, x = a2 - b2)

  # a1 <- C[1] - D[1] * s
  # b1 <- C[2] + D[2] * s
  # kappa1 <- sqrt(a1^2 + b1^2)
  # alpha1 <- atan2(y = b1, x = a1)

  (kappa2 * cos(2 * theta - alpha2)) #+ (kappa1 * cos(theta - alpha1))
}

lden_num_matrixbingham22_4param_part1 <- function(theta = 0.5, s = 1, A, B, C, D) {
  # a2 <- -(A[1, 1] + B[2, 2])
  # b2 <- -(A[2, 2] + B[1, 1])
  # c2 <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])
  # r2 <- sqrt((a2-b2)^2 + c2^2)
  # # alpha <- atan2(y = c, x = a-b)
  # kappa2 <- r2/2
  # alpha2 <- atan2(y = c2, x = a2-b2)

  a1 <- C[1] - D[1] * s
  b1 <- C[2] + D[2] * s
  kappa1 <- sqrt(a1^2 + b1^2)
  alpha1 <- atan2(y = b1, x = a1)

  # (kappa2 * cos(2*theta - alpha2)) #+
  (kappa1 * cos(theta - alpha1))
}



sigma2_from_kappa <- function(kappa) {
  # out <- if (kappa < 30) {
  #   2 * (
  #     log(besselI(kappa, 0, expon.scaled = TRUE)) -
  #       log(besselI(kappa, 1, expon.scaled = TRUE))
  #   )
  # # }
  # # else if (kappa < 50) {
  # #   2 * (
  # #     Bessel::besselIasym(kappa, 0, expon.scaled = TRUE, log = TRUE) -
  # #       Bessel::besselIasym(kappa, 1, expon.scaled = TRUE, log = TRUE)
  # #   )
  # # }
  # } else {
  #   # asymptotic normal approxiamtion
  #   1 / kappa
  # }

  out <- 1 / kappa

  out
}

calc_normal_params <- function(A,
                               B,
                               C = numeric(2),
                               D = numeric(2), ...) {
  # factor 1, depends on s
  fact1_s <- c(1, -1) %>%
    setNames(., .) %>%
    lapply(
      function(s) {
        a1 <- C[1] - D[1] * s
        b1 <- C[2] + D[2] * s
        kappa1 <- sqrt(a1^2 + b1^2)
        alpha1 <- atan2(y = b1, x = a1)

        list(mu = alpha1, sigma2 = sigma2_from_kappa(kappa1))
      }
    )

  # factor 2, does not depend on s
  a2 <- -(A[1, 1] + B[2, 2])
  b2 <- -(A[2, 2] + B[1, 1])
  c2 <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])
  r2 <- sqrt((a2 - b2)^2 + c2^2)
  # alpha <- atan2(y = c, x = a-b)
  kappa2 <- r2 / 2
  alpha2 <- atan2(y = c2, x = a2 - b2)
  fact2 <- list(mu = alpha2 / 2, sigma2 = sigma2_from_kappa(kappa2) / 4)

  # final normal density parameters

  out <- fact1_s %>%
    lapply(
      function(fact1) {
        sigma2 <- 1 / (1 / fact1$sigma2 + 1 / fact2$sigma2)
        mu <- (fact1$mu / fact1$sigma2 + fact2$mu / fact2$sigma2) * sigma2
        list(mu = mu, sigma2 = sigma2, sd = sqrt(sigma2))
      }
    )
}


calc_mix_norm_log_den <- function(x, mu1, mu2, sigma1, sigma2, pi1, pi2) {
  log_den <- dnorm(c(x, x), mean = c(mu1, mu2), sd = c(sigma1, sigma2), log = TRUE)
  log_pi <- c(log(pi1), log(pi2))
  log_den_plus_log_pi <- log_den + log_pi
  log_sum_exp(log_den_plus_log_pi)
}



# using normal approximation as suggested in mardia
rgmatrixbingham22_4param_an <- function(A,
                                        B,
                                        C = numeric(2),
                                        D = numeric(2),
                                        log_den = TRUE, ...) {
  normal_params <- calc_normal_params(A, B, C, D)
  # mixture weights proportional to the sd's
  s_1_prop <- normal_params$`1`$sd / (normal_params$`1`$sd + normal_params$`-1`$sd)
  s <- s_draw <- ifelse(runif(1) <= s_1_prop, 1, -1)

  log_dens_s <- if (s == 1) log(s_1_prop) else log(1 - s_1_prop)

  cond_params <- normal_params[[as.character(s)]]
  theta <- rnorm(
    n = 1,
    mean = cond_params$mu,
    sd = cond_params$sd
  )
  # log_den_theta_given_s <- dnorm(
  #   theta,
  #   mean = cond_params$mu,
  #   sd = cond_params$sd,
  #   log = TRUE
  # )
  # out_den <- log_dens_s + log_den_theta_given_s

  Z <- matrix(
    c(
      cos(theta), s * sin(theta),
      sin(theta), -s * cos(theta)
    ),
    ncol = 2,
    byrow = TRUE
  )


  # if  (!log_den) out_den <- exp(out_den)

  list(
    Z = Z,
    log_den = NA # out_den
  )
}


dgmatrixbingham22_4param_an <- function(Z, A, B,
                                        C = numeric(2),
                                        D = numeric(2),
                                        log_den = TRUE,
                                        nbreak = 5, ...) {
  theta <- atan2(y = Z[1, 2], x = Z[1, 1])
  s <- s_obs <- ifelse(Z[1, 2] * Z[2, 1] >= 0, 1, -1)

  # see random generation for more details
  normal_params <- calc_normal_params(A, B, C, D)
  s_1_prop <- normal_params$`1`$sd / (normal_params$`1`$sd + normal_params$`-1`$sd)
  log_dens_s <- if (s == 1) log(s_1_prop) else log(1 - s_1_prop)
  cond_params <- normal_params[[as.character(s)]]

  log_den_theta_given_s <- dnorm(
    theta,
    mean = cond_params$mu,
    sd = cond_params$sd,
    log = TRUE
  )
  out <- log_dens_s + log_den_theta_given_s

  if (!log_den) out <- exp(out)
  #
  out
}

dgmatrixbingham22_4param <- function(Z, A, B,
                                     C = numeric(2),
                                     D = numeric(2),
                                     log_den = TRUE,
                                     nbreak = 5) {
  # browser()
  #
  theta <- atan2(y = Z[1, 2], x = Z[1, 1])
  s <- ifelse(Z[1, 2] * Z[2, 1] >= 0, 1, -1)

  log_den_this_theta_s <- lden_num_matrixbingham22_4param(
    theta, s,
    A = A, B = B,
    C = C, D = D
  )

  theta_all <- seq(0, 2 * pi, length.out = nbreak)

  tmp1 <- lden_num_matrixbingham22_4param(
    theta_all,
    s = 1,
    A = A, B = B, C = C, D = D
  )
  tmp2 <- if (all(D == 0)) {
    tmp1
  } else {
    lden_num_matrixbingham22_4param(
      theta_all,
      s = 1,
      A = A, B = B, C = C, D = D
    )
  }

  log_den_table <- cbind(tmp1, tmp2)
  max_entry <- max(log_den_table, log_den_this_theta_s)

  log_norm_const <- exp(log_den_table - max_entry) %>%
    sum() %>%
    log()


  out <- log_den_this_theta_s - max_entry - log_norm_const



  # den_table <- log_den_table %>%
  #   {. - max(., log_den_this_theta_s)} %>%
  #   exp() %>%
  #   {./sum(.)}


  #   idx <- ifelse(s == 1, 1, 2)
  #
  #   consts_all <- normalizing_matrixbingham22_4param(
  #     A = A, B = B,
  #     C = C, D = D,
  #     nbreak = nbreak
  #   )
  #
  #   # browser()
  #
  #   params <- consts_all$fun_log_den_part_s(s)
  #
  #   den <- if (params$do_gauss_appox) {
  #     params$fun_approx_full(theta)
  #     # dnorm(theta, mean = params$mu_gauss_comb, sd = params$sigma_gauss_comb, log = TRUE)
  #   } else {
  #     params$fun_exact_num(theta) -
  #       consts_all$log_int_den_parts[idx]
  #   }
  #
  if (!log_den) out <- exp(out)
  #
  out
}


# random generation fromm density proportional to
# exp(kappa1 * cos(2*theta - alpha1) + kappa2 * cos(theta - alpha2))
# for kappa1, kappa2 > 0

# rvm_two_terms_proposal <- function(kappa1, alpha1, kappa2, alpha2) {
#
# }
#
# rvm_two_terms <- function(kappa1, alpha1, kappa2, alpha2) {
#   # browser()
#   log_2_pi <- log(2*pi)
#
#   log_const_1 <- -log(2) + log_2_pi + kappa1 + besselI(kappa1, 0, expon.scaled = TRUE)
#   log_const_2 <- log_2_pi + kappa2 + besselI(kappa2, 0, expon.scaled = TRUE)
#
#   if (log_const_1 <= log_const_2) {
#     mix_prop_1 <- 1 - 1/(1 + exp(log_const_1 - log_const_2))
#   } else {
#     mix_prop_1 <- 1/(1 + exp(log_const_2 - log_const_1))
#   }
#
#   mix_prop_2 <- 1 - mix_prop_1
#
#   log_mix_prop_1 <- log(mix_prop_1)
#   log_mix_prop_2 <- log(mix_prop_2)
#
#   accpt <- 0
#   count <- 1
#
#   while (accpt < 1) {
#
#     if (runif(1) <= mix_prop_1) {
#       # proposal from first component
#       proposal_theta <- BAMBI::rvm(
#         n = 1,
#         kappa = kappa1,
#         mu = alpha1
#       ) / 2
#     } else {
#       proposal_theta <- BAMBI::rvm(
#         n = 1,
#         kappa = kappa2,
#         mu = alpha2
#       )
#     }
#
#     term1 <- kappa1 * cos(2* proposal_theta - alpha1)
#     exponent1 <- term1 - log_const_1 + log_mix_prop_1
#
#     term2 <- kappa2 * cos(proposal_theta - alpha2)
#     exponent2 <- term2 - log_const_2 + log_mix_prop_2
#
#     if (exponent1 >= exponent2) {
#       log_den_bound <-
#         exponent1 + log( 1 + exp(exponent2 - exponent1) )
#     } else {
#       log_den_bound <-
#         exponent2 + log( 1 + exp(exponent1 - exponent2) )
#     }
#
#     log_den_target <- term1 + term2
#
#     log_ratio <- log_den_target - log_den_bound
#
#     if (log(runif(1)) <= log_ratio) {
#       accpt <- 1
#     } else {
#       count <- count + 1
#     }
#
#   }
#
#   proposal_theta
#
# }
#

rinvwish <- function(dim, Phi, nu) {
  p <- dim
  Phi.inv <- (solve_chol(Phi))
  Sigma.inv <- rWishart(1, nu, Phi.inv)[, , 1]
  (solve_chol(Sigma.inv))
}

rmvnorm <- function(n = 1, mu, sigma, dim = NULL) {
  if (is.null(dim)) {
    dim <- ncol(sigma)
  }
  p <- dim
  # eig <- eigen(Sigma)
  # Sigma.half <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Sigma.half <- t(chol(sigma))
  z <- rnorm(p)
  as.vector(mu + as.vector(Sigma.half %*% z))
}

solve_chol <- function(mat) chol2inv(chol(mat))


sqrtmat <- function(mat) {
  e <- eigen(mat, symmetric = TRUE)
  if (length(e$values) == 1) {
    sqrt(mat)
  } else {
    e$vectors %*% diag(c(sqrt(e$values))) %*% t(e$vectors)
  }
  # t(chol(mat))
}

sqrtmatinv <- function(mat) {
  e <- eigen(mat, symmetric = TRUE)
  if (length(e$values) == 1) {
    1 / sqrt(mat)
  } else {
    e$vectors %*% diag(1 / c(sqrt(e$values))) %*% t(e$vectors)
  }
  # t(chol(solve_chol(mat)))
}


mat_svd_sqrt_inv <- function(mat) {
  e <- eigen(mat, symmetric = TRUE)
  d_mat <- dim(mat)
  d1 <- d_mat[1]
  d2 <- d_mat[2]
  Dmat_list <- list(
    # orig = diag(e$values, d1, d2),
    inv = diag(1 / e$values, d1, d2),
    # sqrt = diag(sqrt(e$values), d1, d2),
    inv_sqrt = diag(1 / sqrt(e$values), d1, d2)
  )
  out <- lapply(
    Dmat_list,
    function(this_D) {
      e$vectors %*% this_D %*% t(e$vectors)
    }
  )
  out
}

# matrix volume
mat_vol <- function(A, logarithm = FALSE) {
  # sqrt(det(t(A) %*% A))
  log_det <- c(
    determinant(
      crossprod(A),
      logarithm = TRUE
    )$modulus
  )

  out <- 0.5 * log_det
  if (!logarithm) {
    out <- exp(out)
  }

  out
}
rMatrixNormal <- function(M, V1, V2) {
  # M is the mean matrix, and the V1, V2 are symmetric positive definite matrices

  alldims <- dim(M)
  n.rows <- alldims[1]
  n.cols <- alldims[2]
  z <- matrix(rnorm(n.rows * n.cols), nrow = n.rows, ncol = n.cols)

  A <- t(chol(V1)) # V1 = A %*% t(A)
  B <- chol(V2) # V2 = t(B) %*% B
  M + A %*% z %*% B
}

get_mle_resp <- function(X, Y, u) {
  ml <- env(X, Y, u)
  r <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  gamma <- ml$Gamma

  if (u > 0 & u < r) {
    G1 <- as.matrix(gamma[1:u, , drop = FALSE])
    # check if G1 is invertible - else reorganize the responses
    if (abs(det(G1)) < 1e-7) {
      gamma.t <- t(gamma)
      Y.order <- qr(gamma.t, tol = 1e-7)$pivot
      Y <- Y[, Y.order]
      gamma <- gamma[Y.order, , drop = FALSE]
    } else {
      Y.order <- 1:r
    }
  }

  Sig.Y <- ml$Sigma

  beta <- ml$beta
  mu_shift <- mu.start <- Y.bar <- colSums(as.matrix(Y)) / n
  mu <- mu_shift + beta %*% colMeans(X)

  if (u > 0 & u < r) {
    A <- find_A_from_gamma(gamma)
    gamma_gamma0 <- find_gammas_from_A(A)
    gamma <- gamma.start <- gamma_gamma0$gamma
    gamma0 <- gamma0.start <- gamma_gamma0$gamma0
    # Sig.res <- ml$Sigma
    Omega <- Omega.start <- t(gamma.start) %*% Sig.Y %*% gamma.start
    Omega0 <- Omega0.start <- t(gamma0.start) %*% Sig.Y %*% gamma0.start
    Sig <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0)
    Sig.inv <- (solve_chol(Sig))
    eta <- eta.start <- solve_chol(gamma_gamma0$CAtCA) %*%
      crossprod(gamma_gamma0$CA, beta)
  } else if (u == r) {
    A <- NULL
    gamma <- diag(1, r)
    gamma0 <- NULL
    Omega <- Sig.Y
    Omega0 <- NULL
    Sig <- Sig.Y
    eta <- beta
  } else if (u == 0) {
    A <- NULL
    gamma0 <- diag(1, r)
    gamma <- NULL
    Omega0 <- Sig.Y
    Omega <- NULL
    Sig <- Sig.Y
    eta <- NULL
  }


  list(
    mu = mu,
    mu_shift = mu_shift,
    eta = eta,
    beta = beta, A = A,
    gamma = gamma, gamma0 = gamma0,
    Omega = Omega, Omega0 = Omega0,
    Sigma = Sig.Y, Sigma.inv = solve_chol(Sig.Y)
  )
}

calc_gamma_lpd_resp <- function(gamma, gamma0, G.tilde, YctYc, Omega.inv, Omega0.inv) {
  -0.5 * (
    sum(
      Omega.inv * (crossprod(gamma, G.tilde) %*% gamma)
    ) +
      sum(
        Omega0.inv * (crossprod(gamma0, YctYc) %*% gamma0)
      )
  )
}

calc_lpiror_A_resp <- function(A, A0, K.half.inv, L.half.inv) {
  -0.5 * sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
}

lpd_A_resp_marginal <- function(A,
                                nu,
                                G.tilde,
                                Psi,
                                nu0,
                                YctYc,
                                Psi0,
                                n, r, u,
                                K.half.inv,
                                L.half.inv,
                                A0,
                                grad = FALSE,
                                ...) {
  # if (grad) {
  #   tmp <- find_gammas_from_A(
  #     A,
  #     jacobian = FALSE,
  #     return_jacob_mat = TRUE,
  #     jacobian_only_gamma = FALSE
  #   )
  # }


  tmp <- find_gammas_from_A(
    A,
    jacobian = FALSE,
    return_jacob_mat = FALSE
  )

  gamma <- tmp$gamma
  gamma0 <- tmp$gamma0

  gamma.t_G.tilde <- crossprod(gamma, G.tilde)
  gamma.t_G.tilde_gamma <- gamma.t_G.tilde %*% gamma

  gamma0.t_YctYc <- crossprod(gamma0, YctYc)
  gamma0.t_YctYc_gamma0 <- gamma0.t_YctYc %*% gamma0

  mat1 <- Psi + gamma.t_G.tilde_gamma
  term1 <- (n + nu - 1) * determinant(mat1, logarithm = TRUE)$modulus

  mat2 <- Psi0 + gamma0.t_YctYc_gamma0
  term2 <- (n + nu0 - 1) * determinant(mat2, logarithm = TRUE)$modulus

  term3 <- sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)

  # if (grad) {
  #   grad_term1_vec <- (n + nu - 1) *
  #     c(solve_chol(mat1)) %*% (
  #       Iu_kronecker_matM(u, gamma.t_G.tilde) +
  #         matM_times_commutation(
  #           matM_kronecker_Iu(u, gamma.t_G.tilde),
  #           k1 = r, k2 = u
  #         )
  #     ) %*% tmp$jacob_mat_gamma_over_A
  #   grad_term1 <- matrix(grad_term1_vec, r-u, u)
  #
  #   grad_term2_vec <- (n + nu0 - 1) *
  #     c(solve_chol(mat2)) %*% (
  #       Iu_kronecker_matM(r-u, gamma0.t_YctYc) +
  #         matM_times_commutation(
  #           matM_kronecker_Iu(r-u, gamma0.t_YctYc),
  #           k1 = r, k2 = r-u
  #         )
  #     ) %*% tmp$jacob_mat_gamma0_over_A
  #   grad_term2 <- matrix(grad_term2_vec, r-u, u)
  #
  #   grad_term3 <- 2 * mat3
  #
  #   grad_full <- -0.5 * (grad_term1 + grad_term2 + grad_term3)
  # }

  out <- -0.5 * (term1 + term2 + term3)
  attr(out, "gamma_gamma0") <- tmp
  # if (grad) attr(out, "grad") <- grad_full
  out
}

lpd_A_resp_semi_marginal <- function(A,
                                     G.tilde,
                                     YctYc,
                                     Omega.inv,
                                     Omega0.inv,
                                     K.half.inv,
                                     L.half.inv,
                                     A0, ...) {
  tmp <- find_gammas_from_A(A, ...)
  gamma <- tmp$gamma
  gamma0 <- tmp$gamma0

  gamma.t_G.tilde_gamma <- crossprod(gamma, G.tilde) %*% gamma
  gamma0.t_YctYc_gamma0 <- crossprod(gamma0, YctYc) %*% gamma0

  term1 <- sum(Omega.inv * gamma.t_G.tilde_gamma)
  term2 <- sum(Omega0.inv * gamma0.t_YctYc_gamma0)
  term3 <- sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)

  out <- -0.5 * (term1 + term2 + term3)
  attr(out, "gamma_gamma0") <- tmp
  out
}

.symm <- function(M) (M + t(M)) * 0.5
.proj_orth2 <- function(Z){
  # nearest O(2) via SVD
  sv <- svd(Z); sv$u %*% t(sv$v)
}
.safe_rgB22 <- function(A, B, max_retry = 3, jitter0 = 1e-10){
  A <- .symm(A); B <- .symm(B)
  jit <- jitter0
  for(k in 1:max_retry){
    Z <- try(rgmatrixbingham22(A = A + jit*diag(2), B = B + jit*diag(2)), silent = TRUE)
    if (!inherits(Z, "try-error") && all(is.finite(Z))) return(.proj_orth2(as.matrix(Z)))
    jit <- jit * 10
  }
  # fallback: identity rotation
  diag(2)
}
.safe_dgB22 <- function(Z, A, B, adjust = TRUE){
  A <- .symm(A); B <- .symm(B)
  out <- try(dgmatrixbingham22(Z = Z, A = A, B = B, log = TRUE, adjust = adjust), silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) return(-Inf)
  out
}

sample_A_metrop_manifold_proposal2 <- function(u, n, p, r,
                                              Omega.inv,
                                              Omega0.inv,
                                              current_A,
                                              current_gamma_gamma0,
                                              current_det_jacob_A_over_gamma_gamma0,
                                              G.tilde,
                                              YctYc,
                                              col_pairs_updt_set,
                                              A0,
                                              L.half.inv,
                                              K.half.inv,
                                              jacobian_only_gamma = FALSE,
                                              adjust_jacobian = FALSE,
                                              browser = FALSE,            # kept for signature
                                              random_scan = TRUE,
                                              adjust_gb22 = TRUE,
                                              n_update_max = 250,
                                              ...){

  # prior at start
  lprior_A_start <- calc_lpiror_A_resp(
    A = current_A, A0 = A0,
    K.half.inv = K.half.inv, L.half.inv = L.half.inv
  )

  gamma  <- current_gamma_gamma0$gamma
  gamma0 <- current_gamma_gamma0$gamma0

  # eigendecomp (SPD 가정 + 작은 ridge)
  eps <- 1e-12
  e1 <- eigen(.symm(Omega.inv),  symmetric = TRUE)
  e0 <- eigen(.symm(Omega0.inv), symmetric = TRUE)
  vals1 <- pmax(e1$values, eps)
  vals0 <- pmax(e0$values, eps)
  log_jacob_A_over_gamma_gamma0_start <- current_det_jacob_A_over_gamma_gamma0

  # proposal metric scalars
  lambda <- 1 / c(vals1, vals0)

  # base lpd
  lpd_gamma_gamma0_start <- calc_gamma_lpd_resp(
    gamma  = gamma, gamma0 = gamma0,
    G.tilde = G.tilde, YctYc = YctYc,
    Omega.inv = Omega.inv, Omega0.inv = Omega0.inv
  )
  lpd_A_start <- lpd_gamma_gamma0_start + log_jacob_A_over_gamma_gamma0_start + lprior_A_start

  # build Omat
  Omat <- cbind(
    gamma  %*% e1$vectors,
    gamma0 %*% e0$vectors
  )

  H_i <- function(i) if (i <= u) G.tilde else YctYc

  transition_adj <- 0
  count <- 0L

  update_set <- seq_len(length(col_pairs_updt_set))
  if (random_scan) update_set <- sample(update_set, replace = FALSE)

  for (idx in update_set) {
    count <- count + 1L
    if (count > n_update_max) break

    i <- col_pairs_updt_set[[idx]][1]
    j <- col_pairs_updt_set[[idx]][2]

    Nij <- Omat[, c(i, j), drop = FALSE]  # r x 2
    H_over_i <- .symm(H_i(i)) / lambda[i]
    H_over_j <- .symm(H_i(j)) / lambda[j]

    # forward kernels
    tmp_A_f <- 0.5 * crossprod(Nij, H_over_i) %*% Nij  # 2x2
    tmp_B_f <- 0.5 * crossprod(Nij, H_over_j) %*% Nij  # 2x2

    Z_GB22  <- .safe_rgB22(tmp_A_f, tmp_B_f)
    ld_new_given_old <- .safe_dgB22(Z_GB22, tmp_A_f, tmp_B_f, adjust = adjust_gb22)

    # propose Oij_new
    Oij_new <- Nij %*% Z_GB22
    Omat_old <- Omat
    Omat[, c(i, j)] <- Oij_new

    # backward kernels (using new block)
    tmp_A_b <- 0.5 * crossprod(Oij_new, H_over_i) %*% Oij_new
    tmp_B_b <- 0.5 * crossprod(Oij_new, H_over_j) %*% Oij_new
    ld_old_given_new <- .safe_dgB22(t(Z_GB22), tmp_A_b, tmp_B_b, adjust = adjust_gb22)

    if (!is.finite(ld_new_given_old) || !is.finite(ld_old_given_new)) {
      # rollback and skip this pair
      Omat <- Omat_old
      next
    }

    transition_adj <- transition_adj + (ld_old_given_new - ld_new_given_old)
  }

  # map back to (gamma, gamma0)
  gamma_end  <- tcrossprod(Omat[, 1:u,      drop = FALSE], e1$vectors)
  gamma0_end <- tcrossprod(Omat[, -(1:u),   drop = FALSE], e0$vectors)

  # A from gamma
  A_end <- find_A_from_gamma(gamma = gamma_end)

  lprior_A_end <- calc_lpiror_A_resp(
    A = A_end, A0 = A0,
    K.half.inv = K.half.inv, L.half.inv = L.half.inv
  )

  # two formulations → pick larger lpd
  gamma_gamma0_A_end <- find_gammas_from_A(
    A_end, jacobian = adjust_jacobian, log = TRUE,
    jacobian_only_gamma = jacobian_only_gamma
  )
  if (adjust_jacobian) {
    log_jacob_A_over_gamma_gamma0_end <- -gamma_gamma0_A_end$det_Jacobian_Omat_over_A
  } else {
    log_jacob_A_over_gamma_gamma0_end <- 0
  }

  lpd_pair <- sapply(
    list(list(gamma = gamma_end,  gamma0 = gamma0_end),
         gamma_gamma0_A_end),
    function(xx){
      calc_gamma_lpd_resp(
        gamma = xx$gamma, gamma0 = xx$gamma0,
        G.tilde = G.tilde, YctYc = YctYc,
        Omega.inv = Omega.inv, Omega0.inv = Omega0.inv
      )
    }
  )
  lpd_gamma_gamma0_end <- max(lpd_pair)

  lpd_A_end <- lpd_gamma_gamma0_end + log_jacob_A_over_gamma_gamma0_end + lprior_A_end

  # MH ratio (finite guard)
  metrop_log_ratio <- lpd_A_end - lpd_A_start + transition_adj
  accept <- is.finite(metrop_log_ratio) && (log(runif(1)) < metrop_log_ratio)

  if (accept) {
    A_out         <- A_end
    lpd_out       <- lpd_A_end
    accpt_out     <- 1L
    gamma_gamma0_out <- gamma_gamma0_A_end
    det_jacob_out <- log_jacob_A_over_gamma_gamma0_end
  } else {
    A_out         <- current_A
    lpd_out       <- lpd_A_start
    accpt_out     <- 0L
    gamma_gamma0_out <- current_gamma_gamma0
    det_jacob_out <- log_jacob_A_over_gamma_gamma0_start
  }

  attr(lpd_out, "gamma_gamma0") <- gamma_gamma0_out

  list(
    A = A_out,
    lpd = lpd_out,
    accpt = accpt_out,
    det_jacob_A_over_gamma_gamma0 = det_jacob_out
  )
}

sample_A_metrop_manifold_proposal <- function(u, n, p, r,
                                              Omega.inv,
                                              Omega0.inv,
                                              current_A,
                                              current_gamma_gamma0,
                                              current_det_jacob_A_over_gamma_gamma0,
                                              G.tilde,
                                              YctYc,
                                              col_pairs_updt_set,
                                              A0,
                                              L.half.inv,
                                              K.half.inv,
                                              jacobian_only_gamma = FALSE,
                                              adjust_jacobian = FALSE,
                                              browser = FALSE,
                                              random_scan = TRUE,
                                              adjust_gb22 = TRUE,
                                              n_update_max = 250,
                                              ...) {
  lprior_A_start <- calc_lpiror_A_resp(
    A = current_A,
    A0 = A0,
    K.half.inv = K.half.inv,
    L.half.inv = L.half.inv
  )

  gamma <- current_gamma_gamma0$gamma
  gamma0 <- current_gamma_gamma0$gamma0

  # convert from A to Omat (track the jacobian)
  Omega.inv_eig <- eigen(Omega.inv, symmetric = TRUE)
  Omega0.inv_eig <- eigen(Omega0.inv, symmetric = TRUE)
  log_jacob_A_over_gamma_gamma0_start <- current_det_jacob_A_over_gamma_gamma0

  # needed for the GB22 distribution
  lambda <- 1 / c(Omega.inv_eig$values, Omega0.inv_eig$values)


  lpd_gamma_gamma0_start <- calc_gamma_lpd_resp(
    gamma = gamma,
    gamma0 = gamma0,
    G.tilde = G.tilde,
    YctYc = YctYc,
    Omega.inv = Omega.inv,
    Omega0.inv = Omega0.inv
  )

  lpd_A_start <- lpd_gamma_gamma0_start +
    # Jacobian A -> (gamma, gamma0) -> A
    log_jacob_A_over_gamma_gamma0_start +
    lprior_A_start

  Omat <- cbind(
    gamma %*% Omega.inv_eig$vectors,
    gamma0 %*% Omega0.inv_eig$vectors
  )

  H_i <- function(i) {
    if (i <= u) {
      return(G.tilde)
    } else {
      return(YctYc)
    }
  }

  transition_adj <- 0
  count <- 0
  # all_transition <- c()

  counter <- 1
  # Omat_new <- Omat

  update_set <- seq_len(length(col_pairs_updt_set))
  if (random_scan) update_set <- sample(update_set, replace = FALSE)
  for (idx in update_set) {
    count <- count + 1
    if (count > n_update_max) break

    i <- col_pairs_updt_set[[idx]][1]
    j <- col_pairs_updt_set[[idx]][2]
    Nij <- Omat[, c(i, j)]
    H_over_i <- H_i(i) / lambda[i]
    H_over_j <- H_i(j) / lambda[j]


    # browser()

    tmp_A_forward <- 0.5 * crossprod(Nij, H_over_i) %*% Nij
    tmp_B_forward <- 0.5 * crossprod(Nij, H_over_j) %*% Nij

    log_dens_diff <- -Inf
    # log_dens_diff_thresh <- ifelse(counter == 1, -Inf, -5)


    # log_den_new_given_old <- Inf

    Z_GB22 <- tryCatch(
      rgmatrixbingham22(
        A = tmp_A_forward,
        B = tmp_B_forward
      ),
      error = function(e) e
    )

    log_den_new_given_old <- dgmatrixbingham22(
      Z = Z_GB22,
      A = tmp_A_forward,
      B = tmp_B_forward,
      log = TRUE,
      adjust = adjust_gb22
    )

    # }
    Oij_new <- Nij %*% Z_GB22
    Omat_old <- Omat
    Omat[, c(i, j)] <- Oij_new

    tmp_A_backward <- 0.5 * crossprod(Oij_new, H_over_i) %*% Oij_new
    tmp_B_backward <- 0.5 * crossprod(Oij_new, H_over_j) %*% Oij_new

    log_den_old_given_new <- dgmatrixbingham22(
      Z = t(Z_GB22),
      A = tmp_A_backward, # 0.5 * crossprod(Oij_new, H_over_i) %*% Oij_new,
      B = tmp_B_backward, # 0.5 * crossprod(Oij_new, H_over_j) %*% Oij_new,
      log = TRUE,
      adjust = adjust_gb22
    )


    # log_den_old_given_new <- dgmatrixbingham22(
    #   Z = t(Z_GB22),
    #   A = crossprod(Z_GB22, tmp_A_forward) %*% Z_GB22, #0.5 * crossprod(Oij_new, H_over_i) %*% Oij_new,
    #   B = crossprod(Z_GB22, tmp_B_forward) %*% Z_GB22, #0.5 * crossprod(Oij_new, H_over_j) %*% Oij_new,
    #   log = TRUE
    # )

    if (is(Z_GB22, "error")) browser()

    log_dens_diff <- log_den_old_given_new - log_den_new_given_old

    if (is.na(transition_adj) | is.na(log_dens_diff) | is.infinite(log_dens_diff)) browser()

    if (is.na(log_dens_diff)) {
      Omat <- Omat_old
      next
    } else {
      transition_adj <- transition_adj + log_dens_diff
    } # counter <- counter + 1
    # all_transition <- c(all_transition, log_dens_diff)
  }


  gamma_end <- tryCatch(
    tcrossprod(
      Omat[, 1:u, drop = FALSE],
      Omega.inv_eig$vectors
    ),
    error = function(e) e
  )

  if (is(gamma_end, "error")) browser()

  gamma0_end <- tcrossprod(
    Omat[, -(1:u), drop = FALSE],
    Omega0.inv_eig$vectors
  )


  # if (nrow(current_A) >= ncol(current_A)) {
  # r-u >= A
  A_end <- find_A_from_gamma(gamma = gamma_end)
  # } else {
  # A_end <- find_A_from_gamma0(gamma0 = gamma0_end)
  # }

  lprior_A_end <- calc_lpiror_A_resp(
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
  tmp_lpd_gamma_gamma0_vec <- sapply(
    list(
      list(
        gamma = gamma_end,
        gamma0 = gamma0_end
      ),
      gamma_gamma0_A_end
    ),
    function(this_gamma_gamma0) {
      calc_gamma_lpd_resp(
        gamma = this_gamma_gamma0$gamma,
        gamma0 = this_gamma_gamma0$gamma0,
        G.tilde = G.tilde,
        YctYc = YctYc,
        Omega.inv = Omega.inv,
        Omega0.inv = Omega0.inv
      )
    }
  )
  # take the maximum
  lpd_gamma_gamma0_end <- max(tmp_lpd_gamma_gamma0_vec)

  # browser()

  lpd_A_end <- lpd_gamma_gamma0_end +
    # Jacobian (gamma, gamma0) -> A
    log_jacob_A_over_gamma_gamma0_end +
    lprior_A_end

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

  attr(lpd_A, "gamma_gamma0") <- gamma_gamma0


  list(
    A = A,
    lpd = lpd_A,
    accpt = accpt,
    det_jacob_A_over_gamma_gamma0 = det_jacob
  )
}

find_A_from_gamma <- function(gamma, ...) {
  m <- ncol(gamma)
  G1 <- as.matrix(gamma[1:m, , drop = FALSE])
  # check if G1 is invertible - else reorganize the predictors
  if (abs(det(G1)) < 1e-7) {
    gamma.t <- t(gamma)
    X.order <- qr(gamma.t, tol = 1e-7)$pivot
    gamma <- gamma[X.order, , drop = FALSE]
  }
  G1 <- as.matrix(gamma[1:m, , drop = FALSE])
  G2 <- gamma[-(1:m), , drop = FALSE]
  G2 %*% solve(G1)
}

find_A_from_gamma0 <- function(gamma0, ...) {
  r <- nrow(gamma0)
  r_minus_u <- ncol(gamma0)
  u <- r - r_minus_u
  G2 <- as.matrix(gamma0[-(1:u), ])
  # check if G1 is invertible - else reorganize the predictors
  if (abs(det(G2)) < 1e-7) {
    gamma0.t <- t(gamma0)
    X.order <- qr(gamma0.t, tol = 1e-7)$pivot
    gamma0 <- gamma0[X.order, ]
  }
  G2 <- as.matrix(gamma0[-(1:u), ])
  G1 <- gamma0[1:u, ]
  -t(G1 %*% solve(G2))
}


find_A_from_gamma_gamma0 <- function(gamma, gamma0, ...) {
  0.5 * (
    find_A_from_gamma(gamma) +
      find_A_from_gamma0(gamma0)
  )
}


find_gammas_from_A <- function(A,
                               jacobian = FALSE,
                               log = TRUE,
                               jacobian_only_gamma = TRUE,
                               return_jacob_mat = FALSE) {
  dims <- dim(A)
  u <- dims[2]
  r <- sum(dims)

  CA <- matrix(0, nrow = r, ncol = u)
  CA[(u + 1):r, ] <- A
  CA[1:u, 1:u] <- diag(1, u)
  CAtCA <- crossprod(CA)
  svd_CAtCA <- mat_svd_sqrt_inv(CAtCA)
  CAtCA_minus_half <- svd_CAtCA$inv_sqrt
  CAtCA_inv <- svd_CAtCA$inv
  # CAtCA_half <- svd_CAtCA$sqrt
  gamma <- CA %*% CAtCA_minus_half


  DA <- matrix(0, nrow = r, ncol = r - u)
  DA[1:u, ] <- -t(A)
  DA[-(1:u), ] <- diag(1, r - u)
  DAtDA <- crossprod(DA)
  svd_DAtDA <- mat_svd_sqrt_inv(DAtDA)
  DAtDA_minus_half <- svd_DAtDA$inv_sqrt
  DAtDA_inv <- svd_DAtDA$inv
  # DAtDA_half <- svd_DAtDA$sqrt
  gamma0 <- DA %*% DAtDA_minus_half


  if (jacobian | return_jacob_mat) {
    Jmat_gamma_A <- jacobian_gamma_A(
      u = u, r = r, A = A,
      CA = CA, CAtCA = CAtCA,
      # CAtCA_half,
      CAtCA_inv = CAtCA_inv,
      CAtCA_minus_half = CAtCA_minus_half
    )


    if (jacobian_only_gamma) {
      # returns jacobian from only the gamma part
      Jmat_gamma0_A <- NULL
    } else {
      Jmat_gamma0_A <- jacobian_gamma0_A(
        u = u, r = r, A = A,
        DA = DA, DAtDA = DAtDA,
        # DAtDA_half,
        DAtDA_inv = DAtDA_inv,
        DAtDA_minus_half = DAtDA_minus_half
      )
    }

    if (jacobian) {
      Jmat_A <- rbind(Jmat_gamma_A, Jmat_gamma0_A)
      det_Jacobian_Omat_A <- mat_vol(Jmat_A, logarithm = log)
    } else {
      Jmat_A <- det_Jacobian_Omat_A <- NULL
    }
  } else {
    Jmat_gamma_A <- Jmat_gamma0_A <- det_Jacobian_Omat_A <- NULL
  }

  out <- list(
    gamma = gamma,
    gamma0 = gamma0,
    CA = CA, CAtCA = CAtCA,
    DA = DA, DAtDA = DAtDA,
    det_Jacobian_Omat_over_A = det_Jacobian_Omat_A
  )

  if (return_jacob_mat) {
    out$jacob_mat_gamma_over_A <- Jmat_gamma_A
    out$jacob_mat_gamma0_over_A <- Jmat_gamma0_A
  }

  out
}


jacobian_gamma_A_numeric <- function(A,
                                     jacobian_only_gamma = FALSE,
                                     log = TRUE, ...) {
  m <- nrow(A)
  n <- ncol(A)

  jaco <- numDeriv::jacobian(
    func = function(vecA) {
      matA <- matrix(vecA, m, n)
      tmp <- find_gammas_from_A(matA, jacobian = FALSE)
      if (jacobian_only_gamma) tmp$gamma0 <- NULL
      Omat <- cbind(tmp$gamma, tmp$gamma0)
      c(Omat)
    },
    x = c(A),
    ...
  )

  out <- mat_vol(jaco, logarithm = log)

  out
}

jacobian_gamma_A <- function(u, r, A,
                             CA, CAtCA,
                             CAtCA_inv,
                             CAtCA_minus_half) {
  # for d_CA.d_A
  # browser()
  M <- matrix(0, r, r - u)
  # I_r_minus_u <- diag(1, r - u, r - u)
  # M[-(1:u), ] <- I_r_minus_u
  diag(M[-(1:u), ]) <- 1


  d_CA.d_A <- Iu_kronecker_matM(u = u, M = M)

  # for d_CAtCA_minus_half.d_CA
  tmp11 <- Iu_kronecker_matM(u = u, M = CAtCA_minus_half)
  tmp12 <- commutation_times_matM(
    matM_times_commutation(tmp11, u, u),
    u, u
  ) # matM_kronecker_Iu(u, CAtCA_minus_half)
  tmp1 <- tmp11 + tmp12

  CAtCA_inv_CAt <- tcrossprod(CAtCA_inv, CA)

  tmp21 <- kronecker(CAtCA_inv, CAtCA_inv_CAt)
  tmp22 <- matM_times_commutation(
    kronecker(CAtCA_inv_CAt, CAtCA_inv),
    r, u
  )
  tmp2 <- tmp21 + tmp22

  d_CAtCA_minus_half.d_CA <- -solve_chol(tmp1) %*% tmp2

  # for d_CA_CAtCA_minus_half.d_CA
  term1 <- Iu_kronecker_matM(u, CA) %*% d_CAtCA_minus_half.d_CA
  term2 <- matM_kronecker_Iu(r, CAtCA_minus_half)

  d_CA_CAtCA_minus_half.d_CA <- term1 + term2

  # final result
  out <- d_CA_CAtCA_minus_half.d_CA %*% d_CA.d_A

  out
}


jacobian_gamma0_A <- function(u, r, A, DA, DAtDA,
                              DAtDA_inv,
                              DAtDA_minus_half) {
  # for d_DA.d_A
  I_r_minus_u_times_u <- diag(1, u * (r - u), u * (r - u))
  # M <- matrix(0, u * r, (r - u) * u)
  M <- matrix(0, r * (r - u), u * (r - u))
  # M[1:((r - u) * u), ] <- -I_r_minus_u_times_u
  diag(M[1:((r - u) * u), 1:((r - u) * u)]) <- -1

  # I_r_minus_u_times_u <- diag(1, u * (r - u), u * (r - u))
  # M <- matrix(0, r * (r - u), u * (r - u))
  # M[1:((r - u) * u), ] <- -I_r_minus_u_times_u

  d_DA.d_A <- commutation_times_matM(M, r - u, r)

  # for d_DAtDA_minus_half.d_DA
  tmp11 <- Iu_kronecker_matM(r - u, DAtDA_minus_half)
  tmp12 <- matM_kronecker_Iu(r - u, DAtDA_minus_half)
  tmp1 <- tmp11 + tmp12

  DAtDA_inv_DAt <- tcrossprod(DAtDA_inv, DA)
  tmp21 <- kronecker(DAtDA_inv, DAtDA_inv_DAt)
  tmp22 <- matM_times_commutation(
    kronecker(DAtDA_inv_DAt, DAtDA_inv),
    r, r - u
  )
  tmp2 <- tmp21 + tmp22
  d_DAtDA_minus_half.d_DA <- -solve_chol(tmp1) %*% tmp2

  term1 <- Iu_kronecker_matM(r - u, DA) %*% d_DAtDA_minus_half.d_DA
  term2 <- matM_kronecker_Iu(r, DAtDA_minus_half)
  d_DA_DAtDA_minus_half.d_DA <- term1 + term2

  # final result

  out <- d_DA_DAtDA_minus_half.d_DA %*% d_DA.d_A

  out
}

#' Fit Bayesian response envelope regression
#' @param X design (predictor) matrix. Must have the same number of rows as \code{Y}.
#' @param Y response matrix. Must have the same number of rows as \code{X}.
#' @param u response envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' @param n.iter Number of Markov chain iterations to run in each chains. *Includes burn-in*.
#' @param n.chains Number of independent chains to run.
#' @param autotune logical. Should the Metropolis tuning parameter
#'  be tuned during burn-in via adaptation?
#' @param tau,tune.accpt.prop.lower,tune.accpt.prop.upper,tune.burin.prop,tune.incr
#' \code{tau} is the tuning parameter in Random-walk Metropolis step.
#' @param init initialization. Available options are "mle" or "map".
#'
#'
#' @details `Benvlp_MC_resp` and `Benvlp_MC_resp_gibbs` are aliases.
#'
#' @examples
#' \dontrun{
#' library(future)
#' plan(multiprocess)
#'
#' set.seed(1)
#' MC <- Benvlp_MC_resp(
#'   X = wheatprotein[, 4:7],
#'   Y = wheatprotein[, 1:3],
#'   u = 1, n.iter = 1e4
#' )
#' }
#'
#' @md
#'
#' @export
Benvlp_MC_resp_gibbs <- function(X, Y, u,
                                 n.iter = 1000,
                                 n.chains = 1,
                                 tau = 1,
                                 autotune = TRUE,
                                 tune.accpt.prop.lower = 0.3,
                                 tune.accpt.prop.upper = 0.6,
                                 tune.incr = 0.05,
                                 burnin.prop = 0.3,
                                 tune.burnin.prop = 0.8,
                                 tune_nterm = 50,
                                 tune_scale_nterm = 50,
                                 checkstep = FALSE,
                                 autotune_size = "all",
                                 show_progress = TRUE,
                                 init = "map",
                                 chains_parallel = FALSE,
                                 rand_start_sd = 10,
                                 scan = "random",
                                 A_proposal = "stiefel",
                                 mh_type = "rw",
                                 nsteps_hmc = 1,
                                 jacobian_only_gamma = TRUE,
                                 posterior_gamma_space = TRUE,
                                 adjust_gb22 = TRUE,
                                 compute_llik = TRUE,
                                 rwmh_elementwise = TRUE,
                                 rwmh_elementwise_fraction = 1,
                                 n_A_update_max = Inf,
                                 ...) {

  tt1 <- Sys.time()

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  stopifnot(A_proposal %in% c("stiefel", "rwmh", "rwmh_marginal"))

  n <- nrow(X); p <- ncol(X); r <- ncol(Y)

  # center
  Yc   <- scale(Y, center = TRUE, scale = FALSE)
  YctYc <- crossprod(Yc)
  Y.bar <- attr(Yc, "scaled:center")

  Xc   <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X.bar <- attr(Xc, "scaled:center")

  # burn-in / tune
  n.burnin <- ceiling(n.iter * burnin.prop)
  n.iter.final <- n.iter - n.burnin
  n.tune <- ceiling(n.burnin * tune.burnin.prop)
  final_iter <- (n.burnin + 1L):n.iter

  dots <- list(...)
  prior_params <- dots$prior_params
  if (is.null(prior_params)) prior_params <- list()
  pp <- prior_params

  # ------- prior defaults (치수 0도 안전) -------
  if (is.null(pp$e))    pp$e    <- matrix(0, r, p)
  if (is.null(pp$M))    pp$M    <- 1e-6 * diag(1, nrow = p)

  if (is.null(pp$Psi))  pp$Psi  <- if (u>0) 1e-6 * diag(1, u) else matrix(0,0,0)
  if (is.null(pp$nu))   pp$nu   <- u
  
  if (is.null(pp$Psi0)) pp$Psi0 <- if (r-u>0) 1e-6 * diag(1, r-u) else matrix(0,0,0)
  if (is.null(pp$nu0))  pp$nu0  <- (r-u)

  if (u > 0 && u < r) {
    if (is.null(pp$A0)) pp$A0 <- matrix(0, r - u, u)
    if (is.null(pp$K))  pp$K  <- 1e6 * diag(1, r - u)
    if (is.null(pp$L))  pp$L  <- 1e6 * diag(1, u)
    pp$K.half.inv <- sqrtmatinv(pp$K)
    pp$L.half.inv <- sqrtmatinv(pp$L)
    if (length(tau) == 1) tau <- matrix(tau, r - u, u)
  }

  # common terms
  e  <- pp$e
  M  <- pp$M
  eM <- e %*% M
  XctXc_plus_M <- XctXc + M
  XctXc_plus_M.inv <- solve_chol(XctXc_plus_M)
  e.tilde <- (crossprod(Yc, Xc) + eM) %*% solve_chol(XctXc_plus_M)
  G.tilde <- YctYc + e %*% tcrossprod(M, e) - e.tilde %*% tcrossprod(XctXc_plus_M, e.tilde)

  Psi  <- pp$Psi
  Psi0 <- pp$Psi0
  nu   <- pp$nu
  nu0  <- pp$nu0

  if (u < 0 || u > r || u != as.integer(u)) {
    stop(sprintf("'u' must be an integer in [0, %d].", r))
  }

  # ---------- initializer ----------
  get_init <- function(...) {
    if (init == "mle") {
      get_mle_resp(...)
    } else if (init == "map") {
      Benvlp_resp_map(
        ...,
        A0,
        K = K, L = L,
        Psi = Psi, Psi0 = Psi0,
        e = e, M = M,
        nu = nu, nu0 = nu0,
        maxiter = 100
      )
    } else if (init == "random") {
      get_mle_rand_resp(..., sd = rand_start_sd)
    } else {
      stop("init must be one of 'mle'|'map'|'random'.")
    }
  }

  starting <- get_init(X = X, Y = Y, u = u)

  # gamma 열선택 안정화 (u in (0,r)일 때만)
  if (u > 0 && u < r) {
    det_gamma <- det(starting$gamma[1:u, , drop = FALSE])
    if (abs(det_gamma) < 1e-7) {
      Y.order <- qr(t(starting$gamma), tol = 1e-7)$pivot
      Y <- Y[, Y.order, drop = FALSE]
      # 재시작
      starting <- get_init(X = X, Y = Y, u = u)
    } else {
      Y.order <- seq_len(r)
    }
  } else {
    Y.order <- seq_len(r)
  }

  log_sqrt_2pi <- log(sqrt(2*pi))

  runMC <- function(starting = NULL, chain_no) {
    # progress bar
    pb <- NULL
    if (show_progress) {
      pb <- tryCatch(
        tcltk::tkProgressBar(title = paste("Chain", chain_no),
                             label = "Progress: 0%", min = 0, max = n.iter),
        error   = function(e) NULL,
        warning = function(w) NULL
      )
    }

    one_n <- rep(1, n)
    A_exists <- (u > 0 && u < r)

    mu_shift <- starting$mu
    eta  <- starting$eta
    beta <- starting$beta
    A    <- starting$A

    if (A_exists) {
      gamma_gamma0 <- find_gammas_from_A(
        A,
        jacobian = !posterior_gamma_space,
        log = TRUE,
        jacobian_only_gamma = jacobian_only_gamma
      )
    } else if (u == 0) {
      gamma_gamma0 <- list(gamma = NULL, gamma0 = diag(1, r))
    } else { # u == r
      gamma_gamma0 <- list(gamma = diag(1, r), gamma0 = NULL)
    }
    gamma  <- gamma_gamma0$gamma
    gamma0 <- gamma_gamma0$gamma0
    det_jacob_A_over_gamma_gamma0 <- 0

    Omega  <- starting$Omega
    Omega0 <- starting$Omega0

    if (u > 0)  Omega.inv  <- solve_chol(Omega)  else Omega.inv  <- NULL
    if (u < r)  Omega0.inv <- solve_chol(Omega0) else Omega0.inv <- NULL

    # Sigma
    if (A_exists) {
      Sigma <- gamma %*% tcrossprod(Omega, gamma) +
        gamma0 %*% tcrossprod(Omega0, gamma0)
    } else if (u == 0) {
      Sigma <- Omega0
    } else {
      Sigma <- Omega
    }
    Sigma.inv <- solve_chol(Sigma)

    mu <- mu_shift - beta %*% X.bar
    Y_minus_mu_shift <- Y - tcrossprod(one_n, mu_shift)

    # 사전 계산된 A 로그우도 (A_exists일 때만)
    if (A_exists) {
      lpd_A <- lpd_A_resp_semi_marginal(
        A = A,
        G.tilde = G.tilde,
        YctYc = YctYc,
        Omega.inv = Omega.inv,
        Omega0.inv = Omega0.inv,
        K.half.inv = pp$K.half.inv,
        L.half.inv = pp$L.half.inv,
        A0 = pp$A0,
        jacobian = (A_proposal == "stiefel") & !posterior_gamma_space,
        log = TRUE
      )
      gamma_gamma0 <- attr(lpd_A, "gamma_gamma0")
      gamma  <- gamma_gamma0$gamma
      gamma0 <- gamma_gamma0$gamma0
    }

    mu_shift.list <- mu.list <- eta.list <- beta.list <-
      Omega.list <- Omega0.list <- A.list <-
      gamma.list <- gamma0.list <-
      Sigma.inv.list <- Sigma.list <-
      llik.contri.list <- accpt.A.all <-
      vector("list", n.iter)

    lpd.A.all <- lpd.full.all <- llik.all <- accpt.A.ave <- rep(0, n.iter)

    # 열쌍 업데이트 세트 (A_exists일 때만 사용)
    if (A_exists) {
      col_pairs_all <- combn(1:r, 2)
      col_pairs_list <- lapply(seq_len(ncol(col_pairs_all)), function(jj) col_pairs_all[, jj])
      if (u <= r - u) {
        col_pairs_list <- col_pairs_list[sapply(col_pairs_list, function(x) any(x %in% 1:u))]
      } else {
        col_pairs_list <- col_pairs_list[sapply(col_pairs_list, function(x) any(x %in% (u + 1):r))]
      }
    }

    for (iter in 1:n.iter) {

      # ---- BLOCK: sample A (if 0<u<r) ----
      if (A_exists) {
        if (A_proposal == "stiefel") {
          #set.seed(100)
          mh_step <- sample_A_metrop_manifold_proposal2(
            u = u, n = n, p = p, r = r,
            Omega.inv = Omega.inv, Omega0.inv = Omega0.inv,
            current_A = A,
            current_gamma_gamma0 = list(gamma = gamma, gamma0 = gamma0),
            current_det_jacob_A_over_gamma_gamma0 = det_jacob_A_over_gamma_gamma0,
            G.tilde = G.tilde, YctYc = YctYc,
            col_pairs_updt_set = col_pairs_list,
            A0 = pp$A0,
            L.half.inv = pp$L.half.inv,
            K.half.inv = pp$K.half.inv,
            jacobian_only_gamma = jacobian_only_gamma,
            adjust_jacobian = !posterior_gamma_space,
            adjust_gb22 = adjust_gb22,
            n_update_max = n_A_update_max
          )
          det_jacob_A_over_gamma_gamma0 <- mh_step$det_jacob_A_over_gamma_gamma0

        } else if (A_proposal == "rwmh") {
          mh_step <- rwmh_colwise_A(
            A_start = A,
            tau = tau,
            random_scan = scan %in% c("random", "mixed"),
            n_update_max = n_A_update_max,
            elementwise = rwmh_elementwise,
            random_update_fraction = rwmh_elementwise_fraction,
            lpd_func = lpd_A_resp_semi_marginal,
            G.tilde = G.tilde, YctYc = YctYc,
            Omega.inv = Omega.inv, Omega0.inv = Omega0.inv,
            K.half.inv = pp$K.half.inv, L.half.inv = pp$L.half.inv,
            A0 = pp$A0
          )

        } else if (A_proposal == "rwmh_marginal") {
          mh_step <- rwmh_colwise_A(
            A_start = A,
            tau = tau,
            random_scan = scan %in% c("random", "mixed"),
            n_update_max = n_A_update_max,
            elementwise = rwmh_elementwise,
            random_update_fraction = rwmh_elementwise_fraction,
            lpd_func = lpd_A_resp_marginal,
            nu = nu, G.tilde = G.tilde, Psi = Psi,
            nu0 = nu0, YctYc = YctYc, Psi0 = Psi0,
            n = n, r = r, u = u,
            K.half.inv = pp$K.half.inv, L.half.inv = pp$L.half.inv,
            A0 = pp$A0, grad = FALSE
          )
        }

        # autotune (rwmh류만)
        if (autotune &&
            iter >= tune_nterm &&
            iter %% 5 == 0 &&
            iter <= n.burnin * tune.burnin.prop &&
            grepl("rwmh", A_proposal)) {
          nterm <- tune_nterm
          accpt.mean <- mean(accpt.A.ave[pmax(1, iter - nterm + 1):iter])
          if (accpt.mean > tune.accpt.prop.upper)      tau <- tau * (1 + tune.incr)
          else if (accpt.mean < tune.accpt.prop.lower) tau <- tau * (1 - tune.incr)
        }

        accpt.A <- mh_step$accpt
        lpd_A   <- mh_step$lpd
        accpt.A.all[[iter]] <- accpt.A
        A.list[[iter]] <- A <- mh_step$A
        lpd.A.all[iter] <- c(lpd_A)
        gg <- attr(lpd_A, "gamma_gamma0")
        gamma <- gg$gamma; gamma0 <- gg$gamma0

      } else if (u == r) {
        gamma.list[[iter]] <- diag(1, r)
      } else if (u == 0) {
        gamma0.list[[iter]] <- diag(1, r)
      }

      # ---- BLOCK: sample Omega / Omega0 ----
      if (A_exists) {
        Psi.tilde <- Psi + crossprod(gamma,  G.tilde) %*% gamma
        Omega.list[[iter]] <- Omega <- rinvwish(dim = u, Phi = Psi.tilde, nu = nu + n - 1)
        Omega.inv <- solve_chol(Omega)
      } else if (u == r) {
        Psi.tilde <- Psi + G.tilde
        Omega.list[[iter]] <- Omega <- rinvwish(dim = u, Phi = Psi.tilde, nu = nu + n - 1)
        Omega.inv <- solve_chol(Omega)
      } else { # u == 0
        Omega <- NULL; Omega.inv <- NULL
      }

      if (A_exists) {
        Psi.tilde.0 <- Psi0 + crossprod(gamma0, YctYc) %*% gamma0
        Omega0.list[[iter]] <- Omega0 <- rinvwish(dim = r - u, Phi = Psi.tilde.0, nu = nu0 + n - 1)
        Omega0.inv <- solve_chol(Omega0)
      } else if (u == 0) {
        Psi.tilde.0 <- Psi0 + YctYc
        Omega0.list[[iter]] <- Omega0 <- rinvwish(dim = r - u, Phi = Psi.tilde.0, nu = nu0 + n - 1)
        Omega0.inv <- solve_chol(Omega0)
      } else { # u == r
        Omega0 <- NULL; Omega0.inv <- NULL
      }

      # Sigma / Sigma.inv
      if (u > 0) {
        Sigma.inv.1 <- gamma  %*% tcrossprod(Omega.inv,  gamma)
      } else {
        Sigma.inv.1 <- matrix(0, r, r)
      }
      if (u < r) {
        Sigma.inv.2 <- gamma0 %*% tcrossprod(Omega0.inv, gamma0)
      } else {
        Sigma.inv.2 <- matrix(0, r, r)
      }
      Sigma.inv.list[[iter]] <- Sigma.inv <- Sigma.inv.1 + Sigma.inv.2

      if (A_exists) {
        Sigma.list[[iter]] <- Sigma <- solve_chol(Sigma.inv)
      } else if (u == 0) {
        Sigma.list[[iter]] <- Sigma <- Omega0
      } else { # u == r
        Sigma.list[[iter]] <- Sigma <- Omega
      }

      # ---- BLOCK: sample mu_shift, eta, then beta ----
      mu_shift.list[[iter]] <- mu_shift <- mvtnorm::rmvnorm(1, mean = c(Y.bar), sigma = Sigma / n)[1, ]
      mu_minus_Ybar <- mu_shift - Y.bar
      Y_minus_mu_shift <- Y - tcrossprod(one_n, mu_shift)

      if (u > 0) {
        eta.list[[iter]] <- eta <- rMatrixNormal(
          M = crossprod(gamma, e.tilde),
          V1 = Omega,
          V2 = XctXc_plus_M.inv
        )
      }

      if (A_exists) {
        beta.list[[iter]] <- beta <- gamma %*% eta
      } else if (u == r) {
        beta.list[[iter]] <- beta <- eta
      } else { # u == 0
        beta.list[[iter]] <- beta <- matrix(0, r, p)
      }

      mu.list[[iter]] <- mu <- mu_shift - beta %*% X.bar

      # ---- BLOCK: likelihood / posterior (optional) ----
      if (compute_llik || A_proposal == "rw") {
        log_det_Omega  <- if (u > 0)  log(det(Omega))  else 0
        log_det_Omega0 <- if (u < r)  log(det(Omega0)) else 0

        resi <- Y_minus_mu_shift - tcrossprod(Xc, beta)

        llik.contri.list[[iter]] <-
          -0.5 * r * log_sqrt_2pi +
          -0.5 * (log_det_Omega + log_det_Omega0 +
                    rowSums((resi %*% Sigma.inv.1) * resi) +
                    rowSums((Y_minus_mu_shift %*% Sigma.inv.2) * Y_minus_mu_shift))
        llik.all[iter] <- llik <- sum(llik.contri.list[[iter]])

        if (A_exists) {
          lpd.full.all[iter] <- llik - 0.5 * (
            (nu + u + p + 1) * log_det_Omega +
              sum((Psi + tcrossprod((eta - crossprod(gamma, e)) %*% sqrtmat(M))) * Omega.inv) +
              (nu0 + r - u + 1) * log_det_Omega0 +
              sum(Psi0 * Omega0.inv) +
              sum((pp$K.half.inv %*% (A - pp$A0) %*% pp$L.half.inv)^2)
          )
        } else if (u == r) {
          lpd.full.all[iter] <- llik - 0.5 * (
            (nu + r + p + 1) * log_det_Omega +
              sum((Psi + tcrossprod((eta - e) %*% sqrtmat(M))) * Omega.inv)
          )
        } else { # u == 0
          lpd.full.all[iter] <- llik - 0.5 * (
            (nu0 + r + 1) * log_det_Omega0 +
              sum(Psi0 * Omega0.inv)
          )
        }
      }

      if (!is.null(pb)) {
        tcltk::setTkProgressBar(
          pb, iter,
          label = sprintf("Progress: %d%%  %s",
                          round(iter / n.iter * 100),
                          ifelse(iter <= n.burnin, "(Burn-in)", "(Sampling)"))
        )
      }
    } # end for(iter)

    MC <- list(
      mu_shift = mu_shift.list[final_iter],
      mu       = mu.list[final_iter],
      eta      = eta.list[final_iter],
      beta     = beta.list[final_iter],
      Omega    = Omega.list[final_iter],
      Omega0   = Omega0.list[final_iter],
      Sigma    = Sigma.list[final_iter],
      Sigma.inv= Sigma.inv.list[final_iter],
      A        = A.list[final_iter],
      gamma    = gamma.list[final_iter],
      gamma0   = gamma0.list[final_iter],
      accpt.A.all = accpt.A.all[final_iter],
      lpd.A    = lpd.A.all[final_iter],
      accpt.A.ave = accpt.A.ave[final_iter],
      lpd.full = lpd.full.all[final_iter],
      llik     = llik.all[final_iter],
      llik.contri = llik.contri.list[final_iter],
      tau = tau,
      starting = starting,
      prior_param = pp
    )

    if (!is.null(pb) && inherits(pb, "tkProgressBar")) close(pb)
    MC
  }

  # 멀티체인 실행
  if (n.chains == 1 || !chains_parallel) {
    all_MCs <- lapply(seq_len(n.chains), function(j) runMC(starting, chain_no = j))
  } else {
    all_MCs <- future.apply::future_lapply(seq_len(n.chains), function(j) runMC(starting, chain_no = j),
                                           future.seed = TRUE)
  }

  vars <- setdiff(names(all_MCs[[1]]), "prior_param")
  out_vars <- lapply(vars, function(x) {
    out <- lapply(all_MCs, "[[", x)
    names(out) <- paste0("Chain_", seq_len(n.chains))
    out
  })
  names(out_vars) <- vars

  total_time <- as.numeric(difftime(Sys.time(), tt1, units = "secs"))

  out_res <- c(
    list(
      n.iter = n.iter,
      n.iter.final = n.iter.final,
      X = X, Y = Y,
      Y.order = Y.order,
      u = u,
      burnin = n.burnin,
      tune = n.tune,
      n.chains = n.chains,
      n_par = r * (r + 1) / 2 + r + p * u,
      prior_param = all_MCs[[1]]$prior_param,
      A_proposal = A_proposal,
      total_time = total_time
    ),
    out_vars
  )
  class(out_res) <- c("Benvlp", "Benvlp_resp")
  out_res
}

#' @rdname Benvlp_MC_resp_gibbs
#' @export
