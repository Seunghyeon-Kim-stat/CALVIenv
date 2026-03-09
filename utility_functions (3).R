
library(Matrix)
library(CholWishart)
library(matrixNormal)

sqrtmatinv <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*% diag(1 / sqrt(eig$values), nrow = nrow(mat)) %*% t(eig$vec)
}

sqrtmat <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*%
    diag(sqrt(eig$values), nrow = nrow(mat)) %*%
    t(eig$vec)
}

bmatrix <- function(x){

  if(!is.list(x)) stop("x not a list")

  sqrt(length(unlist(x)))

  n <- (-1 + sqrt(8*length(x)+1))/2
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]

  if(n==1) return(matrix(unlist(x), rr, cc))

  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
    stop("Zero-length component in x"))

  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(6, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,] <- rcum[n]
  ind[4,-1] <- ccum[-n]
  ind[5,] <- ccum
  ind[6,] <- ccum[n]
  imat <- array(1:(rsum * csum), c(rsum, csum))

  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[4]+1):y[6]], imat=imat)

  iuset <- apply(ind, 2, function(y, imat) imat[(y[4]+1):y[6],
                                                (y[1]+1):y[2]], imat=imat)

  iused <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                                (y[4]+1):y[5]], imat=imat)

  iuse <- as.vector(unlist(iuse))
  iuset <- as.vector(unlist(iuset))
  iused <- as.vector(unlist(iused))

  out[iuse] <- unlist(x)
  symmetric_iuse <- iuset[!iuset %in% iused]

  out[symmetric_iuse] <- t(out)[symmetric_iuse]
  return(out)
}

divide_matrix <- function(H, P, commutation = FALSE) {

  # library(matrixcalc)

  if (dim(H)[1] == 0 || dim(P)[1] == 0){
    return(0)
  }

  b <- dim(P)[1]
  a <- dim(H)[1] / b

  if(a == 1){
    commutation <- FALSE
  }

  if(b == 1){

    P <- c(P)
    return(P * H)

  }else if(b > 0){

    if(commutation == TRUE){
      H <- commutation.matrix(a,b) %*% H %*% commutation.matrix(b,a)
    }

    tr <- function(x) sum(diag(x))

    compute_trace_matrix <- function(blocks, P) {

      trace_values <- sapply(blocks, function(block) {
        tr(P %*% block)
      })

      matrix(trace_values, nrow = a, ncol = a, byrow = TRUE)
    }

    block_size <- b

    blocks <- list()

    for (i in 0:(a-1)) {
      for (j in 0:(a-1)) {
        row_start <- i * block_size + 1
        row_end <- (i + 1) * block_size
        col_start <- j * block_size + 1
        col_end <- (j + 1) * block_size
        block_name <- paste("H", i*a + j + 1, sep="")
        blocks[[block_name]] <- H[row_start:row_end, col_start:col_end]
      }
    }

    mat <- compute_trace_matrix(blocks, P)

    return(mat)
  }
}

make_positive_definite <- function(mat, tol = 1e-8) {

  nearPDmat <- nearPD(mat)

  mat_positive_definite <- as.matrix(nearPDmat$mat)

  return(mat_positive_definite)
}

logdetfunction <- function(x) {
  tryCatch({
    chol_decomp <- chol(x)
    return(2 * sum(log(diag(chol_decomp))))
  }, error = function(e) {
    # cat("Error in Cholesky decomposition:", conditionMessage(e), "\n")
    return(log(det(x)))
  })
}

HAfunction <- function(result){

  residualss <- lapply(result$A$Chain_1, function(X) X - (Reduce("+",result$A$Chain_1) / result$n.iter.final))
  udim <- dim(result$A$Chain_1[[1]])[1]
  vdim <- dim(result$A$Chain_1[[1]])[2]

  UU <- matrix(0, udim, udim)
  for (res in residualss){
    UU <- UU + res %*% t(res) / vdim
  }

  VV <- matrix(0, vdim, vdim)
  for (res in residualss){
    VV <- VV + t(res) %*% res / udim
  }

  return(kronecker((VV / length(residualss)), (UU / length(residualss))))
}

softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk)))
  }
  val <- exp(par - Lk)
  return(val)
}


normal_entropy <- function(mu, Sigma){

  r <- dim(Sigma)[1]

  if (r >= 1){

    value <-
      r * 0.5 * log(2 * pi * exp(1)) +
      0.5 * logdetfunction(Sigma)

  }else{

    value <- 0

  }
  return(value)

}


matrix_normal_entropy <- function(M, U, V){

  n <- dim(M)[1]
  m <- dim(M)[2]

  if (n >= 1 & m >= 1){

    value <-
      0.5 * n * m * log(2 * pi * exp(1)) +
      0.5 * n * logdetfunction(V) +
      0.5 * m * logdetfunction(U)

  }else{

    value <- 0

  }
  return(value)

}

invWishart_entropy <- function(nu, invPsi){

  p <- dim(invPsi)[1]

  if (p >= 1){

    value <-
      - (p + 1) * 0.5 * logdetfunction(invPsi) +
      nu * p * 0.5 * log(2 * exp(1)) +
      CholWishart::lmvgamma(0.5 * nu, p) -
      (nu + p + 1) * 0.5 * (CholWishart::mvdigamma(0.5 * nu, p) + p * log(2))

  }
  else{

    value <- 0

  }

  return(value)

}

expectation_loginvWishart <- function(nu, invPsi){

  p <- dim(invPsi)[1]

  if (p >= 1){

    value <-
      - CholWishart::mvdigamma(0.5 * nu, p) -
      p * log(2) -
      logdetfunction(invPsi)

  }else{

    value <- 0

  }

  return(value)

}

invWishart_logconst <- function(nu, Psi, p){

  # required package
  # CholWishart

  if (p >= 1){

    if (is.matrix(Psi)){
      logdetPsi <-  logdetfunction(Psi)
    }else{
      logdetPsi <-  logdetfunction(Psi * diag(p))
    }

    value <-
      -0.5 * (nu * p) * log(2) +
      0.5 * nu * logdetPsi -
      CholWishart::lmvgamma(0.5 * nu, p)

  }
  else{

    value <- 0

  }

  return(value)
}

normal_logconst <- function(H){

  n <- dim(H)[1]

  if (n >= 1){

    value <-
      -0.5 * n * log(2 * pi) -
      0.5 * logdetfunction(H)

  }
  else{

    value <- 0

  }

  return(value)
}

matrix_normal_logconst <- function(U, V){

  n <- dim(U)[1]
  p <- dim(V)[1]

  if (n >= 1 & p >= 1){

    value <-
      -0.5 * n * p * log(2 * pi) -
      0.5 * n * logdetfunction(V) -
      0.5 * p * logdetfunction(U)

  }
  else{

    value <- 0

  }

  return(value)

}

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

