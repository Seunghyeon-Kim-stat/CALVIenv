
library(mvtnorm) # multivariate normal distribution
library(matrixNormal) # matrix normal distribution
library(LaplacesDemon) # inverse wishart distribution
library(Renvlp) # envelope model MLE
library(CholWishart)
library(ggplot2)
library(matrixcalc)
library(tcltk)
library(profvis)
library(dplyr)

xenv_CALVI <- function(maxiter = 100,
                     n_iter = 10000,
                     show_progress=TRUE,
                     tol=1e-6,
                     make_positive=TRUE,
                     init,
                     compute_llik = TRUE,
                     X, Y, m, ...){

  progress_type <- NULL
  if (show_progress) {
    pb <- tryCatch({
      tcltk::tkProgressBar(title = paste0("n=", n, ", u=", u),
                           label = "Progress: 0%",
                           min = 0, max = n_iter)
    }, error = function(e) NULL, warning = function(w) NULL)
    progress_type <- "tk"
  } else {
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    progress_type <- "console"
  }

  start_time <- Sys.time()

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  dim_X <- dim(X)
  dim_Y <- dim(Y)

  n <- dim_X[1]
  p <- dim_X[2]
  r <- dim_Y[2]

  Y <- data.matrix(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  Y_bar <- attr(Yc,"scaled:center")

  X <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  X_bar <- attr(Xc,"scaled:center")

  onevector <- matrix(1, n, 1)

  Z <- matrix(0, nrow=(p-m), ncol=m)
  K <- rbind(diag(1, m), Z)
  L <- rbind(t(Z), diag(1, p-m))

  prior_pars <- xenv_get_prior(r, p, m)

  B0 <- prior_pars$B0
  M <- prior_pars$M

  BMB <- B0%*%M%*%t(B0)

  variational_pars <- init_pars <- xenv_get_init_MLE(X, Y, n, r, p, m)

  Lvaluelist <- list()

  best_Lvalue <- -Inf
  i <- 1

  while (i < n_iter) {

    if (progress_type == "tk") {
      tcltk::setTkProgressBar(pb, i,
                              label = paste0("Progress: ", round(i/n_iter * 100), "%"))
    } else if (progress_type == "console") {
      setTxtProgressBar(pb, i)
    }

    GpX <- crossprod( (X - onevector %*% t(variational_pars[["muX_q"]])) ) + n * variational_pars[["SigmaX_q"]]
    GpY <- crossprod( (Y - onevector %*% t(variational_pars[["muY_q"]])) ) + n * variational_pars[["SigmaY_q"]]

    tYmuXmu <- t( (Y - onevector %*% t(variational_pars[["muY_q"]])) ) %*% (X - onevector %*% t(variational_pars[["muX_q"]]))

    Lvalue0 <- do.call(
      xenv_convergence,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    hatA <- do.call(
      xenv_A_LVI,
      c(list(maxiter = maxiter, method = "BFGS"),
        prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )
    ## Update mean of q(A)
    variational_pars[["hatA"]] <- hatA[[1]]

    ## Update variance of q(A)
    tmpHA <- hatA[[2]]

    if (m > 0 & m < p){
      if (make_positive == TRUE & !is.positive.definite(tmpHA)){
        tmpHA <- make_positive_definite(tmpHA)
      }
    }

    variational_pars[["HA"]] <- tmpHA

    variational_pars[["CA"]] <- hatA[[3]]
    variational_pars[["DA"]] <- hatA[[4]]

    # Update tilde_eta_q
    tilde_eta_q <- do.call(
      xenv_eta_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(eta)
    variational_pars[["eta_q"]] <- tilde_eta_q[[1]]

    tmpSigmaeta_q <- tilde_eta_q[[2]]

    if (m > 0 & m < p){
      if (make_positive == TRUE & !is.positive.definite(tmpSigmaeta_q)){
        tmpSigmaeta_q <- make_positive_definite(tmpSigmaeta_q)
      }
    }

    variational_pars[["Sigmaeta_q"]] <- tmpSigmaeta_q


    # #Update SigmaYcX
    SigmaY_q <- do.call(
      xenv_SigmaYcX_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(Omega)
    variational_pars[["nuY_q"]] <- SigmaY_q[[1]]

    tmpinvPsiY_q <- SigmaY_q[[2]]

    if (m > 0 & m < p){
      if (make_positive == TRUE & !is.positive.definite(tmpinvPsiY_q)){
        tmpinvPsiY_q <- make_positive_definite(tmpinvPsiY_q)
      }
    }

    variational_pars[["invPsiY_q"]] <- tmpinvPsiY_q

    # Update tilde_Omega_q
    tilde_OmegaX1_q <- do.call(
      xenv_Omega_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(Omega)
    variational_pars[["nuX1_q"]] <- tilde_OmegaX1_q[[1]]

    tmpinvPsiX1_q <- tilde_OmegaX1_q[[2]]

    if (m > 0 & m < p){
      if (make_positive == TRUE & !is.positive.definite(tmpinvPsiX1_q)){
        tmpinvPsiX1_q <- make_positive_definite(tmpinvPsiX1_q)
      }
    }

    variational_pars[["invPsiX1_q"]] <- tmpinvPsiX1_q

    # Update tilde_Omega0_q
    tilde_OmegaX0_q <- do.call(
      xenv_Omega0_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(Omega0)
    variational_pars[["nuX0_q"]] <- tilde_OmegaX0_q[[1]]

    tmpinvPsiX0_q <- tilde_OmegaX0_q[[2]]
    if (m > 0 & m < p){
      if (make_positive == TRUE & !is.positive.definite(tmpinvPsiX0_q)){
        tmpinvPsiX0_q <- make_positive_definite(tmpinvPsiX0_q)
      }
    }

    variational_pars[["invPsiX0_q"]] <- tmpinvPsiX0_q


    # Update muX_q
    muX_q <- do.call(
      xenv_muX_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(mu)
    variational_pars[["muX_q"]] <- muX_q[[1]]
    variational_pars[["SigmaX_q"]] <- muX_q[[2]]

    # Update muY_q
    muY_q <- do.call(
      xenv_muY_LVI,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    ## Update parameters of q(mu)
    variational_pars[["muY_q"]] <- muY_q[[1]]
    variational_pars[["SigmaY_q"]] <- muY_q[[2]]

    GpX <- crossprod( (X - onevector %*% t(variational_pars[["muX_q"]])) ) + n * variational_pars[["SigmaX_q"]]
    GpY <- crossprod( (Y - onevector %*% t(variational_pars[["muY_q"]])) ) + n * variational_pars[["SigmaY_q"]]

    tYmuXmu <- t( (Y - onevector %*% t(variational_pars[["muY_q"]])) ) %*% (X - onevector %*% t(variational_pars[["muX_q"]]))

    Lvalue <- do.call(
      xenv_convergence,
      c(prior_pars,
        variational_pars,
        list(X=X, Y=Y, X_bar=X_bar, Y_bar=Y_bar, m=m, n=n, r=r, p=p,
             GpX=GpX, GpY=GpY, BMB=BMB, tYmuXmu=tYmuXmu, K=K, L=L))
    )

    Lvaluelist <- c(Lvaluelist, Lvalue)


    # Check for convergence

    if (abs(Lvalue0 - Lvalue) < tol * abs(Lvalue0)){
      break
    }else{
      Lvalue0 <- Lvalue
      i <- i + 1
    }

  }

  last_Lvalue <- Lvalue
  end_time <- Sys.time()  # End timing

  variational_pars[["beta"]] <- t(variational_pars[["CA"]] %*% variational_pars[["eta_q"]])

  if (compute_llik) {
    llik_vec <- xenv_llik(
      X, Y, N = 1,
      muX_q       = variational_pars$muX_q,
      SigmaX_q    = variational_pars$SigmaX_q,
      muY_q       = variational_pars$muY_q,
      SigmaY_q    = variational_pars$SigmaY_q,
      eta_q      = variational_pars$eta_q,
      Sigmaeta_q = variational_pars$Sigmaeta_q,
      nuX1_q       = variational_pars$nuX1_q,
      invPsiX1_q   = variational_pars$invPsiX1_q,
      nuX0_q      = variational_pars$nuX0_q,
      invPsiX0_q  = variational_pars$invPsiX0_q,
      nuY_q      = variational_pars$nuY_q,
      invPsiY_q  = variational_pars$invPsiY_q,
      hatA       = variational_pars$hatA,
      HA         = variational_pars$HA
    )
    llik_mc <- llik_vec
    LVI_BIC <- -2 * llik_mc + log(n) * (r + r*(r+1)/2 + p + p*(p+1) + r*m)
  } else {
    llik_mc <- NA
    LVI_BIC <- NA
  }

  total_run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  close(pb)

  return(list(variational_pars,
              init_pars,
              total_run_time,
              last_Lvalue,
              Lvaluelist,
              llik_mc,
              LVI_BIC
  ))

}
