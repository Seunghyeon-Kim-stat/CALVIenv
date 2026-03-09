library(rstan)
library(magrittr)

# rstan_options(auto_write = FALSE)
# stan_file <- list.files(recursive = TRUE, pattern = "\\.stan$") %>%
#   grep("response_envelope", ., value = TRUE, ignore.case = TRUE)
# response_envelope_stan <- rstan::stan_model(file = "response_envelope.stan")

# expose_stan_functions(response_envelope_stan)

lpd_A_resp_new <- function(
    A,
    Y_minus_mu_shift,
    Xc,
    eta,
    Omega.inv,
    Omega0.inv,
    e,
    M.half,
    K.half.inv,
    L.half.inv, 
    A0) {
  
  tmp <- find_gammas_from_A(A)
  gamma <- tmp$gamma
  gamma0 <- tmp$gamma0
  
  Y_minus_mu_shift_minus_Xbeta <- Y_minus_mu_shift - tcrossprod(Xc, gamma %*% eta)
  
  term1 <- sum(Omega.inv * crossprod(Y_minus_mu_shift_minus_Xbeta %*% gamma))
  
  term2 <- sum(Omega0.inv * crossprod(Y_minus_mu_shift %*% gamma0))
  
  term3 <- sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
  
  out <-  - 0.5 * (term1 + term2 + term3)
  attr(out, "gamma_gamma0") <- tmp
  out
}

Benvlp_resp_vb <- function(X, Y, u,
                           n.iter = 1000,
                           n.chains = 1,
                           tau = 1,
                           autotune = TRUE,
                           # tune.accpt.prop.lower = 0.4,
                           # tune.accpt.prop.upper = 0.5,
                           # tune.burnin.prop = 0.5,
                           # tune.incr = 0.05,
                           # burnin.prop = 0.3,
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
                           # tune_nterm = 50,
                           init = "map",
                           chains_parallel = FALSE,
                           rand_start_sd = 10,
                           scan = "random",
                           mh_type = "rw",
                           nsteps_hmc = 1,
                           do_mcmc = FALSE,
                           ...)
{
  tt1 <- Sys.time()

  X <- data.matrix(X)
  Y <- data.matrix(Y)


  dim_X <- dim(X)
  dim_Y <- dim(Y)

  n <- dim_X[1]
  p <- dim_X[2]
  r <- dim_Y[2]


  Y <- Y.orig <- data.matrix(Y)
  Y.bar <- colMeans(Y)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  YctYc <- crossprod(Yc)
  Y.bar <- attr(Yc,"scaled:center")


  X <- X.orig <- data.matrix(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  XctXc <- crossprod(Xc)
  X.bar <- attr(Xc,"scaled:center")

  dim_U <- dim(Y)
  r <- dim_U[2]
  n <- dim_U[1]
  p <- ncol(X)


  n.burnin <- ceiling(n.iter * burnin.prop)
  burnin_iter <- seq_len(n.burnin)
  n.tune <- ceiling(n.burnin * tune.burnin.prop)
  tune_iter <- seq_len(n.tune)

  dots <- list(...)
  prior_params <- dots$prior_params


  if (is.null(prior_params)) {
    prior_params = list()
  }

  pp <- prior_params

  if (is.null(pp$e)) {
    e <- matrix(0, r, p)
  }

  if (is.null(pp$M)) {
    M <- 1/1e6 * diag(1, nrow = p)
  }


  M.half <- sqrtmat(M)
  M.inv <- solve_chol(M)
  M.half.inv <- sqrtmat(M.inv)
  XctXc_plus_M <- XctXc + M
  XctXc_plus_M.half <- sqrtmat(XctXc_plus_M)
  XctXc_plus_M.inv <- solve_chol(XctXc_plus_M)
  # YtXc <- crossprod(Y, Xc)
  eM <- e %*% M
  # eM.inv <- e %*% solve_chol(M)


  e.tilde <- (crossprod(Yc, Xc) + eM) %*% solve_chol(XctXc_plus_M)
  # M.tilde <- XctXc_plus_M
  # M.tilde.inv <- XctXc_plus_M.inv
  G.tilde <- YctYc + e %*% tcrossprod(M, e) -
    e.tilde %*% tcrossprod(XctXc_plus_M, e.tilde)
  # G.tilde.half <- sqrtmat(G.tilde)
  # tmp <- tryCatch(chol(G.tilde), error = function(e) e)
  # if (is(tmp, "error")) browser()


  # browser()

  if (is.null(pp$Psi)) {
    Psi <- u *  diag(1, nrow = u) / 1e6
  }

  if (is.null(pp$nu)) {
    nu <- u
  }

  if (is.null(pp$Psi0)) {
    Psi0 <- (r-u) *  diag(1, nrow = r-u) / 1e6
  }

  if (is.null(pp$nu0)) {
    nu0 <- r-u
  }

  if (u > 0 & u < r) {

    if (is.null(pp$A0)) {
      A0 <- matrix(0, r-u, u)
    }

    if (is.null(pp$K)) {
      K <- 1e6 * diag(1, nrow = r-u)
    }

    K.inv <- solve_chol(K)
    K.half.inv <- sqrtmatinv(K)

    if (is.null(pp$L)) {
      L <- 1e6 * diag(1, nrow = u)
    }

    L.inv <- solve_chol(L)
    L.half.inv <- sqrtmatinv(L)

    if (length(tau) == 1)
      tau <- matrix(tau, r-u, u)

  }

  n.burnin <- ceiling(n.iter * burnin.prop)
  burnin_iter <- seq_len(n.burnin)
  n.tune <- ceiling(n.burnin * tune.burnin.prop)
  tune_iter <- seq_len(n.tune)

  n.final <- n.iter - n.burnin

  if (u < 0 || u > r || u != as.integer(u))
    stop(paste0("\'u\' must be an integer between 0 and ", r, " (r)") )



  get_init <- function(...) {
    if (init == "mle") {
      get_mle_resp(...)
    } else if (init == "map") {
      Benvlp_resp_map(
        ...,
        A0, K = K, L = L,
        Psi = Psi, Psi0 = Psi0,
        e = e, M = M,
        nu = nu, nu0 = nu0,
        maxiter = 100
      )
    } else if (init == "random") {
      get_mle_rand_resp(..., sd = rand_start_sd)
    }
  }

  starting <- get_init(X = X, Y = Y, u = u)

  if (u > 0 & u < r) {
    det_gamma <- det(starting$gamma[1:u, , drop = FALSE])
  } else {
    det_gamma <- 1
  }
  if (u > 0 & abs(det_gamma) < 1e-7) {
    # gamma.t <- t(gamma)
    Y.order <- qr(t(gamma), tol = 1e-7)$pivot
    Y <- Y[, Y.order]
    starting <- get_init(X = X, Y = Y, u = u)
  } else {
    Y.order <- 1:r
  }

  log_sqrt_2pi <- log(sqrt(2*pi))

  # browser()
  #
  # # rejection sampling -- for testing in small dimension
  # tmpfun <-  function(A) {
  #   A <- matrix(c(A), nrow(A0), ncol(A0))
  #   -c(lpd_A_resp_marginal(
  #     A,
  #     nu,
  #     G.tilde,
  #     Psi,
  #     nu0,
  #     YctYc,
  #     Psi0,
  #     n, r, u,
  #     K.half.inv,
  #     L.half.inv,
  #     A0
  #   )) - 0.5 * sum((K.half.inv %*% (A - A0) %*% L.half.inv)^2)
  # }
  #
  # tmp_opt <- optim(
  #   par = starting$A,
  #   fn = tmpfun,
  #   method = "OBFGS"
  # )


  runMC <- function(starting = NULL, chain_no)
  {
    pb <- tryCatch(tkProgressBar(title = paste("Chain", chain_no),
                                 label = "Progress: 0%",
                                 min = 0, max = n.iter),
                   error = function(e) e,
                   warning = function(w) w)

    if (is(pb, "error") | is(pb, "warning")) {
      no_pb <- TRUE
      show_progress <- FALSE
      warning(paste("\'tcltk\' could not be loaded. \'show_progress\' set to FALSE.\n\n"))
    }else{
      no_pb <- FALSE
    }


    runif(chain_no)

    one_n <- rep(1, n)

    mu_shift <- starting$mu
    eta <- starting$eta
    beta <- starting$beta
    A <- starting$A
    gamma_gamma0 <- find_gammas_from_A(A)
    gamma <- gamma_gamma0$gamma
    gamma0 <- gamma_gamma0$gamma0

    # during iterations, will be computed via "find_gamma_gamma0"
    Omega <- starting$Omega
    if (u > 0)  {
      Omega.inv <- solve_chol(Omega)
    } else {
      Omega.inv <- NULL
    }

    Omega0 <- starting$Omega0
    if (u < r) {
      Omega0.inv <- solve_chol(Omega0)
    } else {
      Omega0.inv <- NULL
    }

    A_exists <- u > 0 & u < r

    if (A_exists) {
      Sigma <- gamma %*% tcrossprod(Omega, gamma) +
        gamma0 %*% tcrossprod(Omega0, gamma0)
    } else if (u == 0) {
      Sigma <- Omega0
    } else if (u == r) {
      Sigma <- Omega
    }
    Sigma.inv <- solve_chol(Sigma)

    mu <- mu_shift - beta %*% X.bar
    Y_minus_mu_shift <- Y - tcrossprod(one_n, mu_shift)

    # if (A_exists) {
    #   # lpd_A <- lpd_A_resp_new(A = A,
    #   #                         Y_minus_mu_shift = Y_minus_mu_shift,
    #   #                         Xc = Xc,
    #   #                         eta = eta,
    #   #                         Omega.inv = Omega.inv,
    #   #                         Omega0.inv = Omega0.inv,
    #   #                         e = e,
    #   #                         M.half = M.half,
    #   #                         K.half.inv = K.half.inv,
    #   #                         L.half.inv = L.half.inv,
    #   #                         A0 = A0)
    #   #
    #   # gamma_gamma0 <- attr(lpd_A, "gamma_gamma0")
    #   # gamma <- gamma_gamma0$gamma
    #   # gamma0 <- gamma_gamma0$gamma0
    # }

    mu_shift.list <- mu.list <-
      eta.list <- beta.list <-
      Omega.list <- Omega0.list <- A.list <-
      gamma.list <- gamma0.list <-
      Sigma.inv.list <- Sigma.list <-
      llik.contri.list <-
      accpt.A.all <- vector("list", n.iter)

    lpd.full.all <- lpd.A.all <- llik.all <-  accpt.A.ave <- rep(0, n.iter)


    # if (A_exists) {
    #   browser()
    #   K_r_r.minus.u <- matrixcalc::commutation.matrix(r, r-u)
    #   Del_vec_CA <- 1
    # }


    # if (mh_type == "rw") {
    #   mh_fun <- rwmh_colwise_A
    # } else if (mh_type == "langevin") {
    #   # mh_fun <-
    #   #   function(...) langevin_A(..., nsteps_hmc = nsteps_hmc)
    # }

    # browser()

    Diag1 <- diag(1e-3, u)
    Diag0 <- diag(1e-3, (r-u))

    if (A_exists) {
      stan_dt <- list(
        n = n, #int<lower=0> n;
        r = r, #int<lower=0> r;
        u = u, #int<lower=0> u;
        r_minus_u = r-u, #int<lower=0> r_minus_u;
        nu = nu, #int<lower=0> nu;
        nu0 = nu0, #int<lower=0> nu0;
        Psi = Psi, #matrix[u,u] Psi;
        Diag1 = Diag1,
        Psi0 = Psi0, #matrix[r_minus_u,r_minus_u] Psi0;
        Diag0 = Diag0,
        G_tilde = G.tilde, #matrix[r,r] G_tilde;
        YctYc = YctYc, #matrix[r,r] YctYc;
        A0 = A0, #matrix[r_minus_u,u] A0;
        K_inv = K.inv, #matrix[r_minus_u,r_minus_u] K_inv;
        L_inv = L.inv #matrix[u,u] L_inv;
      )


      # browser()

      attempt <- 1
      stan_fit <- NULL

      while (is.null(stan_fit) && attempt <= 10) {
        stan_fit <- tryCatch(
          {if (!do_mcmc) {
            rstan::vb(
              response_envelope_stan,
              data = stan_dt,
              output_samples = n.final,
              init = list(A = A),
              tol_rel_obj = 1e-4,
              importance_resampling = FALSE
            )
          } else {
            rstan::sampling(
              response_envelope_stan,
              data = stan_dt
            )
          }}
          , error = function(e) {
            message(sprintf("Error in (Attempt %d): %s", attempt, e$message))
            return(NULL)  # Return NULL in case of error
          })
        attempt <- attempt + 1
      }

      all_samps <- rstan::extract(stan_fit)

      A.list[-burnin_iter] <- lapply(
        seq_len(n.final),
        function(ii) matrix(all_samps$A[ii, , ], r-u, u)
      )

      gamma.list[-burnin_iter] <- lapply(
        seq_len(n.final),
        function(ii) matrix(all_samps$gamma_A[ii, , ], r, u)
      )

      gamma0.list[-burnin_iter] <- lapply(
        seq_len(n.final),
        function(ii) matrix(all_samps$gamma0_A[ii, , ], r, r-u)
      )


      # all_lpd <- all_samps$lp__

    }

    # stan_fn <- rstan::vb

    # stan_fit <- stan_fn(
    #   response_envelope_stan,
    #   data = stan_dt
    # )
    #
    #
    #     stan_fit <- rstan::vb(
    #       response_envelope_stan,
    #       data = stan_dt,
    #       output_samples =
    #     )
    #
    #     all_A_samps <- rstan::extract(stan_fit)


    check_count <- 0

    for(iter in 1:n.iter) {

      # generate A, if 0 < u < r from the *marginal* posterior of A
      if (A_exists) {

        if (!is.null(A.list[[iter]])) {
          A <- A.list[[iter]]
          lpd.A.all[iter] <- lpd_A_resp_new(
            A = A,
            Y_minus_mu_shift = Y_minus_mu_shift,
            Xc = Xc,
            eta = eta,
            Omega.inv = Omega.inv,
            Omega0.inv = Omega0.inv,
            e = e,
            M.half = M.half,
            K.half.inv = K.half.inv,
            L.half.inv = L.half.inv,
            A0 = A0
          )
          gamma <- gamma.list[[iter]]
          gamma0 <- gamma0.list[[iter]]
        }
        # lpd.A.all[iter] <- c(lpd_A)
        # gamma_gamma0 <- attr(lpd_A, "gamma_gamma0")
        # gamma.list[[iter]] <- gamma <- gamma_gamma0$gamma
        # gamma0.list[[iter]] <- gamma0 <- gamma_gamma0$gamma0
      } else if (u == r) {
        gamma.list[[iter]] <- diag(1, r)
      } else if (u == 0) {
        gamma0.list[[iter]] <- diag(1, r)
      }



      # # tune A if 0 < u < r & iter <= n.tune
      # if (A_exists & autotune & iter >= tune_nterm &
      #     iter %% 5 == 0 & iter <= n.tune) {
      #
      #   if (iter <= n.tune) {
      #     tau <- autotune_param(
      #       draw_param_list = A.list[1:iter],
      #       tune_param = tau,
      #       accpt_list = accpt.A.all[1:iter],
      #       tune_nterm = tune_nterm,
      #       tune.incr = tune.incr,
      #       tune.accpt.prop.lower = tune.accpt.prop.lower,
      #       tune.accpt.prop.upper = tune.accpt.prop.upper,
      #       adjust_scale = (iter %% tune_scale_nterm) == 0
      #     )
      #   } else {
      #     tau <- autotune_SCAM(
      #       draw_param_list = A.list[1:iter],
      #       epsilon = 0.01,
      #       prev_var = attr(tau, "curr_var"),
      #       prev_mean = attr(tau, "curr_mean")
      #     )
      #
      #   }
      # }

      # generate all other parameters given A
      # these are exact samples, so drawn only after burnin
      if (iter >= n.burnin) {

        check_count <- check_count + 1

        # Generate Omega, Omega0 given A


        # if (check_count > 1) browser()

        # browser()

        # generate Omega, if u > 0
        if (A_exists) {
          Psi.tilde <- Psi + crossprod(gamma, G.tilde) %*% gamma
          Omega.list[[iter]] <- Omega <-
            rinvwish(dim = u, Phi = Psi.tilde, nu = nu + n - 1)
          Omega.inv <- solve_chol(Omega)

        } else if (u == r) {
          Psi.tilde <- Psi + G.tilde
          Omega.list[[iter]] <- Omega <-
            rinvwish(dim = u, Phi = Psi.tilde, nu = nu + n - 1)
          Omega.inv <- solve_chol(Omega)
        } else {
          Omega <- Omega.inv <- 0
        }


        # generate Omega0, if u < r
        if (A_exists) {
          Psi.tilde.0 <- Psi0 + crossprod(gamma0, YctYc) %*% gamma0
          Omega0.list[[iter]] <- Omega0 <-
            rinvwish(dim = r - u, Phi = Psi.tilde.0, nu = nu0 + n - 1)
          Omega0.inv <- solve_chol(Omega0)
        } else if (u == 0) {
          Psi.tilde.0 <- Psi0 + YctYc
          Omega0.list[[iter]] <- Omega0 <-
            rinvwish(dim = r - u, Phi = Psi.tilde.0, nu = nu0 + n - 1)
          Omega0.inv <- solve_chol(Omega0)
        } else {
          Omega0 <- Omega0.inv <- 0
        }


        # calculate Sigma and  Sigma.inv
        if (u > 0) {
          Sigma.inv.1 <- gamma %*% tcrossprod(Omega.inv, gamma)
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
        } else if (u == r) {
          Sigma.list[[iter]] <- Sigma <- Omega
        }


        # generate mu_shift, eta given Omega, Omega0 and A

        # generate mu_shift
        mu_shift.list[[iter]] <- mu_shift <- rmvnorm(mu = Y.bar, sigma = Sigma / n)
        mu_minus_Ybar <- mu_shift - Y.bar
        Y_minus_mu_shift <- Y - tcrossprod(one_n, mu_shift)


        # generate eta, if u > 0
        # e.tilde <- (crossprod(Y_minus_mu_shift, Xc) + eM) %*% XctXc_plus_M.inv
        if (u > 0) {
          eta.list[[iter]] <- eta <- rMatrixNormal(
            M = crossprod(gamma, e.tilde),
            V1 = Omega,
            V2 = XctXc_plus_M.inv
          )
        }


        # calculate beta
        if (A_exists) {
          beta.list[[iter]] <- beta <- gamma %*% eta
        } else if (u == r) {
          beta.list[[iter]] <- beta <- eta
        } else if (u == 0) {
          beta.list[[iter]] <- beta <- matrix(0, r, p)
        }



        # calculate mu
        mu.list[[iter]] <- mu <- mu_shift - beta %*% X.bar





        # calculate llik and lpd
        if (u > 0) {
          log_det_Omega <- log(det(Omega))
        } else {
          log_det_Omega <- 0
        }

        if (u < r) {
          log_det_Omega0 <- log(det(Omega0))
        } else {
          log_det_Omega0 <- 0
        }



        ## final residual after calculating beta
        # Y_cent <- Y - tcrossprod(one_n, mu_shift)
        resi <- Y_minus_mu_shift - tcrossprod(X, beta)

        llik.contri.list[[iter]] <-
          llik.contri <-
          - r * log_sqrt_2pi +
          - 0.5 * (log_det_Omega + log_det_Omega0 +
                     rowSums((resi %*% Sigma.inv.1) * resi) +
                     rowSums((Y_minus_mu_shift %*% Sigma.inv.2) * Y_minus_mu_shift))
        llik.all[iter] <- llik <- sum(llik.contri)



        if (A_exists) {

          # lpd_full <- - 0.5 * ( (n + p + nu + u + 1) * log_det_Omega +
          #                         (n + nu0 + r - u + 1) * log_det_Omega0 +
          #                         sum(Psi * Omega.inv) +
          #                         sum(Psi0 * Omega0.inv)
          # ) + c(lpd_A)
          lpd.full.all[iter] <- lpd_full <- llik - 0.5 * (
            (nu + u + p + 1) * log_det_Omega +
              sum(
                (Psi + tcrossprod((eta - crossprod(gamma, e)) %*% M.half)) *
                  Omega.inv
              ) +
              (nu0 + r - u + 1) * log_det_Omega0 +
              sum(Psi0 * Omega0.inv) +
              sum((K.half.inv %*% (A-A0) %*% L.half.inv)^2)
          )

        } else if (u == r) {
          lpd.full.all[iter] <- lpd_full <-
            llik - 0.5 * ((nu + r + p + 1) * log_det_Omega
                          + sum((Psi + tcrossprod((eta - e) %*% M.half)) * Omega.inv))
        } else if (u == 0) {
          lpd.full.all[iter] <- lpd_full <-
            llik - 0.5 * ((nu0 + r + 1) * log_det_Omega0
                          + sum(Psi0 * Omega0.inv))
        }



        # if (check_count == 1) browser()


      }


      if (show_progress & !no_pb) {
        setTkProgressBar(pb, iter,
                         label = paste0("Progress: ",
                                        round(iter/n.iter*100),
                                        "%  ",
                                        ifelse(iter <= n.burnin,
                                               "(Burn-in)",
                                               "(Sampling)")))
      }

    }


    MC <- list(mu_shift = mu_shift.list[-burnin_iter],
               mu = mu.list[-burnin_iter],
               eta = eta.list[-burnin_iter],
               beta = beta.list[-burnin_iter],
               Omega = Omega.list[-burnin_iter],
               Omega0 = Omega0.list[-burnin_iter],
               Sigma = Sigma.list[-burnin_iter],
               Sigma.inv = Sigma.inv.list[-burnin_iter],
               A = A.list[-burnin_iter],
               gamma = gamma.list[-burnin_iter],
               gamma0 = gamma0.list[-burnin_iter],
               accpt.A.all = accpt.A.all[-burnin_iter],
               lpd.A = lpd.A.all[-burnin_iter],
               accpt.A.ave = accpt.A.ave[-burnin_iter],
               lpd.full = lpd.full.all[-burnin_iter],
               llik = llik.all[-burnin_iter],
               llik.contri = llik.contri.list[-burnin_iter],
               tau = tau,
               prior_param = list())


    if (!no_pb) close(pb)

    MC

  }

  lapply_ <- function(...) {
    if (n.chains == 1 | !chains_parallel) {
      lapply(...)
    } else {
      future.apply::future_apply(...,
                                 future.seed = TRUE)
    }
  }

  all_MCs <- lapply_(1:n.chains,
                     function(j)
                       runMC(starting, chain_no = j))


  # out_scalar <- c("n.iter", "X", "Y", "Y.order", "u", "burnin")
  #

  # reshape the final list
  # each parameter, e.g. Omega is an
  # n.chain long list of n.iter.final long lists.

  vars <- setdiff(names(all_MCs[[1]]), "prior_param")

  out_vars <- lapply(vars,
                     function(x)
                     {
                       out <- lapply(all_MCs, "[[", x)
                       names(out) <- paste0("Chain_", 1:n.chains)
                       out
                     }
  )
  names(out_vars) <- vars

  tt2 <- Sys.time()

  total_time <- difftime(tt2, tt1, units = "secs") %>% as.numeric()

  out_res <- c(list("n.iter" = n.iter,
                    "Y.order" = Y.order,
                    "u" = u, "burnin" = n.burnin,
                    tune = n.tune,
                    n.chains = n.chains,
                    prior_param = all_MCs[[1]]$prior_param,
                    total_time = total_time),
               out_vars)

  class(out_res) <- c("Benvlp", "Benvlp_resp")

  out_res

}
