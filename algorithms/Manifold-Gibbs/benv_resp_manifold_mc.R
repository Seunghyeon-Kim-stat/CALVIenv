#' Main function to generate the Markov chain
#' @param Y: n X r Numeric matrix containing  all the observed valued of the response variables. n: number of observations, r: number of response variables.
#' @param X: n X p Numeric matrix containing  the observed values of the predictors/covariates. n: number of observations, p: number of covariates/ predictor variables.
#' @param u: Integer variable u is the dimension  envelope subspace. 0< u< r where  r: number of response variables (number of columns in Y).
#' @param mc_length: Number of samples to be sampled using the MCMC algorithm, length of the markov chain.
#' @param prior: An appropriate customized Benvlp_prior object or character specifying the standard automatic prior choices. Two charcter input options are 1) "Empirical" which stands for emperical prior that uses data to construct a proper prior (Informative prior) and 2) "Uniform" which specifies the uniform improper flat prior (non inforative or weakly  informative prior).
#' @param name: Name of the list variable containing all the variables of the Markov chain.
#' @param backupMc : logical input. If set to True, the generated Markov chain will be saved in a folder name output
#' @param central_space_sampling_type :  takes values in c("ALLPOSSIBLEPAIR","RANDOMPAIR","ALLCOLUMN1", "ALLCOLUMN") or an integer between 1 to r choose 2. r: number of responses. Default is "ALLPOSSIBLEPAIR"
#'                                        ALLPOSSIBLEPAIR tkes maximum time while the corrosponding MCMC algorithm have beter mixing.
#'
#' @examples
#' \dontrun{
#' library(future)
#' plan(multiprocess)
#'
#' set.seed(1)
#' MC <- Benvlp_resp_manifold(
#'   X = wheatprotein[, 4:7],
#'   Y = wheatprotein[, 1:3],
#'   u = 1, mc_length = 1000
#' )
#' }
#'

#' @export
Benvlp_resp_manifold <- function(X, Y, u,
                                 n.iter = 100,
                                 save_summary_only = TRUE,
                                 n.chains = 1,
                                 prior = "Uniform",
                                 burnin.prop = 0.3,
                                 central_space_sampling_type = "allpossiblepair",
                                 Max_rejection = 20,
                                 partition = 1000,
                                 name = "DefaultName",
                                 backupMc = FALSE,
                                 show_progress = TRUE,
                                 chains_parallel = FALSE,
                                 ...) {

  tt1 <- Sys.time()
  mc_length <- n.iter
  no_of_chains <- n.chains


  ######################################################################################## #####################################
  ##########################                                              ################ #####################################
  ########################### validation_of_response_covariate_input(X,Y)  ################ #####################################
  ##########################                                              ################ #####################################
  ######################################################################################## #####################################

  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.numeric(Y)) {
    stop(" Y:The matrix containing all the responses, must be a numeric. ")
  }
  if (!is.numeric(X)) {
    stop(" X:The matrix containing all the covriates, must be a numeric. ")
  }

  ########################## -------------------##########################################################################
  ########################### duplication Check ##########################################################################
  ########################## -------------------##########################################################################

  if (sum(duplicated(cbind(Y), MARGIN = 2)) > 0) stop("There maybe duplicated columns in  Y, the matrix containing the observed values of all the response variables. Same variable/column seems to be used more than once in Y")
  if (sum(duplicated(cbind(X), MARGIN = 2)) > 0) stop("There maybe duplicated columns in  X, the matrix containing the observed values of all the covariates/predictors. Same variable/column seems to be used more than once in X. ")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) stop("Some responses also appear in the predictors. i.e. one multiple columns in X is also columns of Y. Y:matrix containing the observed values of all the response variables, X:matrix containing the observed values of all the covariates/predictors. ")

  ########################## ----------------##########################################################################
  ############################ dimension check ########################################################################
  ########################## ----------------##########################################################################
  if (is.vector(X)) {
    X <- as.matrix(X)
  }
  if (length(dim(Y)) != 2) {
    stop("please provide Two dimensional array/ matrix input for Y: resposne.")
  }
  if (length(dim(X)) != 2) {
    stop("please provide Two dimensional array/ matrix input for X: the covariates.")
  }

  ########################## ----------------##########################################################################
  ########################## continuty check ##########################################################################
  ########################## ----------------##########################################################################

  check_continuty(Y, cut_off = .1)

  U <- center_data_matrix(as.matrix(Y))

  Y.orig <- Y
  X.orig <- X
  X <- center_data_matrix(as.matrix(X))
  X.bar <- attr(X, "scaled:center")


  n.burnin <- ceiling(n.iter * burnin.prop)
  burnin_iter <- seq_len(n.burnin)



  if (dim(U)[2] < 2) {
    print("There is only One response variable. The Envelope regression is equivalant the standard linear regression.")
  }
  if (dim(X)[1] != dim(Y)[1]) {
    stop(" Number of obervations in response Y and covariates X are not same.  ")
  }


  ######################################################################################## #####################################
  ########################### --------------------------------------####################### #####################################
  ########################### validation_of_other integer variables ###################### #####################################
  ########################## --------------------------------------####################### #####################################

  if (!is.numeric(no_of_chains)) {
    stop(" no_of_chains must be positive integer.")
  }
  if (length(no_of_chains) > 1) {
    stop("no_of_chains must be integer.")
  }
  no_of_chains <- as.integer(no_of_chains)
  if (no_of_chains < 1) {
    stop("no_of_chains must be poitive integer.")
  }

  ###########################################################################

  if (!is.numeric(mc_length)) {
    stop(" mc_length must be positive integer.")
  }
  if (length(mc_length) > 1) {
    stop("mc_length must be integer.")
  }
  mc_length <- as.integer(mc_length)
  if (mc_length < 1) {
    stop("mc_length must be poitive integer.")
  }

  #############################################################################

  if (!is.logical(show_progress)) {
    stop("show_progress must be logical TRUE or FALSE.")
  }


  ########################## -----------------------------------############################## #####################################
  ########################### Identification of A few variables ################ ############# #####################################
  ########################## ----------------------------------############################## #####################################


  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- dim(U)[2]
  Y.bar <- colMeans(Y)
  # P_X=X%*%(solve(t(X)%*%X))%*%t(X); #n=dim(P_X)[1] ;  #I=diag(n)

  ########################################################################################### ######################################
  ############################ -----------------------------------------------################ ######################################
  ###########################  A few more validation that depends on n, p r  ################ ######################################
  ########################### ------------------------------------------------################ ######################################
  ########################################################################################### ######################################
  central_space_sampling_type <- validator_central_space_sampling_type(central_space_sampling_type = central_space_sampling_type, number_of_response = r)

  ######################################## validatin on u ########################################################

  if (!is.numeric(u)) {
    stop(" u must be positive integer.")
  }
  if (length(u) > 1) {
    stop("u must be integer.")
  }
  u <- as.integer(u)
  if (u > r) {
    stop("u: the dimension of the response envelope space can not exceed the number of response variables. ")
  }
  if (u < 0) {
    stop("u: the dimension of the response envelope space can not be negative. ")
  }


  ########################## ----------------########################################################## ##############################
  ########################## Starting value ########################################################## ##############################
  ########################## ----------------########################################################## ##############################

  lst <- startingValue_genU(U, X, u)

  gamma <- lst[["gamma"]]
  gamma0 <- lst[["gamma0"]]
  omega <- lst[["omega"]]
  omega0 <- lst[["omega0"]]
  O <- cbind(gamma, gamma0)


  nontrivial <- (u < r) & (u > 1)

  ########################## -------------------------------------------################### ############################################
  ########################## prior selection tools (depends o starting)################### ############################################
  ########################## -------------------------------------------################### ############################################

  if (!is.element(class(prior), c("Benvlp_prior", "character"))) {
    prior <- "Uniform"
    print(" The argument for Prior is ignored as it is invalid. The default argument prior='Uniform' is used to run the sampler.")
  }


  if (is.character(prior)) {
    prior <- toupper(prior)
    if (!is.element((prior), c("EMPIRICAL", "UNIFORM"))) {
      prior <- "UNIFORM"
      print(" The argument for Prior is ignored as it is invalid. The default argument prior='Uniform' is used to run the sampler.")
    }
    if (prior == "EB") {
      prior_obj <- Benvlp_prior_empirical(U, X, u)
      Prior_name <- "Empirical"
    }
    if (prior == "UNIFORM") {
      prior_obj <- Benvlp_prior_uniform(number_of_responses = r, number_of_covariates = p, response_envelope_dim = u)
      Prior_name <- "Uniform"
    }
  }

  if (class(prior) == "Benvlp_prior") {
    prior_obj <- validation_of_Benvlp_prior_object(prior)
    Prior_name <- "Customized"
  }

  e_prior <- prior_obj$eta_hyper_e
  C_prior <- prior_obj$eta_hyper_C
  G_prior <- prior_obj$O_hyper_G
  D_prior <- diag(prior_obj$O_hyper_D) #    c(omega,omega0) # actually value of D_prior doesnot matter as G=0 and all we need is G/D[i]
  lambda_prior <- prior_obj$omega_hyper_rate
  lambda0_prior <- prior_obj$omega0_hyper_rate
  alpha_DF_prior <- prior_obj$omega_hyper_shape
  alpha0_DF_prior <- prior_obj$omega0_hyper_shape


  ############################################################################################################### #######################
  ############################################### ---------------------######################## ################## #######################
  ##########################                       preparation for MC                                   ######### #######################
  ##########################  computation of statistic that would remain constant through out the chain ######### #######################
  ##########################                                                                            ######### #######################
  ############################################################################################################### #######################


  # U=U-t(replicate(n=dim(U)[1],expr=(apply(U,2,"mean") ) ))
  eta_Sigma2 <- solve(t(X) %*% X + C_prior)
  e_tilde <- (t(U) %*% X + e_prior %*% C_prior) %*% eta_Sigma2 # solve(t(X)%*% X+C_prior) eta_Sigma2 is varianc for the posterior of Eta
  G_tilde <- t(U) %*% U + e_prior %*% C_prior %*% t(e_prior) - e_tilde %*% (t(X) %*% X + C_prior) %*% t(e_tilde)


  burnin <- ceiling(mc_length * burnin.prop)
  burnin_iter <- seq_len(burnin)
  one_n <- rep(1, n)


  ############### ---------------------------######################################################################## #######################
  ############## Storing of the MC objects ######################################################################## ########################
  ############### --------------------------######################################################################### #######################


  eta <- beta <- list()
  omega_all <- omega
  omega0_all <- omega0
  O_all <- list(O)
  mu.list <- mu_unshift.list <- vector("list", n.iter)
  Sigma_inv_list <- vector("list", n.iter)
  llik.contri.list <- vector("list", n.iter)
  gamma.list <- gamma0.list <- vector("list", n.iter)
  Omega_mat_all <- Omega0_mat_all <- vector("list", n.iter)
  lpd.full.all <- lpd.A.all <- llik.all <- rep(0, mc_length)

  #################################################################### ########################################### ######################
  #################################################################### ########################################### ######################
  ##########################                  ####################### ########################################### #######################
  ########################## Mc main function ######################## ########################################### ######################
  ##########################                  ######################## ########################################### ######################
  #################################################################### ########################################### ######################



  n.burnin <- ceiling(n.iter * burnin.prop)
  burnin_iter <- seq_len(n.burnin)



  runMC <- function(starting = NULL, chain_no = 1) {
    #####################################################################################################
    ########################## progressbar start ########################################################
    #####################################################################################################
    
    if (show_progress) {
      pb <- tryCatch(
        tkProgressBar(
          title = paste("Chain", chain_no),
          label = "Progress: 0%",
          min = 0, max = mc_length
        ),
        error = function(e) e,
        warning = function(w) w
      )

      if (is(pb, "error") | is(pb, "warning")) {
        show_progress <- FALSE
        cat(paste("\'tcltk\' could not be loaded. \'show_progress\' set to FALSE."))
      }
    }


    ############################################################################################# ############
    ############################################################################################# ############
    ############################################################################################# ############
    ############################################################################################# ############
    runif(chain_no)
    
    # ---- summary accumulators ----
    sum_beta <- matrix(0, r, p)
    sum_beta2 <- matrix(0, r, p)
    
    sum_mu <- rep(0, r)
    sum_mu2 <- rep(0, r)
    
    sign_fix <- (u == 1)   # u=1일 때만 자동 적용 (원하면 옵션으로 빼도 됨)
    ref_beta_vec <- NULL
    
    # eta는 u가 바뀌므로 저장할지 선택 (필요 없으면 빼도 됨)
    if (u > 0) {
      sum_eta <- matrix(0, u, p)
      sum_eta2 <- matrix(0, u, p)
    }
    
    # loglik max (BIC용)
    llik_max <- -Inf
    
    # sample count (burn-in 제외)
    s_count <- 0L
    
    
    i <- 1
    while (i <= mc_length) {
      ################################################################### ######################
      if (u == 0) {
        # gamma=as.matrix(O[,1:u])
        gamma0.list[[i]] <- gamma0 <- as.matrix(O[, (u + 1):r])
        # A=t(gamma)%*% t(U)%*%(I-P_X)%*%U %*%gamma
        # B=t(gamma0)%*%t(U)%*%U%*%gamma0
        B <- t(gamma0) %*% (t(U) %*% U) %*% gamma0 + 2 * lambda0_prior

        omega <- NULL
        Sigma_inv_1 <- matrix(0, r, r)

        omega0 <- UpdateOmega(omega0, diag(B) / 2, alpha = (n + 2 * alpha0_DF_prior - 1) / 2, Max_rejection = Max_rejection)

        Omega0 <- diag(omega0, nrow = (r - u))
        Omega0_mat_all[[i]] <- Matrix(Omega0, sparse = TRUE)
        Omega0_inv <- diag(1 / omega0, nrow = (r - u))
        Sigma_inv_0 <- gamma0 %*% tcrossprod(Omega0_inv, gamma0)

        Sigma_inv_list[[i]] <- Sigma_inv <- Sigma_inv_0

        Sigma <- gamma0 %*% tcrossprod(Omega0, gamma0)
      }


      ################################################################# ########################
      if ((0 < u) * (u < r)) {
        gamma.list[[i]] <- gamma <- as.matrix(O[, 1:u])
        gamma0.list[[i]] <- gamma0 <- as.matrix(O[, (u + 1):r])
        # A=t(gamma)%*% t(U)%*%(I-P_X)%*%U %*%gamma
        A <- t(gamma) %*% G_tilde %*% gamma #+ 2*lambda_prior
        B <- t(gamma0) %*% (t(U) %*% U) %*% gamma0 #+ 2*lambda0_prior

        omega <- UpdateOmega(omega, (diag(A) / 2 + lambda_prior), alpha = (n + 2 * alpha_DF_prior - 1) / 2, Max_rejection = Max_rejection) ### pass the diagonal
        omega0 <- UpdateOmega(omega0, (diag(B) / 2 + lambda0_prior), alpha = (n + 2 * alpha0_DF_prior - 1) / 2, Max_rejection = Max_rejection)



        Omega <- diag(omega, nrow = (u))
        Omega_mat_all[[i]] <- Matrix(Omega, sparse = TRUE)

        Omega_inv <- diag(1 / omega, nrow = (u))
        Sigma_inv_1 <- gamma %*% tcrossprod(Omega_inv, gamma)

        Omega0 <- diag(omega0, nrow = (r - u))
        Omega0_mat_all[[i]] <- Matrix(Omega0, sparse = TRUE)

        Omega0_inv <- diag(1 / omega0, nrow = (r - u))
        Sigma_inv_0 <- gamma0 %*% tcrossprod(Omega0_inv, gamma0)


        Sigma <- gamma %*% tcrossprod(Omega, gamma) + gamma0 %*% tcrossprod(Omega0, gamma0)
        Sigma_inv_list[[i]] <- Sigma_inv <- Sigma_inv_1 + Sigma_inv_0
      }

      #################################################################### #####################
      if (u == r) {
        gamma.list[[i]] <- gamma <- as.matrix(O[, 1:u])
        # gamma0=as.matrix(O[,(u+1):r]) ###A=t(gamma)%*% t(U)%*%(I-P_X)%*%U %*%gamma
        A <- t(gamma) %*% G_tilde %*% gamma + 2 * lambda_prior
        # B=t(gamma0)%*%t(U)%*%U%*%gamma0

        omega <- UpdateOmega(omega, diag(A) / 2, alpha = (n + 2 * alpha_DF_prior - 1) / 2, Max_rejection = Max_rejection) ### pass the diagonal

        Omega <- diag(omega, nrow = (u))
        Omega_inv <- diag(1 / omega, nrow = (u))
        Omega0_mat_all[[i]] <- Matrix(Omega, sparse = TRUE)

        Sigma_inv_1 <- gamma %*% tcrossprod(Omega_inv, gamma)


        omega0 <- NULL
        Sigma_inv_0 <- matrix(0, r, r)

        Sigma <- gamma %*% tcrossprod(Omega, gamma)
        Sigma_inv_list[[i]] <- Sigma_inv <- Sigma_inv_1
      }

      ############################################################################################### #################################
      ############################################################################################### #################################
      O <- UpdateO(O = O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U = U, omega = omega, omega0 = omega0, allstep = central_space_sampling_type, discrete = F, MAX_rejection = Max_rejection, partition = partition)
      O_all[[i]] <- O



      ## gamma=O[ , 1:u];  ##gamma0= O[ , (u+1):r]
      ############################################################################ #####################
      if (u == 0) {
        beta[[i]] <- matrix(0, r, p)
      }

      ########################################################################### #######################
      if (u > 0) {
        gamma <- O[, 1:u]
        # e_tilde=(  t(U) %*% X+ e_prior%*% C_prior  ) %*% Sigma2
        eta_Mean <- t(gamma) %*% e_tilde # check this
        eta_Sigma1 <- Omega #   diag(omega,nrow=length(omega));
        eta[[i]] <- rMatrixNormal(eta_Mean, eta_Sigma1, eta_Sigma2)

        # eta[[i]]=SampleEta(X, U, omega, O[,1:u],C_prior=C_prior,e_prior=e_prior)

        beta[[i]] <- as.matrix(O[, 1:u]) %*% eta[[i]]
      }
      ########################################################################### ########################
      ########################################################################### ########################


      mu.list[[i]] <- mu <- (rmvnorm(n = 1, mu = Y.bar, sigma = Sigma / n))


      ######################################################################################## ######################### #################
      ##################                ####################################################### ######################### #################
      ##################  storing Omega ###################################################### ######################### #################
      ##################                ###################################################### ######################### #################
      ######################################################################################## ######################### #################
      if (i == 1) {
        omega_all <- omega
        omega0_all <- omega0_all
      }

      if (i > 1) {
        omega_all <- cbind(omega_all, omega)
        omega0_all <- cbind(omega0_all, omega0)
      }




      # calculate mu_unshift
      mu_unshift.list[[i]] <- mu_unshift <- mu - beta[[i]] %*% X.bar




      ######################################################################################## ########################## ################
      ###################                                                     ################ ########################## ################
      ##################  Calculating likelihood and posterior contributions ################# ########################## ################
      ###################                                                     ################ ########################## ################
      ######################################################################################## ########################## ################


      log_det_Omega <- ifelse(u > 0, sum(log(omega)), 0)
      log_det_Omega0 <- ifelse(u < r, sum(log(omega0)), 0)


      ######################################################## ########################
      ## final residual after calculating beta
      Y_cent <- Y - tcrossprod(one_n, (mu)) #
      # Y_cent=U
      resi <- Y_cent - tcrossprod(X, beta[[i]])
      ######################################################## ########################



      ################# --------------------------############# ########################
      ################# Likelihood contributions ############# ########################
      ################# --------------------------############# ########################
      llik.contri.list[[i]] <-
        llik.contri <- -0.5 * (log_det_Omega + log_det_Omega0 +
          rowSums((resi %*% Sigma_inv_1) * resi) +
          rowSums((Y_cent %*% Sigma_inv_0) * Y_cent))

      llik.all[i] <- llik <- sum(llik.contri)
      # llik.all=-.5*( tr(   (U%*%(gamma)-X%*%t(eta))  %*%  Omega_inv  %*%   t(U%*%(gamma)-X%*%t(eta))    )
      #  +   tr( U%*%(gamma0)  %*% Omega0_inv  %*% t(U%*% gamma0) )
      #   +n*(   sum(log(theta$omega))   + sum(log(theta$omega0))))


      # ---- after computing beta[[i]], mu, llik ----
      if (i > burnin) {
        s_count <- s_count + 1L
        
        # beta
        b <- beta[[i]]
        
        if (sign_fix) {
          bv <- c(b)
          if (is.null(ref_beta_vec)) {
            ref_beta_vec <- bv
          } else {
            if (sum(bv * ref_beta_vec) < 0) {
              b  <- -b
              bv <- -bv
            }
            # reference를 “running average”로 업데이트(안정적)
            ref_beta_vec <- ref_beta_vec + bv
          }
        }
        
        sum_beta  <- sum_beta  + b
        sum_beta2 <- sum_beta2 + b^2
        
        # mu (unshifted mu: mu_unshift)
        m <- mu_unshift
        sum_mu  <- sum_mu  + m
        sum_mu2 <- sum_mu2 + m^2
        
        # eta (optional)
        if (u > 0) {
          e <- eta[[i]]
          sum_eta  <- sum_eta  + e
          sum_eta2 <- sum_eta2 + e^2
        }
        
        # loglik max
        if (is.finite(llik) && llik > llik_max) llik_max <- llik
      }
      


      ########## -----------------------------######################## ########################
      ######### posterior density evaluation ######################## ########################
      ########## -----------------------------######################## ########################

      if (Prior_name == "Uniform") {
        total_log_prior_contribution <- 0
      }

      if (Prior_name != "Uniform") {
        log_prior_contribution_omega <- ifelse(u > 0, -sum((alpha_DF_prior + p / 2 + 1) * omega) - sum(lambda_prior * (1 / omega)), 0)
        log_prior_contribution_omega0 <- ifelse(u < r, -sum((alpha0_DF_prior + 1) * omega0) - sum(lambda0_prior * (1 / omega0)), 0)

        log_prior_contribution_O <- -.5 * (tr(diag(1 / D_prior, nrow = r) %*% t(O) %*% G_prior %*% O))

        if (u == 0) {
          log_prior_contribution_eta <- 0
        }
        if (u > 0) {
          exp_resid_eta <- eta[[i]] - t(as.matrix(gamma)) %*% e_prior
          log_prior_contribution_eta <- -.5 * (tr(Omega_inv %*% (exp_resid_eta) %*% C_prior %*% t(exp_resid_eta)))
        }
        total_log_prior_contribution <- log_prior_contribution_omega +
          log_prior_contribution_omega0 +
          log_prior_contribution_O +
          log_prior_contribution_eta
      }

      lpd.full.all[i] <- log_posterior <- llik + total_log_prior_contribution







      ###################################################################### ###################
      ###########################                        ################### ###################
      ##########################   progress bar updation ################### ###################
      ##########################                         ################### ###################
      ###################################################################### ###################
      i <- i + 1
      if (show_progress) {
        setTkProgressBar(pb, i,
          label = paste0(
            "Progress: ",
            round(i / mc_length * 100), "%  ",
            ifelse(i <= burnin, "(Burn-in)", "(Sampling)")
          )
        )
      }
    } #### end of the Markov chain While loop
    ###################################################################################### ########################
    ##########################                ############################################ ########################
    ##########################   End of While ############################################ ########################
    ##########################                ############################################ ########################
    ###################################################################################### ########################

    if (show_progress) {
      close(pb)
    }

    ###################################################################################### ########################
    ######################## preparation for the final output ############################ ########################
    ###################################################################################### ########################

# 
# 
# 
#     allList <- list()
#     allList$O_all <- O_all
#     allList$omega_all <- omega_all
#     allList$Omega <- Omega_mat_all
#     allList$omega0_all <- omega0_all
#     allList$Omega0 <- Omega0_mat_all
#     allList$Sigma.inv <- Sigma_inv_list
#     allList$u <- u
#     allList$eta <- eta
#     allList$beta <- beta
#     allList$prior <- prior
#     allList$X <- X
#     allList$U <- U
#     allList$mu <- mu_unshift.list
#     allList$mu_shift <- mu.list
#     allList$iterationLength <- mc_length
#     allList$burnIn <- 0
#     allList$llik_contri_list <- llik.contri.list
#     allList$llik <- llik.all
#     allList$lpd <- lpd.full.all
# 
#     if (nontrivial) {
#       allList$gamma <- 1
#     }
# 
#     # Sigma
# 
# 
#     if (backupMc) {
#       if (!file.exists("output")) {
#         dir.create("output/")
#       }
#       save(allList, file = paste("output/run_", Prior_name, name, "u_equals_", u, ".RData", sep = ""))
#     }
# 
# 
#     return(allList)
    
    # ---- posterior summaries (burn-in 제외) ----
    beta_mean <- sum_beta / max(1L, s_count)
    beta_sd   <- sqrt(pmax(0, (sum_beta2 / max(1L, s_count)) - beta_mean^2))
    
    mu_mean <- sum_mu / max(1L, s_count)
    mu_sd   <- sqrt(pmax(0, (sum_mu2 / max(1L, s_count)) - mu_mean^2))
    
    if (u > 0) {
      eta_mean <- sum_eta / max(1L, s_count)
      eta_sd   <- sqrt(pmax(0, (sum_eta2 / max(1L, s_count)) - eta_mean^2))
    }
    
    allList <- list()
    allList$u <- u
    allList$iterationLength <- mc_length
    allList$burnIn <- burnin
    allList$prior <- prior
    
    # BIC용: max loglik 저장
    allList$llik_max <- llik_max
    
    # 시간 비교용: 필요하면 아래에서 total_time에 넣고, 여기서는 생략 가능
    
    if (save_summary_only) {
      # summary only
      allList$beta_mean <- beta_mean
      allList$beta_sd   <- beta_sd
      allList$mu_mean   <- mu_mean
      allList$mu_sd     <- mu_sd
      
      if (u > 0) {
        allList$eta_mean <- eta_mean
        allList$eta_sd   <- eta_sd
      }
      
      # (선택) O 마지막 값만 남기기
      allList$O_last <- O
      
    } else {
      # 기존 full 저장 (필요하면)
      allList$O_all <- O_all
      allList$omega_all <- omega_all
      allList$Omega <- Omega_mat_all
      allList$omega0_all <- omega0_all
      allList$Omega0 <- Omega0_mat_all
      allList$Sigma.inv <- Sigma_inv_list
      allList$eta <- eta
      allList$beta <- beta
      allList$mu <- mu_unshift.list
      allList$mu_shift <- mu.list
      allList$llik_contri_list <- llik.contri.list
      allList$llik <- llik.all
      allList$lpd <- lpd.full.all
    }
    
    return(allList)
    
  }

  lapply_ <- function(...) {
    if (no_of_chains == 1 | !chains_parallel) {
      lapply(...)
    } else {
      future.apply::future_apply(...,
        future.seed = TRUE
      )
    }
  }



  all_MCs <- lapply_(
    1:no_of_chains,
    function(j) {
      runMC(starting, chain_no = j)
    }
  )



  vars <- setdiff(names(all_MCs[[1]]), c("X", "U", "u", "bunIn"))

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
  total_time <- difftime(tt2, tt1, units = "secs") %>% as.numeric()
  out_res <- c(
    list(
      "n.iter" = n.iter,
      "X" = X.orig,
      "Y" = Y.orig,
      "Y.order" = 1:ncol(Y),
      "u" = u,
      "burnin" = n.burnin,
      tune = 0,
      n.chains = n.chains,
      prior_param = all_MCs[[1]]$prior,
      total_time = total_time
    ),
    out_vars
  )

  #  all_MCs <- future.apply::future_lapply(1:no_of_chains,
  #                                     function(j){runMC(chain_no = j)},
  #                                    future.seed = TRUE)


  class(out_res) <- c("Benvlp", "Benvlp_resp")
  return(out_res)
}


##################################################################################################### ############
##################################################################################################### ############
##########################                                             ############################### ###########
##########################   End of main function Benvlp_resp_manifold() ############################## ############
##########################                                             ############################## ############
##################################################################################################### ############
##################################################################################################### ############






# "beta" = beta.all, "SigYcX" = SigYcX.all,


# MC <- list("n.iter" = n.iter, "X" = X, "Y" = Y, "X.order" = X.order, "m" = m,
#  "muX" = muX.all, "muY" = muY.all, "eta" = eta.all, "SigYcX" = SigYcX.all,
#  "Omega" = Omega.all, "Omega0" = Omega0.all, "beta" = beta.all,
# "lpd.full" = lpd.full.all, "llik" = llik.all)



##################################################################################################### ############
##################################################################################################### ############
#############################                               ######################################### ############
############################# Individual Validator function ############### ######################### ############
#############################                               ######################################### ############
##################################################################################################### ############
##################################################################################################### ############
# validator
validator_central_space_sampling_type <- function(central_space_sampling_type, number_of_response) {
  if (!is.element(class(central_space_sampling_type), c("numeric", "character"))) {
    central_space_sampling_type <- "ALLPOSSIBLEPAIR"
    print(" The argument for central_space_sampling_type is ignored as it is invalid. The default argument allstep='ALLPOSSIBLEPAIR' is used to run the sampler.")
  }


  if (is.numeric(central_space_sampling_type)) {
    if (length(central_space_sampling_type) > 1) {
      stop("central_space_sampling_type: Must be a scaler integer greater than equal to 1. Input type is not scaler.")
    }
    central_space_sampling_type <- as.integer(central_space_sampling_type)
    if (central_space_sampling_type < 1) {
      print("central_space_sampling_type: Must be a scaler integer greater than equal to 1. The algorithm will choose a random pair of columns to update. ")
    }
    central_space_sampling_type <- max(central_space_sampling_type, 1)
    total_number_of_pairs <- number_of_response * (number_of_response - 1) / 2
    if (central_space_sampling_type > total_number_of_pairs) {
      central_space_sampling_type <- total_number_of_pairs
      print("the provided value for central_space_sampling_type is greater than the total number of all possible pairs of columns, r(r-1)/2, where r=number_of_response. The algorithm will use central_space_sampling_type=ALLPOSSIBLEPAIR, which is more efficient.  ")
    }
  }
  if (!is.numeric(central_space_sampling_type)) {
    if (is.character(central_space_sampling_type)) {
      if (!is.element(toupper(central_space_sampling_type), c("ALLPOSSIBLEPAIR", "RANDOMPAIR", "ALLCOLUMN1", "ALLCOLUMN"))) {
        central_space_sampling_type <- "ALLPOSSIBLEPAIR"
        print(" The argument for allstep is ignored as it is invalid. The default argument allstep='ALLPOSSIBLEPAIR' is used to run the sampler.")
      }
    }
    if (!is.character(central_space_sampling_type)) {
      central_space_sampling_type <- "ALLPOSSIBLEPAIR"
      print(" The argument for central_space_sampling_type is ignored as it is invalid. The default argument allstep='ALLPOSSIBLEPAIR' is used to run the sampler.")
    }
  }

  return(central_space_sampling_type)
}




#######

check_continuty <- function(Y, cut_off = .1) {
  unique_fractions <- apply(X = Y, FUN = function(x) {
    return(length(unique(x)) / length(x))
  }, MARGIN = 2)
  suspected_columns_indices <- which(unique_fractions <= cut_off) ########
  if (length(suspected_columns_indices) == 0) {
    check_cont <- "Continuty check Passed."
  }
  if (length(suspected_columns_indices) > 0) {
    print(paste(suspected_columns_indices, "  th column of the response variable may not be continuous."))
    print(paste(names(Y)[suspected_columns_indices], "  may not be  continuous."))
    check_cont <- "Continuty check Failed."
  }
  return(check_cont)
}
