
# 1. MLE

r <- 3
p <- 15
m_true <- 2
set.seed(100)
all_pars_m2 <- predictor_generate_par(r, p, 2)
set.seed(100)
all_pars_m5 <- predictor_generate_par(r, p, 5)

MLE_result <-

  pblapply(c(100, 200, 500, 1000), function(n){

    datasets <- readRDS(file = paste0("/Users/ksh/Library/CloudStorage/Dropbox/Project/LaplaceVI/Simulation/predictor datasets/datasets_n",n,"_m",m_true,".rds"))

    Post_prob_u <- c()
    tmp_MSE <- c()

    for(i in c(1:100)){

      tmp_MLE <- sapply(0:15, function(m){

        tmp_MLE <- xenv(datasets[[i]]$X, datasets[[i]]$Y, u = m)
        tmp_MSE <- sum((t(tmp_MLE$beta) - all_pars$beta_tru)^2)
        tmp_BIC <- -2 * tmp_MLE$loglik + log(n) * (r+r*(r+1)/2+p+p*(p+1)/2+r*m)

        return(c(tmp_MSE,
                 tmp_BIC))
      })

      Post_prob_u <- rbind(Post_prob_u, softmax(-0.5 * tmp_MLE[2,]))
      tmp_MSE <- rbind(tmp_MSE, tmp_MLE[1,])

    }

    list("p(u|D)" = round(apply(Post_prob_u, 2, mean), 3),
         "MA_mse" = mean(apply(Post_prob_u * tmp_MSE, 1, sum)),
         "MA_sd" = mean(apply(Post_prob_u * tmp_MSE, 1, sd)))

  })

# 2. MCMC-stiefel

# parallel

library(doFuture)
library(foreach)
library(progressr)

xenv_file1  <- normalizePath("xenv_MH_Gibbs/xenv_MC_pred_gibbs.R", winslash = "/", mustWork = TRUE)
xenv_file2  <- normalizePath("xenv_CALVI/xenv_prior_init.R", winslash = "/", mustWork = TRUE)
xenv_file3  <- normalizePath("xenv_CALVI/xenv_LVI.R", winslash = "/", mustWork = TRUE)
xenv_file4  <- normalizePath("xenv_CALVI/xenv_LVI_update.R", winslash = "/", mustWork = TRUE)
xenv_file5  <- normalizePath("xenv_CALVI/xenv_generate_sim.R", winslash = "/", mustWork = TRUE)
util_file  <- normalizePath("xenv_CALVI/utility_functions.R", winslash = "/", mustWork = TRUE)

doFuture::registerDoFuture()
plan(multisession, workers = 40)

softmax_masked <- function(x, mask) {
  # mask: TRUE인 위치만 사용
  if (!any(mask)) return(rep(NA_real_, length(x)))
  z <- x[mask]
  z <- z - max(z, na.rm = TRUE)
  ez <- exp(z)
  p  <- ez / sum(ez)
  out <- numeric(length(x))
  out[mask] <- p
  out[!mask] <- 0       # 결측/무한 BIC에는 가중치 0
  out
}

n_values <- c(100, 200, 500, 1000)
m_values <- 2:4
n_datasets <- 50

handlers(global = TRUE)
handlers("progress")

for (ni in seq_along(n_values)) {
  n <- n_values[ni]
  prefix <- "stiefel"
  datasets <- readRDS(
    paste0("/home/jupyter-tmdgus4970/2025_CALVI/simulation/datasets_n",
           n, "_m", 2, ".rds")
  )
  
  progressr::with_progress({
    
    p <- progressor(along = 11:n_datasets)
    
    res_list <- foreach(
      i = 11:n_datasets,
      .options.future = list(seed = TRUE),
      .export = c(
        "m_values","softmax_masked","xenv_file1","xenv_file2","xenv_file3",
        "xenv_file4","xenv_file5","util_file","n","all_pars_m2"
      ),
      .packages = c("purrr", "progressr")  # 여기에 의존 패키지 추가 (mvtnorm, matrixStats 등)
    ) %dopar% {
      # 워커 초기화
      source(xenv_file1)
      source(xenv_file2)
      source(xenv_file3)
      source(xenv_file4)
      source(xenv_file5)
      source(util_file)
      
      X <- datasets[[i]]$X
      Y <- datasets[[i]]$Y
      p <- ncol(X); r <- ncol(Y)
      
      bic_list  <- rep(NA_real_, length(m_values))
      mse_list  <- rep(NA_real_, length(m_values))
      time_list <- rep(NA_real_, length(m_values))
      
      for (j in seq_along(m_values)) {
        m <- m_values[j]
        fit <- try(Benvlp_MC_pred_gibbs_devel(
          X = X, Y = Y, m = m,
          n.iter = 10000, burnin.prop = 0.5, n.chains = 1,
          A_proposal = "rwmh",
          jacobian_only_gamma = FALSE,
          compute_llik = TRUE
        ), silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          bic_list[j]  <- NA_real_
          mse_list[j]  <- NA_real_
          time_list[j] <- NA_real_
          next
        }
        
        out_dir <- "/home/jupyter-tmdgus4970/2025 CALVI/simulation result/xenv stiefel m2 0225"
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(fit, file.path(out_dir, sprintf("%s_n%d_i%d_m%d.rds", prefix, n, i, m)))
       
      } 
    }
  })
}

n_values <- c(100, 200, 500, 1000)
m_values <- 5:7
n_datasets <- 50

for (ni in seq_along(n_values)) {
  n <- n_values[ni]
  prefix <- "stiefel"
  datasets <- readRDS(
    paste0("/home/jupyter-tmdgus4970/2025_CALVI/simulation/datasets_n",
           n, "_m", 5, ".rds")
  )
  
  progressr::with_progress({
    
    p <- progressor(along = 11:n_datasets)
    
    res_list <- foreach(
      i = 11:n_datasets,
      .options.future = list(seed = TRUE),
      .export = c(
        "m_values","softmax_masked","xenv_file1","xenv_file2","xenv_file3",
        "xenv_file4","xenv_file5","util_file","n","all_pars_m5"
      ),
      .packages = c("purrr", "progressr")  # 여기에 의존 패키지 추가 (mvtnorm, matrixStats 등)
    ) %dopar% {
      # 워커 초기화
      source(xenv_file1)
      source(xenv_file2)
      source(xenv_file3)
      source(xenv_file4)
      source(xenv_file5)
      source(util_file)
      
      X <- datasets[[i]]$X
      Y <- datasets[[i]]$Y
      p <- ncol(X); r <- ncol(Y)
      
      bic_list  <- rep(NA_real_, length(m_values))
      mse_list  <- rep(NA_real_, length(m_values))
      time_list <- rep(NA_real_, length(m_values))
      
      for (j in seq_along(m_values)) {
        m <- m_values[j]
        fit <- try(Benvlp_MC_pred_gibbs_devel(
          X = X, Y = Y, m = m,
          n.iter = 10000, burnin.prop = 0.5, n.chains = 1,
          A_proposal = "rwmh",
          jacobian_only_gamma = FALSE,
          compute_llik = TRUE
        ), silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          bic_list[j]  <- NA_real_
          mse_list[j]  <- NA_real_
          time_list[j] <- NA_real_
          next
        }
        
        out_dir <- "/home/jupyter-tmdgus4970/2025 CALVI/simulation result/xenv stiefel m5 0225"
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(fit, file.path(out_dir, sprintf("%s_n%d_i%d_m%d.rds", prefix, n, i, m)))
        
      } 
    }
  })
}

## ---- 결과 집계 -----

# fname -> (n,i,m)
.parse_nim <- function(bn) {
  m <- regexec("^stiefel_n(\\d+)_i(\\d+)_m(\\d+)\\.rds$", bn)
  g <- regmatches(bn, m)[[1]]
  if (length(g) == 4) c(n = as.integer(g[2]), i = as.integer(g[3]), m = as.integer(g[4])) else NULL
}

# r,p 유도: fit에 p/r 있으면 사용, 없으면 beta에서 유추
.get_rp <- function(fit, beta_est) {
  if (!is.null(fit$r) && !is.null(fit$p)) return(c(r = as.integer(fit$r), p = as.integer(fit$p)))
  if (is.matrix(beta_est)) return(c(r = nrow(beta_est), p = ncol(beta_est)))
  stop("Cannot infer (r,p)")
}

# loglik 최대 안전 추출
.max_loglik <- function(fit) {
  tryCatch({
    v <- unlist(fit$llik$Chain_1)
    if (length(v) == 0) NA_real_ else suppressWarnings(max(v, na.rm = TRUE))
  }, error = function(e) NA_real_)
}

# beta 평균
.beta_mean <- function(fit) {
  tryCatch({
    be <- Reduce("+", fit$beta$Chain_1) / length(fit$beta$Chain_1)
    if (!is.matrix(be)) stop("beta_est not matrix"); be
  }, error = function(e) NULL)
}

# ---- 메인: 저장된 RDS로부터 집계 ----
aggregate_stiefel_results <- function(
    dir,
    n_values,
    m_values = 1:14,
    n_datasets = 100,
    m_true,
    util_files = NULL,        # 예: c(xenv_file5, util_file)  (predictor_generate_par 필요시)
    seed_for_truth = 100L     # 기존 코드와 동일하게 고정 시드 사용
) {
  stopifnot(dir.exists(dir))
  if (!is.null(util_files)) {
    for (f in util_files) if (file.exists(f)) source(f)
  }
  
  files <- list.files(dir, pattern = "^stiefel_n\\d+_i\\d+_m\\d+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No RDS files matched in: ", dir)
  
  # 인덱스 테이블
  idx <- do.call(rbind, lapply(basename(files), function(b) .parse_nim(b)))
  idx <- as.data.frame(idx, stringsAsFactors = FALSE)
  idx$path <- files
  
  # 누락 맵 기록
  missing_map <- vector("list", length(n_values)); names(missing_map) <- as.character(n_values)
  
  BMA_prob_table <- setNames(vector("list", length(n_values)), as.character(n_values))
  BMA_summary_table <- data.frame(
    n = n_values, BMA_mse_mean = NA_real_, BMA_mse_sd = NA_real_, total_time = NA_real_
  )
  
  for (ni in seq_along(n_values)) {
    n <- n_values[ni]
    sub <- idx[idx$n == n, , drop = FALSE]
    
    # (i,m) 존재 여부 맵
    present <- with(sub, tapply(path, list(i, m), function(x) if (length(x)) x[1] else NA_character__))
    all_i <- as.character(1:n_datasets)
    all_m <- as.character(m_values)
    present <- present[all_i, all_m, drop = FALSE]  # 정렬/채우기
    
    # 누락 보고
    miss_where <- which(is.na(present), arr.ind = TRUE)
    if (nrow(miss_where)) {
      missing_map[[as.character(n)]] <- data.frame(
        i = as.integer(rownames(present)[miss_where[, 1]]),
        m = as.integer(colnames(present)[miss_where[, 2]])
      )
    } else {
      missing_map[[as.character(n)]] <- data.frame(i = integer(0), m = integer(0))
    }
    
    # 집계 컨테이너
    probs_mat <- matrix(NA_real_, nrow = n_datasets, ncol = length(m_values))
    colnames(probs_mat) <- paste0("m", m_values)
    bma_mse    <- rep(NA_real_, n_datasets)
    total_time <- rep(NA_real_, n_datasets)
    
    # 각 dataset i 집계
    for (ii in seq_len(n_datasets)) {
      bic_list  <- rep(NA_real_, length(m_values))
      mse_list  <- rep(NA_real_, length(m_values))
      time_list <- rep(NA_real_, length(m_values))
      
      for (jj in seq_along(m_values)) {
        m <- m_values[jj]
        fpath <- present[ii, as.character(m)]
        if (is.na(fpath) || !nzchar(fpath)) next
        
        fit <- try(readRDS(fpath), silent = TRUE)
        if (inherits(fit, "try-error")) next
        
        loglik <- .max_loglik(fit)
        beta_est <- .beta_mean(fit)
        
        if (is.null(beta_est)) next
        
        rp <- try(.get_rp(fit, beta_est), silent = TRUE)
        if (inherits(rp, "try-error")) next
        p <- rp[["r"]]; r <- rp[["p"]] # r,p 바껴야함
        
        if (is.finite(loglik)) {
          k_param <- r + r*(r+1)/2 + p + p*(p+1)/2 + r*m
          bic_list[jj] <- -2 * loglik + log(n) * k_param
        }
        
        # MSE: predictor_generate_par 사용 가능하면 동일 방식으로 truth 생성(주의: 데이터별 truth가 다르면 여길 교체)
        if (exists("predictor_generate_par")) {
          set.seed(seed_for_truth)
          truth <- predictor_generate_par(r, p, m_true)
          if (is.list(truth) && is.matrix(truth$beta_tru)) {
            mse_list[jj] <- sum((t(truth$beta_tru) - beta_est)^2)
          }
        } else {
          mse_list[jj] <- NA_real_  # truth 없으면 계산 불가
        }
        
        time_list[jj] <- suppressWarnings(if (is.numeric(fit$total_time)) fit$total_time else NA_real_)
      }
      
      valid <- is.finite(bic_list) & is.finite(mse_list)
      probs  <- softmax_masked(-0.5 * bic_list, valid)
      probs_mat[ii, ] <- probs
      bma_mse[ii]     <- if (all(is.na(probs))) NA_real_ else sum(probs * replace(mse_list, !is.finite(mse_list), 0))
      total_time[ii]  <- sum(time_list, na.rm = TRUE)
    }
    
    keep <- rowSums(is.finite(probs_mat)) > 0 & is.finite(bma_mse)
    
    if (any(keep)) {
      BMA_prob_table[[as.character(n)]] <- colMeans(probs_mat[keep, , drop = FALSE], na.rm = TRUE)
      BMA_summary_table[ni, "BMA_mse_mean"] <- mean(bma_mse[keep], na.rm = TRUE)
      BMA_summary_table[ni, "BMA_mse_sd"]   <- sd(bma_mse[keep],   na.rm = TRUE)
      BMA_summary_table[ni, "total_time"]   <- mean(total_time[keep], na.rm = TRUE)
    } else {
      BMA_prob_table[[as.character(n)]] <- rep(NA_real_, length(m_values))
    }
  }
  
  list(
    BMA_prob_table = BMA_prob_table,
    BMA_summary_table = BMA_summary_table,
    missing = missing_map
  )
}

out_dir_m2 <- "/home/jupyter-tmdgus4970/2025 CALVI/simulation result/xenv stiefel m2 0225"
res_m2 <- aggregate_stiefel_results(
  dir = out_dir_m2,
  n_values = c(100, 200, 500, 1000),
  m_values = 1:14,
  n_datasets = 10,
  m_true = 2,
  util_files = c(xenv_file5, util_file)  # predictor_generate_par가 여기 들어있다면
)

out_dir_m5 <- "/home/jupyter-tmdgus4970/2025 CALVI/simulation result/xenv stiefel m5 0225"
res_m5 <- aggregate_stiefel_results(
  dir = out_dir_m5,
  n_values = c(100, 200, 500, 1000),
  m_values = 1:14,
  n_datasets = 10,
  m_true = 5,
  util_files = c(xenv_file5, util_file)  # predictor_generate_par가 여기 들어있다면
)

BMA_prob_table    <- res_m2$BMA_prob_table
BMA_summary_table <- res_m2$BMA_summary_table
missing_grid      <- res_m2$missing

BMA_prob_table$`100` %>% round(3)
BMA_prob_table$`200` %>% round(3)
BMA_prob_table$`500` %>% round(3)
BMA_prob_table$`1000` %>% round(3)

BMA_summary_table

BMA_prob_table    <- res_m5$BMA_prob_table
BMA_summary_table <- res_m5$BMA_summary_table
missing_grid      <- res_m5$missing

BMA_prob_table$`100` %>% round(3)
BMA_prob_table$`200` %>% round(3)
BMA_prob_table$`500` %>% round(3)
BMA_prob_table$`1000` %>% round(3)

BMA_summary_table