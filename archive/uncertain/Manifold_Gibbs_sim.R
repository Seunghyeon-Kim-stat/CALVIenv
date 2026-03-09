
library(future)
library(future.apply)

## 초기화
plan(sequential)
gc()

base_dir <- "/home/jupyter-tmdgus4970/2025_CALVI"
code_dir <- file.path(base_dir, "yenv_Manifold_Gibbs")
data_dir <- file.path(base_dir, "simulation")
out_dir  <- file.path(base_dir, "simulation result", "yenv_manifold_gibbs_u2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 필요한 코드 로드
rfiles <- list.files(code_dir, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(rfiles, source))

init_worker_once <- function() {
  if (!exists(".yenv_manifold_inited", envir = .GlobalEnv, inherits = FALSE)) {
    for (f in rfiles) source(f, local = .GlobalEnv)
    assign(".yenv_manifold_inited", TRUE, envir = .GlobalEnv)
  }
}

## 설정
u_grid    <- 5
n_grid    <- c(100, 200, 500, 1000)
n_rep     <- 100
n_workers <- 30

plan(multisession, workers = n_workers)

## BLAS 과병렬 방지
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

for (n in n_grid) {
  
  datasets_path <- file.path(data_dir, sprintf("datasets_n%d_u5.rds", n))
  datasets <- readRDS(datasets_path)
  
  message("n = ", n, " | running Manifold Gibbs for u = 0:20 with ", n_workers, " workers")
  
  grid <- expand.grid(
    i = seq_len(n_rep),
    u = u_grid
  )
  
  done_files <- list.files(out_dir, pattern = sprintf("^ManifoldGibbs_n%d_", n))
  if (length(done_files) > 0) {
    done_key <- sub("^ManifoldGibbs_n[0-9]+_i([0-9]+)_u([0-9]+)\\.rds$", "\\1-\\2", done_files)
    grid_key <- paste(grid$i, grid$u, sep = "-")
    grid <- grid[!(grid_key %in% done_key), , drop = FALSE]
  }
  
  if (nrow(grid) == 0L) {
    message("  -> all (i,u) already completed for n = ", n)
    next
  }
  
  future_lapply(seq_len(nrow(grid)), function(k) {
    
    init_worker_once()
    
    i <- grid$i[k]
    u <- grid$u[k]
    
    out_path <- file.path(out_dir,
                          sprintf("ManifoldGibbs_n%d_i%d_u%d.rds", n, i, u))
    err_path <- file.path(out_dir,
                          sprintf("ManifoldGibbs_n%d_i%d_u%d.ERROR.rds", n, i, u))
    
    Sys.setenv(OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1")
    
    X_i <- datasets[[i]]$X
    Y_i <- datasets[[i]]$Y
    
    max_try <- 3L
    
    for (attempt in seq_len(max_try)) {
      
      res <- tryCatch(
        {
          Benvlp_resp_manifold(
            X = X_i, Y = Y_i, u = u,
            n.iter = 10000,
            n.chains = 1,
            prior = "Uniform",
            burnin.prop = 0.5,
            central_space_sampling_type = "allpossiblepair",
            Max_rejection = 20,
            partition = 1000,
            name = sprintf("n%d_i%d_u%d", n, i, u),
            backupMc = FALSE,
            show_progress = FALSE,
            chains_parallel = FALSE
          )
        },
        error = function(e) e
      )
      
      if (!inherits(res, "error")) {
        saveRDS(res, out_path)
        if (file.exists(err_path)) file.remove(err_path)
        return(NULL)
      }
      
      if (attempt < max_try) next
      
      saveRDS(list(
        n = n, i = i, u = u,
        message = conditionMessage(res),
        time = as.character(Sys.time())
      ), err_path)
      
      return(NULL)
    }
    
    NULL
  }, future.seed = TRUE)
  
  rm(datasets)
  gc()
}

plan(sequential)


library(dplyr)
library(purrr)

frobenius_sq <- function(A) sum(A^2)

safe_readRDS <- function(f) {
  tryCatch(readRDS(f), error = function(e) NULL)
}

## true beta 미리 정의
set.seed(100)
r <- 20
p <- 7
u_true <- 2
n_rep <- 100
u_grid <- 1:20

all_pars <- response_generate_par(r, p, u_true)
beta_true <- all_pars$beta_tru
r <- ncol(beta_true)
p <- nrow(beta_true)

summary_BMA <- purrr::map_dfr(n_grid, function(n) {
  
  message("Summarizing BMA for n = ", n)
  
  res_i <- purrr::map(seq_len(n_rep), function(i) {
    
    rows <- purrr::map(u_grid, function(u) {
      
      f <- file.path(out_dir, sprintf("ManifoldGibbs_n%d_i%d_u%d.rds", n, i, u))
      res <- safe_readRDS(f)
      if (is.null(res)) return(NULL)
      
      beta_hat <- res$beta_mean$Chain_1
      mse_u <- frobenius_sq(beta_hat - beta_true)
      
      loglik_max <- res$llik_max$Chain_1
      k_u <- r + p*u + r*(r+1)/2
      bic_u <- -2*loglik_max + log(n)*k_u
      
      tibble::tibble(
        i = i,
        u = u,
        mse = mse_u,
        bic = bic_u,
        time_u = res$total_time
      )
    })
    
    df <- dplyr::bind_rows(rows)
    if (nrow(df) == 0L) return(NULL)
    
    bic_min <- min(df$bic)
    w <- exp(-0.5 * (df$bic - bic_min))
    w <- w / sum(w)
    
    # (i,u)별 weight 테이블
    df_w <- dplyr::mutate(df, w = w)
    
    # replicate 요약(1행)
    df_rep <- tibble::tibble(
      i = i,
      BMA_mse = sum(w * df$mse),
      u_hat = df$u[which.min(df$bic)],
      time = sum(df$time_u)
    )
    
    list(df_w = df_w, df_rep = df_rep)
  })
  
  # flatten
  df_w_all  <- dplyr::bind_rows(purrr::map(res_i, "df_w"))
  df_rep_all <- dplyr::bind_rows(purrr::map(res_i, "df_rep"))
  
  if (nrow(df_rep_all) == 0L) {
    return(tibble::tibble(
      n = n,
      n_rep = 0L,
      BMA_mse_mean = NA_real_,
      BMA_mse_var  = NA_real_,
      time_mean = NA_real_,
      time_var  = NA_real_,
      P_hat_u = list(NULL),
      P_sel_u = list(NULL)
    ))
  }
  
  # u별 평균 weight: E_i[w_u]
  df_wmean <- df_w_all %>%
    dplyr::group_by(u) %>%
    dplyr::summarise(w_mean = mean(w), .groups = "drop")
  
  # (선택) selection frequency도 같이 저장 가능
  df_sel <- df_rep_all %>%
    dplyr::count(u_hat) %>%
    dplyr::mutate(p = n / sum(n))
  
  tibble::tibble(
    n = n,
    n_rep = nrow(df_rep_all),                 # <- 진짜 replicate 수
    BMA_mse_mean = mean(df_rep_all$BMA_mse),
    BMA_mse_var  = sd(df_rep_all$BMA_mse),
    time_mean = mean(df_rep_all$time),
    time_var  = sd(df_rep_all$time),
    P_hat_u = list(stats::setNames(df_wmean$w_mean, df_wmean$u)),   # 평균 weight
    P_sel_u = list(stats::setNames(df_sel$p, df_sel$u_hat))         # 선택비율(옵션)
  )
})

summary_BMA$BMA_mse_mean
summary_BMA$BMA_mse_var
summary_BMA$time_mean
summary_BMA$P_hat_u
summary_BMA
