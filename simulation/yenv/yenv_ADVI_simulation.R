library(future)
library(doFuture)
library(foreach)
library(progressr)
library(rstan)

Sys.setenv(
  OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
  VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1"
)
yenv_mhg_file <- normalizePath("yenv_MH_Gibbs/benv_MC_resp_gibbs.R", winslash = "/", mustWork = TRUE) 
yenv_advi_file <- normalizePath("yenv_ADVI/Benvlp_resp_vb.R", winslash = "/", mustWork = TRUE) 
utility_file <- normalizePath("yenv_ADVI/utility_functions.R", winslash = "/", mustWork = TRUE) 
source(yenv_mhg_file) 
source(yenv_advi_file) 
source(utility_file)

# rstan 캐시/자동저장 옵션 (권장)
rstan_options(auto_write = TRUE)
options(mc.cores = 50)

handlers(global = TRUE)
handlers("progress")

# ---- settings
n_values   <- c(100, 200, 500, 1000)
u_true     <- 5
n_datasets <- 100

ds_path <- function(n, u_true) sprintf(
  "/home/jupyter-tmdgus4970/2025_CALVI/simulation/datasets_n%d_u%d.rds",
  n, u_true
)

set.seed(100)
all_pars  <- response_generate_par(r=20, p=7, u=u_true)
beta_true <- all_pars$beta_tru

# ---- (중요) 마스터에서 먼저 컴파일 (여기서 한 번만!)
sm_master <- rstan::stan_model(
  file="/home/jupyter-tmdgus4970/2025_CALVI/yenv_ADVI/response_envelope_safe.stan"
)
assign("response_envelope_stan", sm_master, envir=.GlobalEnv)

# ---- 이제 병렬 설정
doFuture::registerDoFuture()
future::plan(multisession, workers = 40)

run_one_rep <- function(ds_i, beta_true, u_true) {
  # 워커에서도 stan_model을 "다시 호출"해서 캐시된 DSO를 로드하게 함(컴파일 X 기대)
  sm <- rstan::stan_model(
    file="/home/jupyter-tmdgus4970/2025_CALVI/yenv_ADVI/response_envelope_safe.stan"
  )
  assign("response_envelope_stan", sm, envir=.GlobalEnv)
  
  fit_advi <- Benvlp_resp_vb(
    X = ds_i$X, Y = ds_i$Y,
    u = u_true,
    n.iter = 1e4,
    burnin.prop = 1/2,
    n.chains = 1,
    scan = "full",
    init = "map",
    jacobian_only_gamma = FALSE,
    compute_llik = FALSE
  )
  
  beta_advi <- Reduce("+", fit_advi$beta$Chain_1) / length(fit_advi$beta$Chain_1)
  sum((beta_advi - beta_true)^2)
}

# ---- datasets master load
datasets_by_n <- lapply(n_values, function(n) readRDS(ds_path(n, u_true)))
names(datasets_by_n) <- as.character(n_values)

grid <- expand.grid(n = n_values, i = seq_len(n_datasets))

with_progress({
  p <- progressor(nrow(grid))
  
  ans <- foreach(
    k = seq_len(nrow(grid)),
    .combine  = rbind,
    .packages = c("rstan")   # Benvlp_resp_vb가 있는 패키지도 필요하면 여기에 추가
  ) %dopar% {
    
    n <- grid$n[k]
    i <- grid$i[k]
    
    ds_i <- datasets_by_n[[as.character(n)]][[i]]
    mse  <- run_one_rep(ds_i, beta_true, u_true)
    
    p(sprintf("n=%d i=%d", n, i))
    data.frame(n = n, i = i, mse_advi = mse)
  }
  
  ans
})

summary_df <- aggregate(mse_advi ~ n, data = ans,
                        FUN = function(x) c(mean=mean(x), sd=sd(x)))

summary_out <- data.frame(
  n        = summary_df$n,
  mean_mse = summary_df$mse_advi[, "mean"],
  sd_mse   = summary_df$mse_advi[, "sd"]
)

saveRDS(list(raw = ans, summary = summary_out),
        file = sprintf("yenv_ADVI_u%d_results.rds", u_true))

