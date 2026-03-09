library(doFuture)
library(foreach)
library(progressr)

doFuture::registerDoFuture()
parallel::detectCores()

plan(multisession, workers = 40)   # 워커 수는 네 환경에 맞게
print(plan())

rstan_options(auto_write = TRUE)
handlers(global = TRUE)
handlers("progress")

## ---- 부트스트랩 함수 ----

setwd("/home/jupyter-tmdgus4970/2025_CALVI/")

`%||%` <- function(a,b) if (!is.null(a)) a else b

source_files_CALVI <- list.files("yenv_CALVI", pattern = "\\.[Rr]$",
                                 full.names = TRUE, recursive = TRUE)

source_files_MH_Gibbs <- list.files("yenv_MH_Gibbs", pattern = "\\.[Rr]$",
                                    full.names = TRUE, recursive = TRUE)

source_files_Manifold_Gibbs <- list.files("yenv_Manifold_Gibbs", pattern = "\\.[Rr]$",
                                          full.names = TRUE, recursive = TRUE)

source_files <- c(source_files_CALVI, source_files_MH_Gibbs, source_files_Manifold_Gibbs)

invisible(lapply(source_files, sys.source, envir = .GlobalEnv))

worker_init <- function() {
  for (f in source_files) sys.source(f, envir = .GlobalEnv)
  for (sym in c("Benvlp_resp_manifold","Benvlp_MC_resp_gibbs","response_LVI")) {
    if (!exists(sym, mode = "function", inherits = TRUE)) {
      stop(sprintf("worker missing function: %s", sym))
    }
  }
  invisible(TRUE)
}

make_beta_ref_and_resid <- function(
    X, Y, u,
    methods = c("Manifold_Gibbs","MH_Gibbs","CALVI"),
    center_resid = TRUE,
    seeds = NULL,
    Manifold_Gibbs_args = list(n.iter = 10000,
                               n.chains = 1,
                               prior = "Uniform",
                               burnin.prop = 0.5,
                               central_space_sampling_type = "allpossiblepair",
                               Max_rejection = 20,
                               partition = 1000,
                               backupMc = FALSE,
                               show_progress = FALSE,
                               chains_parallel = FALSE),
    MH_Gibbs_args = list(n.iter = 1500, burnin.prop = 0.3, show_progress = FALSE),
    CALVI_args  = list(maxiter = 200, n_iter = 400, show_progress = FALSE, init = "map")
){
  stopifnot(is.matrix(X), is.matrix(Y), nrow(X) == nrow(Y))
  n <- nrow(X)
  get_seed <- local({ i <- 0L; function(){ if (is.null(seeds) || !length(seeds)) return(NULL); i <<- i+1L; seeds[(i-1L) %% length(seeds) + 1L] } })
  
  beta_ref    <- list(); resid_input <- list(); kept <- character(0)
  .quiet_try <- function(expr) tryCatch(expr, error=function(e) structure(list(error=e), class="try-error"))
  mk_resid <- function(beta, mu) {
    mu <- as.numeric(mu)
    Yf <- matrix(1, n, 1) %*% t(mu) + X %*% t(beta)
    R  <- Y - Yf
    if (center_resid) R <- sweep(R, 2, colMeans(R), `-`)
    list(R=R, Yfit=Yf)
  }
  
  if ("Manifold_Gibbs" %in% methods) {
    set.seed(get_seed() %||% NULL)
    fitm <- .quiet_try(do.call(Benvlp_resp_manifold, c(list(X=X, Y=Y, u=u), Manifold_Gibbs_args)))
    if (!inherits(fitm, "try-error")) {
      beta_ref$Manifold_Gibbs <- fitm$beta_mean$Chain_1; resid_input$Manifold_Gibbs <- mk_resid(fitm$beta_mean$Chain_1, fitm$mu_mean$Chain_1); kept <- c(kept, "Manifold_Gibbs")
    } else message(sprintf("[init] MLE 실패: %s", fitm$error$message))
  }
  if ("MH_Gibbs" %in% methods) {
    set.seed(get_seed() %||% NULL)
    fitm <- .quiet_try(do.call(Benvlp_MC_resp_gibbs, c(list(X=X, Y=Y, u=u), MH_Gibbs_args)))
    if (!inherits(fitm, "try-error")) {
      betas   <- unlist(fitm$beta, recursive = FALSE)
      betaBar <- Reduce(`+`, betas) / length(betas)
      muBar   <- Reduce(`+`, lapply(unlist(fitm$mu, recursive = FALSE), as.numeric)) / length(fitm$mu$Chain_1)
      beta_ref$MH_Gibbs <- betaBar; resid_input$MH_Gibbs <- mk_resid(betaBar, muBar); kept <- c(kept, "MH_Gibbs")
    } else message(sprintf("[init] MCMC 실패: %s", fitm$error$message))
  }
  if ("CALVI" %in% methods) {
    set.seed(get_seed() %||% NULL)
    fitm <- .quiet_try(do.call(response_LVI, c(list(X=X, Y=Y, u=u), CALVI_args)))
    if (!inherits(fitm, "try-error")) {
      vp <- fitm[[1]]; beta_ref$CALVI <- vp$beta; resid_input$CALVI <- mk_resid(vp$beta, vp$mu_q); kept <- c(kept, "CALVI")
    } else message(sprintf("[init] CALVI 실패: %s", fitm$error$message))
  }
  if (!length(kept)) stop("초기 적합이 모두 실패했습니다. (methods=", paste(methods, collapse=","), ")")
  
  list(methods = kept, beta_ref = beta_ref[kept], resid_input = resid_input[kept])
}

bootstrap_beta_rmse_foreach <- function(
    X, Y, u,
    B = 200,
    n.iter = 5e5,
    methods = c("Manifold_Gibbs","MH_Gibbs","CALVI"),
    center_resid     = TRUE,
    wild             = FALSE,
    common_indices   = FALSE,
    seeds            = NULL,
    bootstrap_baseline = NULL,
    skip_on_error    = TRUE,
    
    # 심플 파라미터
    beta_ref,                       # named list(method -> matrix) 또는 단일 matrix
    fitters   = NULL,               # named list(method -> function(Xb,Yb,u) -> beta_hat matrix)
    resid_input = NULL,             # named list(method 또는 baseline -> list(R=..., Yfit=...))
    packages = character(0),        # 워커에서 필요 패키지
    out_dir  = "boot_out",          # 결과 폴더
    prefix   = "boot",              # 결과 파일 접두사
    overwrite = FALSE,              # TRUE면 덮어씀
    save_beta = FALSE,              # 각 반복 beta 저장 여부
    progress  = c("print","none","progressr")  # ★ 기본은 print만 사용
){
  progress <- match.arg(progress)
  stopifnot(is.matrix(X), is.matrix(Y), nrow(X) == nrow(Y))
  n <- nrow(X)
  
  # beta_ref 정규화
  if (is.matrix(beta_ref)) {
    beta_ref <- setNames(rep(list(beta_ref), length(methods)), methods)
  } else {
    stopifnot(is.list(beta_ref), !is.null(names(beta_ref)))
    if (!all(methods %in% names(beta_ref))) stop("beta_ref 는 methods 이름을 키로 갖는 named list 여야 합니다.")
  }
  
  # 기본 fitters 준비
  if (is.null(fitters)) {
    fitters <- list()
    if ("Manifold_Gibbs" %in% methods && exists("Benvlp_resp_manifold", mode="function")) {
      fitters$Manifold_Gibbs <- function(Xb, Yb, u){
        fit <- do.call(Benvlp_resp_manifold, c(list(X=Xb, Y=Yb, u=u, 
                                                    n.iter = n.iter,
                                                    n.chains = 1,
                                                    prior = "Uniform",
                                                    burnin.prop = 0.5,
                                                    central_space_sampling_type = "allpossiblepair",
                                                    Max_rejection = 20,
                                                    partition = 1000,
                                                    backupMc = FALSE,
                                                    show_progress = FALSE,
                                                    chains_parallel = FALSE
        )))
        fit$beta_mean$Chain_1
      }
    }
    if ("MH_Gibbs" %in% methods && exists("Benvlp_MC_resp_gibbs", mode="function")) {
      fitters$MH_Gibbs <- function(Xb, Yb, u) {
        fit <- do.call(Benvlp_MC_resp_gibbs, c(list(X=Xb, Y=Yb, u=u,
                                                    n.iter =  n.iter, burnin.prop = 1/2, show_progress = FALSE, init = "map")))
        mlist <- unlist(fit$beta, recursive = FALSE)
        mlist_aligned <- lapply(mlist, function(B) {
          if (sum(B * beta_ref[["MH_Gibbs"]]) < 0) -B else B
        })
        
        Reduce(`+`, mlist_aligned) / length(mlist_aligned)
      }
    }
    if ("CALVI" %in% methods && exists("response_LVI", mode="function")) {
      fitters$CALVI <- function(Xb, Yb, u) {
        fit <- do.call(response_LVI, c(list(X=Xb, Y=Yb, u=u,
                                            maxiter =  100, n_iter = n.iter, show_progress = FALSE, init = "map")))
        fit[[1]]$beta
      }
    }
  }
  if (!all(methods %in% names(fitters))) stop(sprintf("fitters 누락: %s", paste(setdiff(methods, names(fitters)), collapse=", ")))
  
  # residual 모드 검증
  use_resid <- !is.null(resid_input)
  if (use_resid) {
    if (!is.null(bootstrap_baseline)) {
      stopifnot(!is.null(resid_input[[bootstrap_baseline]]))
      R_base  <- resid_input[[bootstrap_baseline]]$R
      Y0_base <- resid_input[[bootstrap_baseline]]$Yfit
      stopifnot(nrow(R_base)==n, nrow(Y0_base)==n)
    } else {
      for (m in methods) {
        stopifnot(!is.null(resid_input[[m]]))
        stopifnot(nrow(resid_input[[m]]$R)==n, nrow(resid_input[[m]]$Yfit)==n)
      }
    }
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 메인에서만 진행 출력(워커에서는 어떤 progressr 신호도 발생시키지 않음)
  if (progress == "print") {
    cat(sprintf("[bootstrap] B = %d  (%s)\n", B, if (use_resid) "residual" else "pairs"))
  }
  
  p <- progressor(along = 1:B)
  
  res_files <- foreach::foreach(
    b = seq_len(B),
    .combine = "c",
    .inorder = FALSE,
    .options.future = list(
      seed = TRUE,
      packages = character(0),   # 필요 패키지 있으면 넣기
      globals  = structure(TRUE, add = TRUE)  # ←★ 여기! 리스트 금지
    )
  ) %dofuture% {
    worker_init()  # 워커에서 함수 로드 확인
        
    fpath <- file.path(out_dir, sprintf("%s_%05d.rds", prefix, b))
    if (!overwrite && file.exists(fpath)) return(fpath)
    
    out <- tryCatch({
      if (!is.null(seeds)) set.seed(seeds[b])
      Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
      
      rmse_b <- setNames(rep(NA_real_, length(methods)), methods)
      errs_b <- setNames(rep(NA_character_, length(methods)), methods)
      beta_b <- if (save_beta) setNames(vector("list", length(methods)), methods) else NULL
      
      # 공통 인덱스/부호
      if (!is.null(resid_input) && common_indices) {
        if (wild) signs <- sample(c(-1,1), n, TRUE) else idx <- sample.int(n, n, TRUE)
      } else if (is.null(resid_input)) {
        idx_pairs <- sample.int(n, n, TRUE)
      }
      
      for (m in methods) {
        if (!is.null(resid_input)) {
          if (!is.null(bootstrap_baseline)) {
            Rm <- resid_input[[bootstrap_baseline]]$R; Y0m <- resid_input[[bootstrap_baseline]]$Yfit
          } else {
            Rm <- resid_input[[m]]$R;                    Y0m <- resid_input[[m]]$Yfit
          }
          Yb <- if (!is.null(resid_input) && common_indices) {
            if (wild) Y0m + Rm * matrix(signs, n, ncol(Rm)) else Y0m + Rm[idx, , drop = FALSE]
          } else {
            if (wild) { s <- sample(c(-1,1), n, TRUE); Y0m + Rm * matrix(s, n, ncol(Rm)) }
            else      { id <- sample.int(n, n, TRUE);  Y0m + Rm[id, , drop = FALSE] }
          }
          Xb <- X
        } else {
          Yb <- Y[idx_pairs, , drop = FALSE]
          Xb <- X[idx_pairs, , drop = FALSE]
        }
        
        ans <- try(fitters[[m]](Xb, Yb, u), silent = TRUE)
        if (inherits(ans, "try-error")) {
          errs_b[m] <- as.character(ans)
          if (!skip_on_error) stop(sprintf("b=%d, method=%s 실패: %s", b, m, errs_b[m]))
          next
        }
        beta_hat <- as.matrix(ans)
        
        br <- beta_ref[[m]]
        if (!all(dim(beta_hat) == dim(br))) {
          errs_b[m] <- sprintf("dim mismatch: beta_hat %s vs beta_ref %s",
                               paste(dim(beta_hat), collapse="x"), paste(dim(br), collapse="x"))
          next
        }
        rmse_b[m] <- sqrt(mean((beta_hat - br)^2))
        
        if (save_beta) beta_b[[m]] <- beta_hat
      }
      
      list(
        b = b, methods = methods, rmse = rmse_b,
        errors = errs_b,
        beta = if (save_beta) beta_b else NULL,
        mode = if (!is.null(resid_input)) "residual" else "pairs",
        opts = list(center_resid=center_resid, wild=wild,
                    common_indices=common_indices,
                    bootstrap_baseline=bootstrap_baseline)
      )
    }, error = function(e) {
      list(
        b = b, methods = methods,
        rmse = setNames(rep(NA_real_, length(methods)), methods),
        errors = setNames(rep(conditionMessage(e), length(methods)), methods),
        beta = NULL,
        mode = if (!is.null(resid_input)) "residual" else "pairs",
        opts = list(center_resid=center_resid, wild=wild,
                    common_indices=common_indices,
                    bootstrap_baseline=bootstrap_baseline)
      )
    })
    
    saveRDS(out, fpath)
    p(sprintf("b=%g", b))
    fpath
  }
  
  if (progress == "print") {
    # 간단 진행 출력(파일 개수 기준)
    done <- sum(file.exists(file.path(out_dir, sprintf("%s_%05d.rds", prefix, 1:B))))
    cat(sprintf("[done] %d / %d\n", done, B))
  }
  
  invisible(res_files)
}

## ---- 결과 확인용 함수 ----
collect_bootstrap <- function(out_dir = "boot_out", prefix = "boot") {
  files <- Sys.glob(file.path(out_dir, sprintf("%s_*.rds", prefix)))
  if (!length(files)) return(list(wide = data.frame(), long = data.frame()))
  ord <- order(as.integer(sub("^.*_(\\d+)\\.rds$", "\\1", files)))
  files <- files[ord]
  
  rows <- lapply(files, function(f) {
    z <- readRDS(f)
    # rmse는 named numeric: c(MLE=..., MCMC=..., LVI=...)
    data.frame(
      b = z$b,
      t(z$rmse),                 # wide: rmse_MLE, rmse_MCMC, ...
      file = basename(f),
      check.names = FALSE
    )
  })
  wide <- do.call(rbind, rows)
  method_cols <- setdiff(names(wide), c("b","file"))
  names(wide)[match(method_cols, names(wide))] <- paste0("rmse_", method_cols)
  
  # long형도 생성
  long <- reshape(
    wide,
    varying = grep("^rmse_", names(wide), value = TRUE),
    v.names = "rmse",
    timevar = "method",
    times   = sub("^rmse_", "", grep("^rmse_", names(wide), value = TRUE)),
    direction = "long"
  )
  long <- long[order(long$b, long$method), c("b","method","rmse","file")]
  rownames(long) <- NULL
  
  list(wide = wide, long = long)
}

summarize_bootstrap <- function(df_long) {
  if (!nrow(df_long)) return(data.frame())
  agg <- aggregate(df_long$rmse, by = list(method = df_long$method), FUN = function(x) {
    c(B_eff = sum(is.finite(x)),
      B_miss = sum(!is.finite(x)),
      mean = mean(x, na.rm = TRUE),
      sd   = stats::sd(x, na.rm = TRUE))
  })
  out <- do.call(data.frame, agg)
  names(out) <- c("method", "B_eff", "B_miss", "mean", "sd")
  out
}

ci_table <- function(df_long) {
  if (!nrow(df_long)) return(data.frame())
  methods <- unique(df_long$method)
  out <- lapply(methods, function(m){
    x <- df_long$rmse[df_long$method == m]
    x <- x[is.finite(x)]
    if (!length(x)) return(data.frame(method=m, B_eff=0, mean=NA, sd=NA, 
                                      ci_norm_low=NA, ci_norm_high=NA,
                                      q2.5=NA, q97.5=NA))
    m0 <- mean(x); s0 <- sd(x); B <- length(x)
    data.frame(method = m,
               B_eff = B,
               mean = m0,
               sd   = s0,
               ci_norm_low  = m0 - 1.96 * s0 / sqrt(B),
               ci_norm_high = m0 + 1.96 * s0 / sqrt(B),
               q2.5  = as.numeric(quantile(x, 0.025, names = FALSE)),
               q97.5 = as.numeric(quantile(x, 0.975, names = FALSE)))
  })
  do.call(rbind, out)
}

inspect_bootstrap_errors <- function(out_dir = "boot_out", prefix = "boot") {
  files <- Sys.glob(file.path(out_dir, sprintf("%s_*.rds", prefix)))
  if (!length(files)) return(data.frame())
  ord <- order(as.integer(sub("^.*_(\\d+)\\.rds$", "\\1", files)))
  files <- files[ord]
  rows <- lapply(files, function(f){
    z <- readRDS(f)
    ez <- z$errors
    if (is.null(ez)) {
      data.frame(b = z$b, file = basename(f))
    } else {
      as.data.frame(c(list(b = z$b), setNames(as.list(ez), paste0("err_", names(ez))), list(file = basename(f))),
                    check.names = FALSE)
    }
  })
  do.call(rbind, rows)
}

plot_bootstrap <- function(df_long) {
  if (!nrow(df_long)) { message("No data"); return(invisible()) }
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1,4))
  # 히스토그램: 방법별 겹치지 않게 별도 패널로
  methods <- unique(df_long$method)
  for (m in methods) {
    x <- df_long$rmse[df_long$method == m]
    hist(x, main = paste("Histogram -", m), xlab = "RMSE")
  }
  # 박스플롯(한 화면)
  boxplot(rmse ~ method, data = df_long, main = "Bootstrap RMSE by method", ylab = "RMSE")
}



## ---- 적합 ----
library(purrr)
library(pheatmap)
library(RColorBrewer)
library(profvis)
library(htmltools)

my_palette <- rev(brewer.pal(11, "RdYlBu"))

data <- read.csv(file="real data/20201211update.csv", header=TRUE, sep=",")
head(data)
names(data)
dim(data)

data1 <- data
factor(data1$Education.Level) #2 (High school graduate or GED),3,4
factor(data1$Household.Income) #1,2,3,4 (Under 20,000/20,001 to 50,000/50,001 to 100,000/More than 100,000)
factor(data1$Marital.Status) #1,2,3 (Divorced),4(Widowed)

attach(data1)

College=Graduate <- rep(0, length(Education.Level))
College[which(Education.Level==3)] <- 1
Graduate[which(Education.Level==4)] <- 1

Single = Married <- rep(0, length(Marital.Status))
Single[which(Marital.Status==1)] <- 1
Married[which(Marital.Status==2)] <- 1

Race[which(Race=="WHITE")] <- 1
Race[which(Race=="AFRNAMER")] <- 0

X <- cbind(Gail.Score, Age, BMI, Race,
           College, Graduate,
           Household.Income,
           Single, Married)

name <- colnames(X)
X <- matrix(as.numeric(X), dim(X))
colnames(X) <- name

Y <- data1[,12:21]
names(Y) <- sub("^X", "", names(Y))

detach(data1)
X_ww <- as.matrix(scale(X[1:30,-4]))
Y_ww <- as.matrix(scale(log(Y[1:30,])))
X_bw <- as.matrix(scale(X[31:60,-4]))
Y_bw <- as.matrix(scale(log(Y[31:60,])))

init_ww <- make_beta_ref_and_resid(
  X = X_ww, Y = Y_ww, u = 1,
  methods = c("Manifold_Gibbs", "MH_Gibbs", "CALVI"),
  center_resid = TRUE,
  seeds = 123,  # 재현성 원하면 아무 정수(또는 벡터)
  Manifold_Gibbs_args = list(n.iter = 5e4,
                             n.chains = 1,
                             prior = "Uniform",
                             burnin.prop = 0.5,
                             central_space_sampling_type = "allpossiblepair",
                             Max_rejection = 20,
                             partition = 1000,
                             backupMc = FALSE,
                             show_progress = FALSE,
                             chains_parallel = FALSE),
  MH_Gibbs_args = list(n.iter=5e4, burnin.prop=1/2, show_progress=FALSE, init = "map"),
  CALVI_args  = list(maxiter=100, n_iter=5e4, show_progress=FALSE, init="map")
)

bootstrap_beta_rmse_foreach(
  X = X_ww, Y = Y_ww, u = 1,
  B = 1000,
  n.iter = 5e4,
  methods = c("MH_Gibbs", "CALVI"),
  seeds = 1:1000,
  bootstrap_baseline = "CALVI",     # residual 모드에서 공통 잔차/적합치 기반으로
  resid_input = init_ww$resid_input,      # 없으면 pairs bootstrap로 동작
  beta_ref = init_ww$beta_ref,
  out_dir  = "boot_env_ww_Manifold",
  prefix   = "env",
  overwrite = TRUE,
  save_beta = TRUE,
  packages = c("coda", "posterior")
)

aaa <- readRDS("boot_env_ww_Manifold/env_00002.rds")
aaa$beta

init_bw <- make_beta_ref_and_resid(
  X = X_bw, Y = log(Y_bw), u = 1,
  methods = c("MLE","MCMC","LVI"),
  center_resid = TRUE,
  seeds = 123,  # 재현성 원하면 아무 정수(또는 벡터)
  Manifold_Gibbs_args = list(n.iter = 5e4,
                             n.chains = 1,
                             prior = "Uniform",
                             burnin.prop = 0.5,
                             central_space_sampling_type = "allpossiblepair",
                             Max_rejection = 20,
                             partition = 1000,
                             backupMc = FALSE,
                             show_progress = FALSE,
                             chains_parallel = FALSE)
  # mcmc_args = list(n.iter=5e4, burnin.prop=1/2, show_progress=FALSE, init = "map"),
  # lvi_args  = list(maxiter=100, n_iter=5e4, show_progress=FALSE, init="map")
)


bootstrap_beta_rmse_foreach(
  X = X1, Y = Y1, u = 1,
  B = 1000,
  methods = c("MLE","MCMC","LVI"),
  seeds = 1:1000,
  bootstrap_baseline = NULL,     # residual 모드에서 공통 잔차/적합치 기반으로
  resid_input = init_bw$resid_input,      # 없으면 pairs bootstrap로 동작
  beta_ref = init_bw$beta_ref,
  out_dir  = "boot_env_bw_NULL",
  prefix   = "env",
  overwrite = TRUE,
  save_beta = TRUE,
  packages = c()                  # 워커에서 필요한 패키지명 있으면 추가
)

## ---- 결과 확인 ----

dat <- collect_bootstrap(out_dir = "boot_env_ww_Manifold", prefix = "env")
dat$wide
dat$long

sumtab <- summarize_bootstrap(dat$long)
sumtab

ci_tab <- ci_table(dat$long)
ci_tab

errtab <- inspect_bootstrap_errors(out_dir = "boot_env_ww", prefix = "env")
head(errtab, 10)

plot_bootstrap(dat$long)


## ---- real data mean plot ---- 

bw_real_CALVI <- response_LVI(
  maxiter = 100, n_iter = 50000, tol=1e-6,
  show_progress = FALSE, make_positive = TRUE, init = "map", compute_llik = TRUE,
  X = X_bw, Y = Y_bw,
  nuvalue = 0, nu0value = 0, Avalue = 1e6, Mvalue = 1e6, u = 1
)
bw_real_CALVI[[3]]

ww_real_CALVI <- response_LVI(
  maxiter = 100, n_iter = 50000, tol=1e-6,
  show_progress = FALSE, make_positive = TRUE, init = "map", compute_llik = TRUE,
  X = X_ww, Y = Y_ww,
  nuvalue = 0, nu0value = 0, Avalue = 1e6, Mvalue = 1e6, u = 1
)
ww_real_CALVI[[3]]

bw_real_MH_Gibbs <- Benvlp_MC_resp_gibbs(
  X = X_bw, Y = Y_bw,
  u = 1,
  n.iter = 5e4,
  burnin.prop = 1 / 2,
  n.chains = 1,
  scan = "full",
  init = "map",
  A_proposal = "stiefel",
  jacobian_only_gamma = FALSE,
  compute_llik = TRUE
)
bw_real_MH_Gibbs$total_time

ww_real_MH_Gibbs <- Benvlp_MC_resp_gibbs(
  X = X_ww, Y = Y_ww,
  u = 1,
  n.iter = 5e4,
  burnin.prop = 1 / 2,
  n.chains = 1,
  scan = "full",
  init = "map",
  A_proposal = "stiefel",
  jacobian_only_gamma = FALSE,
  compute_llik = TRUE
)
ww_real_MH_Gibbs$total_time

bw_real_Manifold_Gibbs <- Benvlp_resp_manifold(
  X = X_bw, Y = Y_bw, u = 1,
  n.iter = 50000,
  n.chains = 1,
  prior = "Uniform",
  burnin.prop = 0.5,
  central_space_sampling_type = "allpossiblepair",
  Max_rejection = 20,
  partition = 1000,
  backupMc = FALSE,
  show_progress = FALSE,
  chains_parallel = FALSE
)
bw_real_Manifold_Gibbs$total_time

ww_real_Manifold_Gibbs <- Benvlp_resp_manifold(
  X = X_ww, Y = Y_ww, u = 1,
  n.iter = 50000,
  n.chains = 1,
  prior = "Uniform",
  burnin.prop = 0.5,
  central_space_sampling_type = "allpossiblepair",
  Max_rejection = 20,
  partition = 1000,
  backupMc = FALSE,
  show_progress = FALSE,
  chains_parallel = FALSE
)
ww_real_Manifold_Gibbs$total_time

library(pheatmap)
library(grid)
library(gridExtra)

## ---- 1) beta estimates
beta_ww_mh <- Reduce(`+`, ww_real_MH_Gibbs$beta$Chain_1) / length(ww_real_MH_Gibbs$beta$Chain_1)
beta_bw_mh <- Reduce(`+`, bw_real_MH_Gibbs$beta$Chain_1) / length(bw_real_MH_Gibbs$beta$Chain_1)

beta_ww_calvi <- ww_real_CALVI[[1]]$beta
beta_bw_calvi <- bw_real_CALVI[[1]]$beta

beta_ww_mani <- ww_real_Manifold_Gibbs$beta_mean$Chain_1
beta_bw_mani <- bw_real_Manifold_Gibbs$beta_mean$Chain_1

## ---- 2) labels from data
x_labs <- colnames(X[,-4])
y_labs <- colnames(Y)

set_dimnames <- function(B) { dimnames(B) <- list(y_labs, x_labs); B }

beta_ww_mh    <- set_dimnames(beta_ww_mh)
beta_bw_mh    <- set_dimnames(beta_bw_mh)
beta_ww_calvi <- set_dimnames(beta_ww_calvi)
beta_bw_calvi <- set_dimnames(beta_bw_calvi)
beta_ww_mani  <- set_dimnames(beta_ww_mani)
beta_bw_mani  <- set_dimnames(beta_bw_mani)

## ---- 3) common style (shared color scale)
my_palette <- rev(brewer.pal(11, "RdYlBu"))
bk <- seq(-0.1, 0.1, length.out = 51)

hm_grob <- function(mat, main, show_legend = FALSE,
                    pad_top_pt = 0.5, pad_bot_pt = 0.5, title_pad_pt = 10) {
  gt <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(my_palette)(50),
    breaks = bk,
    border_color = NA,
    cellwidth = 20,
    cellheight = 8,
    angle_col = 270,
    legend = show_legend,
    main = main,
    silent = TRUE
  )$gtable
  
  ## 핵심: pheatmap이 잡아둔 위/아래 여백을 강제로 줄임
  gt$heights[1] <- unit(pad_top_pt, "pt")
  #gt$heights[length(gt$heights)] <- unit(pad_bot_pt, "pt")
  
  title_row <- grep("main", gt$layout$name)
  if (length(title_row) > 0) {
    gt$heights[title_row] <- unit(title_pad_pt, "pt")
  }
  
  gt
}

## ---- build rows
row_block <- function(method_name, mat_ww, mat_bw, legend = FALSE) {
  
  text_grob <- textGrob(
    method_name,
    rot = 90,
    x = 0.85,   # 0=왼쪽, 1=오른쪽 (그림쪽)
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  g_ww <- hm_grob(mat_ww, "White Women", show_legend = FALSE)
  g_bw <- hm_grob(mat_bw, "Black Women", show_legend = legend)
  
  arrangeGrob(
    text_grob, g_ww, g_bw,
    ncol = 3,
    widths = c(0.05, 0.95, 0.95)
  )
}

r1 <- row_block("CALVI", beta_ww_calvi, beta_bw_calvi, legend = TRUE)
r2 <- row_block("MH-Gibbs", beta_ww_mh, beta_bw_mh, legend = TRUE)
r3 <- row_block("Manifold-Gibbs", beta_ww_mani, beta_bw_mani, legend = TRUE)

full_plot <- arrangeGrob(r1, r2, r3, nrow = 3)  # 또는 grid.arrange 대신 arrangeGrob 추천

png("/home/jupyter-tmdgus4970/2025_CALVI/figures/beta_compare.png",
    width = 2400,
    height = 3000, 
    res = 300)

grid.newpage()
grid.draw(full_plot)
dev.off()

## ---- real data sd plot ---- 

library(pheatmap)
library(grid)
library(gridExtra)

## -------------------------
## Helper: build C_A from A
## -------------------------
build_CA <- function(A, u) {
  # C_A = [I_u; A]  (r x u), where A is (r-u) x u
  rbind(diag(u), A)
}

# vec(A) indexing: A is (ru x u), column-major
idx_A <- function(a, k, ru) a + (k - 1L) * ru

## -------------------------
## MH-Gibbs: entrywise quantiles
## -------------------------
mh_summaries <- function(beta_list, probs = c(0.025, 0.975)) {
  S <- length(beta_list)
  r <- nrow(beta_list[[1]])
  p <- ncol(beta_list[[1]])
  arr <- array(unlist(beta_list), dim = c(S, r, p))
  mean_mat  <- apply(arr, c(2,3), mean)
  lower_mat <- apply(arr, c(2,3), quantile, probs = probs[1], type = 8)
  upper_mat <- apply(arr, c(2,3), quantile, probs = probs[2], type = 8)
  list(lower = lower_mat, mean = mean_mat, upper = upper_mat)
}

## ------------------------- 
## CALVI: delta-method intervals for beta = C_A %*% solve(crossprod(C_A)) %*% eta 
## Assumption: A ⟂ eta (MF), Cov(vec(A)) = HA, Cov(vec(eta)) = V_q ⊗ U_q 
## -------------------------
calvi_summaries <- function(A_mean, eta_mean, HA, U_q, V_q, u) {
  ru <- nrow(A_mean)
  p  <- ncol(eta_mean)
  r  <- u + ru
  
  # C, G, invG, P at the mean
  C    <- build_CA(A_mean, u)              # r x u
  G    <- crossprod(C)                     # u x u
  invG <- solve(G)                         # u x u
  Pmat <- C %*% invG                        # r x u
  
  beta_mean <- Pmat %*% eta_mean            # r x p
  
  # eta contribution common term
  PUUt <- Pmat %*% U_q %*% t(Pmat)          # r x r
  
  lower <- matrix(NA_real_, r, p)
  upper <- matrix(NA_real_, r, p)
  sdmat <- matrix(0, r, p)
  
  for (j in seq_len(p)) {
    eta_j <- eta_mean[, j, drop = FALSE]    # u x 1
    vjj   <- V_q[j, j]                      # scalar (marginal)
    
    for (i in seq_len(r)) {
      mu <- beta_mean[i, j]
      
      # 1) Var from eta | A
      var_eta <- vjj * PUUt[i, i]
      
      # 2) Var from A via delta method
      # gA component for beta_{i,j}: length ru*u
      gA <- numeric(ru * u)
      
      for (k in seq_len(u)) {
        for (a in seq_len(ru)) {
          # D = dC/dA_{a,k}: only entry at row (u+a), col k is 1
          # We'll form effect analytically without building full D dense:
          # But r,u are small; easiest is to build D as sparse-like matrix
          D <- matrix(0, r, u)
          D[u + a, k] <- 1
          
          dG  <- crossprod(D, C) + crossprod(C, D)     # u x u
          dP  <- D %*% invG - C %*% invG %*% dG %*% invG  # r x u
          
          dbeta_ij <- as.numeric(dP %*% eta_j)[i]
          gA[idx_A(a, k, ru)] <- dbeta_ij
        }
      }
      
      var_A <- as.numeric(t(gA) %*% HA %*% gA)
      
      var <- var_eta + var_A
      sd  <- sqrt(max(0, var))
      
      sdmat[i, j] <- sd
      lower[i, j] <- mu - 1.96 * sd
      upper[i, j] <- mu + 1.96 * sd
    }
  }
  
  list(lower = lower, mean = beta_mean, upper = upper, sd = sdmat)
}

## -------------------------
## Manifold-Gibbs: assume mean & sd available (Gaussian approx)
## -------------------------
mani_summaries <- function(beta_mean, beta_sd) {
  list(
    lower = beta_mean - 1.96 * beta_sd,
    mean  = beta_mean,
    upper = beta_mean + 1.96 * beta_sd
  )
}

## =========================================================
## 1) Build summaries (WW only)
## =========================================================
# MH
mh <- mh_summaries(ww_real_MH_Gibbs$beta$Chain_1)

# Manifold
if (is.null(ww_real_Manifold_Gibbs$beta_mean) || is.null(ww_real_Manifold_Gibbs$beta_sd)) {
  stop("ww_real_Manifold_Gibbs needs both $beta_mean and $beta_sd for this plot.")
}
mani <- mani_summaries(ww_real_Manifold_Gibbs$beta_mean$Chain_1, ww_real_Manifold_Gibbs$beta_sd$Chain_1)

# CALVI
A_mean   <- ww_real_CALVI[[1]]$hatA
eta_mean <- ww_real_CALVI[[1]]$eta_q
HA       <- ww_real_CALVI[[1]]$HA
U_q      <- ww_real_CALVI[[1]]$U_q
V_q      <- ww_real_CALVI[[1]]$V_q

# infer u if not in env
if (!exists("u")) u <- nrow(eta_mean)

calvi <- calvi_summaries(A_mean, eta_mean, HA, U_q, V_q, u = 1)

## =========================================================
## 2) Optional: add dimnames from X/Y
## =========================================================
x_labs <- colnames(X[,-4])
y_labs <- colnames(Y)

set_names <- function(M) { dimnames(M) <- list(y_labs, x_labs); M }
for (nm in c("lower","mean","upper")) {
  mh[[nm]]    <- set_names(mh[[nm]])
  mani[[nm]]  <- set_names(mani[[nm]])
  calvi[[nm]] <- set_names(calvi[[nm]])
}

## =========================================================
## 3) Plot 3x3
## rows: CALVI / MH-Gibbs / Manifold-Gibbs
## cols: Lower95 / Mean / Upper95
## =========================================================
bk <- seq(-0.5, 0.5, length.out = 51)

hm_grob <- function(mat, show_legend = FALSE, pad_pt = 0.5, title_pad_pt = 10) {
  gt <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(my_palette)(50),
    breaks = bk,
    border_color = NA,
    cellwidth = 20,
    cellheight = 8,
    angle_col = 270,
    legend = show_legend,
    silent = TRUE
  )$gtable
  
  # tighten outer padding
  gt$heights[1] <- unit(pad_pt, "pt")
  #gt$heights[length(gt$heights)] <- unit(pad_pt, "pt")
  
  # give a bit more air above heatmap for column titles (we add titles outside)
  main_row <- grep("main", gt$layout$name)
  if (length(main_row) > 0) gt$heights[main_row] <- unit(title_pad_pt, "pt")
  
  gt
}

row_block3 <- function(method_name, mats3, legend_col = 3L) {
  # mats3: list(lower, mean, upper)
  text_grob <- textGrob(method_name, rot = 90,
                        x = 0.9,
                        gp = gpar(fontsize = 14, fontface = "bold"))
  
  g1 <- hm_grob(mats3[[1]], show_legend = (legend_col == 1L))
  g2 <- hm_grob(mats3[[2]], show_legend = (legend_col == 2L))
  g3 <- hm_grob(mats3[[3]], show_legend = (legend_col == 3L))
  
  arrangeGrob(
    text_grob, g1, g2, g3,
    ncol = 4,
    widths = c(0.03, 1, 1, 1)
  )
}

# column header
col_header <- arrangeGrob(
  nullGrob(),
  textGrob("Lower 95%", gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob("Mean",      gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob("Upper 95%", gp = gpar(fontsize = 16, fontface = "bold")),
  ncol = 4,
  widths = c(0.03, 1, 1, 1)
)

r_calvi <- row_block3("CALVI",          list(calvi$lower, calvi$mean, calvi$upper), legend_col = 3L)
r_mh    <- row_block3("MH-Gibbs",       list(mh$lower,    mh$mean,    mh$upper),    legend_col = 3L)
r_mani  <- row_block3("Manifold-Gibbs", list(mani$lower,  mani$mean,  mani$upper),  legend_col = 3L)

full_plot <- arrangeGrob(
  col_header,
  r_calvi,
  r_mh,
  r_mani,
  nrow = 4,
  heights = c(0.08, 1, 1, 1)
)

png("/home/jupyter-tmdgus4970/2025_CALVI/figures/ww_beta_coverage_3x3.png",
    width = 3500, height = 2400, res = 300)
grid.newpage()
grid.draw(full_plot)
dev.off()

