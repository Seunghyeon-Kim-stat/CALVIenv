#' rTrunkInvGamma
#' generates a truncated invarse gamma random variable by rejection sampling method. If there is too many rejections, then it generates approximate sample according to two steps. initially it tries to get approximate numerical cdf of the truncated invarse gamma
#' random variable using the functions "Gammad" and "Truncate" defined in the r package "distr". If the truncation is too acuate to use the approximation by the package "distr" then a
#' discrete approximation is used to sample from the random variable.
#' @param alpha = shape parameter
#' @param beta = scale parameter
#' @param l= lower truncation
#' @param u=upper truncation
#' @param Max_rejection= once the Max_rejection number of rejections occurs it stops the rejection sampler
#' and calls the approximate sampler based on the package "distr"
#' @param partition : once the approximate sampler based on the package "distr" fails then a discrete approximate sampler
#' is applied. log of the partion ( in base 10) is the precession of the discrete approximation. If precission upto 3 decimal places
#' is desired then use partion=100. Default value of partition is 10
# @export

library(LearnBayes)

rTrunkInvGamma <- function(alpha, beta, l = 0, u = Inf, Max_rejection = 1, partition = 1000, method = 0) {
  # a=nimble::pinvgamma(q=l,shape=alpha, rate=beta)
  # b=nimble::pinvgamma(q=u,shape=alpha, rate=beta)
  # u <- runif(1,a,b )
  # y= nimble::qinvgamma(p=u,shape=alpha, rate=beta)


  Max_rejection <- 5

  y <- -1
  i <- 1
  while ((y < 0) * (i < Max_rejection)) {
    r <- rigamma(n = 1, alpha, beta)
    y <- ifelse((r >= l) * (r <= u), r, -1)
    i <- i + 1
  }
  # print(c(alpha,beta,l,u))
  if (y < 0) {
    # browser()

    if (method == 0) {
      y <- tryCatch(rtruncated_inv_gamma(n = 1, alpha = alpha, beta = beta, l = l, u = u), error = function(e) {
        (rDiscreteTInvGamma(n = 1, alpha = alpha, beta = beta, l = l, u = u, partition))
      })
      if (is.infinite(y)) {
        y <- rDiscreteTInvGamma(n = 1, alpha = alpha, beta = beta, l = l, u = u, partition)
      }
    } else if (method == 1) {
      y <- tryCatch(ars_tinv_gamma(n = 1, alpha = alpha, beta = beta, l = l, u = u), error = function(e) {
        (rDiscreteTInvGamma(n = 1, alpha = alpha, beta = beta, l = l, u = u, partition))
      })
    }
    # y = try(  rTIgammaInvCdf(list(alpha=alpha,beta=beta,l=l,u=u)),  1 )
    # if(!is.numeric(y)){
    else if (method == 2) {
      y <- rDiscreteTInvGamma(n = 1, alpha = alpha, beta = beta, l = l, u = u, partition)
    }
    # }
  }
  return(y)
}




ars_tinv_gamma <- function(n = 1, alpha = alpha, beta = beta, l = l, u = u, epsilon = .00001) {
  if (u == Inf) {
    u <- l * exp(1) * ((1 + exp(1)) * (1 - epsilon) / epsilon)^(1 / alpha)
    # print(paste("truncated at upper=",u,"instead of upper=Inf"))
  }
  if (u <= 0) {
    print("Warning: Upper bound for truncated inverse gamma is  Negative")
    return(NULL)
  }
  if (l < 0) {
    print(l)
    print("Warning: Lower bound for truncated inverse gamma is  Negative")
    l <- 0
  }
  if (l == u) {
    print("Warning: Lower bound and Upper bound  for truncated inverse gamma is SAME")
    return(l)
  }
  if (l > u) {
    print("Serious Warning: Lower bound is greater than the Upper bound  for truncated inverse gamma.")
    return(l)
  }
  # browser()
  lf <- function(z) {
    return(-beta * z + (alpha - 1) * log(z))
  }
  lf_prime <- function(z) {
    return(-beta + (alpha - 1) / (z))
  }
  l1 <- 1 / u
  u1 <- ifelse(l == 0, 1000 * l1, 1 / l)
  startingpoints <- c(.99 * l1 + .01 * u1, .5 * l1 + .5 * u1, .01 * l1 + .99 * u1)
  # y=tryCatch(ars(n=1,f=lf,fprima=lf_prime,x=startingpoints,ns=100,m=3,emax=64,lb=TRUE,ub=TRUE,xlb=Inf,xub=u1), error=function(e){ return(3) })
  y <- ars(n = 1, f = lf, fprima = lf_prime, x = startingpoints, ns = 100, m = 3, emax = 64, lb = TRUE, ub = TRUE, xlb = l1, xub = u1)
  return(1 / y)
}



####################



rtruncated_inv_gamma <- function(n = 1, alpha = alpha, beta = beta, l = l, u = u, epsilon = .00001) {
  if (u == Inf) {
    u <- l * exp(1) * ((1 + exp(1)) * (1 - epsilon) / epsilon)^(1 / alpha)
    # print(paste("truncated at upper=",u,"instead of upper=Inf"))
  }
  if (u <= 0) {
    print("Warning: Upper bound for truncated inverse gamma is  Negative")
    return(NULL)
  }
  if (l < 0) {
    print(l)
    print("Warning: Lower bound for truncated inverse gamma is  Negative")
    l <- 0
  }
  if (l == u) {
    print("Warning: Lower bound and Upper bound  for truncated inverse gamma is SAME")
    return(l)
  }
  if (l > u) {
    print("Serious Warning: Lower bound is greater than the Upper bound  for truncated inverse gamma.")
    return(l)
  }
  # browser()
  l1 <- 1 / u
  u1 <- 1 / l
  # ifelse(l==0, 1000*l1 , 1/l)

  y <- rtrunc(1, "gamma", l = l1, u = u1, shape = alpha, rate = beta)
  return(1 / y)
}

qtrunc <- function(p, spec, l = -Inf, u = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(l, ...) + p * (G(u, ...) - G(l, ...)), ...)
  return(tt)
}
rtrunc <- function(n, spec, l = -Inf, u = Inf, ...) {
  u_rnd <- runif(n, min = 0, max = 1)
  x <- qtrunc(u_rnd, spec, l = l, u = u, ...)
  return(x)
}
