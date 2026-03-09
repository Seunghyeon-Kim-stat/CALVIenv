#' returns approximate sample from truncated invarse gamma distribution
#' @param lis: a list containing lis$alpha: Shape parameter, lis$beta: rate parameter, lis$l = lower bound , lis$u = upper bound
#'
rTIgammaInvCdf <- function(lis) {
  # print("inside rTIgammaInvCdf")
  alpha <- lis$alpha
  beta <- lis$beta
  Lower <- lis$l
  Upper <- lis$u
  G0 <- Gammad(scale = (1 / beta), shape = alpha) ## generates a Gamma distribution
  G <- 1 / G0 ## the corresponding inverse Gamma
  TG <- Truncate(G, lower = Lower, upper = Upper)
  return((q(TG)(runif(1, 0, 1))))
}
