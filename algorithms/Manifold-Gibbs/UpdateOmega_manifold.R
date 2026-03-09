#' Updates one step of the Omega vector. Omega is  diagonal elements of the a variance covariance matrix.
#' usses the rTrunkInvGamma for sampling of the diagonal elements of Omega matrix
#' @param omega = omega's are positive elements arranged in decreasing order.
#' @param A= Beta rate parameters for the omega components. vector of length u.
#' @param alpha= shape parameter for all the omega components.shape here is a scaler quantity.
# @export
UpdateOmega <- function(omega, A, alpha, Max_rejection = 10) {
  # make sure that the omega components are in gereasing order.
  u <- length(omega)
  omega_tmp <- c(Inf, omega, 0)
  for (i in 2:(u + 1)) {
    # print(paste("OMEGA=",omega_tmp))
    omega_tmp[i] <- rTrunkInvGamma(alpha = alpha[i - 1], beta = A[i - 1], l = omega_tmp[i + 1], u = omega_tmp[i - 1], Max_rejection = Max_rejection)
  }
  # print(omega_tmp)
  return(omega_tmp[2:(u + 1)])
}
