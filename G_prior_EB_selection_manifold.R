#' Returns G_prior matrix for the Emperical Bayes Prior selection method
# @export
G_prior_EB_selection_new <- function(O, omega, omega0, sample_size = NULL, comfidence_level = "LOW") {
  if (is.null(sample_size)) {
    exponent_factor <- 0.01
  }



  if (comfidence_level == "LOW") {
    exponent_factor <- sample_size * 0.001
  }
  if (comfidence_level == "MEDIUM") {
    exponent_factor <- sample_size * 0.05
  }
  if (comfidence_level == "HIGH") {
    exponent_factor <- sample_size * 0.08
  }
  if (exponent_factor > 5) {
    exponent_factor <- 5
  }

  r <- dim(O)[2]
  # ordered_chisq=sort(rchisq(r,20,0)) #### the variability 20 can be designed in a systemetic way.
  ordered_chisq <- (1:r)^exponent_factor
  rnk <- rank(c(omega, omega0))

  D <- diag(ordered_chisq[r - rnk + 1])
  G <- O %*% D %*% t(O)
  return(G)
}

G_prior_EB_selection <- function(O, omega, omega0) {
  r <- dim(O)[2]
  ordered_chisq <- sort(rchisq(r, 20, 0)) #### the variability 20 can be designed in a systemetic way.
  rnk <- rank(c(omega, omega0))

  D <- diag(ordered_chisq[r - rnk + 1])
  G <- O %*% D %*% t(O)
  return(G)
}
