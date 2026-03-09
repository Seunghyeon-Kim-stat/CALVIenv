#' Centering a data matrix
#' Takes data matrix where columns represents variabes and center each columns
#' @param M is a real matrix with variables stored in columns
#' @return Centered Data matrix
# @export
center_data_matrix <- function(M) {
  # return(M-t(replicate(dim(M)[1],(apply(M,2,"mean")))))
  # return(t(t(M)-(apply(M,2,"mean"))))
  scale(M, scale = FALSE)
}
