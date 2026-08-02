#' Logistic Function
#'
#' Compute the logistic transformation of a numeric input.
#'
#' @param x A numeric value or vector.
#'
#' @return A numeric value or vector with elements in (0,1).
#'
#' @examples
#' logistic(0)
#' logistic(c(-1, 0, 1))
#'
#' @export
logistic<-function(x){
  return(1/(1+exp(-x)))
  } # define lositic function
