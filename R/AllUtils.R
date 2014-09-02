#' Check whether all elements of a numeric vector are equal within machine precision
#' @param x a numeric vector.
#' @param tol tolerance.
#' 
#' @return TRUE if all elements of the vector are identical (within machine 
#' precision). FALSE in all other cases, including if the vector contains any 
#' NAs.
#' 
#' @export
#' 
#' @note This function is based on Hadley and John's answer to 
#' http://stackoverflow.com/q/4752275
almostEqual <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) {
    val <- TRUE
  } 
  if (any(is.na(x)) & any(is.na(x))){
    val <- FALSE
  } else{
    val <- (abs(max(x) - min(x)) < tol)
  }
  return(val)
}
