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
allAlmostEqual <- function(x, tol = .Machine$double.eps ^ 0.5) {
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



mx = matrix(1:100, ncol=10)
rg = GRanges(seqnames="chr1", IRanges(start=1:10, end=11:20), strand="+")
df = DataFrame(blah=paste0("d",LETTERS[1:10]), row.names=LETTERS[1:10])
## colnames(mx) = paste0("col", c(1:10))
## rownames(mx) = paste0("row", c(1:10))
se = SummarizedExperiment(mx,
                     rowData=rg,
                     colData=df
                     )
