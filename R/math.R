dsig <- function(a, q) {
    mq <- quantile(a, q)
    alpha_l <- mq[2] - mq[1]
    alpha_r <- mq[3] - mq[2]
    a <- a - mq[2]
    lsel <- (a < 0.)
    rsel <- (a >= 0.)
    a[lsel] <- a[lsel] / (-0.5 * (alpha_l + .Machine$double.eps))
    a[rsel] <- a[rsel] / (-0.5 * (alpha_r - .Machine$double.eps))
    a <- (exp(a) + 1)**(-1)
    return(a)
}

qthresh <- function(a, q) {
    mq <- quantile(a, q)
    a <- pmin(a/(mq + .Machine$double.eps), 1)
    return(a)
}
