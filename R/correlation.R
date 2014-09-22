.correlationScoreInner <- function(ax, ay, method) {
    ## this all is too slow
    if (method=="cosine") {
        cc <- sapply(1:nrow(ax), function(i) {
            crossprod(ax[i,], ay[1,]) / sqrt(crossprod(ax[1,]) * crossprod(ay[1,]))
        })
    } else {
        cc <- sapply(1:nrow(ax), function(i, ...) {cor(ax[i,], ay[i,], ...)}, method=method)
    }
    return(list(cc))
}

.correlationScore <- function(mx1, mx2, links, method, step) {
    cs <- bplapply(seq(1, nrow(links), step), function(idx) {
        sel <- idx:min(nrow(links),(idx+step-1))
        lsel <- links[sel]
        .correlationScoreInner(
          ax=mx1[lsel[[1]],],
          ay=mx2[lsel[[2]],],
          method=method
        )})
    ## 3. enumerate and convert to data.table
    cs <- data.table(
      raw=unlist(unname(sapply(cs, "[", 1)))
    )
    return(cs)
}

correlationLinks <- function(sites, coverages, links, step=5e5, method="spearman", absolute=TRUE,
                             norm.method="sigmoid",
                             norm.quantiles=c(0.0, 0.0, 0.95),
                             renorm.method="sigmoid",
                             renorm.quantiles=c(0.05, 0.5, 0.95)
                             ) {
    mx <- .numericAssay(coverages)
    mx[is.na(mx)] <- 0.0
    ## 0. normalize the input data
    if (norm.method == "sigmoid") {
        mx <- 2 * (apply(mx, 2, dsig, norm.quantiles) - 0.5)
    } else if (norm.method == "quantile") {
        normalize.quantiles(mx, copy=FALSE) #preprocessCore
    } else if ("qthresh") {
        mx <- apply(mx, 2, qthresh, normalize.quantiles[3])
    }
    sites1 <- mx[assay(sites)[,colnames(links)[1]],]
    sites2 <- mx[assay(sites)[,colnames(links)[2]],]
    ## 1-3. calculate raw co-occurence scores
    cs <- .correlationScore(sites1, sites2, links, method, step)
    ## 4-6. re-normalize the score
    if (renorm.method == "sigmoid") {
        cs[,score:=2*(dsig(raw, renorm.quantiles) - 0.5)]
    } else if (renorm.method == "qthresh") {
        cs[,score:=qthresh(raw, renorm.quantiles[[3]])]
    }
    if (absolute) {
        cs[,score:=abs(score),]
    }
    clinks <- cbind(
      links[,c(1,2), with=FALSE],
      cs[,c("raw", "score"),with=FALSE]
    )
    return(clinks)
}
