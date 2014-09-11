dsig <- function(a, q) {
    mq = quantile(a, q)
    alpha_l = mq[2] - mq[1]
    alpha_r = mq[3] - mq[2]
    a = a - mq[2]
    lsel = (a < 0.)
    rsel = (a >= 0.)
    a[lsel] = a[lsel] / (-0.5 * alpha_l)
    a[rsel] = a[rsel] / (-0.5 * alpha_r)
    a = (exp(a) + 1)**(-1)
}

.computeBackground <- function(se) {
    mx <- assay(se)
    rs <- rowSums(mx)
    tmp <- lapply(seq(1, ncol(mx)),function(i) {
        imx = mx[rs == i,,drop=FALSE]
        pmin((colSums(imx) + 1) / (nrow(imx) + 1), 1)
    })
    background <- do.call(cbind, tmp)
    return(background)
}

.cooccurenceScoreInner <- function(x, y, xbg, ybg) {
    xrs <- rowSums(assay(x))
    yrs <- rowSums(assay(y))
    ## x and y have to be the same length
    xis <- yis <- probs <- vector(mode="numeric", length=nrow(x))
    for (xi in 1:ncol(xbg)) {
        xbgi <- xbg[,xi,drop=TRUE]
        for (yi in 1:ncol(ybg)) {
            ##
            ybgi <- ybg[,yi,drop=TRUE]
            sel <- ((xrs == xi) & (yrs == yi))
            xs <- assay(x)[sel,,drop=FALSE]
            ys <- assay(y)[sel,,drop=FALSE]
            ##
            p <- matrix(1, nrow=nrow(xs), ncol=ncol(xs))
            p11 <- matrix(rep(xbgi * ybgi, nrow(xs)), nrow=nrow(xs), byrow=TRUE)
            p[xs & ys] <- p11[xs & ys]
            probs[sel] <- rowSums(log10(p))
            ## 
            xis[sel] <- xi
            yis[sel] <- yi
        }
    }
    return(list(xis, yis, probs))
}

.cooccurenceScore <- function(regions, dlinks, step) {
    ## 1. compute background distributions
    pbg <- .computeBackground(regions$promoters)
    ebg <- .computeBackground(regions$enhancers)
    ## 2. compute raw scores
    links <- dlinks$pe
    cs <- bplapply(seq(1, nrow(links), step), function(idx) {
        sel <- idx:min(nrow(links),(idx+step-1))
        .cooccurenceScoreInner(
          x=regions$promoters[  links$queryHits[sel]],
          y=regions$enhancers[links$subjectHits[sel]],
          xbg=pbg, ybg=ebg
        )})
    ## 3. enumerate and convert to data.table
    cs <- data.frame(
      xi=unlist(unname(sapply(cs, "[", 1))),
      yi=unlist(unname(sapply(cs, "[", 2))),
      raw=unlist(unname(sapply(cs, "[", 3)))
    )
    return(cs)
}

.limitByCount <- function(cs, min_count) {
    tmp <- cs[ # tabulate by xi and yi
             ,.N,list(xi, yi)
              ][ # filter cells with less than min_count counts
                N<min_count
                ][ # for each offending cell figure out which margin has a higher count
                 ,list(xi, yi, top=pmax(xi, yi))
                  ][ # ... 
                   ,list(xory=ifelse(xi==top, "xi", "yi"),
                         i=ifelse(xi==top, xi, yi)),
                    ][ # determine the lowest offending cut-off for both margins
                     ,list(i_cut=min(i)), xory
                      ]
    cs[,xi:=pmin(xi, tmp[xory == "xi", i_cut]),]
    cs[,yi:=pmin(yi, tmp[xory == "yi", i_cut]),]
}

.checkQuantile <- function(cs, upper_quantile) {
    min_uq <- max(cs[,list(mu=mean(raw==0)),by=list(xi,yi)]$mu)
    if (min_uq >= upper_quantile) {
        warning("upper quantile to low ... adjusting")
        upper_quantile <- min_uq
    }
    return(upper_quantile)
}

.scaleBySigmoid <- function(cs, upper_quantile, min_count) {
    cs <- data.table(cs)
    cs$idx <- 1:nrow(cs)
    ## 4. determine xi and yi cutoffs
    .limitByCount(cs, min_count)
    ## 5. group by xis yis pairs within each group calculate the
    ## normalized score~sigmoid(-log10(p))
    cs[,score:=2*(dsig(-raw, c(0, 0, upper_quantile))-0.5), by=list(xi, yi)]
    ## 6. sort scores by initial order
    return(cs[order(idx)]$score)
}

.scaleByQuantile <- function(cs, upper_quantile, min_count) {
    cs <- data.table(cs)
    cs$idx <- 1:nrow(cs)
    ## 4. determine xi and yi cutoffs
    .limitByCount(cs, min_count)
    ## 5. upper quantile scaling
    xyuq <- cs[,list(uq=quantile(-raw, upper_quantile)),by=list(xi, yi)]
    setkey(cs, "xi", "yi")
    setkey(xyuq, "xi", "yi")
    cs[xyuq, score:=pmin(-raw/uq, 1)]
    ## 6. sort scores by initial order
    return(cs[order(idx)]$score)
}

cooccurenceLinks <- function(regions, dlinks, step=5e5, method="sigmoid", upper_quantile=0.95, min_count=1000) {
    ## 1-3. calculate raw co-occurence scores
    cs <- .cooccurenceScore(regions, dlinks, step)
    if (method == "sigmoid") {
        cs$score <- .scaleBySigmoid(cs, upper_quantile, min_count)
    } else {
        cs$score <- .scaleByQuantile(cs, upper_quantile, min_count)
    }
    clinks <- cbind(
      dlinks$pe[,c(1,2)],
      cs[,c("raw", "score"),with=FALSE]
    )
    return(cs)
}
