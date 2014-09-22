.computeBackground <- function(mx) {
    rs <- rowSums(mx)
    tmp <- lapply(seq(1, ncol(mx)),function(i) {
        imx = mx[rs == i,,drop=FALSE]
        pmin((colSums(imx) + 1) / (nrow(imx) + 1), 1)
    })
    background <- do.call(cbind, tmp)
    return(background)
}

.cooccurenceScoreInner <- function(ax, ay, xbg, ybg) {
    xrs <- rowSums(ax)
    yrs <- rowSums(ay)
    ## x and y have to be the same length
    xis <- yis <- probs <- vector(mode="numeric", length=nrow(x))
    for (xi in 1:ncol(xbg)) {
        xbgi <- xbg[,xi,drop=TRUE]
        for (yi in 1:ncol(ybg)) {
            ##
            ybgi <- ybg[,yi,drop=TRUE]
            sel <- ((xrs == xi) & (yrs == yi))
            xs <- ax[sel,,drop=FALSE]
            ys <- ay[sel,,drop=FALSE]
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

.cooccurenceScore <- function(mx1, mx2, links, step) {
    ## 1. compute background distributions
    bg1 <- .computeBackground(mx1)
    bg2 <- .computeBackground(mx2)
    ## 2. compute raw scores
    cs <- bplapply(seq(1, nrow(links), step), function(idx) {
        sel <- idx:min(nrow(links),(idx+step-1))
        lsel <- links[sel]
        .cooccurenceScoreInner(
          ax=mx1[lsel[[1]],],
          ay=mx2[lsel[[2]],],
          xbg=bg1, ybg=bg2
        )})
    ## 3. enumerate and convert to data.table
    cs <- data.table(
      xi=unlist(unname(sapply(cs, "[", 1))),
      yi=unlist(unname(sapply(cs, "[", 2))),
      raw=unlist(unname(sapply(cs, "[", 3)))
    )
    return(cs)
}

.limitByCount <- function(cs, min.count) {
    tmp <- cs[ # tabulate by xi and yi
             ,.N,list(xi, yi)
              ][ # filter cells with less than min.count counts
                N<min.count
                ][ # for each offending cell figure out which margin has a higher count
                 ,list(xi, yi, top=pmax(xi, yi))
                  ][ # ... 
                   ,list(xory=ifelse(xi==top, "xi", "yi"),
                         i=ifelse(xi==top, xi, yi)),
                    ][ # determine the lowest offending cut-off for both margins
                     ,list(icut=min(i)), xory
                      ]
    cs[,xi:=pmin(xi, tmp[xory == "xi", icut]),]
    cs[,yi:=pmin(yi, tmp[xory == "yi", icut]),]
}

.scaleBySigmoid <- function(cs, upper.quantile, min.count) {
    cs <- data.table(cs)
    cs$idx <- 1:nrow(cs)
    ## 4. determine xi and yi cutoffs
    .limitByCount(cs, min.count)
    ## 5. group by xis yis pairs within each group calculate the
    ## normalized score~sigmoid(-log10(p))
    cs[,score:=2*(dsig(-raw, c(0, 0, upper.quantile))-0.5), by=list(xi, yi)]
    ## 6. sort scores by initial order
    return(cs[order(idx)]$score)
}

.scaleByQuantile <- function(cs, upper.quantile, min.count) {
    cs <- data.table(cs)
    cs$idx <- 1:nrow(cs)
    ## 4. determine xi and yi cutoffs
    .limitByCount(cs, min.count)
    ## 5. upper quantile scaling
    xyuq <- cs[,score:=qthresh(-raw, upper.quantile), by=list(xi, yi)]
    ## 6. sort scores by initial order
    return(cs[order(idx)]$score)
}

cooccurenceLinks <- function(sites, overlaps, links, step=5e5, method="sigmoid", upper.quantile=0.95,
                             min.count=1000) {
    
    mx <- .numericAssay(overlaps)
    mx1 <- mx[assay(sites)[,colnames(links)[1]],]
    mx2 <- mx[assay(sites)[,colnames(links)[2]],]
    
    ## 1-3. calculate raw co-occurence scores
    cs <- .cooccurenceScore(mx1, mx2, links, step)
    ## 4-6. normalize the score
    if (method == "sigmoid") {
        cs$score <- .scaleBySigmoid(cs, upper.quantile, min.count)
    } else {
        cs$score <- .scaleByQuantile(cs, upper.quantile, min.count)
    }
    clinks <- cbind(
      links[,c(1,2), with=FALSE],
      cs[,c("raw", "score"),with=FALSE]
    )
    return(clinks)
}
