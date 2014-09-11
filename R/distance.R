.medianTss <- function(gs, tss) {
    tssByGene <- split(tss, as.character(tss$gene_id))
    medianTss <- unname(quantile(start(tssByGene), probs=0.5, type=3))
    start(gs) <- medianTss
    end(gs) <- medianTss
    return(gs)
}


.distanceLinks <- function(x, y, probable_distance, alpha) {
    radius <- qcauchy(alpha/2, location=0, scale=probable_distance, lower.tail=FALSE)
    if (missing(y)) {
        y = x
        hits = findOverlaps(x, ignoreSelf=TRUE, ignoreRedundant=TRUE, maxgap=radius)
    } else {
        hits = findOverlaps(x, y, maxgap=radius)
    }
    links <- as.data.frame(hits)
    links$raw <- unlist(
        bplapply(seq(1, nrow(links), 1e6), function(idx) {
            distance(
              x[links$queryHits[idx:min(nrow(links),(idx+1e6-1))]],
              y[links$subjectHits[idx:min(nrow(links),(idx+1e6-1))]]
            )}))
    links$score=2*pcauchy(-links$raw, 0, probable_distance)
    return(links)
}

geneTranscriptDistanceLinks <- function(models, probable_distance=50000) {
    gs <- models$genes
    ts <- models$transcripts
    tss <- promoters(ts, 0, 1)
    mTss <- .medianTss(gs, tss)
    dist <- distance(mTss[unlist(tss$gene_id)], tss)
    links <- data.frame(
      queryHits=match(unlist(tss$gene_id), names(gs)),
      subjectHits=1:length(tss),
      raw=dist,
      score=2*pcauchy(-dist, 0, probable_distance))
    return(links)
}

transcriptPromoterDistanceLinks <- function(models, promoters, probable_distance, alpha) {
    x <- promoters(models$transcripts, 0, 1)
    y <- granges(promoters)
    links <- .distanceLinks(x, y, probable_distance=probable_distance, alpha=alpha)
    return(links)
}

promoterEnhacerDistanceLinks <- function(promoters, enhancers, probable_distance, alpha) {
    x <- granges(promoters)
    y <- granges(enhancers)
    links <- .distanceLinks(x, y, probable_distance=probable_distance, alpha=alpha)
    return(links)
}

distanceLinks <- function(models, regions,
                          gt.dist=50000,
                          tp.dist=375,
                          tp.alpha=0.05,
                          pe.dist=50000,
                          pe.alpha=0.05) {
    promoters <- regions$promoters
    enhancers <- regions$enhancers
    gtdl <- geneTranscriptDistanceLinks(models, probable_distance=gt.dist)
    tpdl <- transcriptPromoterDistanceLinks(models, promoters,
                                            probable_distance=tp.dist, alpha=tp.alpha)
    pedl <- promoterEnhacerDistanceLinks(promoters, enhancers,
                                            probable_distance=pe.dist, alpha=pe.alpha)
    links <- list(gt=gtdl, tp=tpdl, pe=pedl)
    return(links)
}
