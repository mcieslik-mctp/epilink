.distanceLinks <- function(x, y, probable.distance, alpha, col.names) {
    radius <- qcauchy(alpha/2, location=0, scale=probable.distance, lower.tail=FALSE)
    if (missing(y)) {
        y <- x
        hits <- findOverlaps(x, ignoreSelf=TRUE, ignoreRedundant=TRUE, maxgap=radius)
    } else {
        hits <- findOverlaps(x, y, maxgap=radius)
    }
    links <- as.data.frame(hits)
    links$raw <- unlist(
        bplapply(seq(1, nrow(links), 1e6), function(idx) {
            distance(
              x[links$queryHits[idx:min(nrow(links),(idx+1e6-1))]],
              y[links$subjectHits[idx:min(nrow(links),(idx+1e6-1))]]
            )}))
    links$score <- 2*pcauchy(-links$raw, 0, probable.distance)
    links <- data.table(links)
    setnames(links, c(col.names, "raw", "score"))
    return(links)
}

geneTranscriptDistanceLinks <- function(models, probable.distance=50000) {
    gs <- models$genes
    ts <- models$transcripts
    tss <- promoters(ts, 0, 1)
    mTss <- .medianTss(models)
    dist <- distance(mTss[unlist(tss$gene_id)], tss)
    links <- data.frame(
      queryHits=match(unlist(tss$gene_id), names(gs)),
      subjectHits=1:length(tss),
      raw=dist,
      score=2*pcauchy(-dist, 0, probable.distance))
    links <- data.table(links)
    setnames(links, c("gene", "transcript", "raw", "score"))
    return(links)
}

transcriptPromoterDistanceLinks <- function(models, sites, probable.distance=375, alpha=0.05) {
    col.names = c("transcript", "promoter")
    x <- promoters(models$transcripts, 0, 1)
    y <- rowData(sites[assay(sites)[,col.names[2]]])
    links <- .distanceLinks(x, y, probable.distance=probable.distance, alpha=alpha, col.names=col.names)
    return(links)
}

promoterEnhacerDistanceLinks <- function(sites, probable.distance=50000, alpha=0.05) {
    col.names = c("promoter", "enhancer")
    x <- rowData(sites[assay(sites)[,col.names[1]]])
    y <- rowData(sites[assay(sites)[,col.names[2]]])
    links <- .distanceLinks(x, y, probable.distance=probable.distance, alpha=alpha, col.names=col.names)
    return(links)
}
