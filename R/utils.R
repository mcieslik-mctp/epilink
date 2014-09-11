.sampleDataFromFileNames <- function(fs) {
    colData <- DataFrame(path=fs,
      row.names=sapply(strsplit(basename(fs), "\\."), "[", 1))
    return(colData)
}

.combineBedScores <- function(individual, merged, scoreColumn=NULL, decreasing=TRUE) {
    hitsScores <- bplapply(individual, function(x) {
        hits <- as.data.frame(findOverlaps(merged, x))
        if (!is.null(scoreColumn)) {
            scores <- mcols(x)[hits[,"subjectHits"], scoreColumn]
            ## hits sorted from high to low scores by default
            hits <- hits[order(scores, decreasing=decreasing),]
            ## transfer score from individual
            hits$score <- as.integer(mcols(x)[hits[,"subjectHits"], scoreColumn])
        } else {
            hits$score = 1L
        }
        ## keep only one hit per merged region (the one with the highest score)
        hitsBst <- hits[!duplicated(hits[,"queryHits"]),]
        return(hitsBst)
    })
    scores = matrix(0L, ncol=length(individual), nrow=length(merged))
    for (i in seq_along(hitsScores)) {
        tmp = hitsScores[[i]]
        scores[tmp[,"queryHits"],i] = tmp[,"score"]
    }
    return(scores)
}

importSites <- function(fileNames, sampleData, ...) {
    if (missing(sampleData)) {
        sampleData = .sampleDataFromFileNames(fileNames)
    }
    individual = GRangesList(lapply(fileNames, import))
    merged = reduce(unlist(individual))
    scores = .combineBedScores(individual, merged, ...)
    sites = SummarizedExperiment(scores, rowData=merged, colData=sampleData)
    return(sites)
}

defineRegions <- function(models, sites, tp.dist=375, tp.alpha=0.05) {
    radius <- qcauchy(tp.alpha/2, location=0, scale=tp.dist, lower.tail=FALSE)
    tss <- promoters(models$transcripts, 0, 1)
    hits <- findOverlaps(tss, sites, maxgap=radius)
    regions <- list(
         promoters=sites[unique(hits@subjectHits)],
         enhancers=sites[setdiff(1:length(granges(sites)), unique(hits@subjectHits))])
    return(regions)
}

prepareGeneModels <- function(txdb) {
    ## we loose genes with no gene_id
    gs <- genes(txdb)
    ts <- transcripts(txdb, columns=c("tx_id", "gene_id"))
    tsByGs <- split(ts, as.character(ts$gene_id))
    ## the sort is required ...
    ## fortunately this is TRUE
    ## all(unlist(split(unlist(ts), as.character(unlist(ts)$gene_id))) == unlist(ts))
    ts <- unlist(tsByGs[names(gs)])
    names(ts) <- ts$tx_id
    return(list(genes=gs, transcripts=ts))
}
