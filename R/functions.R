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

.medianTss <- function(gs, tss) {
    tssByGene <- split(tss, as.character(tss$gene_id))
    medianTss <- unname(quantile(start(tssByGene), probs=0.5, type=3))
    start(gs) <- medianTss
    end(gs) <- medianTss
    return(gs)
}

.definePromoters <- function(models, regions, probable_distance, alpha) {
    radius <- qcauchy(alpha/2, location=0, scale=probable_distance, lower.tail=FALSE)
    tss <- promoters(models$transcripts, 0, 1)
    hits <- findOverlaps(tss, regions, maxgap=radius)
    tmp <- list(promoters=regions[unique(hits@subjectHits)],
         enhancers=regions[setdiff(1:length(granges(regions)), unique(hits@subjectHits))])
    return(tmp)
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

importRegions <- function(fileNames, sampleData, ...) {
    if (missing(sampleData)) {
        sampleData = .sampleDataFromFileNames(fileNames)
    }
    individual = GRangesList(lapply(fileNames, import))
    merged = reduce(unlist(individual))
    scores = .combineBedScores(individual, merged, ...)
    regions = SummarizedExperiment(scores, rowData=merged, colData=sampleData)
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
    gtdl <- geneTranscriptDistanceLinks(models, probable_distance=gt.dist)
    tmp <- .definePromoters(models, regions, probable_distance=tp.dist, alpha=tp.alpha)
    promoters <- tmp$promoters
    enhancers <- tmp$enhancers
    tpdl <- transcriptPromoterDistanceLinks(models, promoters,
                                            probable_distance=tp.dist, alpha=tp.alpha)
    pedl <- promoterEnhacerDistanceLinks(promoters, enhancers,
                                            probable_distance=pe.dist, alpha=pe.alpha)
    links <- list(gt=gtdl, tp=tpdl, pe=pedl)
    return(links)
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
fs = list.files("inst/extdata/chip", "*bed.gz", full.names=TRUE)
models <- prepareGeneModels(txdb)
regions <- importRegions(fs)
links <- distanceLinks(models, regions)






