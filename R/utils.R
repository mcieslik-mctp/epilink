.medianTss <- function(models) {
    gs <- models$genes
    ts <- models$transcripts
    tss <- promoters(ts, 0, 1)
    tssByGene <- split(tss, as.character(tss$gene_id))
    medianTss <- unname(quantile(start(tssByGene), probs=0.5, type=3))
    start(gs[names(tssByGene)]) <- medianTss
    end(gs[names(tssByGene)]) <- medianTss
    return(gs)
}

.numericAssay <- function(se) {
    tmp <- assay(se)
    mx <- matrix(unlist(tmp), nrow=nrow(tmp), ncol=ncol(tmp),
           dimnames=list(rownames(tmp), colnames(tmp)))
    return(mx)
}

.sampleDataFromFileNames <- function(fs) {
    colData <- DataFrame(path=fs,
                         row.names=sapply(strsplit(basename(fs), "\\."), "[", 1))
    return(colData)
}

.combineBedScores <- function(sites, peaks, scoreColumn=NULL, decreasing=TRUE) {
    hitsScores <- bplapply(peaks, function(x) {
        hits <- as.data.frame(findOverlaps(sites, x))
        if (!is.null(scoreColumn)) {
            scores <- mcols(x)[hits[,"subjectHits"], scoreColumn]
            ## hits sorted from high to low scores by default
            hits <- hits[order(scores, decreasing=decreasing),]
            ## transfer score from x peak
            hits$score <- as.integer(mcols(x)[hits[,"subjectHits"], scoreColumn])
        } else {
            hits$score = rep(1L, nrow(hits))
        }
        ## keep only one hit per site (the one with the highest score)
        ## one peak can overlap multiple sites...
        hitsBst <- hits[!duplicated(hits[,"queryHits"]),]
        return(hitsBst)
    })
    scores = matrix(0L, ncol=length(peaks), nrow=nrow(sites))
    for (i in seq_along(hitsScores)) {
        tmp = hitsScores[[i]]
        scores[tmp[,"queryHits"],i] = tmp[,"score"]
    }
    colnames(scores) <- names(peaks)
    return(scores)
}

.importPeaks <- function(fileNames, sites, peakData) {
    peaks <- GRangesList(lapply(fileNames, import))
    seqinfo(peaks) <- merge(seqinfo(peaks), seqinfo(sites))
    peaks <- sortSeqlevels(peaks)
    if (missing(peakData)) {
        peakData <- .sampleDataFromFileNames(fileNames)
    }
    names(peaks) <- row.names(peakData)
    return(peaks)
}

.defineRegions <- function(sites, models, tp.dist=375, tp.alpha=0.05) {
    radius <- qcauchy(tp.alpha/2, location=0, scale=tp.dist, lower.tail=FALSE)
    tss <- promoters(models$transcripts, 0, 1)
    hits <- findOverlaps(tss, sites, maxgap=radius)
    is.promoter = 1:length(sites) %in% unique(hits@subjectHits)
    is.enhancer = !is.promoter
    type <- do.call(cbind, list(
      promoter=is.promoter,
      enhancer=is.enhancer))
    sites <- SummarizedExperiment(SimpleList(type=type), rowData=sites)
    return(sites)
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
    models <- list(genes=sort(sortSeqlevels(gs)), transcripts=sort(sortSeqlevels(ts)))
    return(models)
}

importSites <- function(bedFiles, models){
    tmp <- GRangesList(lapply(bedFiles, import))
    sites <- reduce(unlist(tmp))
    seqinfo(sites) <- merge(seqinfo(sites), seqinfo(models$genes))
    sites <- .defineRegions(sortSeqlevels(sites), models)
    return(sites)
}

sitesOverlap <- function(sites, bedFiles, bedData, ...) {
    peaks <- .importPeaks(bedFiles, sites)
    scores <- .combineBedScores(sites, peaks, ...)
    if (!missing(bedData)) {
        overlaps <- SummarizedExperiment(scores, rowData=rowData(sites), exptData=exptData(sites), colData=bedData)
    } else {
        overlaps <- SummarizedExperiment(scores, rowData=rowData(sites), exptData=exptData(sites))
    }
    return(overlaps)
}

sitesCoverage <- function(sites, wigFiles, wigData, type="mean", ...) {
    levels <- summary(BigWigFileViews(wigFiles, fileRange=rowData(sites)), type=type, ...)
    exptData(levels) <- exptData(sites)
    if (!missing(wigData)) {
        colData(levels) <- wigData
    } else {
        colData(levels) <- .sampleDataFromFileNames(wigFiles)
    }
    return(levels)
}

sitesCounts <- function(sites, bamFiles, bamData, ...) {
    counts <- summarizeOverlaps(rowData(sites), BamFileViews(bamFiles), ...)
    exptData(counts) <- exptData(sites)
    if (!missing(bamData)) {
        colData(counts) <- bamData
    } else {
        colData(counts) <- .sampleDataFromFileNames(bamFiles)
    }
    return(counts)
}
