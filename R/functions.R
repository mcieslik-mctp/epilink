combineBedScores <- function(individual, merged, scoreColumn=NULL, decreasing=TRUE) {
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
}

.prepareGenesAndTranscripts <- function(txdb) {
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

.medianTss <- function(gs, tss) {
    tssByGene <- split(tss, as.character(tss$gene_id))
    medianTss <- unname(quantile(start(tssByGene), probs=0.5, type=3))
    start(gs) <- medianTss
    end(gs) <- medianTss
    return(gs)
}

tssLinks <- function(txdb, probable_distance=50000) {
    tmp <- .prepareGenesAndTranscripts(txdb)
    gs <- tmp$genes
    ts <- tmp$transcripts
    tss <- promoters(ts, 0, 1)
    mTss <- .medianTss(gs, tss)
    dist <- distance(mTss[unlist(tss$gene_id)], tss)
    links = data.frame(
      queryHits=match(unlist(tss$gene_id), names(gs)),
      subjectHits=1:length(tss),
      raw=dist,
      score=2*pcauchy(-dist, 0, probable_distance)
      )
    return(links)
}





distanceLinks <- function(x, y, probable_distance, alpha) {
    radius <- qcauchy(alpha/2, location=0, scale=probable_distance, lower.tail=FALSE)
    if (missing(y)) {
        y = x
        ## unnecessary, bug in bioC
        hits = findOverlaps(granges(x), ignoreSelf=TRUE, ignoreRedundant=TRUE,
          maxgap=radius)
    } else {
        hits = findOverlaps(x, y, maxgap=radius)
    }
    links <- as.data.frame(hits)
    links$dist <- unlist(
        bplapply(seq(1, nrow(links), 1e6), function(idx) {
            distance(
              x[links$queryHits[idx:min(nrow(links),(idx+1e6-1))]],
              y[links$subjectHits[idx:min(nrow(links),(idx+1e6-1))]]
            )}))
    links$p=2*pcauchy(-links$dist, 0, probable_distance)
    return(links)
}


promoterDistanceLinks <- function(regions, txdb, probable_distance=375, alpha=0.05) {
    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
    
    dl = distanceLinks(regions, tss, probable_distance=probable_distance, alpha=alpha)
    return(dl)
}

regionDistanceLinks <- function(regions, probable_distance=50000, alpha=0.05) {
    dl = distanceLinks(regions, probable_distance=probable_distance, alpha=alpha)
    return(dl)
}


                                                                                          
sampleDataFromFileNames <- function(fs) {
    colData <- DataFrame(path=fs,
      row.names=sapply(strsplit(basename(fs), "\\."), "[", 1))
    return(colData)
}


fs = list.files("inst/extdata/chip", "*bed.gz", full.names=TRUE)

importBed = function(fileNames, sampleData, ...) {
    if (missing(sampleData)) {
        sampleData = sampleDataFromFileNames(fileNames)
    }
    individual = GRangesList(lapply(fileNames, import))
    merged = reduce(unlist(individual))
    scores = combineBedScores(individual, merged, ...)
    regions = SummarizedExperiment(scores, rowData=merged, colData=sampleData)
    return(regions)
}

    if (missing(txdb)) {
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    }


distanceLinks = function(q, s, probable_distance, alpha) {
    radius = qcauchy(alpha/2, location=0, scale=probable_distance, lower.tail=FALSE)
    hits = findOverlaps(q, s, maxgap=radius)
    nhits = length(hits)
    step = 1000000
    dl = data.table(
        q_idx = hits@queryHits,
        s_idx = hits@subjectHits,
        dist=dists)
    setkey(dl, "q_idx", "s_idx")
    dl[,score:=,]
    return(dl)
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
