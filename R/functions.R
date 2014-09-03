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


medianTSS = function(annotation) {
    ## gene data
    g = annotation[annotation$type == "gene"]
    gd = data.table(data.frame(mcols(g)))
    gd$chr = as.character(seqnames(g))
    gd$strand = as.character(strand(g))
    ## transcript by gene
    tx = annotation[annotation$type == "transcript"]
    #tx = tx[order(tx$gene_id)]# sorted tx
    tss = promoters(tx, 0, 1)[,"gene_id"]
    tbyg = split(tss, tss$gene_id)
    ## 
    tmp = data.table(
        gene_id = names(tbyg),
        start=quantile(start(tbyg), 0.5, type=3))
    tmp = merge(tmp, gd, by="gene_id")
    GRanges(seqnames=tmp$chr, ranges=IRanges(start=tmp$start, end=tmp$start),
            strand=tmp$strand, gene_id=tmp$gene_id, type="transcript")
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
    txs = transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, columns=c("tx_id", "gene_id"))
    tss = promoters(txs, 0, 1)
    
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
