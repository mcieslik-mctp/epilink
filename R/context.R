contextLinks <- function(sites, models, links, p=0.5) {
    ## 1. define 'folds' i.e. regions spanning between
    ## two linked sites
    sites1 <- granges(sites[assay(sites)[,colnames(links)[1]]])
    sites2 <- granges(sites[assay(sites)[,colnames(links)[2]]])
    ls1 <- sites1[links[[colnames(links)[1]]]]
    ls2 <- sites2[links[[colnames(links)[2]]]]
    folds <- punion(ls1, ls2, fill.gap=TRUE)
    ## 2. for each gene define the canonical (median) TSS
    mTss <- .medianTss(models)
    ## 3. Find overlaps between cononical TSS sites and each fold
    hits <- suppressWarnings(data.table(data.frame(findOverlaps(mTss, folds))))
    ## 4. For each fold count the number of covered canonical TSS sites
    ## if the cover is 0 this means that the TSS is outside the promoter
    ## region; we add 1 although this technically might be incorrect.
    tmp <- hits[,.N,by=list(subjectHits)]
    setkey(tmp, "subjectHits")
    counts <- tmp[J(1:nrow(links))]$N
    counts[is.na(counts)] <- 1
    ## 5. compute the normalized score
    clinks <- data.table(
      links[[colnames(links)[1]]],
      links[[colnames(links)[2]]],
      counts,
      2*dgeom(counts-1, p)
      )
    setnames(clinks, names(links))
    return(clinks)
}
