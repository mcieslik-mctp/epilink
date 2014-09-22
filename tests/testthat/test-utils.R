context("utils tests")

GT.DIST <- 50000

test_that("import list of BED files", {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    models <- prepareGeneModels(TxDb.Hsapiens.UCSC.hg19.knownGene)
    ##
    siteFNs <- list.files(system.file("inst", "extdata", "dhs", package="epilink"), "*bed.gz", full.names=TRUE)
    sites <- importSites(siteFNs, models)
    sites21 <- sites[seqnames(sites) == "chr21",]
    ##
    peakFNs <- list.files(system.file("inst", "extdata", "chip", package="epilink"), "*bed.gz", full.names=TRUE)
    overlaps <- sitesOverlap(sites21, peakFNs)
    ##
    wigFNs <- list.files(system.file("inst", "extdata", "hist", package="epilink"), "*bigWig", full.names=TRUE)
    coverages_sim <- sitesCoverage(sites21, wigFNs)
    ##
    bamFNs <- list.files(system.file("inst", "extdata", "hist", package="epilink"), "*bam$", full.names=TRUE)
    counts <- sitesCounts(sites21, bamFNs)
    ##
    wigFNs <- list.files(system.file("inst", "extdata", "H3K27ac", package="epilink"), "*bigWig", full.names=TRUE)
    coverages_cor <- sitesCoverage(sites21, wigFNs)

    corr
    ##
    gt_links <- geneTranscriptDistanceLinks(models)
    links <- promoterEnhacerDistanceLinks(sites21)

    
    links <- pe_links
    sites <- sites21
    
    
    
    
})


tbl <- data.frame(
  size=width(sites21),
  cov=unlist(assay(coverages)[,"H3K9ac"]),
  cnt=unlist(assay(counts)[,"H3K9ac"]))


plt <- ggplot(tbl) + aes(x=cov, y=cnt/size) + geom_point() + theme_bw()
ggsave("t.png", plt)
