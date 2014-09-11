context("utils tests")

test_that("import list of BED files", {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    fs <- list.files(system.file("inst", "extdata", "chip", package="epilink"), "*bed.gz", full.names=TRUE)
    ##
    models <- prepareGeneModels(txdb)
    sites <- importSites(fs)
    regions <- defineRegions(models, sites)
    ## 
    dlinks <- distanceLinks(models, regions)
    
})
