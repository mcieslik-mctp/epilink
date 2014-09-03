#' epilink
#' Infer epigenetic links between genes and their cis-regulatory regions from epigenetic data.
#' 
#' @name epilink
#' @docType package
#' @import methods GenomicRanges rtracklayer
library(methods)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
