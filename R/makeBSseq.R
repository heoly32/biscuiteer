#' make an in-core BSseq object from a biscuit BED
#'
#' @param tbl         a tibble (from read_tsv) or a data.table (from fread())
#' @param params      parameters (from checkBiscuitBED)
#' @param simplify    simplify sample names by dropping .foo.bar.hg19 & similar
#'
#' @return an in-core BSseq object
#' 
#' @import GenomicRanges
#' @import bsseq 
#'
#' @seealso makeBSseq_HDF5
#'
#' @export 
makeBSseq <- function(tbl, params, hdf5=FALSE, simplify=FALSE) {

  gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1) 
  if (params$how == "data.table") { 
    betaIdx <- match(params$betaCols, names(tbl))
    covgIdx <- match(params$covgCols, names(tbl))
    M <- as.matrix(tbl[, ..betaIdx])
    Cov <- as.matrix(tbl[, ..covgIdx])
  } else { 
    M <- as.matrix(tbl[,betaCols])
    Cov <- as.matrix(tbl[, covgCols])
  }
  
  if(!params$sparse) {
    loci_without_zero_cov <- rowSums(Cov == 0) == 0
    
    M <- M[loci_without_zero_cov, ]
    Cov <- Cov[loci_without_zero_cov, ]
    gr <- gr[loci_without_zero_cov]
  } else {
    M[is.na(M)] <- 0
  }
  
  colnames(M) <- NULL
  colnames(Cov) <- NULL
  
  if(hdf5) {
    hdf5M <- writeHDF5Array(M)
    hdf5Cov <- writeHDF5Array(Cov)
    res <- BSseq(gr=gr, M=hdf5M, Cov=hdf5Cov, sampleNames=rownames(params$pData))
  } else {
    res <- BSseq(gr=gr, M=M, Cov=Cov, sampleNames=rownames(params$pData)) 
  }
  
  if (simplify) res <- simplifySampleNames(res)
  return(res)

}
