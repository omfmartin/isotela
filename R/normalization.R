#' @title Standard normalization factors
#'
#' @description Estimate standard normalization factors using \link[scran]{computeSumFactors}
#'
#' @inheritParams tissueSpecificRegression
#' @param ... Furthur arguments passed to \link[scran]{computeSumFactors}
#'
#' @return \code{numeric} \code{vector} of normalization factors
#'
#' @export
normalizationFactor <- function(counts, ...) {
  nf <- scran::calculateSumFactors(counts, ...)
  names(nf) <- colnames(counts)
  nf
}

#' @title Tissue-specific normalization factors
#'
#' @description Estimate normalization factors using \link[scran]{computeSumFactors}
#' separately for tissue-unique genes
#'
#' @inheritParams tissueSpecificRegression
#' @param ... Furthur arguments passed to \link[scran]{computeSumFactors}
#'
#' @return \code{list} of \code{numeric} \code{vector} of tissue-specific normalization factors
#'
#' @export
normalizationFactorTissueSpecific <- function(counts, tissues, ...) {
  uniq_tissues <- sort(unique(tissues))
  nf <- lapply(uniq_tissues, function(tissue) {
    genes <- names(tissues)[tissues == tissue]
    tissue_counts <- counts[rownames(counts) %in% genes, ]
    normalizationFactor(tissue_counts, ...)
  })
  names(nf) <- uniq_tissues
  nf
}

#' @title Standard normalization
#'
#' @description Normalize counts using standard normalization
#'
#' @inheritParams tissueSpecificRegression
#' @param nf Standard normalization factors (output of \link[isotela]{normalizationFactor})
#'
#' @return Normalized read counts
#'
#' @export
normalizeCounts <- function(counts, nf) {
  sweep(counts, 2, nf, "/")
}

#' @title Tissue-specific normalization
#'
#' @description Normalize counts using tissue-specific normalization
#'
#' @inheritParams tissueSpecificRegression
#' @param nf_list Standard normalization factors (output of \link[isotela]{normalizationFactorTissueSpecific})
#'
#' @return Tissue-normalized read counts of tissue-unique genes
#'
#' @export
normalizeCountsTissueSpecific <- function(counts, nf_list, tissues) {
  uniq_tissues <- names(nf_list)
  out <- lapply(uniq_tissues, function(tissue) {
    genes <- names(tissues)[tissues == tissue]
    tissue_counts <- counts[rownames(counts) %in% genes, ]
    normalizeCounts(counts = tissue_counts, nf = nf_list[[tissue]])
  })
  out <- Reduce(rbind, out)
  out
}
