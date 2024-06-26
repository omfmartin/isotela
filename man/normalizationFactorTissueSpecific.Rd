% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalizationFactorTissueSpecific}
\alias{normalizationFactorTissueSpecific}
\title{Tissue-specific normalization factors}
\usage{
normalizationFactorTissueSpecific(counts, tissues, ...)
}
\arguments{
\item{counts}{RNA-seq biological count \code{matrix} (genes as rows, samples as columns)}

\item{tissues}{\code{vector} with tissues as values and genes as names (have to match \code{rownames(counts)})}

\item{...}{Furthur arguments passed to \link[scran]{computeSumFactors}}
}
\value{
\code{list} of \code{numeric} \code{vector} of tissue-specific normalization factors
}
\description{
Estimate normalization factors using \link[scran]{computeSumFactors}
separately for tissue-unique genes
}
