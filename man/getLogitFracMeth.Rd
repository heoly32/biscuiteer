% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLogitFracMeth.R
\name{getLogitFracMeth}
\alias{getLogitFracMeth}
\alias{getMvals}
\title{helper function for compartment inference (shrink-by-smoothing logit frac5mC)}
\usage{
getLogitFracMeth(x, minCov = 3, minSamp = 2, k = 0.1, r = NULL)
}
\arguments{
\item{x}{a BSseq object with methylated and total reads}

\item{minCov}{minimum read coverage for landmarking samples (3)}

\item{minSamp}{minimum landmark samples with >= minCov (2)}

\item{k}{pseudoreads for smoothing (0.1)}

\item{r}{regions to collapse over (default is NULL, do it by CpG)}
}
\value{
\preformatted{       smoothed logit(M/Cov) matrix with coordinates as row names
}
}
\description{
We want something with nominally Gaussian error for compartment inference, so
this function grabs suitable (default >= 3 reads in >=2 sample) measurements
and turns them into lightly moderated, logit-transformed methylated-fraction
estimates (also known, unfortunately, as M-values) for compartment calling,
by performing Dirichlet smoothing (adding \code{k} reads to M and U support).
}
