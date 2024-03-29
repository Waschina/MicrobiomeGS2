% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fba_tools.R
\name{get_flux_distribution}
\alias{get_flux_distribution}
\title{Get flux distribution}
\usage{
get_flux_distribution(
  mod,
  exclude.unused = T,
  algorithm = "mtf",
  multi.thread = TRUE,
  ncores = NULL
)
}
\arguments{
\item{mod}{Model of class `modelorg` or a named list of these objects.}

\item{exclude.unused}{Exclude reactions with zero fluxes from output table.}

\item{algorithm}{character. Algorithm to use to calculate flux distribution.
Parameter is passed on to \link[sybil]{optimizeProb}.}

\item{multi.thread}{logical. Indicating if parallel processing of models is
used.}

\item{ncores}{integer. Number of CPUs that are used in case of parallel
processing. If NULL, the number of available CPUs is detected.}
}
\value{
A data.table with columns: `rxn` - reaction ID, `flux` - predicted flux, and additional reaction attributes.
}
\description{
Predict all fluxes of the constrained model.
}
