% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/join_mult_models.R
\name{join_mult_models}
\alias{join_mult_models}
\title{Join multiple metabolic models}
\usage{
join_mult_models(model.list, scale.boundaries = 1, merge.lb.method = "none")
}
\arguments{
\item{model.list}{A list of metabolic models in 'sybil' format.}

\item{scale.boundaries}{Scaling factor for all reaction bounds. Might be useful for large communities (>20 species). Default: 1}

\item{merge.lb.method}{Method to use in cases where the individual models
differ in the lower bound values for the same exchange reactions. "minimum"
uses the minimum LB values (i.e. largest maximum uptake rate) found among the
individual models for the final joined model. "maximum" analogously. "none"
will cause the method to stop with an error, if lower boundaries differ
between input models. Further options: "median","mean".}
}
\value{
A list with three objects:
  model.IDs - a data.table, which states the newly assigned internal model IDs (M1, M2, ...) to the corresponding organism name.
  modj - A community model also in sybil format. Each organism has its own compartment ('c' and 'e',), but metabolites can freely move between the 'e' compartments of different species.
  ex.rxns - a data.table of all exchange reaction of the new community model.
}
\description{
The function merges multiple metabolic models into one model
where each organisms is within an own compartment. Organism compartments are
connected by a shared extracellular space (i.e. environment), from where the
organisms can take up nutrients an release metabolic by-products. Exchange
reactions are attached to to shared extracelluar space.
}
