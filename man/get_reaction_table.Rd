% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_analysis_tools.R
\name{get_reaction_table}
\alias{get_reaction_table}
\title{Get table of reactions}
\usage{
get_reaction_table(mod)
}
\arguments{
\item{mod}{Model of class `modelorg`}
}
\value{
data.table with the columns: (`react.id`) Reactions IDs, (`react.name`)
Reaction names, (`ec`) the reactions' EC identifiers. If models were built
with gapseq the additional columns `bitscore`, `pident`, `status`, `gs.origin`
are part of the table.
}
\description{
Produces a table with the model's reactions (ID, name,
ec, charge).
}
\examples{
eure <- readRDS(system.file("extdata", "eure.RDS", package = "MicrobiomeGS2"))
dt <- get_reaction_table(eure)
dt

}
