#' Get reaction involving a specific metabolite
#'
#' @description Gets all reactions IDs, whose reactions involve one of the
#' metabolites listed in `met`
#'
#' @param mod Model of class `modelorg`
#' @param met Metabolite ID(s) of interest.
#'
#' @return Vector with corresponding reaction IDs.
#'
#' @export
get_reactions_with_metabolite <- function(mod, met) {
  metind <- which(mod@met_id %in% met)

  rxnind <- c()
  for(mi in metind) {
    rxnind <- c(rxnind, which(mod@S[mi,] != 0))
  }


  return(mod@react_id[rxnind])
}

#' Get table of metabolites
#'
#' @description Produces a table with the model's metabolites (ID, name,
#' chemical formula, charge).
#'
#' @param mod Model of class `modelorg`
#'
#' @return data.table with the columns: (`met.id`) Metabolit IDs, (`met.name`)
#' Metabolite names, (`formula`) Chemical formulae of metabolites, (`charge`)
#' Charge of metabolites.
#'
#' @examples
#' eure <- readRDS(system.file("extdata", "eure.RDS", package = "MicrobiomeGS2"))
#' dt <- get_metabolite_table(eure)
#' dt
#'
#' @export
get_metabolite_table <- function(mod) {

  dt <- data.table(met.id   = mod@met_id,
                   met.name = mod@met_name,
                   formula  = mod@met_attr$chemicalFormula,
                   charge   = mod@met_attr$charge)

  return(dt)

}


#' Get table of reactions
#'
#' @description Produces a table with the model's reactions (ID, name,
#' ec, charge).
#'
#' @param mod Model of class `modelorg`
#'
#' @return data.table with the columns: (`react.id`) Reactions IDs, (`react.name`)
#' Reaction names, (`ec`) the reactions' EC identifiers. If models were built
#' with gapseq the additional columns `bitscore`, `pident`, `status`, `gs.origin`
#' are part of the table.
#'
#' @examples
#' eure <- readRDS(system.file("extdata", "eure.RDS", package = "MicrobiomeGS2"))
#' dt <- get_reaction_table(eure)
#' dt
#'
#' @export
get_reaction_table <- function(mod) {

  dt <- data.table(react.id   = mod@react_id,
                   react.name = mod@react_name,
                   ec         = mod@react_attr$ec)

  # in case of gapseq models
  if("gs.origin" %in% colnames(mod@react_attr)) {
    dt$bitscore  <- mod@react_attr$bitscore
    dt$pident    <- mod@react_attr$pident
    dt$status    <- mod@react_attr$status
    dt$gs.origin <- mod@react_attr$gs.origin
  }

  return(dt)

}
