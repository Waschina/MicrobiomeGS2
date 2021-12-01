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
