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

#' @title Pathway coverage
#'
#' @description Analyses a pathway by stating which reactions are part of a
#' model and which not. Works only for gapseq reconstructions.
#'
#' @param models Object of class `modelorg` or list of `modelorg` objects.
#' @param reactions data.table or list of data.table for gapseq predictions
#' (files: *-all-Reactions.tbl). See example.
#' @param pathways data.table or list of data.table for gapseq predictions
#' (files: *-all-Pathways.tbl). See example.
#' @param pathways.of.interest character vector of pathway IDs to be tested for
#' coverage. If 'NULL', all pathways in 'pathways' are considered.
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#'
#' @return Returns a data.table. See details
#'
#' @details The output data.table has the columns: pathway (Pathway-ID),
#' rxn.metacyc (Reaction within pathway), rxn.name (Reaction name), spontaneous
#' (TRUE if reaction is non-enzymatic), rxns.models (reaction ids in the model),
#' prediction (TRUE if reaction is present in the model).<br>
#' Please note that he order of objects in models, reactions, and pathways are
#' required to be the same with respect to the models.
#'
#' @examples
#' mods <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", subset = 100)
#' rxns <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/",
#'                                file.type = "reactions", subset = 100)
#' pwys <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/",
#'                                file.type = "pathways", subset = 100)
#'
#' pwys_cov <- get_pathway_coverage(mods, rxns, pwys,
#'                                  pathways.of.interest = c("TRPSYN-PWY",
#'                                                           "TRPSYN-PWY2",
#'                                                           "HISTSYN-PWY"))
#'
#' @import data.table
#' @import sybil
#'
#' @export
get_pathway_coverage <- function(models, reactions, pathways,
                                 pathways.of.interest = NULL,
                                 multi.thread = TRUE) {

  if(class(models) == "modelorg") {
    models    <- list(models)
    reactions <- list(reactions)
    pathways  <- list(pathways)
    names(models) <- models[[1]]@mod_id
  }

  if(length(models) != length(reactions) | length(models) != length(pathways))
    stop("Input objects models/reactions/pathways are of differnt lengths.")

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  n.cores <- min(c(n.cores, length(models)))
  cl <- makeCluster(max(c(1,n.cores)))
  clusterExport(cl, c("pathways.of.interest"), envir=environment())

  # remove leading and trailing pipes in pathway IDs
  reactions <- parLapply(cl, reactions, fun = worker_filter_reactions)
  pathways <- parLapply(cl, pathways, fun = worker_filter_pathways)

  # construct combined object lists
  combObj <- lapply(1:length(models), FUN = function(k) {
    list(mod = models[[k]],
         rxn = reactions[[k]],
         pwy = pathways[[k]])
  })
  names(combObj) <- names(models)

  tmp_cov <- parLapply(cl, combObj, fun = worker_pathway_cov)
  stopCluster(cl)

  tmp_cov <- rbindlist(tmp_cov, idcol = "model")


  return(tmp_cov)
}


#' @import data.table
worker_filter_pathways <- function(x) {
  x[, ID := gsub("^\\||\\|$","", ID)]
  if(!is.null(pathways.of.interest))
    x[ID %in% pathways.of.interest]
  x[]
  return(x)
}

#' @import data.table
worker_filter_reactions <- function(x) {
  x[, pathway := gsub("^\\||\\|$","", pathway)]
  if(!is.null(pathways.of.interest))
    x[pathway %in% pathways.of.interest]
  x[]
  return(x)
}

#' @import sybil
#' @import data.table
worker_pathway_cov <- function(co) {
  mod <- co$mod
  rxn <- co$rxn
  pwy <- co$pwy

  # Pathways of interest
  pwyOI <- pathways.of.interest
  if(is.null(pwyOI))
    pwyOI <- unique(rxn$pathway)
  tmp_pwy <- copy(rxn[pathway %in% pwyOI])
  tmp_pwy <- tmp_pwy[!duplicated(paste(rxn, pathway, sep = "$"))]

  # get which reactions are in the model
  rxns.in.model <- mod@react_id
  rxns.in.model <- rxns.in.model[grepl("^rxn", rxns.in.model)]
  rxns.in.model <- gsub("_.0$","", rxns.in.model)

  # init output table
  out_pwycov <- copy(tmp_pwy[,.(pathway,
                                rxn.metacyc = rxn,
                                rxn.name = name,
                                ec = ec,
                                spontaneous = status == "spontaneous",
                                rxns.model = NA_character_,
                                prediction = FALSE)])

  # checking reaction presence
  for(i in 1:nrow(tmp_pwy)) {
    ms.rxn <- unlist(strsplit(tmp_pwy[i, dbhit], " "))
    ms.rxn <- ms.rxn[ms.rxn %in% rxns.in.model]
    out_pwycov[i, rxns.model := paste(ms.rxn, collapse = " ")]
  }
  out_pwycov[spontaneous == TRUE, prediction := TRUE]
  out_pwycov[spontaneous == FALSE, prediction := ifelse(rxns.model != "" & !is.na(rxns.model),TRUE,FALSE)]
  out_pwycov <- merge(out_pwycov, pwy[, .(pathway = ID, pathway.name = Name)],
                      by = "pathway")


  return(out_pwycov)
}
