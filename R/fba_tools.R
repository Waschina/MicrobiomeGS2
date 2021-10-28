#' @title Predicting auxotrophies
#'
#' @description Performs an FBA on the input model and on the model, where the
#' lower bound (i.e. influx) of the compounds is set each (individually) to 0,
#' and growth rates are compared.
#'
#' @param mod Model of type `modelorg`
#' @param compounds character. IDs of the model compound whose lower bound
#' should be set to 0 in order to test its essentiality. No suffix for the
#' external compartment is required.
#' @param min.growth double. Minimum required growth of the original model.
#' @param min.growth.fraction double. Minimum growth fraction of the model
#' without the compound relative to the original model that is required to call
#' the model prototrophic for the compound.
#'
#' @return Named numeric vector with the same length and order as `compounds`.
#' Entry of 1 indicates prototrophy and 0 auxotrophy.
#'
#' @import sybil
#'
#' @export
predict_auxotrohies <- function(mod, compounds = NULL, min.growth = 0.005,
                                min.growth.fraction = 0.05) {
  if("cplexAPI" %in% rownames(utils::installed.packages())) {
    lpsolver <- "cplexAPI"
    SYBIL_SETTINGS("SOLVER",lpsolver)
    okcode   <- 1
  }

  if(!any(grepl("cpd00002",mod@met_id))) {
    stop("Model is not in modelseed/gapseq namespace. Predictions currently work only for those types.")
  }

  if(is.null(compounds) || compounds[1] == "amino acids") {
    compounds <- c(Ala = "cpd00035",
                   Val = "cpd00156",
                   Met = "cpd00060",
                   Leu = "cpd00107",
                   Ile = "cpd00322",
                   Pro = "cpd00129",
                   Trp = "cpd00065",
                   Phe = "cpd00066",
                   Lys = "cpd00039",
                   Arg = "cpd00051",
                   His = "cpd00119",
                   Tyr = "cpd00069",
                   Thr = "cpd00161",
                   Glu = "cpd00023",
                   Gln = "cpd00053",
                   Gly = "cpd00033",
                   Ser = "cpd00054",
                   Cys = "cpd00084",
                   Asp = "cpd00041",
                   Asn = "cpd00132",
                   Chor = "cpd00216"
                   )
  }

  compounds <- gsub("^EX_|_.0$", "", compounds)
  if(is.null(names(compounds)))
    names(compounds) <- compounds

  # init output
  auxo_out <- rep(NA_real_, length(compounds))
  names(auxo_out) <- names(compounds)

  # get orig growth rate
  m0_growth <- get_growth(mod)
  if(m0_growth < min.growth) {
    warning(paste0("Model ('",mod@mod_id,"') has a too low or zero growth rate."))
    return(auxo_out)
  }

  # checking auxotrophies
  for(i in 1:length(compounds)) {
    imet <- compounds[i]
    ex_id <- paste0("EX_",imet,"_e0")
    if(ex_id %in% mod@react_id) {
      mod_tmp <- changeBounds(mod, ex_id, lb = 0)
      sol_tmp <- optimizeProb(mod_tmp)
      m1_growth <- sol_tmp@lp_obj
      auxo_out[i] <- m1_growth / m0_growth
    } else {
      auxo_out[i] <- 1
    }
  }

  auxo_out <- ifelse(auxo_out >= min.growth.fraction, 1, 0)

  return(auxo_out)
}

#' Get growth rate / Value of objective function
#'
#' @description Uses a simple FBA to predict the models optimal value for the
#' objective functions (most commonly the growth rate).
#'
#' @param mod Model of type `modelorg`
#'
#' @return Numeric. Value of optimal value of objective function. (e.g. growth
#' rate)
#'
#' @export
get_growth <- function(mod) {
  sol <- optimizeProb(mod)

  sol.feasible <- F
  if(sol@solver == "glpkAPI" & sol@lp_stat != 5)
    warning("LP solution infeasible.")
  if(sol@solver == "cplexAPI" & sol@lp_stat != 1)
    warning("LP solution infeasible.")

  return(sol@lp_obj)
}

#' Get produced metabolites
#'
#' @description Predicts metabolite production using MTF-FBA and FVA.
#'
#' @param mod Model object of class `modelorg`
#'
#' @return data.table with predicted production rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `rxn.name` - exchange reaction name, `l` - lower bound of reaction flux, `u` - upper bound, `mtf.flux` - predicted production flux.
#'
#' @details Ranges (`l` & `u`) are based on a FVA-derived method. This method attempts to prevent cases where a
#' nutrient is taken up from the environment, transformed to another compound that is produced, but without any
#' contribution to the organism's growth rate.
#'
#' @export
get_produced_metabolites <- function(mod) {

  # get MTF solution
  #rxn.coef <- ifelse(mod@react_attr$status %in% c("good_blast",NA,"no_seq_data"),1,1)
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])

  # this following two lines are there to prevent the case that a nutrient (e.g. L-Lactate)
  # from the environment is taken up, and thus enables the production of D-Lactate.
  dt.mtf.tmp[mtf.flux > 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, lb = dt.mtf.tmp$mtf.flux)


  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])

  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-6 & l >= 0) | (u > 1)]



  dt <- merge(dt, dt.mtf, by = "ex")

  return(dt[order(-u)])
}

#' Get utilized metabolites
#'
#' @description Predicts metabolite consumption using MTF-FBA and FVA.
#'
#' @param mod Model object of class `modelorg`
#'
#' @return data.table with predicted consumption rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `rxn.name` - exchange reaction name, `l` - lower bound of reaction flux, `u` - upper bound, `mtf.flux` - predicted utilization flux,
#' `flux.at.limit` - character specifying if the uptake is at it's constraint limit. "*" if at limit and "" if not at limit.
#'
#' @export
get_utilized_metabolites <- function(mod) {

  # get MTF solution
  #rxn.coef <- ifelse(mod@react_attr$status %in% c("good_blast",NA,"no_seq_data"),1,1)
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num],
                        lb = mod@lowbnd)
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])
  dt.mtf.tmp[mtf.flux < 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, ub = dt.mtf.tmp$mtf.flux)
  #model.tmp <- mod

  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])

  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(l < -1e-4 & u <= 0) | (l < -1)]


  dt <- merge(dt, dt.mtf, by = "ex")
  dt[, flux.at.limit := ifelse(mtf.flux <= lb*0.999, "*","")]

  return(dt[order(l)])
}



#' Get flux distribution
#'
#' @description Predict all fluxes of the constrained model.
#'
#' @param mod Model of class `modelorg`
#' @param exclude.unused Exclude reactions with zero fluxes from output table.
#' @param algorithm character. Algorith to use to calculate flux distribution.
#' Parameter is passed on to \link[sybil]{optimizeProb}.
#'
#' @return A data.table with columns: `rxn` - reaction ID, `flux` - predicted flux, and additional reaction attributes.
#'
#' @export
get_flux_distribution <- function(mod, exclude.unused = T, algorithm = "mtf") {
  sol <- sybil::optimizeProb(mod, algorithm = algorithm)

  dt <- data.table(rxn = mod@react_id, flux = sol@fluxdist@fluxes[1:mod@react_num])

  dt <- cbind(dt, data.table(mod@react_attr))
  colnames(dt)[duplicated(colnames(dt))] <- paste0(colnames(dt)[duplicated(colnames(dt))],"_2")

  if(exclude.unused)
    dt <- dt[flux != 0]

  dt$equation <- print_reaction(mod, react = dt$rxn)
  dt[, annotation := NULL]
  return(dt)
}

#' Get reduced costs
#'
#' @description Calculated all reactions' reduced costs. For the interpretation
#' of reduced costs values see for instance (https://dx.doi.org/10.1002%2Fbiot.201200291)
#' section 6.
#'
#' @param mod Model of class `modelorg`
#'
#' @return data.table with first columns being the reaction ID and the second column the reduced cost value.
#'
#' @export
get_reduced_costs <- function(mod) {

  mod_warm <- sysBiolAlg(mod)
  mod_sol  <- optimizeProb(mod_warm)

  rc <- getRedCosts(mod_warm@problem)

  out <- data.table(rxn = mod@react_id, red.costs = rc)

  return(out)
}

