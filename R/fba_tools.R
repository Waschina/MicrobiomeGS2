#' @title Predicting auxotrophies
#'
#' @description Performs a FBA on the input model and on the model, where the
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
predict_auxotrophies <- function(mod, compounds = NULL, min.growth = 0.005,
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
#' @param easyConstraints Additional constraints. See easyConstraints in the sybil
#' documentation.
#'
#' @return data.table with predicted production rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `rxn.name` - exchange reaction name, `l` - lower bound of reaction flux, `u` - upper bound, `mtf.flux` - predicted production flux.
#'
#' @details Ranges (`l` & `u`) are based on a FVA-derived method. This method attempts to prevent cases where a
#' nutrient is taken up from the environment, transformed to another compound that is produced, but without any
#' contribution to the organism's growth rate. When providing additional
#' constraints, please note that these constraints are not (yet) applied for the
#' flux variability analysis.
#'
#' @export
get_produced_metabolites <- function(mod, easyConstraints = NULL) {

  # get MTF solution
  if(!is.null(easyConstraints)) {
    mod_warm <- sysBiolAlg(mod,
                           algorithm = "mtfEasyConstraint",
                           easyConstraint = easyConstraints)
  } else {
    mod_warm <- sysBiolAlg(mod,
                           algorithm = "mtf")
  }

  sol.mtf <- optimizeProb(mod_warm)
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf$fluxes[1:mod@react_num])
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
#' @param easyConstraints Additional constraints. See easyConstraints in the sybil
#' documentation.
#'
#' @return data.table with predicted consumption rate and flux ranges. Columns: `ex` - exchange reaction ID,
#' `rxn.name` - exchange reaction name, `l` - lower bound of reaction flux, `u` - upper bound, `mtf.flux` - predicted utilization flux,
#' `flux.at.limit` - character specifying if the uptake is at it's constraint limit. "*" if at limit and "" if not at limit.
#' When providing additional constraints, please note that these constraints are
#' not (yet) applied for the flux variability analysis.
#' @export
get_utilized_metabolites <- function(mod, easyConstraints = NULL) {

  # get MTF solution
  if(!is.null(easyConstraints)) {
    mod_warm <- sysBiolAlg(mod,
                           algorithm = "mtfEasyConstraint",
                           easyConstraint = easyConstraints)
  } else {
    mod_warm <- sysBiolAlg(mod,
                           algorithm = "mtf")
  }

  sol.mtf <- optimizeProb(mod_warm)
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf$fluxes[1:mod@react_num],
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
#' @param algorithm character. Algorithm to use to calculate flux distribution.
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

#' Prediction of maximum production rates
#'
#' @description Gets the maximum production capacity of metabolite \code{met}.
#' Implementation: A new reaction with the outflow of the metabolite of interest
#' is introduced and it's flux maximized
#' using FBA.
#'
#' @param mod Model of class \code{modelorg}
#' @param met Character vector of metabolite IDs, whose production capacity is predicted.
#'
#' @return A named numeric vector with individual maximum production rates. Names correspond to metabolite IDs.
#'
#' @export
get_metabolite_production_capacity <- function(mod, met) {
  met_tmp <- met[met %in% mod@met_id]

  if(length(met_tmp) == 0) {
    stop("None of the metabolite IDs is part of the model. Returning empty vector.")
    return(numeric(0L))
  }

  if(length(met_tmp) != met) {
    lost_mets <- met[!(met %in% met_tmp)]
    warning(paste("Follwing metabolite IDs are not part of the model:\n",paste(lost_mets, collapse = ", ")))
  }

  met <- met_tmp

  res <- c()
  for(moi in met) {
    mod.tmp <- addReact(mod,
                        id = "OUTFLOW_TMP",
                        met = met,
                        Scoef = -1,
                        reversible = F)

    mod.tmp <- changeObjFunc(mod.tmp, react = "OUTFLOW_TMP")

    sol <- optimizeProb(mod.tmp)
    res <- c(res, sol@lp_obj)
  }

  names(res) <- met

  return(res)
}

#' Get essential nutrients
#'
#' @description Predicts which metabolites are essential for the model's growth.
#'
#' @param mod Model of class `modelorg`
#' @param min.growth.fraction Minimum required predicted growth without metabolite relative to the predicted growth with the focal metabolite, to consider the metabolite as non-essential. Default: 0.01 (i.e. 1\% of original growth)
#'
#' @return data.table with the columns: 1 - ID of exchange reaction for the focal nutrient metabolite, 2 - Name of the exchange reaction, 3 - predicted growth rate w/o metabolite, 4 - Predicted ration of growth w/o metabolite to growth w/ metabolite. 5 - Essentiality prediction.
#'
#' @examples
#' eure <- readRDS(system.file("extdata", "eure.RDS", package = "MicrobiomeGS2"))
#' dt <- get_essential_nutrients(eure)
#' dt
#'
#' @export
get_essential_nutrients <- function(mod, min.growth.fraction = 0.01) {
  orig.growth <- get_growth(mod)

  if(orig.growth == 0)
    stop("Model cannot grow...")

  if(orig.growth < 0.01)
    warning(paste0("Model's growth rate is very small: ",round(orig.growth, digits = 8)))

  out.dt <- data.table(ex.rxn = mod@react_id,
                       name   = mod@react_name,
                       lb     = mod@lowbnd)
  out.dt <- out.dt[grepl("^EX_", ex.rxn) & lb < 0]
  out.dt[, lb := NULL]

  out.dt[, growth.rate.wo.met := NA_real_]
  out.dt[, growth.fraction := NA_real_]

  for(i in 1:nrow(out.dt)) {
    mod_tmp <- changeBounds(mod, react = out.dt[i, ex.rxn], lb = 0)
    sol_tmp <- get_growth(mod_tmp)

    out.dt[i, growth.rate.wo.met := sol_tmp]
    out.dt[i, growth.fraction := sol_tmp/orig.growth]
  }

  out.dt[, essential := F]
  out.dt[growth.fraction < min.growth.fraction, essential := T][]

  return(out.dt)
}

#' Constrain model
#'
#' @description  Constrains a model to a specific set of lower bounds for exchange reactions.
#'
#' @param mod Model file of class `modelorg`
#' @param media.id ID for available sybil.tools medium to which model is constraint to.
#' @param media.file Comma- or Tab-separated table with three columns: First - Compound ID; Second - Compound name; Third - maximum uptake flux.
#' @param nasp Namespace of the model. Default: "auto". See details for supported identifier namespaces. Has only an effect if media is selected via `media.id`
#'
#' @return A model object of class `modelorg`.
#'
#' @details Three metabolite namespaces are supported: `seed` for ModelSEED, `bigg` for BIGG, and `vmh` for Virtual Metabolic Human. If `auto`, the functions takes an educated guess.
#'
#' @export
constrain_mod <- function(mod, media.id = NULL, media.file = NULL, nasp = "auto") {

  # consistency check
  if(!is.null(media.id) & !is.null(media.file))
    stop("Only one of 'media.id' or 'media.file' should be specified.")
  if(is.null(media.id) & is.null(media.file))
    stop("Either 'media.id' or 'media.file' need to be specified.")

  nasp = "auto"
  if(nasp == "auto")
    nasp <- get_mod_namespace(mod)

  # read sybil.tools provided media list table
  if(!is.null(media.id)[1]) {
    mediadb <- system.file("extdata", "media.tsv", package="sybil.tools")
    mediadb <- fread(mediadb)
    mediadb <- mediadb[medium == media.id]
    if(nrow(mediadb) == 0) {
      stop(paste("The medium with the ID",media.id,"is not (yet) supported or part of the Medium-Database."))
    }
  }

  # read custom/user table
  if(!is.null(media.file)) {
    media.tab <- fread(media.file)
    if(ncol(media.tab) != 3)
      stop("Media file should have exactly three columns: First - Compound ID; Second - Compound name; Third - minimum flux (negative for uptake).")

    if(nasp == c("bigg"))
      colnames(media.tab) <- c("compound", "name", "maxFlux")
    if(nasp == c("vmh"))
      colnames(media.tab) <- c("vmh.id", "name", "maxFlux")
    if(nasp == c("seed"))
      colnames(media.tab) <- c("modelseed", "name", "maxFlux")

    mediadb <- media.tab
  }


  # 1. block all inflows
  mod@lowbnd[grep("^EX_", mod@react_id, fixed = F)] <- 0

  # 2. contrain model according to medium
  if(nasp == "bigg")
    mediadb[,ex.rxns := paste0("EX_",compound,"(e)")]
  if(nasp == "seed")
    mediadb[,ex.rxns := paste0("EX_",modelseed,"_e0")]
  if(nasp == "vmh")
    mediadb[,ex.rxns := paste0("EX_",vmh.id,"(e)")]

  mediadb$mod.rxn.id <- match(mediadb$ex.rxns,mod@react_id)

  # Metabolite already present
  media1 <- copy(mediadb[!is.na(mod.rxn.id)])
  mod@lowbnd[media1$mod.rxn.id] <- -media1$maxFlux

  return(mod)
}

#' Guess compound identifier namespace
#'
#' @description Guesses the namespace of metabolite identifiers
#'
#' @param mod Model of class `modelorg`
#'
#' @return Character. Either 'seed', 'bigg', or 'vhm'
#'
#' @examples
#' ecoli <- readRDS(system.file("extdata", "ecoli.RDS", package = "MicrobiomeGS2"))
#' get_mod_namespace(ecoli)
#'
#' bilo <- readRDS(system.file("extdata", "bilo.RDS", package = "MicrobiomeGS2"))
#' get_mod_namespace(bilo)
#'
#' @export
get_mod_namespace <- function(mod) {
  nasp <- NA
  if(any(grepl("cpd00002",mod@met_id))) {
    nasp = "seed"
    cat("Namespace: SEED\n")
  }
  else if(any(grepl("atp",mod@met_id))) {
    if(any(grepl("^ala_L|^arg_L|^asp_L|^cys_L|^glu_L|^glc_D",mod@met_id))) {
      nasp = "vmh"
      cat("Namespace: VMH\n")
    }
    else if(any(grepl("^ala__L|^arg__L|^asp__L|^cys__L|^glu__L|^glc__D",mod@met_id))) {
      nasp = "bigg"
      cat("Namespace: BIGG\n")
    }
  }

  return(nasp)
}


#' Get all exchange flux values
#'
#' @description Performs FBA(-MTF) to predict fluxes of metabolite exchange
#' reactions
#'
#' @param model Model file of class `modelorg`, or a named-list of `modelorg`
#' objects.
#' @param algorithm Algorithm to use to calculate flux distribution. Either 'mtf'
#' or 'fba'.
#' @param combine.compounds list of character vectors. Optional. This option can
#' be used to add up the exchange fluxes of specific compounds. This for instance
#' makes sense for enantiomers such as D- and L-Lactate, which is also set as
#' default and example.
#'
#' @return A data.table, with the columns: 'model' for model ID, 'rxn' for the
#' ID if the exchange reaction, 'flux' for the predicted flux value, 'in.model'
#' specifying if the exchange reaction is part of the model (TRUE) or
#' if absent (FALSE), 'name' for metabolite names.
#'
#' @details
#' In case an exchange reaction is not part of one model but present in another
#' model of the input models, a flux value of 0 is indicated.
#'
#' @export
get_exchanges <- function(model, algorithm = "mtf",
                          combine.compounds = list(`DL-Lactate` = c("EX_cpd00159_e0",
                                                                    "EX_cpd00221_e0"))) {
  if(is.list(model)) {
    if(!all(unlist(lapply(model, class)) == "modelorg"))
      stop("Not all models in list are of type 'modelorg'")
    if(is.null(names(model)))
      stop("Model list is not a named list.")
    if(any(is.na(names(model))))
      stop("Not all items in the list of models are named.")
  } else {
    if(class(model) != "modelorg")
      stop("model is not of class 'modelorg'")
    model <- list(model)
    names(model)[1] <- model[[1]]@mod_id
  }

  flx <- lapply(model, FUN = function(mod) {
    flxtmp <- get_flux_distribution(mod, exclude.unused = F,
                                    algorithm = algorithm)
    flxtmp <- flxtmp[grepl("^EX", rxn), .(rxn, flux)]
    flxtmp$name <- mod@react_name[match(flxtmp$rxn, mod@react_id)]
    flxtmp <- flxtmp[rxn != "EX_cpd11416_c0"] # that's the biomass
    flxtmp[, name := gsub("-e0 Exchange$| Exchange","", name)]

    return(flxtmp)
  })
  flx <- rbindlist(flx, idcol = "model")

  # save metabolite ids and names
  metinfo <- copy(flx[!duplicated(rxn), .(rxn, name)])

  # extent to have all model ~ exchange combinations
  flx <- dcast(flx, formula = model ~ rxn, value.var = "flux")
  flx <- melt(flx, id.vars = "model", variable.name = "rxn",
              value.name = "flux")
  flx[, in.model := TRUE]
  flx[is.na(flux), in.model := FALSE]
  flx[is.na(flux), flux := 0]

  # re-add metabolite names
  flx <- merge(flx, metinfo, by = "rxn")

  # combine compound fluxes (if applicable)
  if(!is.null(combine.compounds)) {
    if(length(combine.compounds) > 0) {
      for(icpd in 1:length(combine.compounds)) {
        # assign only one id
        flx[rxn %in% combine.compounds[[icpd]],
            rxn := combine.compounds[[icpd]][1]]
        # assign new name
        flx[rxn %in% combine.compounds[[icpd]],
            name := names(combine.compounds)[icpd]]
      }

      # sum fluxes for combined compounds
      flx <- flx[, .(flux = sum(flux), in.model = any(in.model)),
                 by = .(model, rxn, name)]
    }
  }

  return(flx)
}

