#' @title Join multiple metabolic models
#'
#' @description The function merges multiple metabolic models into one model
#' where each organisms is within an own compartment. Organism compartments are
#' connected by a shared extracellular space (i.e. environment), from where the
#' organisms can take up nutrients an release metabolic by-products. Exchange
#' reactions are attached to to shared extracelluar space.
#'
#' @param model.list A list of metabolic models in 'sybil' format.
#' @param scale.boundaries Scaling factor for all reaction bounds. Might be useful for large communities (>20 species). Default: 1
#' @param merge.lb.method Method to use in cases where the individual models
#' differ in the lower bound values for the same exchange reactions. "minimum"
#' uses the minimum LB values (i.e. largest maximum uptake rate) found among the
#' individual models for the final joined model. "maximum" analogously. "none"
#' will cause the method to stop with an error, if lower boundaries differ
#' between input models. Further options: "median","mean".
#'
#' @return A list with three objects:
#'   model.IDs - a data.table, which states the newly assigned internal model IDs (M1, M2, ...) to the corresponding organism name.
#'   modj - A community model also in sybil format. Each organism has its own compartment ('c' and 'e',), but metabolites can freely move between the 'e' compartments of different species.
#'   ex.rxns - a data.table of all exchange reaction of the new community model.
#'
#' @export
join_mult_models <- function(model.list, scale.boundaries = 1,
                             merge.lb.method = "none") {
  require(sybil)
  require(data.table)

  if(!(merge.lb.method %in% c("none","minimum","maximum","mean","median")))
    stop("Unknown option for 'merge.lb.method'")

  n <- length(model.list)


  model.IDs <- data.table(model.name = unlist(lapply(model.list,FUN = function(x) x@mod_id)),
                          model.id = as.character(1:n))
  model.IDs$model.id <- paste0("M",model.IDs$model.id)

  # get a table of all exchange reactions and their lower bounds
  #cat("Extracting all metabolite exchange (i.e. \"EX\") reactions...\n")
  ex.rxns <- data.table(model.name = character(0L), react_id = character(0L), lb = double(0L), react_name = character(0L))
  for(i in 1:n) {
    ex.ind <- grep("^EX_",model.list[[i]]@react_id)
    tmp <- data.table(model.name = model.list[[i]]@mod_id,
                      react_id = model.list[[i]]@react_id[ex.ind],
                      lb = model.list[[i]]@lowbnd[ex.ind],
                      react_name = model.list[[i]]@react_name[ex.ind])
    ex.rxns <- rbind(ex.rxns,tmp)
  }
  ex.rxns[,model.name := NULL]

  if(merge.lb.method == "none") {
    ex.rxns <- ex.rxns[!duplicated(paste0(react_id, lb))]
    if(any(duplicated(ex.rxns$react_id)))
      stop(paste0("unequal lower bounds for at least one exchange reaction. "))
  }
  if(merge.lb.method == "minimum") {
    ex.rxns <- ex.rxns[order(react_id, lb)]
    ex.rxns <- ex.rxns[!duplicated(react_id)]
  }
  if(merge.lb.method == "maximum") {
    ex.rxns <- ex.rxns[order(react_id, -lb)]
    ex.rxns <- ex.rxns[!duplicated(react_id)]
  }
  if(merge.lb.method == "mean") {
    ex.rxns[, lb := mean(lb), by = react_id]
    ex.rxns <- ex.rxns[!duplicated(react_id)]
  }
  if(merge.lb.method == "median") {
    ex.rxns[, lb := median(lb), by = react_id]
    ex.rxns <- ex.rxns[!duplicated(react_id)]
  }

  ex.rxns[,met := gsub("^EX_","",react_id)]
  ex.rxns <- ex.rxns[met != "cpd11416_c0"]

  # rename reactions and metabolites
  #cat("Renaming reaction-, metabolite, and compartment IDs...\n")
  for(i in 1:n){
    model.list[[i]]@react_id    <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@react_id)
    model.list[[i]]@met_id      <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@met_id)
    model.list[[i]]@mod_compart <- paste0(model.IDs$model.id[i],"_", model.list[[i]]@mod_compart)
  }

  # initiating joined/community model.
  compNames <- character(0)
  for(i in 1:n)
    compNames <- c(compNames,model.list[[i]]@mod_compart)

  modj <- modelorg(id = "joined.mod",
                   name = paste0("Number of species in community: ",n),
                   compartment = compNames)

  # get final met X rxn dimensions
  model.inds <- data.table(modelID = model.IDs$model.name,
                           n.met = unlist(lapply(model.list, FUN = function(x) x@met_num)),
                           n.rxn = unlist(lapply(model.list, FUN = function(x) x@react_num)),
                           n.gen = unlist(lapply(model.list, FUN = function(x) length(x@allGenes))))
  model.inds[,met.to := cumsum(n.met)]
  model.inds[,met.from := NA_real_]; model.inds[1, met.from := 1]; model.inds[-1, met.from := model.inds[-n,met.to]+1]
  model.inds[,rxn.to := cumsum(n.rxn)]
  model.inds[,rxn.from := NA_real_]; model.inds[1, rxn.from := 1]; model.inds[-1, rxn.from := model.inds[-n,rxn.to]+1]
  model.inds[,gen.to := cumsum(n.gen)]
  model.inds[,gen.from := NA_real_]; model.inds[1, gen.from := 1]; model.inds[-1, gen.from := model.inds[-n,gen.to]+1]

  modj@met_num    <- sum(model.inds$n.met)
  modj@react_num  <- sum(model.inds$n.rxn)
  modj@S          <- Matrix(0,nrow = modj@met_num, ncol = modj@react_num, sparse = T)
  # Metabolites:
  modj@met_attr   <- model.list[[1]]@met_attr
  modj@met_comp   <- as.integer(rep(0,modj@met_num))
  modj@met_de     <- logical(modj@met_num)
  modj@met_id     <- character(modj@met_num)
  modj@met_name   <- character(modj@met_num)
  modj@met_single <- logical(modj@met_num)
  # Reactions:
  modj@react_attr   <- model.list[[1]]@react_attr
  modj@react_de     <- logical(modj@react_num)
  modj@react_id     <- character(modj@react_num)
  modj@react_name   <- character(modj@react_num)
  modj@react_rev    <- logical(modj@react_num)
  modj@react_single <- logical(modj@react_num)
  modj@subSys       <- Matrix(T,ncol = 1, nrow = modj@react_num, sparse = T)
  # Bounds & objectives:
  modj@lowbnd <- numeric(modj@react_num)
  modj@uppbnd <- numeric(modj@react_num)
  modj@obj_coef <- numeric(modj@react_num)
  # Genes:
  modj@gprRules   <- character(modj@react_num)
  modj@gpr        <- character(modj@react_num)
  modj@genes      <- list()
  modj@allGenes   <- character(sum(model.inds$n.gen))
  modj@rxnGeneMat <- Matrix(0, nrow = modj@react_num, ncol = sum(model.inds$n.gen), sparse = T)


  for(i in 1:n) {
    cat(paste0("\r",i,"/",n))
    # Stoichiometric matrix # das geht wharhscheinlich auch noch schneller
    I  <- Matrix::which(model.list[[i]]@S != 0, arr.ind = T)
    Ij <- I; Ij[,1] <- Ij[,1] + model.inds[i,met.from] - 1; Ij[,2] <- Ij[,2] + model.inds[i,rxn.from] - 1
    modj@S[Ij] <- model.list[[i]]@S[I]
    #modj@S[model.inds[i,met.from]:model.inds[i,met.to],
    #       model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@S
    # Metabolites
    if(i > 1)
      modj@met_attr <- rbind(modj@met_attr, model.list[[i]]@met_attr)
    max.c <- max(modj@met_comp)
    modj@met_comp[model.inds[i,met.from]:model.inds[i,met.to]] <- model.list[[i]]@met_comp + max.c
    modj@met_de[model.inds[i,met.from]:model.inds[i,met.to]] <- model.list[[i]]@met_de
    modj@met_id[model.inds[i,met.from]:model.inds[i,met.to]] <- model.list[[i]]@met_id
    modj@met_name[model.inds[i,met.from]:model.inds[i,met.to]] <- model.list[[i]]@met_name
    modj@met_single[model.inds[i,met.from]:model.inds[i,met.to]] <- model.list[[i]]@met_single

    # Reactions:
    if(i > 1)
      modj@react_attr <- rbind(modj@react_attr, model.list[[i]]@react_attr)
    modj@react_de[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@react_de
    modj@react_id[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@react_id
    modj@react_name[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@react_name
    modj@react_rev[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@react_rev
    modj@react_single[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@react_single

    # Bounds & Objective:
    modj@lowbnd[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@lowbnd
    modj@uppbnd[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@uppbnd
    modj@obj_coef[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@obj_coef

    # Genes:
    modj@gprRules[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@gprRules
    modj@gpr[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@gpr
    modj@genes[model.inds[i,rxn.from]:model.inds[i,rxn.to]] <- model.list[[i]]@genes
    modj@allGenes[model.inds[i,gen.from]:model.inds[i,gen.to]] <- model.list[[i]]@allGenes
    modj@rxnGeneMat[model.inds[i,rxn.from]:model.inds[i,rxn.to],
                    model.inds[i,gen.from]:model.inds[i,gen.to]] <- model.list[[i]]@rxnGeneMat
  }

  # add new external exchanges
  #cat("\rAdding new external exchanges...\n")
  mod_compart(modj) <- c(mod_compart(modj), "e0")
  Nex <- nrow(ex.rxns)
  modj <- addMultiReact(model=modj, ids=ex.rxns$react_id, mets=ex.rxns$met, Scoefs = rep(-1, Nex),
                        reversible = TRUE, lb=ex.rxns$lb, metComp = rep("e0", Nex),
                        reactName = ex.rxns$react_name)
  modj@react_rev[is.na(modj@react_rev)] <- TRUE
  mod_name(modj) <- "A GS - joined model."

  # setting up original exchange interactions to interact with new common "e" compartment
  # e.g. "M1_ac(e) <->"  ==> "M1_ac(e) <-> ac(e)" + removing lower bnd (new: -1000)
  #cat("Modifying original exchange reactions to interact with new common \"e\" compartment...\n")
  r.ind <- grep("^M[0-9]+_EX_",modj@react_id)
  tmp.ind <- grep("^M[0-9]+_EX_cpd11416",modj@react_id)
  r.ind <- r.ind[!(r.ind %in% tmp.ind)]
  tmp.mets <- gsub("^M[0-9]+_EX_","",modj@react_id[r.ind])
  m.ind.e  <- match(tmp.mets, modj@met_id)
  modj@lowbnd[r.ind] <- rep(-1000,length(r.ind))

  I <- matrix(0,ncol = 2, nrow = length(r.ind))
  I[,1] <- m.ind.e; I[,2] <- r.ind
  modj@S[I] <- 1

  # Scale boundaries
  modj@lowbnd <- modj@lowbnd * scale.boundaries
  modj@uppbnd <- modj@uppbnd * scale.boundaries

  return(list(model.IDs=model.IDs,
              modj = modj,
              ex.rxns = ex.rxns))
}

# addMultiReact
#
# extension for sybil to allow to add several reaction from one model to another
#
# TODO: add handling of gprAssoc and subSystem

addMultiReact <- function(model,
                          ids,
                          src = NA, # if defined then all other parameters will be taken from this model
                          mets = NA,   # vector needed
                          Scoefs = NA, # vector needed
                          reversible = FALSE,
                          lb = 0,
                          ub = SYBIL_SETTINGS("MAXIMUM"),
                          obj = 0,
                          subSystem = NA,
                          gprAssoc = NA,
                          reactName = NA,
                          metName = NA,
                          metComp = NA) {


  # ------------------------------------------------------------------------ #
  # check arguments
  # ------------------------------------------------------------------------ #

  if (!is(model, "modelorg")) {
    stop("needs an object of class modelorg!")
  }

  stopifnot(checkVersion(model))

  if (!is(src, "modelorg")){
    Nids <- length(ids)
    if (length(mets) != Nids | length(Scoefs) != Nids){
      stop("all arguments have to be provided for each reaction")
    }
  }

  check_mets <- sapply(1:length(ids), function(i){length(mets[i]) == length(Scoefs[i])})
  if (any(is.na(check_mets) | !all(check_mets))) {
    stop("in each reaction arguments 'met' and 'Scoef' must have the same length")
  }


  # ------------------------------------------------------------------------ #
  # get main parameters from source model
  # ------------------------------------------------------------------------ #

  if(is(src, "modelorg")){
    reactInd  <- match(ids, react_id(src))
    reactName <- react_name(src)[reactInd]
    Crev      <- react_rev(src)[reactInd]
    lb        <- lowbnd(src)[reactInd]
    ub        <- uppbnd(src)[reactInd]
    obj       <- obj_coef(src)[reactInd]
    subSystem <- NA #subSys(src)[reactInd,] # subsystems not implemented yet
    gprAssoc  <- NA # gene association not implemented yet
    mets      <- sapply(reactInd, function(r){met_id(src)[which(S(src)[,r]!=0)]})
    Scoefs    <- sapply(reactInd, function(r){S(src)[,r][which(S(src)[,r]!=0)]})
    met       <- unique(unlist(mets))
    metInd    <- match(met, met_id(src))
    metName   <- met_name(src)[metInd]
    metComp   <- met_comp(src)[metInd]
  }else{
    met       <- unique(unlist(mets))
    nR        <- length(ids)
    if(length(obj)==1 & nR>1) obj <- rep(obj, nR)
    if(length(lb)==1 & nR>1)  lb  <- rep(lb,  nR)
    if(length(ub)==1 & nR>1)  ub  <- rep(ub,  nR)
    if(length(reversible)==1 & nR>1) reversible <- rep(reversible, nR)
    Crev      <- ifelse( ( (ub > 0) & (lb < 0) ) & (!isTRUE(reversible)), TRUE, reversible)
  }


  # ------------------------------------------------------------------------ #
  # check, if we need to add columns and/or rows
  # ------------------------------------------------------------------------ #

  # reactions
  colInd    <- match(ids, react_id(model))
  newR      <- which(is.na(colInd))
  nCols     <- react_num(model)
  nNewCols  <- length(newR)
  addCol    <- FALSE

  for(i in seq(along = newR)){
    addCol <- TRUE
    nCols <- nCols + 1
    colInd[newR[i]] <- nCols
  }


  # metabolites

  rowInd  <- match(met, met_id(model))

  newM    <- which(is.na(rowInd))
  nRows   <- met_num(model)          # number of rows in the model
  nNewRows<- length(newM)            # number of new rows
  addRow  <- FALSE

  for (i in seq(along = newM)) {
    addRow   <- TRUE
    nRows    <- nRows + 1
    rowInd[newM[i]] <- nRows
  }


  if ( (isTRUE(addCol)) || (isTRUE(addRow)) ) {

    # -------------------------------------------------------------------- #
    # make a new model
    # -------------------------------------------------------------------- #

    # -------------------------------------------------------------------- #
    # data structures

    newmet_num      <- met_num(model)
    newmet_id       <- met_id(model)
    newmet_name     <- met_name(model)
    newmet_comp     <- met_comp(model)
    newmet_single   <- met_single(model)
    newmet_de       <- met_de(model)

    newreact_num    <- react_num(model)
    newreact_rev    <- react_rev(model)
    newreact_id     <- react_id(model)
    newreact_name   <- react_name(model)
    newreact_single <- react_single(model)
    newreact_de     <- react_de(model)
    newlowbnd       <- lowbnd(model)
    newuppbnd       <- uppbnd(model)
    newobj_coef     <- obj_coef(model)

    newgprRules     <- gprRules(model)
    newgenes        <- genes(model)
    newgpr          <- gpr(model)
    newallGenes     <- allGenes(model)
    newrxnGeneMat   <- rxnGeneMat(model)
    newsubSys       <- subSys(model)

    newS            <- S(model)

    newMetAttr <- met_attr(model)
    newReactAttr <- react_attr(model)
    newCompAttr <- comp_attr(model)
    newModAttr <- mod_attr(model)


    if (isTRUE(addRow)) {

      # new number of metabolites
      newmet_num  <- nRows

      # new metabolite id's
      newmet_id   <- append(met_id(model), met[newM])

      # new metabolite names
      if (any(is.na(metName))) {
        newmet_name <- append(met_name(model), met[newM])
      }
      else {
        newmet_name <- append(met_name(model), metName[newM])
      }

      # new metabolite compartments
      if (any(is.na(metComp))) {
        newmet_comp <- append(met_comp(model), rep(NA, nNewRows))
      }
      else {
        if (is(metComp, "numeric")) {
          newmet_comp <- append(met_comp(model), metComp[newM])
        }
        else {
          newmet_comp <- append(met_comp(model),
                                match(metComp[newM],
                                      mod_compart(model)))
        }
      }

      # singleton and dead end metabolites (not checked!)
      newmet_single <- append(met_single(model), rep(NA, nNewRows))
      newmet_de     <- append(met_de(model),     rep(NA, nNewRows))

      # new rows in stoichiometric matrix
      newRows <- Matrix::Matrix(0,
                                nrow = nNewRows,
                                ncol = react_num(model))
      newS <- rbind(newS, newRows)

      # new met attrs
      if(ncol(newMetAttr) > 0){
        newMetAttr[nrow(newMetAttr)+1:nNewRows, ] <- NA
      }
    }

    if (isTRUE(addCol)) {                        # we add at least one column
      # new number of reactions
      newreact_num  <- nCols

      # new reaction ids
      newreact_id   <- append(react_id(model), ids[newR])

      # new reaction names
      if (any(is.na(reactName))) {
        newreact_name <- append(react_name(model), ids[newR])
      }
      else {
        newreact_name <- append(react_name(model), reactName[newR])
      }

      # reaction contains singleton or dead end metabolites (not checked!)
      newreact_single <- append(react_single(model), rep(NA, nNewCols))
      newreact_de     <- append(react_de(model),     rep(NA, nNewCols))

      # reversibility, lower and upper bounds, objective coefficient
      newreact_rev <- append(react_rev(model), Crev[newR])
      newlowbnd    <- append(lowbnd(model),    lb[newR])
      newuppbnd    <- append(uppbnd(model),    ub[newR])
      newobj_coef  <- append(obj_coef(model),  obj[newR])

      # new columns in stoichiometric matrix
      newColumns <- Matrix::Matrix(0,
                                   ncol = nNewCols,
                                   nrow = nrow(newS))

      newS <- cbind(newS, newColumns)

      # new react Attr
      if(ncol(newReactAttr) > 0){
        newReactAttr[nrow(newReactAttr)+1:nNewCols, ] <- NA
      }

      # subsystems
      if (any(is.na(subSystem))) {
        ss <- subSys(model)
        if(ncol(ss)==0){ # if no subSys defined, rbind (see else) failed
          dim(ss) <- c(nrow(ss)+nNewCols, ncol(ss)) # new reactions are rows in this matrix
          newsubSys <- ss
        }
        else {
          newsubSys <- rbind(ss, rep(rep(FALSE, ncol(subSys(model))), nNewCols))
        }
      }
      # else {
      #   if (is(subSystem, "logical")) {
      #     newsubSys <- rBind(subSys(model), subSystem)
      #   }
      #   else {
      #     nSubsRow  <- colnames(subSys(model)) %in% subSystem
      #     newsubSys <- rBind(subSys(model), nSubsRow)
      #   }
      # }


      # gpr association
      if (ncol(rxnGeneMat(model)) > 0) {
        newrxnGeneMat   <- rbind(rxnGeneMat(model),
                                 rep(rep(FALSE, ncol(rxnGeneMat(model))), nNewCols))
      }
      else { #if (nrow(rxnGeneMat(model)) > 0) {
        newrxnGeneMat <- rxnGeneMat(model)
        dim(newrxnGeneMat) <- c(nrow(newrxnGeneMat)+nNewCols,
                                ncol(newrxnGeneMat))
      }
      # do above else always.

      if ( (is.na(gprAssoc)) || (gprAssoc == "") ) {
        if ((length(gprRules(model)) > 0)) {
          newgprRules     <- append(gprRules(model), rep("", nNewCols))
          newgenes        <- append(genes(model), rep("", nNewCols))
          newgpr          <- append(gpr(model), rep("", nNewCols))
        }
      }
      # else {
      #   gene_rule <- .parseBoolean(gprAssoc)
      #
      #   geneInd <- match(gene_rule$gene, allGenes(model))
      #
      #   # indices of new genes
      #   new_gene <- which(is.na(geneInd))
      #
      #   # if we have new gene(s), add a column in rxnGeneMat and
      #   # gene name(s) to allGenes
      #   if (length(new_gene) > 0) {
      #     newallGenes <- append(allGenes(model),
      #                           gene_rule[["gene"]][new_gene])
      #
      #     # update geneInd
      #     geneInd <- match(gene_rule[["gene"]], newallGenes)
      #
      #     # if we have an empty modelorg object, we need to
      #     # initialize rxnGeneMat
      #     if (ncol(newrxnGeneMat) == 0) {
      #       newrxnGeneMat <- Matrix::Matrix(FALSE,
      #                                       nCols, max(geneInd))
      #     }
      #     else {
      #       for (i in seq(along = gene_rule[["gene"]][new_gene])) {
      #         newrxnGeneMat <- cBind(newrxnGeneMat,
      #                                rep(FALSE, nrow(newrxnGeneMat)))
      #       }
      #     }
      #   }
      #
      #   # rxnGeneMat
      #   newrxnGeneMat[nCols, geneInd] <- TRUE
      #
      #   # new rule
      #   newgpr <- append(gpr(model), gprAssoc)
      #
      #   # genes per reaction
      #   newgenes <- append(genes(model), list(gene_rule$gene))
      #   newrule  <- gene_rule$rule
      #
      #   # not needed for modelorg version 2.0
      #   #                for (j in 1 : length(geneInd)) {
      #   #                    pat  <- paste("x(", j, ")", sep = "")
      #   #                    repl <- paste("x[", geneInd[j], "]", sep = "")
      #   #
      #   #                    newrule <- gsub(pat, repl, newrule, fixed = TRUE)
      #   #                }
      #
      #   newgprRules <- append(gprRules(model), newrule)
      # }
    }

    # values for stoichiometric matrix
    newS[ , colInd] <- sapply(1:length(colInd), function(i){
      newCol <- rep(0, length = nrow(newS))
      curRInd <- colInd[i]
      #curMInd <- metInd[match(mets[[i]], met)]
      curMInd <- rowInd[match(mets[[i]], met)]
      newCol[curMInd] <- Scoefs[[i]]
      newCol
    })


    # -------------------------------------------------------------------- #
    # new model
    # -------------------------------------------------------------------- #

    if (is(model, "modelorg_irrev")) {
      mod_out <- modelorg_irrev(mod_id(model), mod_name(model))
      irrev(mod_out)     <- TRUE
      matchrev(mod_out)  <- append(matchrev(model), 0L)

      revReactId <- as.integer(max(irrev2rev(model))+1)
      irrev2rev(mod_out) <- append(irrev2rev(model), revReactId)
      rev2irrev(mod_out) <- rbind(rev2irrev(model), c(nCols, nCols))
    } else {
      mod_out <- modelorg(mod_id(model), mod_name(model))
    }

    mod_desc(mod_out)    <- mod_desc(model)
    mod_compart(mod_out) <- mod_compart(model)


    met_num(mod_out)      <- as.integer(newmet_num)
    met_id(mod_out)       <- newmet_id
    met_name(mod_out)     <- newmet_name
    met_comp(mod_out)     <- as.integer(newmet_comp)
    met_single(mod_out)   <- newmet_single
    met_de(mod_out)       <- newmet_de

    react_num(mod_out)    <- as.integer(newreact_num)
    react_rev(mod_out)    <- newreact_rev
    react_id(mod_out)     <- newreact_id
    react_name(mod_out)   <- newreact_name
    react_single(mod_out) <- newreact_single
    react_de(mod_out)     <- newreact_de
    lowbnd(mod_out)       <- newlowbnd
    uppbnd(mod_out)       <- newuppbnd
    obj_coef(mod_out)     <- newobj_coef

    gprRules(mod_out)     <- newgprRules
    genes(mod_out)        <- newgenes
    gpr(mod_out)          <- newgpr
    allGenes(mod_out)     <- newallGenes
    rxnGeneMat(mod_out)   <- newrxnGeneMat
    subSys(mod_out)       <- newsubSys

    S(mod_out)            <- newS

    react_attr(mod_out) <- newReactAttr
    met_attr(mod_out) <- newMetAttr
    comp_attr(mod_out) <- newCompAttr
    mod_attr(mod_out) <- newModAttr


  } else{
    stop("Nothing to change")
  }

  check <- validObject(mod_out, test = TRUE)

  if (check != TRUE) {
    msg <- paste("Validity check failed:", check, sep = "\n    ")
    warning(msg)
  }

  return(mod_out)

}

#' @title Couple organisms' internal reactions to their own biomass production
#'
#' @param modj Community model as returned from join_mult_models()$modj
#' @param cpl_c Coupling parameter c
#' @param cpl_u Coupling parameter u
#'
#' @export
get_coupling_constraints_mult <- function(modj, cpl_c = 400, cpl_u = 0.01) {
  coupling <- list(react=list(), x=list(), lb=vector(), ub=vector(), rtype=vector())
  #intern_rea  <- grep("EX_", modj@react_id, invert = T)
  intern_rea  <- 1:modj@react_num
  allObj.ind  <- grep("^M[0-9]+_bio1",modj@react_id)

  m.obj.dt <- data.table(m.id = paste0("M",1:length(allObj.ind)),
                         corr.BM.id  = allObj.ind)
  dt.cpl <- data.table(int.rxn.id = intern_rea,
                       rev = modj@lowbnd[intern_rea] < 0)
  dt.cpl[, m.id := gsub("_.*", "", react_id(modj)[int.rxn.id])]
  dt.cpl <- merge(dt.cpl,m.obj.dt, by="m.id")
  dt.cpl[, dir := "U"]
  dt.cpl.tmp <- copy(dt.cpl[rev==T])
  dt.cpl.tmp[, dir := "L"]
  dt.cpl <- rbind(dt.cpl, dt.cpl.tmp)

  coupling$react <- structure(as.list(as.data.frame(apply(dt.cpl, 1, FUN = function(x) c(as.integer(x["int.rxn.id"]),as.integer(x["corr.BM.id"])))))
                              , names = NULL)
  coupling$x <- structure(as.list(as.data.frame(apply(dt.cpl, 1, FUN = function(x) c(1,ifelse(x["dir"]=="U",-cpl_c,cpl_c)))))
                          , names = NULL)
  coupling$lb <- dt.cpl[,ifelse(dir=="L",-cpl_u,NA_real_)] ## check this: Yes: Definitely "-cpl_u"
  coupling$ub <- dt.cpl[,ifelse(dir=="U",cpl_u,NA_real_)]
  coupling$rtype <- dt.cpl[,dir]

  return(coupling)
}


#' @title Community FBA with fixed biomass ratios
#'
#' @param models List of models of class `modelorg`.
#' @param model.prop Vector of relative proportions of each model's biomass in
#' the total community biomass. Should be of the same length as number of
#' models.
#' @param pFBAcoeff parsimonious FBA coefficient that corresponds to the weight of absolute flux values relative to the biomass mass production.
#' @param lp.method Solver algorith to use. Default: hybbaropt.
#' @param scale.boundaries Factor (x) that scaled the lower bounds of Exchange Reactions. Experimental, thus it is recommended to use the default: 1 (no scaling).
#' @param cpl_c Coupling constraint \code{c}. See Heinken et al (2013) Gut Microbes (doi: 10.4161/gmic.22370). Default: 400
#' @param cpl_u Coupling constraint \code{u}. See Heinken et al (2013) Gut Microbes (doi: 10.4161/gmic.22370). Default: 0.01
#'
#' @export
communityFBA_FB <- function(models, model.prop, scale.boundaries = 1,
                            cpl_c = 400, cpl_u = 0.01, pFBAcoeff = 1e-6,
                            lp.method = "hybbaropt") {

  sybil::SYBIL_SETTINGS("METHOD", lp.method)
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

  model.prop <- model.prop/(sum(model.prop))
  spec.ratio <- data.table(spec  = unlist(lapply(models, function(x) x@mod_id)),
                           ratio = model.prop)
  if(any(duplicated(spec.ratio$model.id)))
    stop("There are models with identical model IDs.")

  # check if models can form BM alone
  org.gr    <- lapply(models, FUN = function(x) return(optimizeProb(x)@lp_obj))
  nogr.orgs <- names(which(org.gr == 0))

  if(length(nogr.orgs) > 0)
    warning(paste0("Following models do not grow alone: ",paste(nogr.orgs, collapse = ", ")))

  # construct community model
  mod.joined <- join_mult_models(models, scale.boundaries = scale.boundaries,
                                 merge.lb.method = "minimum")
  mod.joined$model.IDs <- merge(mod.joined$model.IDs,spec.ratio,
                                by.x="model.name",by.y="spec")
  bm.mets <- paste0(mod.joined$model.IDs$model.id,"_cpd11416[c0]")

  mod.joined$modj <- addReact(model = mod.joined$modj, id = "EX_BIOMASS",
                              met = bm.mets,
                              Scoef = -mod.joined$model.IDs$ratio,
                              reversible = F,
                              reactName = "joined Biomass with fixed ratios")
  mod.joined$modj@obj_coef <- rep(0, length(mod.joined$modj@react_id))
  mod.joined$modj@obj_coef[grep("EX_BIOMASS", mod.joined$modj@react_id)] <- 1
  # block individual Biomass outflow reactions
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11416",mod.joined$modj@react_id)] <- 0

  # block quinon exchanges
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11606",mod.joined$modj@react_id)] <- 0.5
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11451",mod.joined$modj@react_id)] <- 0.5

  # Flux coupling
  coupling <- get_coupling_constraints_mult(mod.joined$modj,
                                            cpl_c = cpl_c, cpl_u = cpl_u)
  # Simulation
  modj_warm <- sysBiolAlg(mod.joined$modj,
                          algorithm = "mtfEasyConstraint2",
                          easyConstraint=coupling,
                          pFBAcoeff = pFBAcoeff,
                          scaling = scale.boundaries)
  mod.joined$solj <- optimizeProb(modj_warm)

  # Get community growth
  mod.joined$community.growth <- mod.joined$solj$fluxes[grep("EX_BIOMASS", mod.joined$modj@react_id)]
  # Get metabolic interchange
  mod.joined$met.interchange <- get_metabolic_interchange(mod.joined$modj, mod.joined$solj)
  # Optimization successful?
  mod.joined$stat <- mod.joined$solj$stat

  return(mod.joined)
}

#
# Get a quantitative measure for the metabolic interchange of individual metabolites
#
get_metabolic_interchange <- function(modj, solj) {
  int.ext.rxns <- grep("^M[0-9]+_EX_",modj@react_id)
  ext.ext.rxns <- grep("^EX_",modj@react_id)

  dt <- data.table(rxn = modj@react_id[int.ext.rxns],
                   flux = solj$fluxes[int.ext.rxns])
  dt[, rxn := gsub("M[0-9]+_","",rxn)]
  dt[, flux := abs(flux)]
  dt <- dt[,.(flux = sum(flux)),by="rxn"]

  # external exchange
  dt2 <- data.table(rxn = modj@react_id[ext.ext.rxns],
                    o.flux = solj$fluxes[ext.ext.rxns])

  dt <- merge(dt,dt2,by="rxn")
  dt[,flux := flux-abs(o.flux)]
  #dt <- dt[flux > 1e-6]

  dt[order(-flux)]
}
