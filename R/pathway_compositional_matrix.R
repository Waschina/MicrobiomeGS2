#' @title Calculate pathway composition matrix (experimental)
#'
#' @description In principle, this function translates the OTU/ASV count table
#' to a pathway count table based on the gapseq predictions of the reference
#' models/genomes. Please note, that here not only one reference model is used
#' but the ensemble of models that have mapped with the same statistics to the
#' corresponding OTU/ASV sequence.
#'
#' @param mic A Microbiome object.
#' @param model.dir Path to the directory, where gapseq model prediction files
#' are stored. Files should be stored in model-named subdirectories.
#' @param abun.ex.thres Fraction of models in an ensemble, that needs to be
#' exceeded in order to coin that the pathway is also present in the respective
#' ASV/OTU.
#' @param count.freq Logical. Indicating whether the relative pathway
#' composition matrix (TRUE) or the absolute count table (FALSE, default) is
#' returned.
#'
#' @return A list. TODO: provide results
#'
#' @export
pathway_compositional_matrix <- function(mic, model.dir, abun.ex.thres = 0.5,
                                         count.freq = FALSE) {

  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  # Ensemble-X-Sample abundance table #
  # E_S                               #
  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

  # identify common ensembles
  ens_map <- mic@model.mapping.ensemble[,.(target.model, query.label)]

  ens_map <- ens_map[, .(ensemble = paste(target.model, collapse = "$")), by = "query.label"]
  ens_map <- ens_map[, .(cl.query = paste(query.label, collapse = "$")), by = "ensemble"]
  ens_map <- ens_map[order(cl.query)]
  ens_map[, ensID := paste0("Ens_",1:.N)]

  E_S <- matrix(0, nrow = nrow(ens_map), ncol = ncol(mic@uniq.table))
  rownames(E_S) <- ens_map$ensID
  colnames(E_S) <- colnames(mic@uniq.table)
  E_S <- as.table(E_S)

  k <- 1
  for(ens in ens_map$ensID) {
    #cat("\r",k,"/",nrow(ens_map))

    iclust <- ens_map[ensID == ens, cl.query]
    iclust <- unlist(str_split(iclust,"\\$"))

    if(length(iclust) > 1) {
      E_S[ens,] <- colSums(mic@uniq.table[iclust,])
    } else {
      E_S[ens,] <- mic@uniq.table[iclust,]
    }

    k <- k + 1
  }
  #cat("\n")

  if(count.freq == T) {
    E_S <- t(t(E_S)/colSums(E_S))
  }

  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  # Genome-X-Pathway table  #
  # G_P                     #
  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

  # get gapseq's pathway predictions
  pwys <- list()
  missing.models <- c()
  for(i in mic@model.mapping.ensemble$target.model) {
    file.tmp <- paste0(model.dir,"/",i,"/",i,"-all-Pathways.tbl.gz")
    if(file.exists(file.tmp)) {
      a <- fread(file.tmp) # reading model's pathway prediction table

      a[, ID := gsub("(^\\|)|(\\|$)","", ID)] # removing starting and trailing pipes in PWY-IDs.
      a <- a[!duplicated(ID)] # removing duplicated entries table
      a <- a[, .(ID,Prediction)] # extracting only the relevant columns

      pwys[[i]] <- a
    } else {
      missing.models <- c(missing.models,i)
    }
  }

  if(length(missing.models) > 0) {
    warning(paste0("Following models are missing:\n    ", paste(missing.models, collapse = "\n    ")))
  }

  pwys <- rbindlist(pwys, idcol = "model")
  pwys$Prediction <- as.numeric(pwys$Prediction) # Predicton column from Boolean to {0,1}
  #pwys # inspecting table... three columns: model, PWY-ID and its prediction. looks good


  # Constructing G x P table
  # this is a long-to-wide transformation - function of choice: dcast()
  pwy_tab <- dcast(pwys, model ~ ID, value.var = "Prediction", fill = 0)

  m_names <- pwy_tab$model # temp storing model names in correct order
  pwy_tab <- as.matrix(pwy_tab[,-1]) # data.table to simple table for later algebra
  rownames(pwy_tab) <- m_names # restoring rownames with modelnames

  #dim(pwy_tab) # inspect dimensions

  # remove pathways that were never predicted to be present
  pwy_tab <- pwy_tab[, colSums(pwy_tab) > 0] # Removing columns/pathways with zero-only entries
  pwy_tab <- pwy_tab[, colSums(pwy_tab) < nrow(pwy_tab)] # Removing columns/pathways with one-only entries

  #dim(pwy_tab) # inspect dimensions . okay
  G_P <- pwy_tab

  # ~ ~ ~ ~ ~ ~ ~ #
  # get E_P table #
  # ~ ~ ~ ~ ~ ~ ~ #
  E_P <- matrix(0, nrow = nrow(E_S), ncol = ncol(G_P))
  colnames(E_P) <- colnames(G_P)
  rownames(E_P) <- rownames(E_S)

  for(i in rownames(E_P)) {
    iG <- ens_map[ensID == i, ensemble]
    iG <- unlist(str_split(iG,"\\$"))

    if(length(iG) > 1) {
      E_P[i, ] <- colSums(G_P[iG,])/length(iG)
    } else {
      E_P[i, ] <- G_P[iG,]
    }
  }
  E_P <- ifelse(E_P > abun.ex.thres, 1, 0) # OR >= ???
  E_P <- E_P[,which(colSums(E_P) > 0)] # only zero
  E_P <- E_P[,which(colSums(E_P) < nrow(E_P))]

  # ~ ~ ~ ~ ~ ~ ~ #
  # get P_S table #
  # ~ ~ ~ ~ ~ ~ ~ #
  P_S <- t(E_P) %*% E_S

  return(list(P_S = P_S,
              E_P = E_P,
              E_S = E_S,
              G_P = G_P,
              ens_map = ens_map))
}
