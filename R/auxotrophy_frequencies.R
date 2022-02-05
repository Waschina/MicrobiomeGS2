#' @title Calculate auxotrophy frequency table
#'
#' @description Bases on the ASV/OTU-Count table and the mapping of OTUs/ASVs to
#' models/genomes, this function calculated the relative abundancies of
#' auxotrophic bacteria in each sample by metabolites (e.g. amino acids).
#'
#' @param mic A Microbiome object.
#' @param model.dir Path to the directory, where gapseq model prediction files
#' are stored. Files should be stored in model-named subdirectories.
#'
#' @return A list. TODO: provide results
#'
#' @import parallel
#'
#' @export
auxotrophy_frequencies <- function(mic,
                                   model.dir,
                                   compounds = NULL,
                                   min.growth = 0.005,
                                   min.growth.fraction = 0.05,
                                   n.cores = NULL) {
  if(length(mic@model.table) == 0)
    stop("Model count table not yet generated.")

  if(is.null(n.cores))
    n.cores <- min(10, detectCores()-1)

  cl <- makeCluster(n.cores)
  clusterExport(cl, "model.dir", envir=environment())

  mod_ids <- rownames(mic@model.table)[-1]

  models <- parLapply(cl, mod_ids, fun = function(x) {
    mod <- readRDS(paste0(model.dir, x,"/",x,".RDS"))
  })
  names(models) <- mod_ids

  model_auxo <- parLapply(cl, models, function(mod) {
    auxtmp <- predict_auxotrophies(mod, compounds, min.growth,
                                   min.growth.fraction)
    return(auxtmp)
  })

  model_auxo <- rbindlist(lapply(model_auxo, FUN = function(x) {
    data.table(compound = names(x),
               prototroph = x)
  }), idcol = "model")
  model_auxo[is.na(prototroph), prototroph := 1]

  auxo_mat <- dcast(model_auxo, formula = model ~ compound, value.var = "prototroph")
  tmp_names <- auxo_mat$model
  auxo_mat <- as.matrix(auxo_mat[,-1])
  rownames(auxo_mat) <- tmp_names

  model_cts_mat <- mic@model.table[-1,]
  n_mapped <- colSums(model_cts_mat)
  auxo_mat <- auxo_mat[rownames(model_cts_mat),]

  auxo_freq <- t(as.matrix(as.data.frame.matrix(model_cts_mat))) %*% auxo_mat
  auxo_freq <- t(auxo_freq/n_mapped)

  dt.auxo <- data.table(as.table(auxo_freq))
  colnames(dt.auxo) <- c("compound", "sample","freq")
  dt.auxo <- merge(dt.auxo, mic@sample.description)

  stopCluster(cl)

  return(list(model.auxotrophies = auxo_mat,
              proto.freq.mat = auxo_freq,
              proto.freq.dt  = dt.auxo))
}
