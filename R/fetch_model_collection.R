#' @title Fetch a collection of gapseq models
#'
#' @description Fetches a collection of genome-scale metabolic models
#' reconstructed by gapseq.
#'
#' @param model.dir character. Path to the directory that contains the models'
#' RDS files.
#' @param IDs character vector that specifies which models should be retrieved
#' from the collection. If `NULL` all models from the collection will be
#' retrieved.
#' @param file.type character. Select which gapseq output file should be read.
#' One of 'model' (gap-filled model), 'draft' (draft network), 'reactions',
#' 'pathways', 'transporters', 'medium' (predicted growth medium).
#' @param subset integer. For testing purposes you can choose the maximum number
#' of model files to be read.
#' @param entries In case 'file.type' is "pathways" or "reactions", the argument
#' can be used to limit the output to specific pathways or reactions. If `NULL`,
#' all entries are returned.
#' @param multi.thread logical. Indicating if parallel processing of models is
#' used.
#' @param ncores integer. Number of CPUs that are used in case of parallel
#' processing. If NULL, the number of available CPUs is detected.
#'
#' @return A named (IDs) list of with elements of class `modelorg` (in case of
#' file.type is 'model' or 'draft') or elements of class `data.table` otherwise.
#'
#' @import sybil
#'
#' @export
fetch_model_collection <- function(model.dir, IDs = NULL, file.type = "model",
                                   subset = NULL, entries = NULL,
                                   multi.thread = TRUE,
                                   ncores = NULL) {
  # subset argument for debugging
  if(is.null(subset) || subset < 1)
    subset <- Inf

  #model.dir <- "/mnt/nuuk/2021/HRGM/models/"

  if(file.type == "model") {
    mod.files <- dir(model.dir, recursive = T, pattern = "\\.RDS$", full.names = T)
    mod.files <- mod.files[!grepl("-rxnWeights\\.RDS$|-rxnXgenes\\.RDS$|-draft\\.RDS$",
                                  mod.files)]

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("\\.RDS$","",mod.names)
  }
  if(file.type == "draft") {
    mod.files <- dir(model.dir, recursive = T, pattern = "-draft\\.RDS$", full.names = T)

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("-draft\\.RDS$","",mod.names)
  }
  if(file.type == "reactions") {
    mod.files <- dir(model.dir, recursive = T, pattern = "-all-Reactions\\.tbl", full.names = T)

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("-all-Reactions\\.tbl.*$","",mod.names)
  }
  if(file.type == "pathways") {
    mod.files <- dir(model.dir, recursive = T, pattern = "-all-Pathways\\.tbl", full.names = T)

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("-all-Pathways\\.tbl.*$","",mod.names)
  }
  if(file.type == "transporters") {
    mod.files <- dir(model.dir, recursive = T, pattern = "-Transporter\\.tbl", full.names = T)

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("-Transporter\\.tbl.*$","",mod.names)
  }
  if(file.type == "medium") {
    mod.files <- dir(model.dir, recursive = T, pattern = "-medium\\.csv", full.names = T)

    mod.names <- gsub("^.*/","",mod.files)
    mod.names <- gsub("-medium\\.csv.*$","",mod.names)
  }

  if(length(mod.names) == 0)
    stop(paste0("No model files('", file.type, "') found."))

  if(any(duplicated(mod.names)))
    stop("Duplicated model names. Please check your model path ('model.dir').")

  names(mod.files) <- mod.names

  if(!is.null(IDs)) {
    tmp.ids <- IDs[IDs %in% mod.names]
    missing.mods <- IDs[!(IDs %in% mod.names)]

    if(length(missing.mods) == length(IDs)) {
      stop(paste0("None of the IDs refer to a model file('",file.type,
                  "') in the collection."))
    }
    if(length(missing.mods) > 0) {
      n.missing <- length(missing.mods)
      tmp.suffix <- ""
      if(n.missing > 25) {
        tmp.suffix <- paste0("... (",n.missing-25," IDs skipped)")
      }
      warning(paste0(n.missing," IDs do not have a model file('",file.type,
                     "') in the collection:\n",
                     paste(missing.mods, collapse = "\n"), "\n",
                     tmp.suffix))
    }

    mod.files <- mod.files[tmp.ids]
  }

  inds <- min(c(length(mod.files), subset))

  # parallel processing?
  n.cores <- ifelse(multi.thread, detectCores()-1, 1)
  if(!is.null(ncores))
    n.cores <- ncores
  n.cores <- min(c(n.cores, inds))
  cl <- makeCluster(max(c(1,n.cores)))
  clusterExport(cl, c("file.type","entries"), envir=environment())

  if(file.type %in% c("draft","model"))
    out <- parLapply(cl, mod.files[1:inds], fun = worker_readRDS)
  if(file.type %in% c("reactions","pathways","transporters","medium"))
    out <- parLapply(cl, mod.files[1:inds], fun = worker_fread)

  stopCluster(cl)
  return(out)
}


worker_readRDS <- function(x) {
  return(readRDS(x))
}

#' @import data.table
worker_fread <- function(x) {
  res <- fread(x)
  if(is.null(entries) | file.type %in% c("transporters","medium"))
    return(res)

  if(file.type == "pathways") {
    entries <- c(entries, paste0("|",entries,"|"))
    res <- res[ID %in% entries]
  }
  if(file.type == "reactions")
    res <- res[rxn %in% entries]
  return(res)
}
