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
#'
#' @return A named (IDs) list of with elements of class `modelorg` (in case of
#' file.type is 'model' or 'draft') or elements of class `data.table` otherwise.
#'
#' @import sybil
#'
#' @export
fetch_model_collection <- function(model.dir, IDs = NULL, file.type = "model",
                                   subset = NULL) {
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

  if(file.type %in% c("draft","model"))
    out <- lapply(mod.files[1:inds], FUN = readRDS)
  if(file.type %in% c("reactions","pathways","transporters","medium"))
    out <- lapply(mod.files[1:inds], FUN = fread)

  return(out)
}
