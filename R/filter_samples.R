#' @title Filter samples
#'
#' @param mic A Microbiome object.
#' @param min.seqs integer. Minimum count of total ASVs/OTUs per sample.
#' @param max.unclassified double. Maximum proportion of counts, that correspond
#' to unclassified ASVs/OTUs, i.e. not mapped to any reference genome/model.
#'
#' @export
filter_samples <- function(object,min.seqs=5000,max.unclassified=0.3) {
  cat(paste0("Excluding samples w/ less than ",min.seqs," sequences and with >=",max.unclassified*100,
             "% unclassified sequences. \n"))
  # proportion of unclassfied seqs
  ex.unclass <- names(which(object@model.table["_unclassified",]/colSums(object@model.table) >= max.unclassified))

  # number of seqs per samples
  ex.nrseqs <- names(which(colSums(object@model.table) < min.seqs))

  # samples to explude
  excl <- c(ex.unclass,ex.nrseqs)
  excl <- excl[!duplicated(excl)]

  print(excl)

  # updating object
  object@model.table <- object@model.table[,!(colnames(object@model.table) %in% excl)]
  object@uniq.table <- object@uniq.table[,!(colnames(object@uniq.table) %in% excl)]
  object@sample.description <- object@sample.description[!(sample %in% excl)]

  # model.table
  # unique.table
  # sample.description

  cat(paste0("Number of samples excluded: ",length(excl),"\n"))
  cat(paste0("Remaining samples: ",nrow(object@sample.description),"\n"))
  return(object)
}
