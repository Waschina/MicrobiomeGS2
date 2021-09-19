#' @title Transform ASV/OTU count table to model count table
#'
#' @param object A Microbiome object.
#'
#' @export
create_model_table <- function(object) {
  rel.ag <- levels(factor(object@model.mapping.single$target.model))
  a <- matrix(0, nrow = length(rel.ag)+1, ncol = ncol(object@uniq.table))
  rownames(a) <- c("_unclassified",rel.ag)
  colnames(a) <- colnames(object@uniq.table)

  for(i in 2:nrow(a)) {
    cat(paste0("\r Creating Model-Table: ",i,"/",nrow(a)))
    ag <- rel.ag[i-1]
    qu <- object@model.mapping.single[target.model==ag,query.label]
    if(length(qu)>1) {
      a[i,] <- colSums(object@uniq.table[qu,])
    } else {
      a[i,] <- object@uniq.table[qu,]
    }
  }
  cat("\n")
  # get number of unclassified OTUs
  a[1,] <- colSums(object@uniq.table) - colSums(a)

  object@model.table <- as.table(a)

  return(object)
}
