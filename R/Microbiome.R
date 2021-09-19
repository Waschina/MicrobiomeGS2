######################################
#                                    #
# Microbiome data set object (slim)  #
#                                    #
######################################

#' An S4 class to represent a Microbiome.
#'
#' @slot uniq.table.file character. Path to OTU/ASV/ESV table.
#' @slot model.mapping.file character. Path to blast m8-output of amplicon
#' sequences to reference 16S database (e.g. HRGM-16S, RefSeq, UHGG, ...).
#' @slot sample.description.file character. Path to sample meta information
#' table. Should be a csv/tsv file with column names, while the first column
#' needs to be names "sample" and should correspond to the column names in
#' 'uniq.table.file'.
#' @slot model.mapping data.table. Data table structure to store the assignment
#' of amplicon sequence identifiers to reference genome/model IDs.
#' @slot model.mapping.single data.table. Same as 'model.mapping', but each
#' amplicpon sequence is mapped only to a single reference genome/model.
#' @slot model.mapping.ensemble data.table TODO
#' @slot uniq.table table. OTU/ASV/ESV count table.
#' @slot model.table table. Aggregated reference model/genome count table.
#' @slot sample.description data.table. Samples' meta information.
#' @slot copy.nr.corrected logical. Indicated, whether count tables are
#' corrected by amplicon gene copy number.
#'
#' @exportClass Microbiome
#' @aliases Microbiome
Microbiome <- setClass(
  # Set the name for the class
  "Microbiome",

  # Define the slots
  slots = c(
    uniq.table.file = "character",
    model.mapping.file = "character",
    sample.description.file = "character",

    model.mapping = "data.table",
    model.mapping.single = "data.table",
    model.mapping.ensemble = "data.table",
    uniq.table = "table",
    model.table = "table",
    sample.description = "data.table",
    copy.nr.corrected = "logical"

  ),

  prototype = list(
    uniq.table.file = NA_character_,
    model.mapping.file = NA_character_,
    sample.description.file = NA_character_,

    model.mapping = data.table(),
    model.mapping.single = data.table(),
    model.mapping.ensemble = data.table(),
    uniq.table = table(NA),
    model.table = table(NA),
    sample.description = data.table(),
    copy.nr.corrected = logical(1)
  )
)

#' @title Initialize a new Microbiome object
#'
#' @description
#' Create a new Microbiome object.
#' @param uniq.table.file Path to OTU/ASV/ESV table.
#' @param model.mapping.file Path to blast m8-output of amplicon sequences to
#' reference 16S database (e.g. HRGM-16S, RefSeq, UHGG, ...).
#' @param sample.description.file Path to sample meta information table. Should
#' be a csv/tsv file with column names, while the first column needs to be names
#' "sample" and should correspond to the column names in 'uniq.table.file'. If
#' 'NULL', sample names are taken from colnames of 'uniq.table.file', without
#' any additional meta info.
#' @param uniq.table.format File format of 'uniq.table.file'. See details for
#' supported formats.
#' @return A new `Microbiome` object.
setMethod(
  f="initialize",
  signature="Microbiome",
  definition=function(.Object,
                      uniq.table.file,
                      model.mapping.file,
                      sample.description.file = NULL,
                      uniq.table.format = "mothur"){
    require(data.table)
    cat("initiating MicrobiomeGS object \n")

    # Get Unique sequence-Model mapping
    .Object@model.mapping <- fread(model.mapping.file,header=F,stringsAsFactors=F,sep = "\t")
    colnames(.Object@model.mapping) <- c('query.label','target.label','identity','alignment.length',
                                         'nr.mismatches','nr.gap.opens','start.pos.query','end.pos.query',
                                         'start.pos.target','end.pos.target','e.value','bit.score')
    rexp <- "([^\\|]*)\\|?(.*)$"
    .Object@model.mapping$target.model <- sub(rexp,"\\1",.Object@model.mapping$target.label) # Model organisms ID
    .Object@model.mapping$target.gene <- sub(rexp,"\\2",.Object@model.mapping$target.label) # Model organisms' 16S gene ID (or position on scaffold/contig)

    # Get Unique sequence table
    if(uniq.table.format == "mothur") {
      tmp <- fread(uniq.table.file,header=T,stringsAsFactors = F)
      .Object@uniq.table <- as.table(as.matrix(tmp[,-(1:2)]))
      colnames(.Object@uniq.table) <- colnames(tmp)[-(1:2)]
      rownames(.Object@uniq.table) <- tmp[[1]]
    }
    if(uniq.table.format == "R_table") {
      .Object@uniq.table <- as.table(as.matrix(read.table(uniq.table.file, check.names = F)))
      tmp <- cbind(data.table(Representative_Sequence = rownames(.Object@uniq.table),
                              total = rowSums(.Object@uniq.table)),
                   data.table(data.table(as.data.frame.matrix(.Object@uniq.table))))
    }


    cat(paste0("Number of sequences:\t\t",sum(tmp$total),"\n"))
    cat(paste0("Number of unique sequences:\t",nrow(tmp),"\n"))

    tmp <- merge(tmp[,1:2],.Object@model.mapping,by.x="Representative_Sequence",by.y="query.label",all.x=T)
    tmp <- tmp[!duplicated(Representative_Sequence)]

    cat(paste0("Per. unique sequences mapped:\t",round(tmp[!is.na(target.model),.N]/nrow(tmp)*100,digits = 1),"%\n"))
    cat(paste0("Per. sequences mapped:\t\t",
               round(tmp[!is.na(target.model),sum(total)]/tmp[,sum(total)]*100,digits = 1),
               "%\n"))

    # get sample descriptions (one column should be named "sample")
    if(!is.null(sample.description.file)) {
      .Object@sample.description <- fread(sample.description.file)
    } else {
      .Object@sample.description <- data.table(sample = colnames(.Object@uniq.table))
    }
    .Object@sample.description <- .Object@sample.description[order(as.character(sample))]

    a <- .Object@sample.description$sample
    b <- colnames(.Object@uniq.table)
    ab <- a[!(a %in% b)]
    ba <- b[!(b %in% a)]

    # samples in meta.data table but without sequences
    if(length(ab)>0) {
      .Object@sample.description <- .Object@sample.description[!(sample %in% ab)]
      warning(paste("Following samples occur in sample description file but have no associated sequences/counts:\n",
                    paste(ab,collapse = ", "),"\nSamples are excluded from sample description table."))
    }

    # samples with sequences but without meta.data
    if(length(ba)>0) {
      .Object@uniq.table <- .Object@uniq.table[,!(colnames(.Object@uniq.table) %in% ba)]
      warning(paste("Following samples do not occur in sample description but have sequences/count information:\n",
                    paste(ba,collapse = ", "),"\nSamples are excluded."))
    }


    cat(paste0("Number of samples:\t\t",nrow(.Object@sample.description),"\n"))

    .Object@copy.nr.corrected <- FALSE

    return(.Object)
  }
)
