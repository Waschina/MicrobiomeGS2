#' @title Filter amplicon seq. to model mappings
#'
#' @param object A Microbiome object.
#' @param allowed.identity.deviation double. Allowed deviation in itentity
#' between best hit and subobtimal hits.
#' @param method.resolve.multiple Method to use, when an ASV/OTU has several
#' optimal hits with the same statistics. Available methods: "user", "random",
#' "first".
#'
#' @export
filter_mapping <- function(object,allowed.identity.deviation=0,method.resolve.multiple='user') {
  a <- object@model.mapping

  # exclude sub-obtimal hits (allow 0.2% deviation)
  n <- nrow(a)
  a[,maxID := max(identity),by=query.label]
  a <- a[identity>=maxID-allowed.identity.deviation]
  a[,maxID := NULL]
  cat(paste0("Sub-optimal hits removed:\t",n-nrow(a),"\n"))

  # Remove muplidple hits of one 16S sequence in one model organism (due to multiple 16S gene copies)
  n <- nrow(a)
  #a <- a[order(query.label,-rank(identity))]
  a <- a[!duplicated(paste(a$query.label,a$target.model)),]
  cat(paste0("Removed mult. hits per org.:\t",n-nrow(a),"\n"))

  # at this stage: save all the mappings for potential essemble analysis
  object@model.mapping.ensemble <- a

  # in case that there are several optimal hits for a single uniq seq choose the model that occurs
  # most often across all samples
  n <- nrow(a)
  b <- data.table(uniq.count=rowSums(object@uniq.table),uniq.label=rownames(object@uniq.table))
  a <- merge(a,b,by.x="query.label",by.y="uniq.label")
  a[,model.count := sum(uniq.count),by=target.model]
  a[,max.model.count := max(model.count),by=query.label]
  a <- a[model.count>max.model.count*0.999]
  a[,uniq.count := NULL]
  a[,model.count := NULL]
  a[,max.model.count := NULL]
  cat(paste0("Removed unusual hits:\t\t",n-nrow(a),"\n"))

  # if OTU -> targe.model is not unique. Test if one of the target.models is unique to another OTU.
  # if yes, discard the other hits
  n <- nrow(a)
  a[,nr.hits := .N, by=query.label]
  unique.models <- a[nr.hits ==1, target.model]
  a[,is.um := target.model %in% unique.models]
  a[, has.um.hit := any(is.um), by=query.label]
  a <- a[is.um == T | (is.um == F & has.um.hit == F)]
  a[,nr.hits := NULL]
  a[,is.um := NULL]
  a[,has.um.hit := NULL]
  cat(paste0("Removed rare clique hits:\t",n-nrow(a),"\n"))

  # Calculating cliques (=strains, which match 16Sreads from microbiome samples with the same identity)
  #zut <- mic@model.mapping.single
  a <- a[order(query.label,target.model)]
  a[,clique := paste(target.model,sep = "$",collapse = "$"),by=query.label]
  a[,clique_org_names := NA]

  a$clique_size <- str_count(a$clique,"\\$")
  #clique.dt <- a[,15:16]
  clique.dt <- a[,.(clique, clique_size, clique_org_names)]
  clique.dt <- clique.dt[!duplicated(clique)]
  cat(paste0("Number of distinct cliques:\t",nrow(clique.dt),"\n"))
  # decide on type strain
  clique.dt$type.strain = NA
  mn <- sum(clique.dt$clique_size > 0)
  mi <- 1
  # if(method.resolve.multiple=="assembly.stats"){
  #   as.stat <- fread("/nuuk/2018/bacref/metadata/assembly_summary.txt", skip = 1)
  #   colnames(as.stat)[1] <- "assembly_accession"
  #   as.stat <- as.stat[assembly_level == "Complete Genome"]
  # }
  for(i in 1:nrow(clique.dt)) {
    if(clique.dt$clique_size[i]==0)
      clique.dt$type.strain[i] <- clique.dt$clique[i]
    else {
      tmp  <- unlist(str_split(clique.dt$clique[i],"\\$"))

      tmp2 <- paste0(1:length(tmp),": ",tmp)

      if(method.resolve.multiple=="user") {
        cat("\n")
        cat(paste0(tmp2,collapse = "\n"))

        # get user input
        k <- readline(prompt=paste("Enter Nr. or type strain for clique ",mi,"/",mn,": "))
        k <- as.integer(k)
      }
      if(method.resolve.multiple=="random") {
        k <- sample(1:length(tmp),size=1)
      }
      if(method.resolve.multiple=="first") {
        k <- 1
      }

      mi <- mi + 1
      clique.dt$type.strain[i] <- tmp[k]
    }
  }
  clique.dt[,clique_size := NULL]

  a <- merge(a,clique.dt,by="clique")

  a <- a[order(-(target.model==type.strain),target.label)]
  a <- a[!duplicated(query.label)]
  a$target.model <- a$type.strain

  object@model.mapping.single <- a


  return(object)
}
