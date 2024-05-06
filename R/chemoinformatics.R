#' @title Get neutral sum formula for charged compound
#'
#' @description Transforms a formula of a charged compound to a neutral formula
#' by removing adding protons.
#'
#' @param formula character (vector) with chemical formulas.
#' @param charge integer (vector) of same length as 'formula' specifying the
#' compounds charge.
#'
#' @export
neutralize_formula <- function(formula, charge) {

  if(length(formula) != length(charge))
    stop("Length of 'formula' and 'charge' vectors not equal.")

  # skip a few formula (returning NA for those)
  skip_inds <- integer(0L)

  skip_formula <- c("Fe","Cr","Co","Mo","Hg","Ag","Ca","Ni","Mg","Cu","Cd","Zn","Mn","H","K","Li","Na","Rb",
                    "Tl","NO2","La","Ti","Sb","Ti","V","Al","W","O2U","Be","Pt","Pb","Sr","Pd","Sn","Ba",
                    "Cs","Z","")

  skip_elements <- "As|X|x"

  skip_inds <- c(skip_inds,which(formula %in% skip_formula)) # anorganic ions
  skip_inds <- c(skip_inds,grep(skip_elements, formula)) # arsenic unspecific formula
  skip_inds <- c(skip_inds,grep("(\\)[0-9]*n)|(R)|(\\.)",formula)) # non-specific sum formula
  skip_inds <- c(skip_inds, which(is.na(charge)))# missing charge

  skip_inds <- unique(skip_inds)

  # uncharge the charged
  out_vec <- rep(NA_character_, length(formula))
  n <- length(formula)
  for(i in setdiff(1:n, skip_inds)) {
    cat("\r",i,"/", n)
    i_charge <- charge[i]

    if(!is.na(formula[i]) & !is.na(i_charge)) {
      if(i_charge == 0) {
        out_vec[i] <- formula[i]
      }
      if(i_charge < 0) {
        add <- paste0("H", -i_charge)
        out_vec[i] <- addsub_formulas(formula[i],add, mode = "add")
      }
      if(i_charge > 0) {
        add <- paste0("H", i_charge)
        out_vec[i] <- addsub_formulas(formula[i],add, mode = "sub")
      }
    }
  }
  cat("\n")
  return(out_vec)
}

#' @title Add/sub a chemical formula to another formula
#'
#' @description Add/sub a chemical formula to another formula
#'
#' @param frml1 character. Chemical formula to which something is
#' added/subtracted.
#' @param frml2 character. Chemical formula that is added/subtracted from
#' 'frml1'
#' @param mode character. Either 'add' or 'sub'
#'
#' @return character. Chemical formula of product.
#'
#' @details If an element is subtracted in a way, that the product has a
#' negative atom count, NA is returned.
#'
#' @import CHNOSZ
#'
#' @export
addsub_formulas <- function(frml1, frml2, mode = "add") {
  multmode <- 1
  if(mode == "sub")
    multmode <- -1

  mu_frml1 <- CHNOSZ::makeup(frml1)
  mu_frml2 <- CHNOSZ::makeup(frml2)

  rel_elements <- unique(c(names(mu_frml1),names(mu_frml2)))

  outvec <- rep(0, length(rel_elements))
  names(outvec) <- rel_elements

  outvec[names(mu_frml1)] <- outvec[names(mu_frml1)] + mu_frml1
  outvec[names(mu_frml2)] <- outvec[names(mu_frml2)] + mu_frml2 * multmode

  if(any(outvec < 0))
    return(NA_character_)

  return(as.chemical.formula(makeup(outvec)))
}


#' @title Calculate PPM error between two mz values
#'
#' @description Calculate PPM error between two mz values (measured vs
#' theoretical)
#'
#' @param theo.mz double. Theoretical MZ values
#' @param meas.mz double. Measured/observed MZ values
#'
#' @return double. Error in 'ppm'.
#'
#' @export
calc_ppmerr <- function(theo.mz, meas.mz) {
  ppmerr <- (meas.mz - theo.mz) / theo.mz * 1000000
  return(ppmerr)
}
