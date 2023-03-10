#  sysBiolAlg_mtfEasyConstraintClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Copyright (C) 2014-2015 Claus Jonathan Fritzemeier, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#           definition of the class sysBiolAlg_mtfEasyConstraint               #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_mtfEasyConstraint2",
         representation(
           maxobj = "numeric",
           easyConstraint = "list"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_mtfEasyConstraint
setMethod(f = "initialize",
          signature = "sysBiolAlg_mtfEasyConstraint2",
          definition = function(.Object,
                                model,
                                wtobj = NULL,
                                react = NULL, lb = NULL, ub = NULL,
                                costcoefbw = NULL,
                                costcoeffw = NULL,
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                easyConstraint = NULL,
                                writeProbToFileName = NULL,
                                pFBAcoeff = 1e-6,
                                ...) {

            if ( ! missing(model) ) {
              stopifnot(is(model, "modelorg"),
                        is(absMAX, "numeric"))

              # If wtobj is longer than 1, mtf algorithm has to run several
              # times. In that case, wtobj is not written in the problem
              # object, it is written separately (maxobj) and used for
              # each iteration.

              #  the problem: minimize:
              #
              #            |      |      |
              #         S  |  0   |  0   |  = b
              #            |      |      |
              #       -------------------------
              #            |      |      |
              #         1  |  1   |  0   | >= 0
              #            |      |      |
              #       -------------------------
              #            |      |      |
              #         -1 |  0   |  1   | >= 0
              #            |      |      |
              #       -------------------------
              #       c_wt |  0   |  0   | >= c^T * v_wt
              #            |      |      |
              #  lb   wt_lb|  0   |  0   |
              #  ub   wt_ub|10000 |10000 |
              #            |      |      |
              #  obj    0  |  1   |  1   |


              # ---------------------------------------------
              # problem dimensions
              # ---------------------------------------------

              nc     <- react_num(model)
              nr     <- met_num(model)

              nCols  <- 3*nc
              nRows  <- nr + 2*nc + 1


              # ---------------------------------------------
              # constraint matrix
              # ---------------------------------------------

              # the initial matrix dimensions
              LHS <- Matrix::Matrix(0,
                                    nrow = nRows,
                                    ncol = nCols,
                                    sparse = TRUE)

              # rows for the mutant strain
              I <- Matrix::which(model@S!=0, arr.ind = T)
              LHS[I] <- S(model)[I]
              #LHS[1:nr,1:nc] <- S(model)

              # location of the mutant strain
              fi <- c(1:nc)

              # rows for the delta match matrix
              I <- matrix(c((nr+1)   :(nr+nc)  ,1       :nc    ),ncol=2); LHS[I] <- 1
              I <- matrix(c((nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)),ncol=2); LHS[I] <- 1
              I <- matrix(c((nr+nc+1):(nr+2*nc),1       :nc    ),ncol=2); LHS[I] <- -1
              I <- matrix(c((nr+nc+1):(nr+2*nc),(2*nc+1):(3*nc)),ncol=2); LHS[I] <- 1
              # diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <-  1
              # diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
              # diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
              # diag(LHS[(nr+nc+1):(nr+2*nc),(2*nc+1):(3*nc)]) <-  1

              # fix the value of the objective function
              ##print(obj_coef(model))
              ##print(maxobj)
              ##print(currmo)
              ##LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)
              LHS[(nr+2*nc+1),1:nc] <- 0


              # ---------------------------------------------
              # columns
              # ---------------------------------------------

              lower  <- c(lowbnd(model), rep(0, 2*nc))
              upper  <- c(uppbnd(model), rep(absMAX, 2*nc))


              # ---------------------------------------------
              # rows
              # ---------------------------------------------

              #rlower <- c(rhs(model), rep(0, 2*nc), currmo)
              #rupper <- c(rhs(model), rep(absMAX, 2*nc + 1))
              ##rlower <- c(rep(0, nr), rep(0, 2*nc), currmo)
              rlower <- c(rep(0, nr), rep(0, 2*nc), 0)
              rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1))
              rtype  <- c(rep("E", nr), rep("L", 2*nc + 1))

              # ---------------------------------------------
              # objective function
              # ---------------------------------------------

              if (is.null(costcoeffw)) {
                fw <- rep(1, nc)
              }
              else {
                stopifnot(is(costcoeffw, "numeric"),
                          (length(costcoeffw) == nc))
                fw <- costcoeffw
              }

              if (is.null(costcoefbw)) {
                bw <- fw
              }
              else {
                stopifnot(is(costcoefbw, "numeric"),
                          (length(costcoefbw) == nc))
                bw <- costcoefbw
              }

              ##cobj <- c(rep(0, nc), bw, fw)
              cobj <- c(obj_coef(model), -bw*pFBAcoeff, -fw*pFBAcoeff)

              # ---------------------------------------------
              # row and column names for the problem object
              # ---------------------------------------------

              if (isTRUE(useNames)) {
                if (is.null(cnames)) {
                  cn <- c(react_id(model),
                          paste("bw", react_id(model), sep = "_"),
                          paste("fw", react_id(model), sep = "_")
                  )
                  colNames <- .makeLPcompatible(cn, prefix = "x")
                }
                else {
                  stopifnot(is(cnames, "character"),
                            length(cnames) == nCols)
                  colNames <- cnames
                }

                if (is.null(rnames)) {
                  rn <- c(met_id(model),
                          paste("bw", 1:nc, sep = "_"),
                          paste("fw", 1:nc, sep = "_"),
                          "obj_wt"
                  )
                  rowNames <- .makeLPcompatible(rn, prefix = "r")
                }
                else {
                  stopifnot(is(rnames, "character"),
                            length(rnames) == nRows)
                  rowNames <- rnames
                }

                if (is.null(pname)) {
                  probName <- .makeLPcompatible(
                    paste("MTF", mod_id(model), sep = "_"),
                    prefix = "")
                }
                else {
                  stopifnot(is(pname, "character"),
                            length(pname) == 1)
                  probName <- pname
                }
              }
              else {
                colNames <- NULL
                rowNames <- NULL
                probName <- NULL
              }

              #add easyConstraints:
              if(!is.null(easyConstraint)){
                if(		length(easyConstraint$react) != length(easyConstraint$x)
                     | 	length(easyConstraint$react) != length(easyConstraint$rtype)
                ){
                  stop("easyConstraint elements have to have equal lengths")
                }
                stopifnot(is.list(easyConstraint$react))
                stopifnot(is.list(easyConstraint$x))
                stopifnot(all(easyConstraint$rtype %in% c("F", "L", "U", "D", "E")))

                # setting and checking rlb
                if(is.null(easyConstraint$lb)){
                  rlower <- c(rlower, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$lb)){
                    stop("easyConstraint$lb length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$lb))
                    rlower <- c(rlower, easyConstraint$lb)
                  }
                }

                # setting and checking rub
                if(is.null(easyConstraint$ub)){
                  rupper <- c(rupper, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$ub)){
                    stop("easyConstraint$ub length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$ub))
                    rupper <- c(rupper, easyConstraint$ub)
                  }
                }

                m <- Matrix(0, ncol=nCols, nrow=length(easyConstraint$react))

                I <- matrix(0, ncol = 2, nrow = 2*length(easyConstraint$react))
                I[,1] <- rep(1:length(easyConstraint$react), each = 2)
                I[,2] <- unlist(easyConstraint$react, use.names = F)

                m[I] <- unlist(easyConstraint$x, use.names = F)

                #for(i in 1:length(easyConstraint$react)){
                #  m[i, easyConstraint$react[[i]]] <- easyConstraint$x[[i]]
                #}


                LHS <- rbind2(LHS, m)
                rtype <- c(rtype, easyConstraint$rtype)
                nRows <- nRows + length(easyConstraint$react)
                if(!is.null(rowNames)){
                  rowNames <- c(rowNames, paste0("easyConstraint", 1:length(easyConstraint$react)))
                }

              }

              # ---------------------------------------------
              # build problem object
              # ---------------------------------------------

              .Object <- callNextMethod(.Object,
                                        sbalg      = "mtf",
                                        pType      = "lp",
                                        scaling    = scaling,
                                        fi         = fi,
                                        nCols      = nCols,
                                        nRows      = nRows,
                                        mat        = LHS,
                                        ub         = upper,
                                        lb         = lower,
                                        obj        = cobj,
                                        rlb        = rlower,
                                        rub        = rupper,
                                        rtype      = rtype,
                                        lpdir      = "max",
                                        ctype      = NULL,
                                        cnames     = colNames,
                                        rnames     = rowNames,
                                        pname      = probName,
                                        algPar     = list("wtobj" = wtobj,
                                                          "costcoefbw" = bw,
                                                          "costcoeffw" = fw),
                                        ...)

              ##.Object@maxobj <- as.numeric(maxobj)
              .Object@easyConstraint <- easyConstraint

              if (!is.null(writeProbToFileName)) {
                writeProb(problem(.Object),
                          fname = as.character(writeProbToFileName))
              }
              #
              #                  # ---------------------------------------------
              #                  # build problem object
              #                  # ---------------------------------------------
              #
              #                  lp <- optObj(solver = solver, method = method)
              #                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
              #
              #                  # ---------------------------------------------
              #                  # set control parameters
              #                  # ---------------------------------------------
              #
              #                  if (!any(is.na(solverParm))) {
              #                      setSolverParm(lp, solverParm)
              #                  }
              #
              #
              #                  loadLPprob(lp,
              #                             nCols = nCols,
              #                             nRows = nRows,
              #                             mat   = LHS,
              #                             ub    = upper,
              #                             lb    = lower,
              #                             obj   = cobj,
              #                             rlb   = rlower,
              #                             rub   = rupper,
              #                             rtype = rtype,
              #                             lpdir = "min"
              #                  )
              #
              #                  if (!is.null(scaling)) {
              #                      scaleProb(lp, scaling)
              #                  }
              #
              #                  .Object@problem   <- lp
              #                  .Object@algorithm <- "mtf"
              #                  .Object@nr        <- as.integer(nRows)
              #                  .Object@nc        <- as.integer(nCols)
              #                  .Object@fldind    <- as.integer(fi)
              #                  validObject(.Object)

            }
            return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                                other methods                                 #
#------------------------------------------------------------------------------#

setMethod("changeMaxObj", signature(object = "sysBiolAlg_mtf"),
          function(object, j) {

            if (!is.null(object@maxobj)) {
              changeRowsBnds(problem(object), i = nr(object),
                             lb = object@maxobj[j], ub = SYBIL_SETTINGS("MAXIMUM"))
            }

            return(invisible(TRUE))
          }
)
