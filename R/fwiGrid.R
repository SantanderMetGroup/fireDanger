#' @title Fire Weather Index applied to multigrids
#' 
#' @description Implementation of the Canadian Fire Weather Index System for multigrids
#' 
#' @param multigrid containing Tm (temperature records in deg. Celsius); H (relative humidity records in \%);
#' r (last 24-h accumulated precipitation in mm); W (wind velocity records in Km/h). See details.
#' @param mask Optional. Binary grid (0 and 1, 0 for sea areas) with \code{dimensions} attribute \code{c("lat", "lon")}.
#' @param what Character string. What component of the FWI system is computed?. Default to \code{"FWI"}.
#' Note that unlike \code{\link{fwi1D}}, only one single component choice is possible in \code{fwiGrid}.
#'   See \code{\link{fwi1D}} for details and possible values.
#' @param nlat.chunks For an efficient memory usage, the computation of FWI can be split into 
#' latitudinal bands (chunks) sepparately. The number of chunks is controlled here. 
#' Default to \code{NULL} (i.e., no chunking applied).
#' @param restart.annual Logical. Should the calculation be restarted at the beginning of every year?
#' If the grid encompasses just one season (e.g. JJA...), this is the recommended option. Default to \code{TRUE}.
#' @param ... Further arguments passed to \code{\link{fwi1D}}.
#' @template templateParallelParams
#' 
#' @return A grid corresponding to the variable defined in argument \code{what}.  
#' 
#' @details 
#' 
#' \strong{Variable names}
#' 
#' The variables composing the input multigrid are expected to have standard names, as defined by the dictionary
#'  (their names are stored in the \code{multigrid$Variable$varName} component).
#' These are: \code{"tas"} for temperature, \code{"tp"} or \code{"pr"} for precipitation, \code{"wss"} or \code{"sfcWind"} for windspeed. In the case of relative humidity,
#' either \code{"hurs"} or \code{"hursmin"} are accepted, the latter in case of FWI calculations according to the \dQuote{proxy} version
#' described in Bedia \emph{et al} 2014.
#' 
#' Note that the order of the variables within the multigrid is not relevant. These are indexed by variable names.
#' 
#' \strong{Landmask definition}
#' 
#' The use of a landsmask is highly recommended when using RCM/GCM data becasue (i) there is no point in calculating
#' FWI over sea areas and (ii) for computational efficiency, as sea gridboxes will be skipped before calculations.
#' 
#' The landmask must be a grid spatially consistent with the input multigrid. You can use 
#' \code{\link[transformeR]{interpGrid}} in combination with the \code{getGrid} method to ensure this condition is fulfilled.  . Its \code{data} component can be either a 2D or 3D array with the \code{dimensions} 
#' attribute \code{c("lat","lon")} or \code{c("time","lat","lon")} respectively. In the latter case, the length of the time 
#' dimension should be 1. Note that values of 0 correspond to sea areas (thus discarded for FWI calculation), being land areas any other 
#' values different from 0 (tipically 1 or 100, corresponding to the land/sea area fraction).   
#' 
#' \strong{Latitudinal chunking}
#' 
#' Splitting the calculation in latitudinal chunks is highly advisable, and absolutely necessary when
#' considering large spatial domains, otherwise running out of memory during the computation. The number
#' of latitudinal chunks need to be estimated on a case-by-case basis, but in general there are no restrictions in the
#' number of chunks that can be used, as long as it does not exceed the number of actual latitudes 
#' in the model grid.
#' 
#' @template templateParallel 
#' 
#' @references
#' \itemize{
#' \item Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' \item van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' \item van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' }
#' 
#' @author J. Bedia \& M. Iturbide
#' 
#' @export
#' 
#' @importFrom abind abind asub
#' @importFrom parallel parLapply splitIndices
#' @importFrom magrittr %>% 
#' @importFrom transformeR redim getDim getShape parallelCheck getYearsAsINDEX subsetGrid array3Dto2Dmat mat2Dto3Darray 

fwiGrid <- function(multigrid,
                    mask = NULL,
                    what = "FWI",
                    nlat.chunks = NULL,
                    restart.annual = TRUE,
                    parallel = FALSE,
                    ncores = NULL,
                    max.ncores = 16,
                    ...) {
      what <- match.arg(what,
                        choices = c("FFMC", "DMC", "DC", "ISI", "BUI", "FWI", "DSR"),
                        several.ok = FALSE)
      fwi1D.opt.args <- list(...)
      # fwi1D.opt.args <- list()
      months <- as.integer(substr(multigrid$Dates[[1]]$start, start = 6, stop = 7))
      fwi1D.opt.args <- c(fwi1D.opt.args, list("what" = what))
      if ("lat" %in% names(fwi1D.opt.args)) {
            message("NOTE: argument 'lat' will be overriden by the actual latitude of gridboxes\n(See help of fwi1D for details).")
            fwi1D.opt.args[-grep("lat", names(names(fwi1D.opt.args)))]
      }
      varnames <- multigrid$Variable$varName
      # Naming adjustments for RCM data
      if (any(grepl("^pr$", varnames))) varnames <- gsub("^pr$", "tp", varnames)
      if ("sfcWind" %in% varnames) varnames <- gsub("sfcWind", "wss", varnames)
      
      proj <- getGridProj(multigrid)
      ycoords <- multigrid$xyCoords$y
      co <- expand.grid(ycoords, multigrid$xyCoords$x)[2:1]
      ## In rotated grids, actual latitude must be extracted from the 'lat' component of xyCoords
      if (!is.null(proj) && proj == "RotatedPole") {
          coords_geo <- data.frame("lon" = as.vector(multigrid$xyCoords$lon),
                                   "lat" = as.vector(multigrid$xyCoords$lat))
          co <- coords_geo[order(co[,1], co[,2]), ]
          coords_geo <- NULL
          gc(verbose = FALSE)
      } 
      
      dimNames.mg <- getDim(multigrid)
      multigrid <- redim(multigrid, member = TRUE)
      n.mem <- getShape(multigrid, "member")
      
      yrsindex <- getYearsAsINDEX(multigrid)
      nyears <- length(unique(yrsindex))
      if (!is.null(mask)) {
            dimNames.mask <- getDim(mask)
      }
      if (is.null(nlat.chunks)) {
            nlat.chunks <- 1L
      }
      if (nlat.chunks <= 0) {
            nlat.chunks <- 1L
            message("Invalid 'nlat.chunks' argument value. It was ignored")
      }
      idx.chunk.list <- parallel::splitIndices(length(ycoords), nlat.chunks)
      if (any(vapply(idx.chunk.list, FUN = "length", FUN.VALUE = numeric(1)) < 2L)) {
            stop("Too many latitudinal chunks. Reduce the value of 'nlat.chunks' to a maximum of ", length(ycoords) %/% 2)
      }
      message("[", Sys.time(), "] Calculating ", what)
      aux.list <- lapply(1:nlat.chunks, function(k) {
            ## Lat chunking
            ind.lat <- idx.chunk.list[[k]]
            dims <- grep("lat", dimNames.mg)
            multigrid_chunk <- multigrid 
            mask_chunk <- mask
            if (nlat.chunks > 1) {
                  aux <- asub(multigrid$Data, idx = ind.lat, dims = dims)
                  attr(aux, "dimensions") <- dimNames.mg
                  multigrid_chunk$Data <- aux
                  multigrid_chunk$xyCoords$y <- multigrid$xyCoords$y[ind.lat]
                  ## Mask chunking
                  if (!is.null(mask)) {
                        aux <- asub(mask$Data, idx = ind.lat, dims = grep("lat", dimNames.mask), drop = FALSE)
                        attr(aux, "dimensions") <- dimNames.mask
                        mask_chunk$Data <- aux
                        mask_chunk$xyCoords$y <- mask_chunk$xyCoords$y[ind.lat]
                  }
                  aux <- NULL
            }
            ## Multigrid subsetting
            Tm1 <- subsetGrid(multigrid_chunk,
                              var = multigrid$Variable$varName[grep("tas", varnames)]) %>% redim(drop = FALSE)
            H1  <- subsetGrid(multigrid_chunk,
                              var = multigrid$Variable$varName[grep("hurs", varnames)]) %>% redim(drop = FALSE)
            r1  <- subsetGrid(multigrid_chunk,
                              var = multigrid$Variable$varName[grep("tp", varnames)]) %>% redim(drop = FALSE)
            W1  <- subsetGrid(multigrid_chunk,
                              var = multigrid$Variable$varName[grep("wss", varnames)]) %>% redim(drop = FALSE)
            multigrid_chunk <- NULL
            gc()
            ## Parallel checks
            parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
            if (n.mem < 2 && isTRUE(parallel.pars$hasparallel)) {
                  parallel.pars$hasparallel <- FALSE
                  message("NOTE: parallel computing only applies to multimember grids. The option was ignored")
            }
            if (parallel.pars$hasparallel) {
                  apply_fun <- function(...) {
                        parallel::parLapply(cl = parallel.pars$cl, ...)
                  }
                  on.exit(parallel::stopCluster(parallel.pars$cl))
            } else {
                  apply_fun <- lapply      
            }
            ## Landmask 
            if (!is.null(mask)) {
                  if (dimNames.mask[1] != "time") {
                        aux <- unname(abind(mask_chunk$Data, along = 0L))
                        attr(aux, "dimensions") <- c("time", dimNames.mask)
                  } else {
                        aux <- mask_chunk$Data
                  }
                  msk <- array3Dto2Dmat(aux)[1,]
                  ind <- which(msk > 0)
                  msk <- NULL
            } else {
                  aux <- suppressWarnings(subsetGrid(Tm1, members = 1))$Data
                  aux <- array3Dto2Dmat(aux)
                  ind <- which(apply(aux, MARGIN = 2, FUN = function(y) !all(is.na(y))))
            }
            aux <- NULL
            gc()
            
            ## FWI calculation
            if (nlat.chunks > 1) {
                message("[", Sys.time(), "] Processing chunk ", k, " out of ", nlat.chunks, "...")
            }
            a <- apply_fun(1:n.mem, function(x) {
                  Tm2 <- array3Dto2Dmat(subsetGrid(Tm1, members = x)$Data)
                  H2 <- array3Dto2Dmat(subsetGrid(H1, members = x)$Data)
                  r2 <- array3Dto2Dmat(subsetGrid(r1, members = x)$Data)
                  W2 <- array3Dto2Dmat(subsetGrid(W1, members = x)$Data)
                  b <- array(dim = dim(Tm2))
                  if (length(ind) > 0) {
                        for (i in 1:length(ind)) {
                              if (isTRUE(restart.annual)) {
                                    ## Iterate over years
                                    annual.list <- lapply(1:nyears, function(j) {
                                          idx <- which(yrsindex == unique(yrsindex)[j])
                                          arg.list2 <- list("months" = months[idx],
                                                            "Tm" = Tm2[idx,ind[i]],
                                                            "H" = H2[idx,ind[i]],
                                                            "r" = r2[idx,ind[i]],
                                                            "W" = W2[idx,ind[i]],
                                                            "lat" = co[ind[i],2])
                                          arg.list <- c(fwi1D.opt.args, arg.list2)
                                          z <- tryCatch({suppressWarnings(drop(do.call("fwi1D",
                                                                                       args = arg.list)))},
                                                        error = function(err) {
                                                            rep(NA, length(idx))
                                                        })
                                          return(z)
                                    })
                                    b[,ind[i]] <- do.call("c", annual.list)
                              } else {
                                    arg.list2 <- list("months" = months,
                                                      "Tm" = Tm2[,ind[i]],
                                                      "H" = H2[,ind[i]],
                                                      "r" = r2[,ind[i]],
                                                      "W" = W2[,ind[i]],
                                                      "lat" = co[ind[i],2])
                                    arg.list <- c(fwi1D.opt.args, arg.list2)
                                    z <- tryCatch({suppressWarnings(drop(do.call("fwi1D",
                                                                                 args = arg.list)))},
                                                  error = function(err) {
                                                      rep(NA, length(months))
                                                  })
                                    b[,ind[i]] <- z
                              }
                        }
                        out <- mat2Dto3Darray(mat2D = b,
                                       x = Tm1$xyCoords$x,
                                       y = Tm1$xyCoords$y)
                        return(out)
                  }
            })
            Tm1 <- r1 <- H1 <- W1 <- NULL
            gc()
            unname(do.call("abind", list(a, along = 0)))
      })
      message("[", Sys.time(), "] Done.")
      ## Final grid and metadata
      fwigrid <- redim(subsetGrid(multigrid, var = varnames[1]), drop = FALSE)
      multigrid <- NULL     
      dimNames <- getDim(fwigrid)
      fwigrid$Data <- unname(do.call("abind", c(aux.list, along = grep("lat", dimNames))))
      aux.list <- NULL
      attr(fwigrid$Data, "dimensions") <- dimNames
      fwigrid$Variable <- list()
      fwigrid$Variable$varName <- what 
      fwigrid$Variable$level <- NA
      desc <- switch(what,
                     "FFMC" = "Fine Fuel Moisture Code",
                     "DMC" = "Duff Moisture Code",
                     "DC" = "Drought Code",
                     "ISI" = "Initial Spread Index",
                     "BUI" = "Builtup Index",
                     "FWI" = "Fire Weather Index",
                     "DSR" = "Daily Severity Rating")
      attr(fwigrid$Variable, "use_dictionary") <- FALSE
      attr(fwigrid$Variable, "description") <- desc
      attr(fwigrid$Variable, "units") <-  "adimensional"
      attr(fwigrid$Variable, "longname") <- paste(desc, "component of the Canadian Fire Weather Index System")
      attr(fwigrid, "calculation") <- "Calculated with the fireDanger package (https://github.com/SantanderMetGroup/fireDanger)"
      return(fwigrid)
}
      
