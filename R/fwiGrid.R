#' @title Fire Weather Index (FWI)
#' 
#' @description Implementation of the Canadian Fire Weather Index System 
#' 
#' @param multigrid A multigrid of the variables needed to compute the FWI.
#' @param mask Optional. Grid of the land mask to be applied to the data.
#' @param return.all Logical. Should all components of the FWI system be returned?. 
#' Default to FALSE, indicating that only FWI is returned.
#' @param init.pars A numeric vector of length 3 with the initialization values for the
#'  FFMC, DMC and DC components, in this order. Default values as proposed by van Wagner (1987).
#' @param parallel Logical. Should parallel execution be used?
#' @param max.ncores Integer. Upper bound for user-defined number of cores.
#' @param ncores Integer number of cores used in parallel computation. Self-selected number of
#'  cores is used when \code{ncpus = NULL} (the default), or when \code{maxcores} exceeds the default \code{ncores} value.
#' 
#' @return A grid, containing the requested components of the FWI system (either all or just FWI). See details.
#' 
#' @section Daylength adjustment factors: 
#' By default, the function applies the original FWI daylength adjustment factors for DC and DMC (van Wagner 1987),
#'  although it is possible to adjust them by as a function of latitude through the argument \code{lat}.
#' The reference values used for each latitudinal range are those indicated in p.71 and Tables A3.1 and A3.2 (p69) in
#' Lawson and Armitage (2008).
#' 
#' @references
#' Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' 
#' @author J. Bedia \& M.Iturbide, partially based on the original FORTRAN code by van Wagner and Pickett (1985)
#' @export
#' @importFrom abind abind
#' @importFrom downscaleR subsetGrid
#' @importFrom downscaleR getYearsAsINDEX



fwiGrid <- function(multigrid,
                     mask = NULL,
                     return.all = FALSE, 
                     init.pars = c(85, 6, 15),
                     parallel = FALSE,
                     max.ncores = 16,
                     ncores = NULL){
      x <-  c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
      
      latLim <- range(multigrid$xyCoords$y)
      lonLim <- range(multigrid$xyCoords$x)
      latind <- findInterval(latLim, x)[1] : findInterval(latLim, x)[2]
      if(x[latind[length(latind)]] < latLim[2]) latind[3] <- latind[2]+1
      x <- x[latind]
      lats <- seq(min(x)+5, max(x)-5, 10) 
      if(x[length(x)] > latLim[2]) x[length(x)] <- latLim[2]
      if(x[1] < latLim[1]) x[1] <- latLim[1]
      
      dataset <- attr(multigrid, "dataset")
      years <- unique(getYearsAsINDEX(multigrid))
      latdim <- which(downscaleR:::getDim(multigrid) == "lat")
      a <- list()
      message("[", Sys.time(), "] Calculating FWI..")
      for(i in 1:(length(x)-1)){
            latLimchunk <- c(x[i], x[i+1])
            lat <- lats[i]
            multigrid_chunk <- subsetGrid(multigrid, latLim = latLimchunk)
            if(is.null(mask) & dataset == "WFDEI"){
                  msk <- subsetGrid(multigrid_chunk, var = "tas")
                  msk$Data <- msk$Data[1,,]
                  msk$Data[which(!is.na(msk$Data))] <- 100
                  msk$Data[which(is.na(msk$Data))] <- 0
                  attr(msk$Data, "dimensions") <- c("lat", "lon")
            }else if(!is.null(mask)){
                  msk <- subsetGrid(mask, latLim = latLimchunk, lonLim = lonLim, outside = T)
            }else{
                  message("The use of a land mask is recommended")
            }
            if(i != (length(x)-1)){
                  xx <- dim(multigrid_chunk$Data)[latdim]-1
            }else{
                  xx <- dim(multigrid_chunk$Data)[latdim]
            }
            o <- lapply(1:length(years), function(x){
                  multigrid.y <- subsetGrid(multigrid_chunk, years = years[x])
                  suppressMessages(
                        fwi(multigrid = multigrid.y, mask = msk, lat = lat, return.all = return.all, 
                        parallel = parallel, init.pars = init.pars, 
                        max.ncores = max.ncores, ncores = ncores)$Data[,,1:xx,,drop=FALSE]
                  )
            })
            o.full <-  unname(do.call("abind", list(o, along = 2)))  
            months <- as.integer(substr(multigrid$Dates[[1]]$start, start = 6, stop = 7))
            multigrid.y <- NULL
            month.ind <- which(months == months[1])
            a[[i]] <- o.full[,-month.ind,,]
      }
      message("[", Sys.time(), "] Done.")
      temp <- subsetGrid(multigrid, var = "tas")
      dimNames <- attr(temp$Data, "dimensions")
      latdim.f <- which(downscaleR:::getDim(temp) == "lat")
      out <- unname(do.call("abind", list(a, along = latdim.f)))
      temp$Data <- out
      attr(temp$Data, "dimensions") <- dimNames
      temp$Variable$varName <- "fwi"
      attr(temp$Variable, "use_dictionary") <- FALSE
      attr(temp$Variable, "description") <- "Fire Weather Index"
      attr(temp$Variable, "units") <-  "none"
      attr(temp$Variable, "longname") <- "Fire Weather Index"
      #       if(length(years) > 1){
      #             yname <- paste(years[1], "_", years[length(years)], sep = "")
      #       }else{
      #             yname <- years
      #       }
      return(temp)
}