#     fffdiGrid.R Finnish Forest Fire Danger Index
#
#     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Finnish Forest Fire Danger Index (FFFDI)
#' @description Implementation of the FFFDI for climate4R grids
#' @param pr A climate4R object containing daily precipitation (in mm)
#' @param pet A climate4R object containing daily (potential) evapotranspiration data (in mm).
#' @param Wvol.init Initialization value for volumetric moisture, in the range 0.1-0.5.
#'  Default to 0.5 (very wet soil), but see Details.
#' @param z reference surface layer thickness (mm). Default to 60.
#' @references Vajda, A., Venalainen, A., Suomi, I., Junila, P. and Makela, H., 2014. Assessment of forest 
#' fire danger in a boreal forest environment: description and evaluation of the operational 
#' system applied in Finland. Meteorol. Appl., 21: 879-887, DOI: 10.1002/met.1425
#' @importFrom transformeR gridArithmetics getDim getShape
#' @importFrom abind asub abind
#' @importFrom transformeR checkDim checkTemporalConsistency redim getShape subsetGrid array3Dto2Dmat getCoordinates
#' @importFrom magrittr %>% %<>% extract2
#' @return A climate4R object containing FFFDI data
#' @details 
#' \strong{Volumetric moisture}
#' The default is 0.5, indicating that the soil is very wet and near field capacity. This is so, assuming that the
#'  index is started in early spring. This value is applied to all locations as a spatially constant initialization value.
#'  However, Vajda \emph{et al.} (2014, Table 1) provide reference values for different soil moisture conditions. This value ranges
#'  from 0.1 (very dry) to 0.5 (very wet).
#' @author J. Bedia (based on the original FORTRAN code from Vajda et al. 2014).
#' @seealso \code{petGrid} function, from package \pkg{drought4R} 
#' is recommended for PET estimation <https://github.com/SantanderMetGroup/drought4R>.
#' @export


fffdiGrid <- function(pr, pet, Wvol.init = 0.5, z = 60) {
    pr %<>% redim()
    pet %<>% redim()
    suppressMessages(checkDim(pr, pet))
    checkTemporalConsistency(pr, pet)
    if ("var" %in% getDim(pet)) stop("Input multigrids are not allowed")
    ntimes <- getShape(pr, "time")
    nmem <- getShape(pr, "member")
    npoints <- getCoordinates(pet) %>% expand.grid() %>% nrow()
    aux.grid <- pet
    fffdi.list <- lapply(1:nmem, function(x) {
        Wvol <- rep(Wvol.init, npoints)
        Wt <- Wvol * z 
        pr.aux <- subsetGrid(pr, members = x) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        pet.aux <- subsetGrid(pet, members = x) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        ref <- pet.aux
        for (i in 1:ntimes) {
            dw <-  (-pet.aux[i,] * 0.757) / (1 + exp(2.74 - 16.67 * (Wvol - 0.1))) + 5.612 * (1 - exp(-pr.aux[i,] / 5.612))
            Wt <- Wt + dw
            Wvol <- Wt/z
            Wvol[which(Wvol > 0.5)] <- 0.5
            Wvol[which(Wvol < 0.1)] <- 0.1
            ref[i,] <- 8.76 - 30.879 * Wvol + 30.712 * Wvol^2
        }
        grd <- mat2Dto3Darray(ref, x = getCoordinates(pet)[["x"]], y = getCoordinates(pet)[["y"]]) 
        aux.grid$Data <- grd
        redim(aux.grid) %>% return()
    })
    out <- suppressMessages(do.call("bindGrid", c(fffdi.list, dimension = "member")))
    out$Variable$varName <- "fffdi"
    descr <- "finnish_forest_fire_danger_index"
    attr(out$Variable, "description") <- descr
    attr(out$Variable, "longname") <- descr
    attr(out$Variable, "daily_agg_cellfun") <- "index"
    return(out)
}
