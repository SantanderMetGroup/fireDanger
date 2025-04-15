## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      fig.align = "center",
                      fig.width = 7.5,
                      fig.asp = 1,
                      fig.path = 'figure/graphics-',
                      message = FALSE,
                      tidy = FALSE,
                      cache = FALSE,
                      cache.path = 'cache/graphics-')

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github(c("SantanderMetGroup/loadeR.java@v1.1-0",
#                   "SantanderMetGroup/loadeR@v1.0-7",
#                   "SantanderMetGroup/loadeR.ECOMS@v1.2.3"))

## ------------------------------------------------------------------------
library(loadeR.ECOMS)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("SantanderMetGroup/transformeR@v0.0.8")

## ------------------------------------------------------------------------
library(transformeR)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("SantanderMetGroup/downscaleR@v2.0-1")

## ------------------------------------------------------------------------
library(downscaleR)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("SantanderMetGroup/visualizeR@v0.2-0")

## ------------------------------------------------------------------------
library(visualizeR)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("SantanderMetGroup/fireDanger@v1.0.3")

## ------------------------------------------------------------------------
library(fireDanger)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("easyVerification")

## ------------------------------------------------------------------------
library(easyVerification)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("RColorBrewer")

## ------------------------------------------------------------------------
library(RColorBrewer)

## ---- eval = FALSE-------------------------------------------------------
#  loginUDG(username = "myuser", password = "mypassword")

## ---- message = FALSE----------------------------------------------------
dataset <- "CFSv2_seasonal"
season <- 5:9
leadMonth <- 0
years <- 1983:2009
latLim <- c(35, 47)
lonLim <- c(-10, 30)
members <- 1:15

## ---- eval=FALSE---------------------------------------------------------
#  ## Load temperature
#  Tm <- loadECOMS(dataset, var = "tas", members = members,
#                  latLim = latLim, lonLim = lonLim,
#                  season = season, years = years, leadMonth = leadMonth,
#                  time = "DD", aggr.d = "mean")
#  ## Load relative humidity
#  H <- loadECOMS(dataset, var = "hurs", members = members,
#                 latLim = latLim, lonLim = lonLim,
#                 season = season, years = years, leadMonth = leadMonth,
#                 time = "DD", aggr.d = "min")
#  ## Load precipitation
#  r <- loadECOMS(dataset, var = "tp", members = members,
#                 latLim = latLim,  lonLim = lonLim,
#                 season = season, years = years, leadMonth = leadMonth,
#                 time = "DD", aggr.d = "sum")
#  ## Load wind speed
#  W <- loadECOMS(dataset, var = "wss", members = members,
#                 latLim = latLim, lonLim = lonLim,
#                 season = season, years = years, leadMonth = leadMonth,
#                 time = "DD", aggr.d = "mean")

## ---- echo = FALSE-------------------------------------------------------
load("../ignore/vignette_data/Tm_hind.rda")

## ---- results='hide', fig.keep='all', fig.align='center', message=FALSE----
plotClimatology(climatology(Tm), backdrop.theme = "coastline",
                main = "CFSv2 T2M climatology (1983-2009)")

## ---- eval=FALSE---------------------------------------------------------
#  ## Convert wss units from m/s to km/h
#  W$Data <- W$Data*3.6
#  ## Update "units" attribute
#  attr(W$Variable, "units") <- "km.h-1"
#  ## make multigrid
#  multigrid_hind <- makeMultiGrid(Tm, H, r, W)

## ---- eval=FALSE---------------------------------------------------------
#  ## Define the target dataset
#  dataset <- "WFDEI"
#  
#  ## Load temperature
#  Tm.obs <- loadECOMS(dataset = dataset, var = "tas",
#                  latLim = latLim, lonLim = lonLim,
#                  season = season, years = years)
#  ## Load relative minimum humidity
#  H.obs <- loadECOMS(dataset = dataset, var = "hursmin",
#                 latLim = latLim, lonLim = lonLim,
#                 season = season, years = years)
#  ## Load precipitation
#  r.obs <- loadECOMS(dataset = dataset, var = "tp",
#                 latLim = latLim,  lonLim = lonLim,
#                 season = season, years = years)
#  ## Load wind speed
#  W.obs <- loadECOMS(dataset = dataset, var = "wss",
#                 latLim = latLim, lonLim = lonLim,
#                 season = season, years = years)

## ---- eval=FALSE---------------------------------------------------------
#  W.obs$Data <- W.obs$Data*3.6
#  attr(W.obs$Variable, "units") <- "km.h-1"

## ---- eval=FALSE---------------------------------------------------------
#  multigrid_obs <- makeMultiGrid(Tm.obs, H.obs, r.obs, W.obs)

## ---- eval=FALSE---------------------------------------------------------
#  load(url("http://www.meteo.unican.es/work/fireDanger/wiki/data/multigrid_hind.rda"))
#  load(url("http://www.meteo.unican.es/work/fireDanger/wiki/data/multigrid_obs.rda"))

## ---- echo= FALSE--------------------------------------------------------
load("../ignore/vignette_data/multigrid_obs.rda")
load("../ignore/vignette_data/multigrid_hind.rda")

## ---- message = FALSE, eval = FALSE--------------------------------------
#  obs <- fwiGrid(multigrid = multigrid_obs)

## ---- message = FALSE, echo = FALSE--------------------------------------
load("../ignore/vignette_data/obs.rda")

## ---- eval = FALSE-------------------------------------------------------
#  cfs_mask <- loadECOMS(dataset, var = "lm", latLim = latLim, lonLim = lonLim)
#  mask <- interpGrid(cfs_mask, getGrid(multigrid_hind))

## ---- echo = FALSE, warning = FALSE, message=FALSE-----------------------
load("../ignore/vignette_data/cfs_mask.rda")
mask <- interpGrid(cfs_mask, getGrid(multigrid_hind))

## ---- message = FALSE, eval = FALSE--------------------------------------
#  hindcast <- fwiGrid(multigrid = multigrid_hind, mask = mask)

## ---- message = FALSE, echo = FALSE--------------------------------------
load("../ignore/vignette_data/hindcast.rda")

## ---- message = FALSE----------------------------------------------------
hindcastJJAS <- subsetGrid(hindcast, season = 6:9)
obsJJAS <- subsetGrid(obs, season = 6:9)


## ---- message = FALSE, fig.align='center'--------------------------------

fwi.colors <- colorRampPalette(c(rev(brewer.pal(9, "YlGnBu")[3:5]), 
                                 brewer.pal(9, "YlOrRd")[3:9]))

plotClimatology(climatology(obsJJAS), 
                backdrop.theme = "coastline", 
                col.regions = fwi.colors,
                main = "WFDEI (observations) - Mean FWI Climatology")

## ---- message = FALSE, fig.align='center'--------------------------------
plotClimatology(climatology(hindcastJJAS), 
                backdrop.theme = "coastline", 
                col.regions = fwi.colors,
                main = "CFSv2 (model) - Mean FWI Climatology")

## ---- message = FALSE, warning = FALSE-----------------------------------
hindcast_aggr <- aggregateGrid(hindcastJJAS, aggr.mem = list(FUN = mean))
bias <- hindcast_aggr
obsJJAS_int <- interpGrid(grid = obsJJAS, new.coordinates = getGrid(hindcast_aggr))
bias$Data <- hindcast_aggr$Data - obsJJAS_int$Data

## ---- message = FALSE, fig.align='center'--------------------------------

bias.colors <- colorRampPalette(c("blue", "white", "red"))
plotClimatology(climatology(bias), 
                backdrop.theme = "coastline", 
                col.regions = bias.colors, 
                at = seq(-20,20,2),
                main = "Mean FWI bias of CFSv2")

## ---- message = FALSE, eval = FALSE--------------------------------------
#  hindcast_bc <- biasCorrection(y = obsJJAS, x = hindcastJJAS,
#                                newdata = hindcastJJAS, method = "eqm")

## ---- echo = FALSE-------------------------------------------------------
load("../ignore/vignette_data/hindcast_bc.rda")

## ---- message = FALSE----------------------------------------------------
fwimean_obs <- aggregateGrid(obsJJAS, aggr.y = list(FUN = "mean", na.rm = TRUE))


## ---- message = FALSE----------------------------------------------------
# raw data
fwimean_hind_raw <- aggregateGrid(hindcast, aggr.y = list(FUN = "mean", na.rm = TRUE))
# bias-corrected data
fwimean_hind <- aggregateGrid(hindcast_bc, aggr.y = list(FUN = "mean", na.rm = TRUE))

## ---- message = FALSE----------------------------------------------------
fwi90_obs <- aggregateGrid(obs, aggr.y = list(FUN = "quantile", probs = .9, na.rm = TRUE))

## ---- message = FALSE----------------------------------------------------
fot30 <- function(x) {
      sum(x > 30, na.rm = TRUE) / length(na.omit(x))
}

fwi30_obs <- aggregateGrid(obs, aggr.y = list(FUN = "fot30"))

## ------------------------------------------------------------------------
fwimean_obs_detrend <- detrendGrid(fwimean_obs)
fwimean_hind_detrend <- detrendGrid(fwimean_hind)

## ---- eval=FALSE---------------------------------------------------------
#  skill <- veriApply(verifun = "EnsRoca",
#                     fcst = fwimean_hind_detrend$Data,
#                     obs = fwimean_obs_detrend$Data,
#                     prob = 2/3,
#                     ensdim = 1,
#                     tdim = 2)

## ----echo=FALSE,message=FALSE--------------------------------------------
load("../ignore/vignette_data/skill.rda")

## ------------------------------------------------------------------------
upper.tercile <- easyVeri2grid(easyVeri.mat = skill$cat2,
                               obs.grid = fwimean_obs_detrend,
                               verifun = "EnsRoca")

## ---- message = FALSE, fig.align='center'--------------------------------
ms.colors <- colorRampPalette(c("blue", "grey90", "brown"))
plotClimatology(upper.tercile, 
                backdrop.theme = "countries",
                main = "ROC AREA (Above-normal FWI predictions)",
                col.regions = ms.colors)

## ---- message=TRUE-------------------------------------------------------
library(visualizeR)

## ---- warning = FALSE, fig.height=4.5, fig.width=4.5---------------------
fwimean_hind_sub <- subsetGrid(fwimean_hind_detrend, 
                               latLim = c(40,43.5), lonLim = c(24,26.5))
fwimean_obs_sub <- subsetGrid(fwimean_obs_detrend, 
                              latLim = c(40,43.5), lonLim = c(24,26.5))
plotClimatology(climatology(fwimean_obs_sub), 
                backdrop.theme = "countries", 
                col.regions = fwi.colors,
                scales = list(draw = TRUE, x = list(alternating = FALSE, tick.number = 6)),
                main = "Obs FWI climatology - Subset Area")

## ---- warning = FALSE----------------------------------------------------
tercilePlot(fwimean_hind_sub, obs = fwimean_obs_sub, color.pal = "ypb")

## ---- warning=FALSE, fig.asp=0.5-----------------------------------------
bubblePlot(hindcast = fwimean_hind, obs = fwimean_obs,
           year.target = 2001,
           bubble.size = 3, 
           size.as.probability = TRUE,
           score = FALSE)

## ------------------------------------------------------------------------
sessionInfo()

