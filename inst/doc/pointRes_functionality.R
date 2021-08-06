## ---- echo = FALSE, message = FALSE-------------------------------------------
library(pointRes)

## -----------------------------------------------------------------------------
library(pointRes)
data(s033) # the result of s033 <- read.rwl('s033.rwl') - a function of the dplR package

## -----------------------------------------------------------------------------
detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
pyc <- pointer.norm(detr_s033, period = NULL, window = 13, method.thresh = "Cropper", C.thresh = 0.75, series.thresh = 75, make.plot = FALSE)
pyn <- pointer.norm(detr_s033, period = NULL, window = 13, method.thresh = "Neuwirth", N.thresh = c(1, 1.28, 1.645), series.thresh = 75, make.plot = FALSE)

## -----------------------------------------------------------------------------
rgc <- pointer.rgc(s033, period = NULL, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, series.thresh = 75, make.plot = FALSE)

## -----------------------------------------------------------------------------
detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
pz <- pointer.zchron(detr_s033, period = NULL, bi.weight = TRUE, z.thresh = 1, t.Test = FALSE, make.plot = FALSE) 

## -----------------------------------------------------------------------------
it <- interval.trend(s033, period = NULL, trend.thresh = 0, IT.thresh = 95, make.plot = FALSE) 

## ---- fig.width = 7, fig.height = 3.5, fig.retina = 3-------------------------
detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
pyn <- pointer.norm(detr_s033, method = "Neuwirth", make.plot = TRUE)

## ---- fig.width = 7, fig.height = 5.3, fig.retina = 3-------------------------
event.plot(pyn, period = c(1950, 2007), x.tick.major = 10, x.tick.minor = 5)

## ---- fig.width = 7, fig.height = 2, fig.retina = 3---------------------------
pointer.plot(list(pyn,pyc,pz,it), sign = "neg", period = c(1950, 2007), labels = c("Neuwirth","Cropper","zChron","IT"))

## -----------------------------------------------------------------------------
res <- res.comp(s033, nb.yrs = c(4,4), max.yrs.rec = 10)

## ---- fig.width = 6, fig.height = 5.3, fig.retina = 3-------------------------
res.plot(res, select.yr = c(1976, 1992, 2003), param = "resist")

## -----------------------------------------------------------------------------
citation()
citation("pointRes")


