#' Calculate resilience indices
#'
#' @description The function calculates resilience indices on a \code{data.frame}, e.g., of tree-ring series, after Lloret et al. (2011; i.e. resistance, recovery, (relative) resilience), Thurm et al. (2016; recovery period, total growth reduction) and Schwarz et al. (2020; average growth reduction, average recovery rate), useful to analyze growth responses of individual trees prior, during and after extreme events / disturbances. The component 'resistance' is conceptually identical to 'abrupt growth changes' as described in Schweingruber et al. (1990; cf. \code{\link{pointer.rgc}}). 'Recovery' is the ability of tree growth to recover after disturbance, whereas 'resilience' reflects the ability of trees to reach pre-disturbance growth levels. Weighting of the resilience by the experienced growth reduction results in 'relative resilience'. 'Recovery period' (or: 'growth recovery time') is the time needed to reach pre-disturbance growth levels again. 'Total growth reduction' reflects the cumulative growth reduction in the year of disturbance as well as the associated years in the recovery period. 'Average growth reduction' is the total growth reduction divided by the length of the recovery period. 'Average recovery rate' is the mean percentual recovery over the recovery period.
#'
#' @usage res.comp(data, nb.yrs = c(4,4), max.yrs.rec = 10)
#'
#' @param data a \code{data.frame} with tree-ring series (raw or detrended) as columns and years as rows (e.g., output of \code{read.rwl}, \code{bai.in} or \code{bai.out} of package dplR)
#' @param nb.yrs a \code{vector} specifying the number of years for pre- and post-disturbance periods to be considered in calculating resilience components after Lloret et al. (2011). Defaults to 4 for both periods.
#' @param max.yrs.rec a \code{numeric} specifying the maximum length of the recovery period to be considered. Defaults to 10.
#' 
#' @details The function calculates various resilience indices after Lloret et al. (2011), Thurm et al. (2016) and Schwarz et al. (2020). The output provides different matrices with resilience indices for individual tree-ring series and years.
#' 
#' In calculating resilience indices the number of pre- and post disturbance years (\code{\var{nb.yrs}}), as well as the maximum number of years to be considered in defining the recovery period (\code{\var{max.yrs.rec}}) can be specified.
#' 
#' @return 
#' The function returns a \code{list} containing the following components:
#' \item{resist}{a \code{matrix} with resistance values (i.e. relative growth changes) for individual tree-ring series}
#' \item{recov}{a \code{matrix} with recovery values for individual tree-ring series}
#' \item{resil}{a \code{matrix} with resilience values for individual tree-ring series}
#' \item{rel.resil}{a \code{matrix} with relative resilience values for individual tree-ring series}
#' \item{rec.period}{a \code{matrix} with recovery periods for individual tree-ring series in years with decimal places (cf. Fig. 2 in Thurm et al. 2016). In case of no growth reduction (and thus no recovery), 0 is given as output. Inf indicates that no recovery occurred within the period as specified by \code{\var{max.yrs.rec}}}
#' \item{avg.rec.rate}{a \code{matrix} with average recovery rates for individual tree-ring series as percentage, calculated as 1 / \code{\var{rec.period}} * 100. In case of no growth reduction (and thus no recovery), NA is given as output. Inf indicates that no recovery occurred within the period as specified by \code{\var{max.yrs.rec}}}
#' \item{tot.abs.grow.red}{a \code{matrix} with total absolute growth reduction values for individual tree-ring series. 0 and Inf are given as output as for \code{\var{rec.period}}}
#' \item{tot.rel.grow.red}{a \code{matrix} with total growth reduction for individual tree-ring series expressed as percentage. 0 and Inf are given as output as for \code{\var{rec.period}}}
#' \item{avg.abs.grow.red}{a \code{matrix} with average absolute growth reduction, i.e. total absolute growth reduction divided by the number of full years needed for recovery. 0 and Inf are given as output as for \code{\var{rec.period}}}
#' \item{avg.rel.grow.red}{a \code{matrix} with average growth reduction expressed as percentage. 0 and Inf are given as output as for \code{\var{rec.period}}}
#' \item{nb.series}{a \code{data.frame} with the number of series for which the diverse indices could be calculated, with years in rows and indices in columns}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#' 
#' @author Marieke van der Maaten-Theunissen, Ernst van der Maaten and Mario Trouillier.
#'
#' @references Lloret, F., Keeling, E.G. and Sala, A. (2011) Components of tree resilience: effects of successive low-growth episodes in old ponderosa pine forests. \emph{Oikos} 120: 1909-1920.
#' @references Schwarz, J., Skiadaresis, G., Kohler, M., Kunz, J., Schnabel, F., Vitali, V. and Bauhus, J. (2020) Quantifying growth responses of trees to drought — a critique of commonly used resilience indices and recommendations for future studies. \emph{Current Forestry Reports} 6: 185-200.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and \enc{Bräker}{Braker}, O.U. (1990) Identification, presentation and interpretation of event years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#' @references Thurm, E.A., Uhl, E. and Pretzsch, H. (2016) Mixture reduces climate sensitivity of Douglas-fir stem growth. \emph{Forest Ecology and Management} 376: 205-220. 
#' 
#' @examples ## Calculate resilience indices on tree-ring series
#' data(s033)
#' res <- res.comp(s033)
#' 
#' @import stats
#' 
#' @export res.comp
#' 
res.comp <- function(data, nb.yrs = c(4,4), max.yrs.rec = 10)
{
  stopifnot(is.numeric(nb.yrs[1]), length(nb.yrs[1]) == 1, is.finite(nb.yrs[1]))
  stopifnot(is.numeric(nb.yrs[2]), length(nb.yrs[2]) == 1, is.finite(nb.yrs[2]))
  if(nb.yrs[1] < 1) {
    stop("'nb.yrs[1]' must be > 0")
  }
  if(nb.yrs[2] < 1) {
    stop("'nb.yrs[2]' must be > 0")
  }
  if(max.yrs.rec < 1) {
    stop("'max.yrs.rec' must be >= 1")
  }
  data2 <- as.matrix(data)
  if (!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  if(ncol(data2) == 1) {
    stop("'data' must contain more than one series")
  }
  data2 <- data2[rowSums(is.na(data2)) != ncol(data2), ]
  data3 <- rbind(data2, matrix(NA, nrow = max.yrs.rec, ncol = ncol(data2)) )
  last.yr <- as.numeric(rownames(data2))[nrow(data2)]
  rownames(data3) <- c(rownames(data2), (last.yr + 1) : (last.yr + max.yrs.rec))
  
  rnames <- rownames(data2)
  if (is.null(rnames)) {
    stop("'data' must have explicit row names")
  }
  yrs <- as.numeric(rnames)
  nyrs <- length(yrs)
  if (nyrs < nb.yrs[1] + nb.yrs[2] + 1) {
    stop("'data' must be longer than the full calculation window")
  }

  start <- nb.yrs[1] + 1
  
  avg.pre <- matrix(nrow = nrow(data2), ncol = ncol(data2), dimnames = list(rownames(data2), colnames(data2)) )
  if(nb.yrs[1] == 1) {
    avg.pre[start:nyrs, ] <- data3[(1:(nyrs - nb.yrs[1])), ]   
  } else {
    for(i in start:nyrs) {
      avg.pre[i, ] <- colMeans(data3[(i - nb.yrs[1]):(i - 1),])
    }
  }
  
  resist <- data2[,, drop = FALSE] / avg.pre[, , drop = FALSE]
  
  avg.post <- matrix(nrow = nrow(data2), ncol = ncol(data2), dimnames = list(rownames(data2), colnames(data2)) )
  if(nb.yrs[2] == 1){
    avg.post[1:(nyrs - nb.yrs[2]), ] <- data3[2:nyrs, ]  
  } else {
    for(i in 1:(nyrs - nb.yrs[2])) {
      avg.post[i,] <- colMeans(data2[(i + 1):(i + nb.yrs[2]),])
    }
  }
  
  recov <- avg.post[,, drop = FALSE] / data2[,, drop = FALSE]
  
  resil <- avg.post / avg.pre
  
  rel.resil <- (avg.post - data2) / avg.pre

  abs.grow.red <- array(data = NA, dim = c(dim(avg.pre), max.yrs.rec + 1),
                        dimnames = list(rownames(avg.pre), colnames(data), c(1:(max.yrs.rec + 1))) )
  rel.grow.red <- abs.grow.red
  row.i <- 1 : nyrs
  for(i in 1:(max.yrs.rec + 1)){
    row.i2 <- row.i + i - 1
    abs.grow.red[ , ,i] <- avg.pre[, , drop = FALSE] - data3[row.i2, , drop = FALSE]
    rel.grow.red[ , ,i] <- 1 - (data3[row.i2, , drop = FALSE] / avg.pre[, , drop = FALSE])
  }
  
  after.rec.set.na <- function(x){
    i <- which(x <= 0 | !is.finite(x))[1]
    if(!is.na(i) && i < length(x)){ x[(i + 1):length(x)] <- NA }
    return(x)
  }
  abs.grow.red <- apply(abs.grow.red, c(1,2), after.rec.set.na)
  abs.grow.red <- aperm(abs.grow.red, c(2,3,1))
  rel.grow.red <- apply(rel.grow.red, c(1,2), after.rec.set.na) 
  rel.grow.red <- aperm(rel.grow.red, c(2,3,1))

  measure.rec.time <- function(x){
    if( !all(is.na(x)) && !(any(is.na(x)) && !any(x[!is.na(x)] <= 0)) ){
      if( !any(x <= 0) ){
        return(Inf)
      }else{
        rec.i <- which(x <= 0)[1]
        if(rec.i == 1){
            return(0)
        } else {
            rec.val <- x[rec.i]
            pre.rec.val <- x[rec.i - 1]
            rec_period <- rec.i - 2 - pre.rec.val / (rec.val - pre.rec.val)
            return(rec_period)
        }
      }
    } else { return(NA) }
  }
  rec.period <- apply(rel.grow.red, c(1,2), measure.rec.time)

  avg.rec.rate <- 1 / rec.period * 100
  avg.rec.rate[avg.rec.rate == Inf] <- NA
  avg.rec.rate[avg.rec.rate == 0] <- Inf
  
  tot.grow.red.fun <- function(x){
    finite.i <- is.finite(x)
    if( !any(finite.i) ){ return(NA)
      } else {
      if( !any(x[finite.i] < 0) ){ return(Inf)
        } else {
        recovery.yr <- which(x < 0)[1]
        return( sum(x[0:(recovery.yr - 1)]) )
        }
      }
  }
  tot.abs.grow.red <- apply(abs.grow.red, c(1,2), tot.grow.red.fun)  
  tot.rel.grow.red <- apply(rel.grow.red, c(1,2), tot.grow.red.fun) * 100

  nb.yrs.below.avg.pre <- rec.period
  any.growth.red.i <- nb.yrs.below.avg.pre > 0 & is.finite(nb.yrs.below.avg.pre)
  nb.yrs.below.avg.pre[any.growth.red.i] <- floor(nb.yrs.below.avg.pre[any.growth.red.i] + 1)
  avg.abs.grow.red <- tot.abs.grow.red / nb.yrs.below.avg.pre
  avg.abs.grow.red[tot.abs.grow.red == Inf] <- Inf
  avg.abs.grow.red[is.nan(avg.abs.grow.red)] <- 0
  
  avg.rel.grow.red <- tot.rel.grow.red / nb.yrs.below.avg.pre
  avg.rel.grow.red[tot.rel.grow.red == Inf] <- Inf
  avg.rel.grow.red[is.nan(avg.rel.grow.red)] <- 0
  
  nb.series <- data.frame(resist = rowSums(is.finite(resist)))
  nb.series$recov <- rowSums(is.finite(recov))
  nb.series$resil <- rowSums(is.finite(resil))
  nb.series$rel.resil <- rowSums(is.finite(rel.resil))
  nb.series$rec.period <- rowSums(is.finite(rec.period))
  nb.series$avg.rec.rate <- rowSums(is.finite(avg.rec.rate))
  nb.series$tot.abs.grow.red <- rowSums(is.finite(tot.abs.grow.red))
  nb.series$tot.rel.grow.red <- rowSums(is.finite(tot.rel.grow.red))
  nb.series$avg.abs.grow.red <- rowSums(is.finite(avg.abs.grow.red))
  nb.series$avg.rel.grow.red <- rowSums(is.finite(avg.rel.grow.red))
  
  spec.param <- data.frame(argument = c("nb.yrs.pre","nb.yrs.post","max.yrs.rec"), 
                           value = c(nb.yrs[1], nb.yrs[2], max.yrs.rec))
  
  output <- list(resist = resist, recov = recov, resil = resil, rel.resil = rel.resil,
                 rec.period = rec.period, avg.rec.rate = avg.rec.rate,
                 tot.abs.grow.red = tot.abs.grow.red, tot.rel.grow.red = tot.rel.grow.red,
                 avg.abs.grow.red = avg.abs.grow.red, avg.rel.grow.red = avg.rel.grow.red,
                 nb.series = nb.series, spec.param = spec.param)
  
  class(output) <- "res.comp"
  return(output)
}
