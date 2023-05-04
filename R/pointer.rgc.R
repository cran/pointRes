#' Calculate pointer years using the relative growth change method
#'
#' @description The function calculates event and pointer years on a \code{data.frame} with tree-ring series using the relative growth change method, described as abrupt growth change method in Schweingruber et al. (1990). This method relates tree growth in year \code{\var{i}} to the average growth of \code{\var{n}} preceding years. Thresholds for event- and pointer-year calculations can be adjusted.
#' 
#' @usage pointer.rgc(data, period = NULL, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, 
#'             series.thresh = 75, make.plot = FALSE)
#'
#' @param data a \code{data.frame} with tree-ring series (raw or detrended) as columns and years as rows (e.g., output of \code{read.rwl} or \code{detrend} of package dplR).
#' @param period a \code{vector} specifying the start and end year of the analysis. Defaults to the full period covered by the data.
#' @param nb.yrs an \code{integer} specifying the number of preceding years to be used in calculating relative growth changes. Defaults to 4.
#' @param rgc.thresh.pos a \code{numeric} specifying the threshold above which a relative growth change (in percentage) for a specific tree and year is considered a positive event year. Defaults to 60.
#' @param rgc.thresh.neg a \code{numeric} specifying the threshold below which a relative growth change (in percentage) for a specific tree and year is considered a negative event year. Defaults to 40.
#' @param series.thresh a \code{numeric} specifying the minimum percentage of trees that should display a positive (or negative) event year for that year to be considered as positive (or negative) pointer year. Defaults to 75.
#' @param make.plot a \code{logical} specifying whether a bar plot, showing mean relative growth changes, should be created. Pointer years are indicated with dark-gray bars. Defaults to FALSE.
#' 
#' @details The function calculates the ratio of tree growth in year \code{\var{t}} and the average growth of \code{\var{n}} preceding years for individual trees. Resulting relative growth changes are used to identify event years for trees, and these event years to define pointer years for the site.
#' 
#' Following Schweingruber et al. (1990), \code{\var{nb.yrs}}, \code{\var{rgc.thresh.pos}}, \code{\var{rgc.thresh.neg}} and \code{\var{series.thresh}} are set to 4, 60, 40 and 75 respectively, meaning that a positive or negative pointer year will be defined when at least 75\% of the tree-ring series display an event year with a growth increase or decrease of at least 60\% or 40\%, respectively, relative to the average growth in the four preceding years. Depending on the sensitivity of the tree-ring series, the thresholds \code{\var{rgc.thresh.pos}} and \code{\var{rgc.thresh.neg}} may be adjusted (e.g., lowered for more complacent series; Jetschke et al. 2019).
#'
#' Note that the resulting time series are truncated by \code{\var{nb.yrs}} at the beginning inherent to the calculation method.
#'
#' @return 
#' The function returns a \code{list} containing the following components:
#' \item{rgc}{a \code{matrix} with relative growth changes for individual tree-ring series}
#' \item{EYvalues}{a \code{matrix} indicating positive (1), negative (-1) and non-event years (0) for individual tree-ring series}
#' \item{out}{a \code{data.frame} containing the following columns: \code{year} - time stamp, \code{nb.series} - number of series considered, \code{perc.pos} - percentage of trees showing a positive event year, \code{perc.neg} - percentage of trees showing a negative event year, \code{nature} - number indicating whether the year is a positive (1), negative (-1) or no pointer year (0), \code{dev_mean} - mean growth deviation in percentage over the available series for each year, and \code{dev_sd} - standard deviation of the growth deviation for each year}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @references Jetschke, G., van der Maaten, E. and van der Maaten-Theunissen, M. (2019) Towards the extremes: A critical analysis of pointer year detection methods. \emph{Dendrochronologia} 53: 55-62.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and Bräker, O.U. (1990) Identification, presentation and interpretation of event years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#' @references In writing the function, the code of the dplR function \code{pointer} (Pierre Mérian) was used as a reference.
#'
#' @examples ## Calculate pointer years on tree-ring series
#' data(s033)
#' py1 <- pointer.rgc(s033)
#' head(py1$out)
#' 
#' ## Calculate pointer years with user-defined arguments
#' data(s033)
#' py2 <- pointer.rgc(s033, period = c(1950,2010), nb.yrs = 5, rgc.thresh.pos = 50,
#'                    rgc.thresh.neg = 50, series.thresh = 50, make.plot = TRUE)
#' head(py2$out)
#' 
#' @import stats
#' @import ggplot2
#' @importFrom plyr round_any
#' 
#' @export pointer.rgc
#' 
pointer.rgc <- function(data, period = NULL, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, series.thresh = 75, make.plot = FALSE)
{
  stopifnot(is.numeric(nb.yrs), length(nb.yrs) == 1, is.finite(nb.yrs))
  if(nb.yrs < 1) {
    stop("'nb.yrs' must be > 0")
  }
  stopifnot(is.numeric(rgc.thresh.pos), is.numeric(rgc.thresh.neg),
            length(rgc.thresh.pos) == 1, length(rgc.thresh.neg) == 1, 
            is.finite(rgc.thresh.pos), is.finite(rgc.thresh.neg))
  if(rgc.thresh.pos < 0 | rgc.thresh.neg < 0) {
    stop("'rgc.thresh.pos' and (or) 'rgc.thresh.neg' must be > 0")
  }
  if(rgc.thresh.pos > 100 | rgc.thresh.neg > 100) {
    warning("'rgc.thresh.pos' and (or) 'rgc.thresh.neg' > 100 is unusual")
  }
  stopifnot(is.numeric(series.thresh), length(series.thresh) == 1, 
            is.finite(series.thresh))
  if(series.thresh < 0 || series.thresh > 100) {
    stop("'series.thresh' must range from 0 to 100")
  }
  if(is.null(period)) {
    start.yr <- min(as.numeric(rownames(data)))
    end.yr <- max(as.numeric(rownames(data)))
  } else{
    if(!is.null(period) && length(period) == 1) {
      stop("'period' needs a start and end year, e.g. c(1950, 2010), or should be NULL")
    }
    if(!is.null(period) && length(period) == 2 && is.numeric(period[1]) && is.numeric(period[2])) {
      start.yr <- period[1] 
      end.yr <- period[2]
    }else{
      stop("in 'period' the start and (or) end year are not numeric")
    }
  }
  if(start.yr > max(rownames(data))) {
    stop("the start year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  if(end.yr > max(rownames(data))) {
    stop("the end year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  data2 <- as.matrix(data[as.character(start.yr:end.yr),])
  if(!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  if(ncol(data2) == 1) {
    stop("'data' must contain more than one series")
  }
  rnames <- rownames(data2)
  if(is.null(rnames)) {
    stop("'data' must have explicit row names")
  }
  
  yrs <- as.numeric(rnames)
  nyrs <- length(yrs)
  if(nyrs < nb.yrs + 1) {
    stop("'data' must be longer than nb.yrs + 1")
  }
  
  if(nb.yrs == 1) {
    rgc <- data2[-1, , drop = FALSE] / data2[-nyrs, , drop = FALSE]
  }  
  else {
    avg.pre <- matrix(nrow = nrow(data2) - nb.yrs, ncol = ncol(data2))
    for(i in (nb.yrs+1):nyrs) {
      avg.pre[i - nb.yrs,] <- colMeans(data2[(i - nb.yrs):(i - 1),])
    }
    rownames(avg.pre) <- yrs[-nb.yrs:-1]
    rgc <- data2[-nb.yrs:-1, , drop = FALSE] / avg.pre[, , drop = FALSE]
  }
  
  pos.thresh <- rgc.thresh.pos/100 + 1
  neg.thresh <- 1 - rgc.thresh.neg/100
  EYvalues <- ifelse(rgc >= pos.thresh, 1, rgc)
  EYvalues <- ifelse(EYvalues <= neg.thresh, -1, EYvalues) 
  EYvalues <- ifelse(EYvalues == 1 | EYvalues == (-1), EYvalues, 0)
  
  year <- yrs[-nb.yrs:-1]
  nb.series <- rowSums(!is.na(rgc))
  perc.pos <- rowSums(rgc >= pos.thresh, na.rm = TRUE)/nb.series * 100
  perc.neg <- rowSums(rgc <= neg.thresh, na.rm = TRUE)/nb.series * 100
  nat.y.1 <- pmax(0, perc.pos - (series.thresh - 1e-07))
  nat.y.2 <- pmax(0, perc.neg - (series.thresh - 1e-07))
  nature <- sign(nat.y.1 - nat.y.2)
  dev_mean <- (rowMeans(rgc, na.rm = TRUE) - 1) * 100
  dev_sd <- apply(rgc, 1, function(x) sd(x, na.rm = TRUE)) * 100
  
  out <- data.frame(year, nb.series, perc.pos, perc.neg, nature, dev_mean, dev_sd, row.names = NULL)
  out[,c(3, 4, 6, 7)] <- round(out[,c(3, 4, 6, 7)], 2)
  
  rgc <- round(rgc, 2)
  
  spec.param <- data.frame(argument = c("nb.yrs", "rgc.thresh.pos", "rgc.thresh.neg", "series.thresh"), 
                           value = c(nb.yrs, rgc.thresh.pos, rgc.thresh.neg, series.thresh))
  output <- list(rgc = rgc, EYvalues = EYvalues, out = out, spec.param = spec.param)
  
  class(output) <- "pointer.rgc"
  
  if(make.plot == TRUE){
    start.yr2 <- round_any(start.yr, 10, f = floor)
    end.yr2 <- round_any(end.yr, 5, f = ceiling)
  
    data2 <- output$out[which(output$out[, "year"] == start.yr+nb.yrs):which(output$out[, "year"] == end.yr),]
    data3 <- as.data.frame(data2)
  
    year <- nature <- dev_mean <- dev_sd <- NULL
  
    nat.levels <- c(-1, 0, 1)
    fill.levels <- c("#636363", "#f0f0f0", "#636363")
  
    pl <- ggplot(data3, aes(x = year, y = dev_mean, fill = factor(nature))) +
      geom_bar(stat = "identity", position = "identity", colour = "black") +
      scale_fill_manual(limits = factor(nat.levels), values = fill.levels) +
      guides(fill = "none") +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, 10), 
                         minor_breaks = seq(start.yr2, end.yr2, 5),
                         limits = c(start.yr2-1, end.yr2+1)) +
      ylab("mean growth deviation (%)") + theme_bw()
    print(pl)
    return(output)
  } else{
    return(output)
    }
}