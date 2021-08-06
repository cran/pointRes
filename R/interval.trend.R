#' Calculate pointer years using the interval trend method
#'
#' @description The function calculates year-to-year growth changes for individual tree-ring series and defines interval trends for (the population of) trees using the pointer interval method according to Schweingruber et al. (1990), which is also extensively described in Jetschke et al. (2019). The minimum percentual growth change and the minimum percentage of trees that should display a negative (or positive) trend for that year to be considered as negative (or positive) pointer year, can be adjusted.
#' 
#' @usage interval.trend(data, period = NULL, trend.thresh = 0, 
#'                IT.thresh = 95, make.plot = FALSE) 
#'              
#' @param data a \code{data.frame} with tree-ring series (raw or detrended) as columns and years as rows (e.g., output of \code{read.rwl} of package dplR).
#' @param period a \code{vector} specifying the start and end year of the analysis. Defaults to the full period covered by the data.
#' @param trend.thresh a \code{numeric} specifying the threshold for a percentual growth change to be considered a trend. Defaults to 0.
#' @param IT.thresh a \code{numeric} specifying the minimum percentage of trees that should display a negative (or positive) trend for that year to be considered as negative (or positive) pointer year. Defaults to 95.
#' @param make.plot a \code{logical} specifying whether a line plot, showing mean annual interval-trend values, should be created. The plot highlights positive and negative pointer years with black triangles pointing up and down, respectively. Defaults to FALSE.
#'  
#' @details The function calculates year-to-year growth changes. For each tree and year, the interval trend is defined as 1 if a positive change exceeds \code{trend.thresh}, as 0 if a negative change falls below minus \code{trend.thresh} and as 0.5 if the absolute change is below \code{trend.thresh}. \code{trend.thresh} defaults to 0\%. The interval trend for a population is defined as the average interval trend of the individual trees. A year is considered a negative (or positive) pointer year if the percentage of trees showing a decreasing (or increasing) trend exceeds \code{IT.thresh} (defaults to 95\%). Hence, in case of a negative pointer year the mean overall interval trend falls below 1 - \code{IT.thresh}, for a positive pointer year the mean overall interval trend exceeds \code{IT.thresh}.
#'
#' @return
#' The function returns a \code{list} containing the following components:
#' \item{perc.diff}{a \code{matrix} with percentual growth changes for individual tree-ring series}
#' \item{ITvalues}{a \code{matrix} indicating positive (1), negative (0) and no interval trends (0.5) for individual tree-ring series}
#' \item{out}{a \code{data.frame} containing the following columns:}
#' \item{}{\code{year} - time stamp}
#' \item{}{\code{nb.series} - number of series considered}
#' \item{}{\code{nature} - number indicating whether the year is a positive (1), negative (-1) or no pointer year (0)}
#' \item{}{\code{IT} - mean overall interval trend}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#'
#' @author Marieke van der Maaten-Theunissen, Ernst van der Maaten and Gottfried Jetschke.
#' 
#' @references Jetschke, G., van der Maaten, E. and van der Maaten-Theunissen, M. (2019) Towards the extremes: a critical analysis of pointer year detection methods. \emph{Dendrochronologia} 53: 55-62.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and Bräker, O.U. (1990) Identification, presentation and interpretation of event years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#' 
#' @examples ## Calculate pointer years using interval.trend
#' ## for a specified period and create a plot
#' data(s033)
#' IT <- interval.trend(s033, period = c(1950,2010), make.plot = TRUE)
#' 
#' ## Calculate pointer years as years with at least 90% of the trees
#' ## showing a positive/negative interval trend
#' data(s033)
#' IT <- interval.trend(s033, IT.thresh = 90)
#' IT$out[which(IT$out$nature == 1),"year"]
#' IT$out[which(IT$out$nature == -1),"year"]
#' 
#' @import ggplot2
#' @import stats
#' @importFrom plyr round_any
#' 
#' @export interval.trend
#'
interval.trend <- function(data, period = NULL, trend.thresh = 0, IT.thresh = 95, make.plot = FALSE) 
{
  stopifnot(is.numeric(IT.thresh), length(IT.thresh) == 
              1, is.finite(IT.thresh))
  if (IT.thresh < 0 || IT.thresh > 100) {
    stop("'IT.thresh' must range from 0 to 100")
  }
  stopifnot(is.numeric(trend.thresh), length(trend.thresh) == 
              1, is.finite(trend.thresh))
  if (trend.thresh < 0 || trend.thresh > 100) {
    stop("'trend.thresh' must range from 0 to 100")
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
  if(start.yr < min(rownames(data))) {
    stop("the start year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  if(end.yr > max(rownames(data))) {
    stop("the end year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  data2 <- as.matrix(data[as.character(start.yr:end.yr),])
  if (!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  if (ncol(data2) == 1) {
    stop("'data' must contain more than one series")
  }
  rnames <- rownames(data2)
  if (is.null(rnames)) {
    stop("'data' must have explicit row names")
  }

  perc.diff <- matrix(nrow = nrow(data2)-1, ncol = ncol(data2))
  for(i in 2:nrow(data2)) {
    perc.diff[i-1,] <- (data2[i,]/data2[i-1,] - 1) * 100
  }
  rownames(perc.diff) <- rownames(data2[2:nrow(data2),])
  colnames(perc.diff) <- colnames(data2)
  
  if(trend.thresh == 0){
    ITvalues <- as.matrix(perc.diff)
    ITvalues[perc.diff > 0] <- 1
    ITvalues[perc.diff = 0] <- 0.5
    ITvalues[perc.diff < 0] <- 0
  }else{
    ITvalues <- as.matrix(perc.diff)
    ITvalues[perc.diff >= trend.thresh] <- 1
    ITvalues[perc.diff > -trend.thresh & perc.diff < trend.thresh] <- 0.5
    ITvalues[perc.diff <= -trend.thresh] <- 0
  }
  rownames(ITvalues) <- rownames(perc.diff)
  
  year <- as.numeric(rnames[-1])
  nb.series <- rowSums(!is.na(ITvalues))
  IT <- rowSums(ITvalues, na.rm = TRUE)/nb.series
  
  nat.y.1 <- pmax(0, IT - (IT.thresh - 1e-07)/100)
  nat.y.2 <- pmax(0, (1 - IT) - (IT.thresh - 1e-07)/100)
  nature <- sign(nat.y.1 - nat.y.2)
  
  out <- data.frame(year, nb.series, nature, IT)
  rownames(out) <- NULL
  
  spec.param <- data.frame(argument = c("trend.thresh", "IT.thresh"), 
                           value = c(trend.thresh, IT.thresh))
  
  output <- list(perc.diff = perc.diff, ITvalues = ITvalues, out = out, spec.param = spec.param)
  class(output) <- c("interval.trend")
  
  if(make.plot == TRUE){
    start.yr2 <- round_any(start.yr, 10, f = floor)
    end.yr2 <- round_any(end.yr, 5, f = ceiling)
    
    data3 <- as.data.frame(out[which(out[, "year"] == start.yr+1):which(out[, "year"] == end.yr),])
    
    year <- nature <- NULL
    
    nat.levels <- c(-1, 0, 1)
    label.levels <- c("negative","none","positive")
    shape.levels <- c(25,21,24)
    fill.levels <- c("black", "gray", "black")
    
    pl <- ggplot(data3, aes(x = year, y = IT)) + geom_line() +
      geom_point(aes(fill = factor(nature), shape = factor(nature))) +
      scale_shape_manual(name = "pointer year class", limits = factor(nat.levels), labels = label.levels, values = shape.levels) + 
      scale_fill_manual(name = "pointer year class", limits = factor(nat.levels), labels = label.levels, values = fill.levels) +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, 10), 
                         minor_breaks = seq(start.yr2, end.yr2, 5),
                         limits = c(start.yr2-1, end.yr2+1)) +
      ylab("mean interval trend") + theme_bw()
    print(pl)
    return(output)
  } else{
    return(output)
  }
}
