#' Calculate pointer years using the normalization in a moving window method
#'
#' @description The function calculates event and pointer years on a \code{data.frame} with tree-ring series using the normalization in a moving window method introduced by Cropper (1979; cf. Schweingruber et al. 1990, Jetschke et al. 2019). This method locally z-transforms tree growth in year \code{\var{i}} within a symmetric moving window of \code{\var{n}} years, thereby providing the number of standard deviations that tree growth deviates in individual years (Cropper values, Ct) from the window average. To identify event years, one absolute threshold on the number of standard deviations can be set (Cropper 1979), or, alternatively, three intensity classes (Neuwirth et al. 2007). Threshold values for defining event and pointer years can be adjusted.
#' 
#' Prior to the calculation of event and pointer years with \code{pointer.norm}, a 13-year weighted low-pass filter, as described by Fritts (1976), may be applied on the tree-ring series using \code{\link{lowpass13}}. According to Cropper (1979), such a filter improves the detection of event and pointer years for complacent series, whereas for sensitive series filtering has little effect. Alternatively, a flexible detrending (e.g., using a cubic-smoothing spline with a 50\% frequency cut-off at 15 years) may be applied prior to \code{pointer.norm}.
#'
#' @usage pointer.norm(data, period = NULL, window = 13, method.thresh = c("Cropper", "Neuwirth"),
#'              C.thresh = 0.75, N.thresh = c(1, 1.28, 1.645), series.thresh = 75, 
#'              make.plot = FALSE) 
#'              
#' @param data a \code{data.frame} with detrended tree-ring series as columns and years as rows (e.g., output of \code{detrend} of package dplR), or a \code{data.frame} with filtered series (output of \code{\link{lowpass13}}).
#' @param period a \code{vector} specifying the start and end year of the analysis. Defaults to the full period covered by the data.
#' @param window an \code{integer} specifying the window size (i.e. number of years) to be used to calculate normalized growth deviations. Must be an odd number (>=3). Defaults to 13.
#' @param method.thresh a \code{character} string of \code{"Cropper"} or \code{"Neuwirth"}, specifying whether one absolute threshold or three intensity classes should be used for defining event years. Argument matching is performed.
#' @param C.thresh a (positive) \code{numeric} specifying the threshold for identification of event years using method \code{"Cropper"}. Defaults to 0.75.
#' @param N.thresh a \code{vector} specifying the (positive) thresholds for identification of weak, strong and extreme event years using method \code{"Neuwirth"}. Defaults to 1, 1.28 and 1.645, respectively.
#' @param series.thresh a \code{numeric} specifying the minimum percentage of trees that should display a positive (or negative) event year for that year to be considered as positive (or negative) pointer year. Defaults to 75.
#' @param make.plot a \code{logical} specifying whether a bar plot, showing mean Cropper values, should be created. Pointer years are indicated with dark-gray bars. If event years were defined using \code{method.thresh "Neuwirth"} (\code{\link{pointer.norm}}), different tones of gray indicate weak, strong and extreme pointer years, based on the most common event year class. Defaults to FALSE.
#' 
#' @details The function locally z-transforms tree growth in year \code{\var{i}} within a symmetric moving window of \code{\var{n}} years. For \code{method.thresh} \code{"Cropper"}, event years are defined as those years having absolute Cropper values (Ct) above a specified threshold (defaults to |Ct| > 0.75). For \code{method.thresh} \code{"Neuwirth"}, three classes of distinct growth deviations can be defined, being 'weak', 'strong' and 'extreme' (defaults to |Ct| > 1, |Ct| > 1.28, and |Ct| > 1.645). The window size can be adjusted (defaults to 13 years), as well as the minimum percentage of trees that should display a positive (or negative) event year for that year to be considered as positive (or negative) pointer year (defaults to 75\%).
#'
#' Note that the resulting time series are truncated by \code{\var{(window-1)/2}} at both ends inherent to the calculation method. 
#'
#' @return
#' The function returns a \code{list} containing the following components:
#'
#' \itemize{\item{for \code{method.thresh} \code{"Cropper"}:}}
#' \item{Cvalues}{a \code{matrix} with Cropper values for individual tree-ring series}
#' \item{EYvalues}{a \code{matrix} indicating positive (1), negative (-1) and non-event years (0) for individual tree-ring series}
#' \item{out}{a \code{data.frame} containing the following columns:}
#' \item{}{\code{year} - time stamp}
#' \item{}{\code{nb.series} - number of series considered}
#' \item{}{\code{perc.pos} - percentage of trees showing a positive event year}
#' \item{}{\code{perc.neg} - percentage of trees showing a negative event year}
#' \item{}{\code{nature} - number indicating whether the year is a positive (1), negative (-1) or no pointer year (0)}
#' \item{}{\code{Cvalues_mean} - mean Cropper value over the available series}
#' \item{}{\code{Cvalues_sd} - standard deviation of Cropper values}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#'
#' \itemize{\item{for \code{method.thresh} \code{"Neuwirth"}:}}
#' \item{Cvalues}{a \code{matrix} with Cropper values for individual tree-ring series}
#' \item{EYvalues}{a \code{matrix} indicating weak (1/-1), strong (2/-2) and extreme (3/-3) positive/negative event years, as well as non-event years (0) for individual tree-ring series}
#' \item{out}{a \code{data.frame} containing the following columns:}
#' \item{}{\code{year} - time stamp}
#' \item{}{\code{nb.series} - number of series considered}
#' \item{}{\code{perc.pos.extreme} - percentage of trees showing a positive extreme event year}
#' \item{}{\code{perc.pos.strong} - percentage of trees showing a positive strong event year}
#' \item{}{\code{perc.pos.weak} - percentage of trees showing a positive weak event year}
#' \item{}{\code{perc.neg.weak} - percentage of trees showing a negative weak event year}
#' \item{}{\code{perc.neg.strong} - percentage of trees showing a negative strong event year}
#' \item{}{\code{perc.neg.extreme} - percentage of trees showing a negative extreme event year}
#' \item{}{\code{nature} - number indicating whether the year is a positive (1), negative (-1) or no pointer year (0)}
#' \item{}{\code{Cvalues_mean} - mean Cropper value over the available series}
#' \item{}{\code{Cvalues_sd} - standard deviation of Cropper values}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#'  
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#' 
#' @references Cropper, J.P. (1979) Tree-ring skeleton plotting by computer. \emph{Tree-Ring Bulletin} 39: 47-59.
#' @references Jetschke, G., van der Maaten, E. and van der Maaten-Theunissen, M. (2019) Towards the extremes: A critical analysis of pointer year detection methods. \emph{Dendrochronologia} 53: 55-62.
#' @references Fritts, H.C. (1976) Tree rings and climate. Academic Press Inc. (London) Ltd.
#' @references Neuwirth, B., Schweingruber, F.H. and Winiger, M. (2007) Spatial patterns of central European pointer years from 1901 to 1971. \emph{Dendrochronologia} 24: 79-89.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and Bräker, O.U. (1990) Identification, presentation and interpretation of event years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#' 
#' @examples ## Calculate pointer years on detrended tree-ring series using method.thresh "Cropper"
#' ## and a user-defined threshold for event-year definition of 1
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
#' pyc <- pointer.norm(detr_s033, method.thresh = "Cropper", C.thresh = 1)
#' head(pyc$out)
#'
#' ## Calculate pointer years on detrended tree-ring series using method.thresh "Neuwirth"
#' ## and a user-defined window of 11 years
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
#' pyn <- pointer.norm(detr_s033, window = 11, method.thresh = "Neuwirth", make.plot = TRUE)
#' head(pyn$out)
#' 
#' @importFrom dplR detrend
#' @import stats
#' @import ggplot2
#' @importFrom plyr round_any
#' 
#' @export pointer.norm
#' 
pointer.norm <- function(data, period = NULL, window = 13, method.thresh = c("Cropper", "Neuwirth"), C.thresh = 0.75, N.thresh = c(1, 1.28, 1.645), series.thresh = 75, make.plot = FALSE) 
{
  stopifnot(is.numeric(window), length(window) == 1, is.finite(window))
  if(window < 3) {
    stop("'window' must be >= 3")
  }
  is.odd <- function(x) x %% 2 != 0
  if(is.odd(window) == FALSE) {
    stop(" a 'window' with an odd number of years is required")
  }
  stopifnot(is.numeric(series.thresh), length(series.thresh) == 
              1, is.finite(series.thresh))
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
  if(start.yr < min(rownames(data))) {
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
  if(nyrs < window) {
    stop("'data' must have more rows than the window length")
  }
  stopifnot(is.numeric(C.thresh), length(C.thresh) == 1, is.finite(C.thresh))
  if(method.thresh == "Cropper" && C.thresh <= 0) {
    stop("'C.thresh' must be larger than 0")
  }
  stopifnot(is.numeric(N.thresh[1]), length(N.thresh[1]) == 1, is.finite(N.thresh[1]))
  stopifnot(is.numeric(N.thresh[2]), length(N.thresh[2]) == 1, is.finite(N.thresh[2]))
  stopifnot(is.numeric(N.thresh[3]), length(N.thresh[3]) == 1, is.finite(N.thresh[3]))
  if(method.thresh == "Neuwirth" && N.thresh[1] <= 0) {
    stop("'N.thresh[1]' must be larger than 0")
  }
  if(method.thresh == "Neuwirth" && (N.thresh[2] <= N.thresh[1] || N.thresh[3] <= N.thresh[2])) {
    stop("the 'N.thresh' thresholds must have increasing values ")
  }
  
  method.thresh2 <- match.arg(method.thresh, c("Cropper", "Neuwirth"))
  
  tail <- (window - 1) / 2
  start <- tail + 1
  
  Cvalues <- matrix(nrow = nrow(data2) - 2*tail, ncol = ncol(data2))
  
  for(i in start:(nrow(data2) - tail)) {
    Cvalues[i - tail,] <- (data2[i,] - colMeans(data2[(i - tail) : (i + tail),])) / 
      apply(data2[(i - tail) : (i + tail),], 2, sd)
  }
  rownames(Cvalues) <- yrs[start : (length(yrs) - tail)]
  colnames(Cvalues) <- colnames(data2)
  
  type.cropper <- method.thresh2 == "Cropper"
  type.neuwirth <- method.thresh2 == "Neuwirth"
  
  if(type.cropper) {
    EYvalues <- as.matrix(Cvalues)
    EYvalues[ EYvalues < C.thresh  & EYvalues > (-C.thresh)] <- 0
    EYvalues[ EYvalues <= (-C.thresh)] <- -1
    EYvalues[ EYvalues >= C.thresh] <- 1
    
            year <- yrs[start : (length(yrs) - tail)]
       nb.series <- rowSums(!is.na(Cvalues))
        perc.pos <- rowSums(EYvalues == 1, na.rm = TRUE)/nb.series * 100
        perc.neg <- rowSums(EYvalues == -1, na.rm = TRUE)/nb.series * 100            
         nat.y.1 <- pmax(0, perc.pos - (series.thresh - 1e-07))
         nat.y.2 <- pmax(0, perc.neg - (series.thresh - 1e-07))
          nature <- sign(nat.y.1 - nat.y.2)
    Cvalues_mean <- rowMeans(Cvalues, na.rm = TRUE)
      Cvalues_sd <- apply(Cvalues, 1, function(x) sd(x, na.rm = TRUE))
    
    out <- data.frame(year, nb.series, perc.pos, perc.neg, nature, Cvalues_mean, Cvalues_sd, row.names = NULL)
    out[,c(3, 4, 6, 7)] <- round(out[,c(3, 4, 6, 7)], 2)

    Cvalues <- round(Cvalues, 2)
    
    spec.param <- data.frame(argument = c("window", "C.thresh", "series.thresh"), 
                             value = c(window, C.thresh, series.thresh))
    
    output <- list(Cvalues = Cvalues, EYvalues = EYvalues, out = out, spec.param = spec.param)
    class(output) <- c("pointer.norm", "Cropper")
  }
  
  if(type.neuwirth) {
    EYvalues <- as.matrix(Cvalues)
    EYvalues[ EYvalues <= (-N.thresh[3]) ] <- -3
    EYvalues[ EYvalues <= (-N.thresh[2]) & EYvalues > (-N.thresh[3])] <- -2
    EYvalues[ EYvalues <= (-N.thresh[1]) & EYvalues > (-N.thresh[2])] <- -1
    EYvalues[ EYvalues >= (N.thresh[3])] <- 3
    EYvalues[ EYvalues >= (N.thresh[2]) & EYvalues < (N.thresh[3])] <- 2
    EYvalues[ EYvalues >= (N.thresh[1]) & EYvalues < (N.thresh[2])] <- 1
    EYvalues[ EYvalues < (N.thresh[1]) & EYvalues > (-N.thresh[1])] <- 0

                year <- yrs[start : (length(yrs) - tail)] 
           nb.series <- rowSums(!is.na(Cvalues))
    perc.pos.extreme <- rowSums(EYvalues == 3, na.rm=TRUE)/nb.series * 100
     perc.pos.strong <- rowSums(EYvalues == 2, na.rm = TRUE)/nb.series * 100
       perc.pos.weak <- rowSums(EYvalues == 1, na.rm = TRUE)/nb.series * 100
       perc.neg.weak <- rowSums(EYvalues == -1, na.rm=TRUE)/nb.series * 100
     perc.neg.strong <- rowSums(EYvalues == -2, na.rm = TRUE)/nb.series * 100
    perc.neg.extreme <- rowSums(EYvalues == -3, na.rm = TRUE)/nb.series * 100
                 pos <- cbind(perc.pos.extreme, perc.pos.strong, perc.pos.weak)
                 neg <- cbind(perc.neg.weak, perc.neg.strong, perc.neg.extreme)
             nat.y.1 <- pmax(0, rowSums(pos) - (series.thresh - 1e-07))
             nat.y.2 <- pmax(0, rowSums(neg) - (series.thresh - 1e-07))
              nature <- sign(nat.y.1 - nat.y.2)
        Cvalues_mean <- rowMeans(Cvalues, na.rm = TRUE)
          Cvalues_sd <- apply(Cvalues, 1, function(x) sd(x, na.rm = TRUE))
    
    out <- data.frame(year, nb.series, perc.pos.extreme, perc.pos.strong, perc.pos.weak,
                      perc.neg.weak, perc.neg.strong, perc.neg.extreme, nature, Cvalues_mean,
                      Cvalues_sd, row.names = NULL)  
    out[, c(3:8, 10, 11)] <- round(out[, c(3:8, 10, 11)], 2)
    
    Cvalues <- round(Cvalues, 2)
  
    spec.param <- data.frame(argument = c("window", "N.thresh[1]", "N.thresh[2]", "N.thresh[3]", "series.thresh"), 
                             value = c(window, N.thresh[1], N.thresh[2], N.thresh[3], series.thresh))
    
    output <- list(Cvalues = Cvalues, EYvalues = EYvalues, out = out, spec.param = spec.param)
    
    class(output) <- c("pointer.norm", "Neuwirth")
  }
  
  if(make.plot == TRUE){
    start.yr2 <- round_any(start.yr, 10, f = floor)
    end.yr2 <- round_any(end.yr, 5, f = ceiling)
    
    data2 <- output$out[which(output$out[, "year"] == start.yr+tail):which(output$out[, "year"] == end.yr-tail),]
    data3 <- as.data.frame(data2)
    
    year <- nature <- Cvalues_mean <- Cvalues_sd <- int.class <- NULL
    limits <- aes(ymax = Cvalues_mean + Cvalues_sd, ymin = Cvalues_mean - Cvalues_sd)
    
    if(colnames(data3)[3] == "perc.pos") {
      nat.levels <- c(-1, 0, 1)
      fill.levels <- c("#636363", "#f0f0f0", "#636363")
      
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(nature))) +
        geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = factor(nat.levels), values = fill.levels) +
        guides(fill = "none") +
        scale_x_continuous(breaks = seq(start.yr2, end.yr2, 10),
                           minor_breaks = seq(start.yr2, end.yr2, 5),
                           limits = c(start.yr2-1, end.yr2+1)) +
        ylab("mean Cropper value") + theme_bw()
    }
    else {
      data3[,12] <- ifelse(data2[, "nature"] == (-1), 
                           max.col(data2[,c(1, 2, 5, 4, 3, 6:11)][, 6:8], ties.method = "first"), 
                           ifelse(data2[, "nature"] == 1,
                                  max.col(data2[,c(1, 2, 5, 4, 3, 6:11)][, 3:5], ties.method = "first"), 0))
      data3[,12] <- ifelse(data2[, "nature"] == (-1), paste("-", data3[, 12], sep = ''), data3[, 12])
      colnames(data3)[12] <- "int.class"
      
      int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
      fill.levels <- c("black", "#636363","#bdbdbd", "#f0f0f0", "#bdbdbd", "#636363", "black")
      
      pl <- ggplot(data3, aes(x = year, y = Cvalues_mean, fill = factor(int.class))) +
        geom_bar(stat = "identity", position = "identity", colour = "black") +
        scale_fill_manual(limits = factor(int.levels), values = fill.levels) +
        guides(fill = "none") +
        scale_x_continuous(breaks = seq(start.yr2, end.yr2, 10), 
                           minor_breaks = seq(start.yr2, end.yr2, 5),
                           limits = c(start.yr2-1, end.yr2+1)) +
        ylab("mean Cropper value") +  theme_bw()
    }
    print(pl)
    return(output)
  } else{
    return(output)}
}
  