#' Calculate pointer years using the normalization in a moving window method
#'
#'
#' The function calculates pointer years on a \code{data.frame} of tree-ring series using the normalization in a moving window method introduced by Cropper (1979; cf. Schweingruber et al. 1990). This method normalizes tree growth in year \code{\var{i}} within a moving window of \code{\var{n}} years, thereby providing the number of standard deviations that tree growth deviates in individual years (Cropper values, C). To identify event years, thresholds on the number of standard deviations can be set according to Cropper (1979) and Neuwirth et al. (2007), optionally with different threshold values. The threshold for defining pointer years can be adjusted.
#'
#'
#' The function normalizes tree growth in year \code{\var{i}} within a moving window of \code{\var{n}} years. For method \code{"Cropper"}, event years are defined as those years having absolute Cropper values above a specified threshold (defaults to |C| > 0.75). For method \code{"Neuwirth"}, three classes of distinct growth deviations can be defined, being 'weak', 'strong' and 'extreme' (defaults to |C| > 1, |C| > 1.28, and |C| > 1.645). The window size can be adjusted, as well as the minimum percentage of trees that should display a positive (or negative) event year for that year to be considered as positive (or negative) pointer year.
#'
#' Note that the resulting time series are truncated by \code{\var{(window-1)/2}} at both ends inherent to the calculation methods. 
#'
#' @usage pointer.norm(data, window = 5, method.thresh = c("Cropper", "Neuwirth"),
#'              C.thresh = 0.75, N.thresh1 = 1, N.thresh2 = 1.28, 
#'              N.thresh3 = 1.645, series.thresh = 40) 
#'              
#' @param data a \code{data.frame} with tree-ring series as columns and years as rows (e.g., output of \code{read.rwl} of package dplR)
#' @param window an \code{integer} specifying the window size (i.e. number of years) to be used to calculate normalized growth deviations. Defaults to 5.
#' @param method.thresh a \code{character} string of \code{"Cropper"} or \code{"Neuwirth"}, specifying the method to be used to calculate event years. Argument matching is performed.
#' @param C.thresh a \code{numeric} specifying the threshold for identification of event years using method \code{"Cropper"}. Defaults to 0.75.
#' @param N.thresh1 a \code{numeric} specifying the threshold for identification of weak event years using method \code{"Neuwirth"}. Defaults to 1.
#' @param N.thresh2 a \code{numeric} specifying the threshold for identification of strong event years using method \code{"Neuwirth"}. Defaults to 1.28.
#' @param N.thresh3 a \code{numeric} specifying the threshold for identification of extreme event years using method \code{"Neuwirth"}. Defaults to 1.645.
#' @param series.thresh a \code{numeric} specifying the minimum percentage of trees that should display a positive (or negative) event year for that year to be considered as positive (or negative) pointer year. Defaults to 40.
#' @return
#'
#' The function returns a \code{list} containing the following components:
#'
#' \itemize{\item{for method \code{"Cropper"}:}}
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
#' \itemize{\item{for method \code{"Neuwirth"}:}}
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
#' @references Neuwirth, B., Schweingruber, F.H. and Winiger, M. (2007) Spatial patterns of central European pointer years from 1901 to 1971. \emph{Dendrochronologia} 24: 79-89.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and Bräker, O.U. (1990) Identification, presentation and interpretation of even years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#' 
#' @examples ## Calculate pointer years on tree-ring series using the method "Cropper"
#' ## and a user-defined threshold for event-year definition of 1:
#' data(s033)
#' py_c <- pointer.norm(s033, window = 5, method.thresh = "Cropper", 
#'                      C.thresh = 1, series.thresh = 40)
#' py_c$out
#'
#' ## Calculate pointer years on tree-ring series using the method "Neuwirth" 
#' ## and a moving window of 13 years:
#' data(s033)
#' py_n <- pointer.norm(s033, window = 13, method.thresh = "Neuwirth", 
#'                      series.thresh = 40)
#' py_n$out
#' @export
#' 
pointer.norm <- function(data, window = 5, method.thresh = c("Cropper", "Neuwirth"), C.thresh = 0.75,
                         N.thresh1 = 1, N.thresh2 = 1.28, N.thresh3 = 1.645, series.thresh = 40) 
{
  stopifnot(is.numeric(window), length(window) == 1, is.finite(window))
  if (window < 3) {
    stop("'window' must be > 3")
  }
  is.odd <- function(x) x %% 2 != 0
  if(is.odd(window) == FALSE) {
    stop(" a 'window' with an odd number of years is required")
  }
  stopifnot(is.numeric(series.thresh), length(series.thresh) == 
              1, is.finite(series.thresh))
  if (series.thresh < 0 || series.thresh > 100) {
    stop("'series.thresh' must range from 0 to 100")
  }
  data2 <- as.matrix(data)
  if (!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  if(ncol(data2) == 1){
    stop("'data' must contain more than one series")
  }
  rnames <- rownames(data2)
  if (is.null(rnames)) {
    stop("'data' must have explicit row names")
  }
  yrs <- as.numeric(rnames)
  nyrs <- length(yrs)
  if (nyrs < window) {
    stop("'data' must have more rows than the window length")
  }
  stopifnot(is.numeric(C.thresh), length(C.thresh) == 1, is.finite(C.thresh))
  if (method.thresh == "Cropper" && C.thresh <= 0) {
    stop("'C.thresh' must be larger than 0")
  }
  stopifnot(is.numeric(N.thresh1), length(N.thresh1) == 1, is.finite(N.thresh1))
  stopifnot(is.numeric(N.thresh2), length(N.thresh2) == 1, is.finite(N.thresh2))
  stopifnot(is.numeric(N.thresh3), length(N.thresh3) == 1, is.finite(N.thresh3))
  if (method.thresh == "Neuwirth" && N.thresh1 <= 0) {
    stop("'N.thresh1' must be larger than 0")
  }
  if (method.thresh == "Neuwirth" && (N.thresh2 <= N.thresh1 || N.thresh3 <= N.thresh2)){
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
  
  if (type.cropper) {
    EYvalues <- as.matrix(Cvalues)
    EYvalues[ EYvalues <= (-C.thresh)] <- -1
    EYvalues[ EYvalues >= C.thresh] <- 1
    EYvalues[ EYvalues < C.thresh  & EYvalues > (-C.thresh)] <- 0
    
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
  
  if (type.neuwirth) {
    EYvalues <- as.matrix(Cvalues)
    EYvalues[ EYvalues <= (-N.thresh3) ] <- -3
    EYvalues[ EYvalues <= (-N.thresh2) & EYvalues > (-N.thresh3)] <- -2
    EYvalues[ EYvalues <= (-N.thresh1) & EYvalues > (-N.thresh2)] <- -1
    EYvalues[ EYvalues >= (N.thresh3)] <- 3
    EYvalues[ EYvalues >= (N.thresh2) & EYvalues < (N.thresh3)] <- 2
    EYvalues[ EYvalues >= (N.thresh1) & EYvalues < (N.thresh2)] <- 1
    EYvalues[ EYvalues < (N.thresh1) & EYvalues > (-N.thresh1)] <- 0

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
  
    spec.param <- data.frame(argument = c("window", "N.thresh1", "N.thresh2", "N.thresh3", "series.thresh"), 
                             value = c(window, N.thresh1, N.thresh2, N.thresh3, series.thresh))
    
    output <- list(Cvalues = Cvalues, EYvalues = EYvalues, out = out, spec.param = spec.param)
    class(output) <- c("pointer.norm", "Neuwirth")
  }
  return(output)
}