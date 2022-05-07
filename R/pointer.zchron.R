#' Calculate pointer years using z-transformation of a site chronology
#' 
#' @description The function calculates pointer years on a \code{data.frame} of tree-ring series using a z-transformation of the site chronology (by default based on a biweight robust mean). The method provides the number of standard deviations that the chronology deviates in individual years. To identify pointer years, one absolute threshold on the number of standard deviations can be set. Optionally, a t-test can be applied to test whether the z-transformed chronology significantly exceeds the selected threshold. The function is intended to use on flexibly detrended data, e.g., with a cubic-smoothing spline with a 50\% frequency cut-off at 15 years (cf. Jetschke et al. 2019).
#'
#' @usage pointer.zchron(data, period = NULL, bi.weight = TRUE,
#'                z.thresh = 1, t.Test = FALSE, make.plot = FALSE) 
#'              
#' @param data a \code{data.frame} with detrended tree-ring series as columns and years as rows (e.g., output of \code{detrend} of package dplR).
#' @param period a \code{vector} specifying the start and end year of the analysis. Defaults to the full period covered by the data.
#' @param bi.weight a \code{logical} flag specifying whether a Tukey's biweight robust mean site chronology should be calculated. Defaults to TRUE.
#' @param z.thresh a \code{numeric} specifying the threshold for identification of pointer years. Defaults to 1.
#' @param t.Test a \code{logical} flag specifying whether a t-test should be performed. Defaults to FALSE.
#' @param make.plot a \code{logical} specifying whether a bar plot indicating pointer years should be created. Defaults to FALSE.
#' 
#' @details The function develops a site chronology, which is z-transformed over its entire length, thereby providing the number of standard deviations that the chronology deviates in individual years. In developing the site chronology, a normal or biweight robust mean can be used. A threshold \code{z.thresh} on the minimum number of standard deviations can be set (cf. Cropper 1979) to define the years to be considered as pointer years. 
#' Optionally, a t-test may be performed to test whether the z-transformed chronology significantly differs from the selected threshold value in a particular year. Therefore, individual tree-ring series are z-transformed as well and compared to the threshold value \code{z.thresh}. In case a biweight robust mean is used in building the site chronology, the t-test is based on the biweight robust estimate of the standard deviation. In all t-tests a significance level of 0.05 is used.
#'
#' @return 
#' The function returns a \code{list} containing the following components:
#' \item{TRIsite}{a \code{data.frame} with the calculated site chronology and corresponding sample depth (output of \code{chron} of package dplR)}
#' \item{out}{a \code{data.frame} containing the following columns: \code{year} - time stamp, \code{nb.series} - number of series considered, \code{nature} - number indicating whether the year is a positive (1), negative (-1) or no pointer year (0), and \code{AVGztrans} - z-transformed site chronology}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#'
#' @author Marieke van der Maaten-Theunissen, Ernst van der Maaten and Gottfried Jetschke.
#' 
#' @references Cropper, J.P. (1979) Tree-ring skeleton plotting by computer. \emph{Tree-Ring Bulletin} 39: 47-59.
#' @references Jetschke, G., van der Maaten, E. and van der Maaten-Theunissen, M. (2019) Towards the extremes: A critical analysis of pointer year detection methods. \emph{Dendrochronologia} 53: 55-62.
#' 
#' @examples ## Calculate pointer years on detrended tree-ring series
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 15)
#' pz1 <- pointer.zchron(detr_s033)
#' head(pz1$out)
#' 
#' ## Calculate pointer years with user-defined arguments
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 15)
#' pz2 <- pointer.zchron(detr_s033, period = c(1950,2010), z.thresh = 1.28,
#'                       make.plot = TRUE)
#' head(pz2$out)
#' 
#' @importFrom dplR chron
#' @importFrom DescTools TukeyBiweight
#' @import stats
#' @import ggplot2
#' @importFrom plyr round_any
#' 
#' @export pointer.zchron
#'
pointer.zchron <- function(data, period = NULL, bi.weight = TRUE, z.thresh = 1, t.Test = FALSE, make.plot = FALSE)
{
  stopifnot(is.numeric(z.thresh), length(z.thresh) == 1, is.finite(z.thresh))
  if(z.thresh <= 0) {
    stop("'z.thresh' must be larger than 0")
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
  if(ncol(data2) == 1) {
    stop("'data' must contain more than one series")
  }
  rnames <- rownames(data2)
  if(is.null(rnames)) {
    stop("'data' must have explicit row names")
  }  

  TRIsite <- chron(data2, prefix = "AVG", biweight = bi.weight, prewhiten = FALSE)
  
  mean.chron <- mean(TRIsite[,1], na.rm = TRUE)
  sd.chron <- sd(TRIsite[,1], na.rm = TRUE)
  
  ztrans <- as.matrix((TRIsite[,1, drop = FALSE] - mean.chron) / sd.chron)
  colnames(ztrans) <- "AVGztrans"
  ztransTrees <- as.matrix((data2 - mean.chron) / sd.chron)
  
  year <- as.numeric(rnames)
  nb.series <- rowSums(!is.na(data2))
  
  AVGztrans <- CI <- NULL
  
  if(t.Test == TRUE) {
    if(bi.weight == TRUE){
      CIz <- data.frame(matrix(ncol = 3, nrow = nrow(ztransTrees), dimnames = list(NULL, c("tbw", "lower", "upper"))))
      seq.yr <- which(nb.series > 1)
      for(i in seq.yr){
        CIz[i,] <- TukeyBiweight(ztransTrees[i,], na.rm = TRUE, conf.level = 0.95,
                                 ci.type = "basic", R = 10000)
      }
      if(1 %in% nb.series) {
        message("For years with only one series, no confidence interval could be calculated. Consequently, these years are excluded as possible pointer years.")
      }
      nature <- as.matrix(ztrans)
      nature[CIz$upper <= -z.thresh] <- -1
      nature[CIz$lower >= z.thresh] <- 1
      nature[rowSums(is.na(CIz)) == 2] <- NA
      nature <- ifelse(nature == 1 | nature == -1 | is.na(nature), nature, 0)
      colnames(nature) <- "nature"
      
      out <- data.frame(year, nb.series, nature, ztrans, row.names = NULL)  
      out[, 3:4] <- round(out[, 3:4], 2)
    } else {
      sd_ztransTrees <- as.data.frame(apply(ztransTrees, 1, sd, na.rm = TRUE))
      
      CIz <- data.frame(matrix(ncol = 2, nrow = nrow(ztransTrees), dimnames = list(NULL, c("upper", "lower"))))
      for(i in 1:length(ztrans)){
        CIz[i,"upper"] <- ztrans[i,] + qt(.975, df = nb.series[i]-1) * sd_ztransTrees[i,]/sqrt(nb.series[i])
        CIz[i,"lower"] <- ztrans[i,] - qt(.975, df = nb.series[i]-1) * sd_ztransTrees[i,]/sqrt(nb.series[i])
      }
      nature <- as.matrix(ztrans)
      nature[CI$upper <= -z.thresh] <- -1
      nature[CI$lower >= z.thresh] <- 1
      nature <- ifelse(nature == 1 | nature == -1, nature, 0)
      colnames(nature) <- "nature"
      
      out <- data.frame(year, nb.series, nature, ztrans, row.names = NULL)  
      out[, 3:4] <- round(out[, 3:4], 2)
    }
  } else {
    nature <- as.matrix(ztrans)
    nature[nature < z.thresh & nature > -z.thresh] <- 0
    nature[nature <= -z.thresh] <- -1
    nature[nature >= z.thresh] <- 1
    nature <- ifelse(nature == 1 | nature == -1, nature, 0)
    colnames(nature) <- "nature" 
    
    out <- data.frame(year, nb.series, nature, ztrans, row.names = NULL)  
    out[, 3:4] <- round(out[, 3:4], 2)
  }
  
  spec.param <- data.frame(argument = c("bi.weight","z.thresh", "t.Test"), value = c(bi.weight, z.thresh, t.Test))
  
  output <- list(TRIsite = TRIsite, out = out, spec.param = spec.param)
  class(output) <- c("pointer.zchron")
  
  if(make.plot == TRUE){
    start.yr2 <- round_any(start.yr, 10, f = floor)
    end.yr2 <- round_any(end.yr, 5, f = ceiling)
    
    data2 <- output$out[which(output$out[, "year"] == start.yr):which(output$out[, "year"] == end.yr),]
    data3 <- as.data.frame(data2)
    
    year <- nature <- NULL
    
    nat.levels <- c(-1, 0, 1)
    fill.levels <- c("#636363", "#f0f0f0", "#636363")
    
    pl <- ggplot(data3, aes(x = year, y = AVGztrans, fill = factor(nature))) +
      geom_bar(stat = "identity", position = "identity", colour = "black") +
      scale_fill_manual(limits = factor(nat.levels), values = fill.levels) +
      guides(fill = FALSE) +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, 10), 
                         minor_breaks = seq(start.yr2, end.yr2, 5),
                         limits = c(start.yr2-1, end.yr2+1)) +
      ylab("z-scores") + theme_bw()
    print(pl)
    return(output)
  } else{
    return(output)
  }
}


