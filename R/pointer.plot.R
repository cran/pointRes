#' Plot pointer years for multiple sites
#'
#' @description The function creates a dot plot showing positive and (or) negative pointer years from \code{lists} of the type as produced by \code{\link{pointer.norm}}, \code{\link{pointer.rgc}}, \code{\link{pointer.zchron}} and (or) \code{\link{interval.trend}}.
#' 
#' @usage pointer.plot(list.sites, sign = c("both", "pos", "neg"),
#'              period = NULL, labels = NULL,
#'              x.tick.major = 10, x.tick.minor = 5) 
#'
#' @param list.sites a \code{list} with \code{lists} as produced by \code{\link{pointer.norm}}, \code{\link{pointer.rgc}}, \code{\link{pointer.zchron}} or \code{\link{interval.trend}} for individual sites (created using list(site1, site2,..)).
#' @param sign a \code{character} string specifying whether both positive and negative (\code{"both"}), or only positive (\code{"pos"}) or negative (\code{"neg"}) pointer years should be displayed. Defaults to \code{"both"}.
#' @param period a \code{vector} specifying the start and end year to be plotted. Defaults to the full period covered by the output of the pointer year analyses.
#' @param labels a \code{character} vector with labels for the sites. Defaults to 'site 1, 2, .., \code{\var{i}}'.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a dot plot showing pointer years for multiple sites. Positive and negative pointer years are indicated with different symbols an (or) colors.
#' 
#' @return 
#' Dot plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot negative pointer years for multiple sites (or different methods)
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
#' py <- pointer.rgc(s033)
#' pyn <- pointer.norm(detr_s033, method = "Neuwirth")
#' it <- interval.trend(s033)
#' zchron <- pointer.zchron(detr_s033)
#' comparison <- list(py, pyn, it, zchron)
#' pointer.plot(comparison, sign = "neg", period = c(1950, 2013),
#'              labels = c("py", "pyn", "it", "zchron")) 
#'
#' ## Plot pointer years for different specifications of pointer.norm (method "Neuwirth")
#' data(s033)
#' w09 <- pointer.norm(detr_s033, window = 9, method.thresh = "Neuwirth")
#' w11 <- pointer.norm(detr_s033, window = 11, method.thresh = "Neuwirth")
#' w13 <- pointer.norm(detr_s033, method.thresh = "Neuwirth")
#' comparison <- list(w09, w11, w13)
#' pointer.plot(comparison, period = c(1950, 2007)) 
#'            
#' @import ggplot2
#' @import stats
#' @importFrom plyr round_any
#' @importFrom TripleR matrix2long
#'       
#' @export pointer.plot
#' 
pointer.plot <- function(list.sites, sign = c("both", "pos", "neg"), period = NULL, labels = NULL, x.tick.major = 10, x.tick.minor = 5) 
{
  stopifnot(is.list(list.sites))
  for(i in 1:length(list.sites)) {
    if(FALSE %in% (class(list.sites[[i]])[1] == "pointer.norm") & FALSE %in% (class(list.sites[[i]])[1] == "pointer.rgc") & FALSE %in% (class(list.sites[[i]])[1] == "pointer.zchron") & FALSE %in% (class(list.sites[[i]])[1] == "interval.trend")){
      stop("'list.sites' contains list(s) that are no output of functions pointer.norm, pointer.rgc, pointer.zchron or interval.trend")
    }
  }
  check2 <- vector()
  for(i in 1:length(list.sites)) {
    check2[i] <- class(list.sites[[i]])[1] == "pointer.norm"
  }
  if((TRUE %in% check2) && (FALSE %in% check2)){
    warning("'list.sites' contains output of multiple functions to calculate pointer years. When performing a method comparison, please ignore this warning")
  }
  if((TRUE %in% check2)) {
    check3 <- vector()
    for(i in 1:length(list.sites)){
      check3[i] <- class(list.sites[[i]])[2] == "Neuwirth"
    }
    check3[is.na(check3)] <- FALSE
    if((TRUE %in% check3) && (FALSE %in% check3)) {
      check3 <- check3
    }
  } else {
      check3 <- vector()
      check3[1] <- FALSE
  }
  if (x.tick.minor > x.tick.major){
    stop("'x.tick.minor' should be smaller then 'x.tick.major'")
  }
  
  sign2 <- match.arg(sign, c("both", "pos", "neg"))
  
  year <- site <- PYvalues <- NULL
  
  if(length(labels) == 0) {
    labels2 <- vector()
    for(i in 1:length(list.sites)) {
      labels2[i] <- paste("site", i, sep = " ")
    }
  } else {
      labels2 <- labels
  }
  
  vec.min <- vec.max <- vector()
  for(i in 1:length(list.sites)) {
    vec.min[i] <- min(list.sites[[i]]$out[, "year"])
    vec.max[i] <- max(list.sites[[i]]$out[, "year"])
  }
  yrs <- seq(min(vec.min), max(vec.max), 1)
  
  if(is.null(period)) {
    start.yr <- min(vec.min)
    end.yr <- max(vec.max)
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
  if(start.yr < min(vec.min)) {
    stop("the start year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  if(end.yr > max(vec.max)) {
    stop("the end year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  
  nature <- matrix(nrow = length(yrs), ncol = length(list.sites))
  rownames(nature) <- yrs
  for(i in 1:length(list.sites)) {
    nature[as.character(list.sites[[i]]$out[, "year"]),i] <- list.sites[[i]]$out[, "nature"]
  }
  colnames(nature) <- labels2
  
  start.yr2 <- round_any(start.yr, 10, f = floor)
  end.yr2 <- round_any(end.yr, 5, f = ceiling)
  
  if(all(check3 == TRUE)) { 
    int.class <- matrix(nrow = length(yrs), ncol = length(list.sites))
    rownames(int.class) <- yrs
    for(i in 1:length(list.sites)) {
      int.class[as.character(list.sites[[i]]$out[, "year"]),i] <- ifelse(list.sites[[i]]$out[, "nature"] == (-1), 
             max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,6:8], ties.method = "first"), 
             max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,3:5], ties.method = "first"))
    }
    colnames(int.class) <- labels2
    
    input <- matrix2long(t(nature), new.ids = FALSE)
    input.int <- matrix2long(t(int.class), new.ids = FALSE)
    input[,4] <- input.int[,3]
    input[,4] <- ifelse(input[,3] == (-1), paste("-", input[,4], sep = ''), input[,4])
    input[,4] <- ifelse(input[,3] == 0, 0, input[, 4])
    input2 <- na.omit(input)
    rownames(input2) <- NULL
    colnames(input2) <-c ("site", "year", "PYvalues", "int.class")
    input3 <- subset(input2, input2[, "year"] >= start.yr & input2[, "year"] <= end.yr)
    rownames(input3) <- NULL
    
    if(sign2 == "both") {
      int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
      label.levels <- c("negative extreme", "negative strong", "negative weak", "none", 
                        "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(25, 25, 25, 95, 24, 24, 24)
      fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
    }
    if(sign2 == "pos") {
      input3[input3$int.class < 0, "int.class"] <- 0
      int.levels <- c(0, 1, 2, 3)
      label.levels <- c("other", "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(95, 24, 24, 24)
      fill.levels <- c("#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
    }
    if(sign2 == "neg") {
      input3[input3$int.class > 0, "int.class"] <- 0
      int.levels <- c(-3, -2, -1, 0)
      label.levels <- c("negative extreme", "negative strong", "negative weak", "other")
      shape.levels <- c(25, 25, 25, 95)
      fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7")
    }
    
    ggplot(input3, aes(x = year, y = site, shape = factor(int.class), 
                             fill = factor(int.class))) +
      geom_point(size = 2, colour = "black") +
      scale_shape_manual(name = "pointer year class", limits = factor(int.levels),
                         labels = label.levels, values = shape.levels) +
      scale_fill_manual(name = "pointer year class", limits = factor(int.levels), 
                        labels = label.levels, values = fill.levels) +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, x.tick.major), 
                         minor_breaks = seq(start.yr2, end.yr2, x.tick.minor),
                         limits = c(start.yr2, end.yr2)) +
      theme_bw() + theme(legend.key = element_blank())
  } else if(((TRUE %in% check3) && (FALSE %in% check3))) { 
      int.class <- matrix(nrow = length(yrs), ncol = length(list.sites))
      rownames(int.class) <- yrs
      seqNeuwirth <- which(check3 == TRUE)
      seqOther <- which(check3 != TRUE | is.na(check3))
      
      for(i in seqNeuwirth) {
        int.class[as.character(list.sites[[i]]$out[, "year"]),i] <- ifelse(list.sites[[i]]$out[, "nature"] == (-1), max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,6:8], ties.method = "first"), max.col(list.sites[[i]]$out[,c(1, 2, 5, 4, 3, 6:11)][,3:5], ties.method = "first"))
      }
      int.class[ int.class == 3 ] <- 4
      int.class[ int.class == 2 ] <- 3
      int.class[ int.class == 1 ] <- 2
      
      for(i in seqOther) {
        int.class[as.character(list.sites[[i]]$out[, "year"]),i] <- abs(list.sites[[i]]$out[,][,"nature"])
      }
      colnames(int.class) <- labels2
      
      input <- matrix2long(t(nature), new.ids = FALSE)
      input.int <- matrix2long(t(int.class), new.ids = FALSE)
      input[,4] <- input.int[,3]
      input[,4] <- ifelse(input[,3] == (-1), paste("-", input[,4], sep = ''), input[,4])
      input[,4] <- ifelse(input[,3] == 0, 0, input[, 4])
      input2 <- na.omit(input)
      rownames(input2) <- NULL
      colnames(input2) <-c ("site", "year", "PYvalues", "int.class")
      input3 <- subset(input2, input2[, "year"] >= start.yr & input2[, "year"] <= end.yr)
      rownames(input3) <- NULL
      
      if(sign2 == "both") {
        int.levels <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
        label.levels <- c("negative extreme", "negative strong", "negative weak", "negative", "none", 
                          "positive","positive weak", "positive strong", "positive extreme")
        shape.levels <- c(25, 25, 25, 25, 95, 24, 24, 24, 24)
        fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "grey", "#f7f7f7", "grey", "#d1e5f0", "#67a9cf", "#2166ac")
      }
      if(sign2 == "pos") {
        input3[input3$int.class < 0, "int.class"] <- 0
        int.levels <- c(0, 1, 2, 3, 4)
        label.levels <- c("other", "positive", "positive weak", "positive strong", "positive extreme")
        shape.levels <- c(95, 24, 24, 24, 24)
        fill.levels <- c("#f7f7f7", "grey", "#d1e5f0", "#67a9cf", "#2166ac")
      }
      if(sign2 == "neg") {
        input3[input3$int.class > 0, "int.class"] <- 0
        int.levels <- c(-4,-3, -2, -1, 0)
        label.levels <- c("negative extreme", "negative strong", "negative weak", "negative", "other")
        shape.levels <- c(25, 25, 25, 25, 95)
        fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "grey", "#f7f7f7")
      }
      
      ggplot(input3, aes(x = year, y = site, shape = factor(int.class), fill = factor(int.class))) +
        geom_point(size = 2, colour = "black") +
        scale_shape_manual(name = "pointer year class", limits = factor(int.levels), labels = label.levels, values = shape.levels) +
        scale_fill_manual(name = "pointer year class", limits = factor(int.levels), labels = label.levels, values = fill.levels) +
        scale_x_continuous(breaks = seq(start.yr2, end.yr2, x.tick.major), 
                           minor_breaks = seq(start.yr2, end.yr2, x.tick.minor),
                           limits = c(start.yr2, end.yr2)) +
        theme_bw() + theme(legend.key = element_blank())
    } else {
    nature2 <- nature[as.character(start.yr:end.yr),]
    input <- matrix2long(t(nature2), new.ids = FALSE)
    input2 <- na.omit(input)
    rownames(input2) <- NULL
    colnames(input2) <-c ("site", "year", "PYvalues")
      
    if(sign2 == "both") {
      int.levels <- c(-1, 0, 1)
      label.levels <- c("negative", "none", "positive")
      shape.levels <- c(25, 95, 24)
    }
    if(sign2 == "pos") {
      input2[input2$PYvalues < 0, "PYvalues"] <- 0
      int.levels <- c(0, 1)
      label.levels <- c("other", "positive")
      shape.levels <- c(95, 24)
    }
    if(sign2 == "neg") {
      input2[input2$PYvalues > 0, "PYvalues"] <- 0
      int.levels <- c(-1, 0)
      label.levels <- c("negative", "other")
      shape.levels <- c(25, 95)
    }
    
    ggplot(input2, aes(x = year, y = site, shape = factor(PYvalues))) +
      geom_point(size = 2, colour = "black", fill = "grey") +
      scale_shape_manual(name = "pointer year", limits = factor(int.levels),
                         labels = label.levels, values = shape.levels) + 
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, x.tick.major), 
                         minor_breaks = seq(start.yr2, end.yr2, x.tick.minor),
                         limits = c(start.yr2, end.yr2)) +
      theme_bw() + theme(legend.key = element_blank())
  }
}






