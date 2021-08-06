#' Plot event years for individual trees
#'
#' @description The function creates a dot plot showing positive and (or) negative event years from a \code{list} of the type as produced by \code{\link{pointer.norm}} or \code{\link{pointer.rgc}}.
#' 
#' @usage event.plot(list.name, sign = c("both", "pos", "neg"),
#'            period = NULL, x.tick.major = 10, x.tick.minor = 5) 
#'
#' @param list.name a \code{list} as produced by \code{\link{pointer.norm}} or \code{\link{pointer.rgc}}
#' @param sign a \code{character} string specifying whether both positive and negative (\code{"both"}), or only positive (\code{"pos"}) or negative (\code{"neg"}) event years should be displayed. Defaults to \code{"both"}.
#' @param period a \code{vector} specifying the start and end year to be plotted. Defaults to the full period covered by the output of the pointer year analysis.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a dot plot showing event years for individual trees. Positive and negative event years are indicated with different symbols and (or) colors.
#' 
#' @return 
#' Dot plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot event years from pointer.rgc output
#' data(s033)
#' py <- pointer.rgc(s033)
#' event.plot(py) 
#'
#' ## Plot negative event years from pointer.norm output (method "Neuwirth") for a specific period
#' data(s033)
#' detr_s033 <- detrend(s033, method = "Spline", nyrs = 30)
#' pyn <- pointer.norm(detr_s033, method.thresh = "Neuwirth")
#' event.plot(pyn, sign = "neg", period = c(1950, 2007)) 
#'            
#' @import ggplot2
#' @import stats
#' @importFrom plyr round_any
#' @importFrom TripleR matrix2long
#' 
#' @export event.plot
#' 
event.plot <- function(list.name, sign = c("both", "pos", "neg"), period = NULL, x.tick.major = 10, x.tick.minor = 5) 
{
  stopifnot(is.list(list.name))
  if(FALSE %in% (class(list.name)[1] == "pointer.norm") & FALSE %in% (class(list.name)[1] == "pointer.rgc")) {
    stop("'list.name' is no list output of function pointer.norm or pointer.rgc or ")
  }
  if(is.matrix(list.name$EYvalues) == FALSE) {
    stop("'list.name' is no list output of function pointer.norm or pointer.rgc")
  }
  if(is.null(period)) {
    start.yr <- min(as.numeric(rownames(list.name$EYvalues)))
    end.yr <- max(as.numeric(rownames(list.name$EYvalues)))
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
  if(start.yr < min(rownames(list.name$EYvalues))) {
    stop("the start year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  if(end.yr > max(rownames(list.name$EYvalues))) {
    stop("the end year in 'period' is out of bounds. By default (period = NULL) the calculations are performed over the whole period covered by the data")
  }
  if(x.tick.minor > x.tick.major) {
    stop("'x.tick.minor' should be smaller then 'x.tick.major'")
  }
 
  sign2 <- match.arg(sign, c("both", "pos", "neg"))
  
  start.yr2 <- round_any(start.yr, 10, f = floor)
  end.yr2 <- round_any(end.yr, 5, f = ceiling)
  
  matrix <- list.name$EYvalues[as.character(start.yr:end.yr),]
  input <- matrix2long(t(matrix), new.ids = FALSE)
  input2 <- na.omit(input) 
  rownames(input2) <- NULL
  colnames(input2) <-c ("tree","year","EYvalues")
  
  year <- tree <- EYvalues <- NULL
  
  col.pos <- colnames(list.name$out)[3] == "perc.pos"
  
  if(col.pos) {
    if(sign2 == "both") {
      int.levels <- c(-1, 0, 1)
      label.levels <- c("negative", "none", "positive")
      shape.levels <- c(25, 95, 24)
    }
    if(sign2 == "pos") {
      input2[input2$EYvalues < 0, "EYvalues"] <- 0
      int.levels <- c(0, 1)
      label.levels <- c("other", "positive")
      shape.levels <- c(95, 24)
    }
    if(sign2 == "neg") {
      input2[input2$EYvalues > 0, "EYvalues"] <- 0
      int.levels <- c(-1, 0)
      label.levels <- c("negative", "other")
      shape.levels <- c(25, 95)
    }
    
    ggplot(input2, aes(x = year, y = tree, shape = factor(EYvalues))) +
      geom_point(size = 2, colour = "black", fill = "black") +
      scale_shape_manual(name = "event year", limits = factor(int.levels),
                         labels = label.levels, values = shape.levels) + 
      theme_bw() + theme(legend.key = element_blank()) +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, x.tick.major), 
                         minor_breaks = seq(start.yr2, end.yr2, x.tick.minor),
                         limits = c(start.yr2, end.yr2))
  }
  else {
    if(sign2 == "both") {
      int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
      label.levels <- c("negative extreme", "negative strong", "negative weak",
                       "none", "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(25, 25, 25, 95, 24, 24, 24)
      fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
    }
    if(sign2 == "pos") {
      input2[input2$EYvalues < 0, "EYvalues"] <- 0
      int.levels <- c(0, 1, 2, 3)
      label.levels <- c("other", "positive weak", "positive strong", "positive extreme")
      shape.levels <- c(95, 24, 24, 24)
      fill.levels <- c("#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
    }
    if(sign2 == "neg") {
      input2[input2$EYvalues > 0, "EYvalues"] <- 0
      int.levels <- c(-3, -2, -1, 0)
      label.levels <- c("negative extreme", "negative strong", "negative weak", "other")
      shape.levels <- c(25, 25, 25, 95)
      fill.levels <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7")
    }
    
    ggplot(input2, aes(x = year, y = tree, shape = factor(EYvalues), fill = factor(EYvalues))) +
      geom_point(size = 2, colour = "black") +
      scale_shape_manual(name = "event year class", limits = factor(int.levels),
                         labels = label.levels, values = shape.levels) +
      scale_fill_manual(name = "event year class", limits = factor(int.levels), 
                        labels = label.levels, values = fill.levels) +
      theme_bw() + theme(legend.key = element_blank()) +
      scale_x_continuous(breaks = seq(start.yr2, end.yr2, x.tick.major), 
                         minor_breaks = seq(start.yr2, end.yr2, x.tick.minor),
                         limits = c(start.yr2, end.yr2))
  }
}