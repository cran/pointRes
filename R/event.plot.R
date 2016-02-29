#' Plot event years for individual trees
#'
#' @description The function creates a dot plot showing positive and negative event year values from a \code{list} of the type as produced by \code{\link{pointer.norm}} or \code{\link{pointer.rgc}}.
#' 
#' @usage event.plot(list.name, start.yr = NULL, end.yr = NULL,
#'            x.tick.major = 10, x.tick.minor = 5) 
#'
#' @param list.name a \code{list} as produced by \code{\link{pointer.norm}} or \code{\link{pointer.rgc}}
#' @param start.yr an \code{integer} specifying the first year to be plotted. Defaults to the first year with data if \code{\var{start.yr}} is \code{NULL}.
#' @param end.yr an \code{integer} specifying the last year to be plotted. Defaults to the last year with data if \code{\var{end.yr}} is \code{NULL}.
#' @param x.tick.major an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.
#' @param x.tick.minor an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.
#' 
#' @details The function makes a dot plot showing event years for individual trees. Positive and negative event years are indicated with different symbols. If event years were defined using \code{method.thresh "Neuwirth"} (\code{\link{pointer.norm}}), different tones of gray indicate weak, strong and extreme event years.
#' 
#' Non-event years are indicated with minus-signs, allowing the assessment of individual series length.
#' 
#' @return 
#' Dot plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot event years for individual trees as generated using pointer.rgc
#' data(s033)
#' py <- pointer.rgc(s033, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, 
#'                   series.thresh = 75)
#' event.plot(py, start.yr = 1950, end.yr = NULL,
#'            x.tick.major = 10, x.tick.minor = 5) 
#'
#' ## Plot event years for individual trees as generated using pointer.norm (method "Neuwirth")
#' data(s033)
#' py_n <- pointer.norm(s033, window = 5, method.thresh = "Neuwirth",
#'                      series.thresh = 75)
#' event.plot(py_n, start.yr = 1950, end.yr = NULL,
#'            x.tick.major = 10, x.tick.minor = 5) 
#'            
#' @import ggplot2
#' @import stats
#' @importFrom plyr round_any
#' @importFrom TripleR matrix2long
#' 
#' @export event.plot
#' 
event.plot <- function(list.name, start.yr = NULL, end.yr = NULL, x.tick.major = 10, x.tick.minor = 5) 
{
  stopifnot(is.list(list.name))
  if(FALSE %in% (class(list.name)[1] == "pointer.rgc") & FALSE %in% (class(list.name)[1] == "pointer.norm")) {
    stop("'list.name' is no list output of function pointer.rgc or pointer.norm")
  }
  if(is.matrix(list.name$EYvalues) == FALSE) {
    stop("'list.name' is no list output of function pointer.rgc or pointer.norm")
  }
  if(!is.null(start.yr) && start.yr < min(list.name$out[, "year"])) {
    stop("'start.yr' is out of bounds. By default (start.yr = NULL) the first year is displayed")
  }
  if(!is.null(end.yr) && end.yr > max(list.name$out[, "year"])) {
    stop("'end.yr' is out of bounds. By default (end.yr = NULL) the last year is displayed")
  }
  if(x.tick.minor > x.tick.major) {
    stop("'x.tick.minor' should be smaller then 'x.tick.major'")
  }
 
  start.yr2 <- ifelse(length(start.yr) != 0, start.yr, min(list.name$out[, "year"]))
  end.yr2 <- ifelse(length(end.yr) != 0, end.yr, max(list.name$out[, "year"]))
  start.yr3 <- round_any(start.yr2, 10, f = floor)
  end.yr3 <- round_any(end.yr2, 5, f = ceiling)
  
  matrix <- list.name$EYvalues[as.character(start.yr2:end.yr2),]
  input <- matrix2long(t(matrix), new.ids=FALSE)
  input2 <- na.omit(input)
  rownames(input2) <- NULL
  colnames(input2) <-c ("tree","year","EYvalues")
  
  year <- tree <- EYvalues <- NULL
  
  col.pos <- colnames(list.name$out)[3] == "perc.pos"
  
  if(col.pos) {
    int.levels <- c(-1, 0, 1)
    label.levels <- c("negative", "none", "positive")
    shape.levels <- c(25, 95, 24)
    
    pl <- ggplot(input2, aes(x = year, y = tree, shape = factor(EYvalues))) 
    pl + geom_point(size = 2, colour = "black", fill = "#bdbdbd") +
      scale_shape_manual(name = "event year", limits = int.levels,
                         labels = label.levels, values = shape.levels) + 
      theme_bw() + theme(legend.key = element_blank()) +
      scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                         minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                         limits = c(start.yr3, end.yr3))
  }
  else {
    int.levels <- c(-3, -2, -1, 0, 1, 2, 3)
    fill.levels <- c( "black", "#bdbdbd", "white" , "#bdbdbd", "white", "#bdbdbd", "black")
    label.levels <- c("negative extreme", "negative strong", "negative weak",
                       "none", "positive weak", "positive strong", "positive extreme")
    shape.levels <- c(25, 25, 25, 95, 24, 24, 24)
    
    pl <- ggplot(input2, aes(x = year, y = tree, shape = factor(EYvalues),
                             fill = factor(EYvalues))) 
    pl + geom_point(size = 2, colour = "black") +
      scale_shape_manual(name = "event year class", limits = int.levels,
                         labels = label.levels, values = shape.levels) +
      scale_fill_manual(name = "event year class", limits = int.levels, 
                        labels = label.levels, values = fill.levels) +
      theme_bw() + theme(legend.key = element_blank()) +
      scale_x_continuous(breaks = seq(start.yr3, end.yr3, x.tick.major), 
                         minor_breaks = seq(start.yr3, end.yr3, x.tick.minor),
                         limits = c(start.yr3, end.yr3))
  }
}






