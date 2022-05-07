#' Plot resilience indices
#'
#' @description The function creates box plots for selected years of the resilience indices as calculated by \code{\link{res.comp}}, and is intended for quick visualization.
#' 
#' @usage res.plot(list.name, select.yr = NULL,
#'          param = c("resist", "recov", "resil", "rel.resil",
#'                    "rec.period", "avg.rec.rate",
#'                    "tot.abs.grow.red", "tot.rel.grow.red",
#'                    "avg.abs.grow.red", "avg.rel.grow.red"))
#'
#' @param list.name a \code{list} as produced by \code{\link{res.comp}}.
#' @param select.yr an \code{integer} or \code{vector} specifying the year(s) to be plotted (e.g., c(1948, 1992)).
#' @param param a \code{character} string specifying the resilience index to be plotted. Argument matching is performed.
#' 
#' @details The function creates a box plot for a selected resilience index showing the full range of variation for individual trees in specific years. Box plots are only created for years for which indices are available for >= 5 series, as this value represents the number of statistics that a box plot represents in its' simplest form.
#'
#' @return 
#' Box plot.
#'
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @examples ## Plot the recovery period for three selected years
#' data(s033)
#' res <- res.comp(s033)
#' res.plot(res, select.yr = c(1976, 1992, 2003), param = "resist")
#' 
#' @import ggplot2
#' @import stats
#' @importFrom TripleR matrix2long
#' 
#' @export res.plot
#' 
res.plot <- function(list.name, select.yr = NULL, param = c("resist", "recov", "resil", "rel.resil", "rec.period", "avg.rec.rate", "tot.abs.grow.red", "tot.rel.grow.red", "avg.abs.grow.red", "avg.rel.grow.red"))
{
  stopifnot(is.list(list.name))
  if(inherits(list.name, "res.comp") == FALSE) {
    stop("'list.name' is no list output of function res.comp")
  }
  if(is.matrix(list.name$resist) == FALSE) {
    stop("'list.name' is no list output of function res.comp")
  }
  if(is.null(select.yr)) {
    stop("'select.yr' is NULL: please specify the years to be plotted")
  }
  
  param2 <- match.arg(param, c("resist", "recov", "resil", "rel.resil", "rec.period", "avg.rec.rate", "tot.abs.grow.red", "tot.rel.grow.red", "avg.abs.grow.red", "avg.rel.grow.red"))
  
  list.param <- as.data.frame(list.name[[noquote(param2)]])
  list.param[list.param == "Inf"] <- NA
  list.param$nb.series <- rowSums(!is.na(list.param))
  
  if(FALSE %in% (select.yr %in% rownames(list.param))) {
    stop("the selection under 'select.yr' contains year(s) that are not in the dataset")
  }
  
  if(length(select.yr) == 1) {
    yrs <- as.character(select.yr)
    list.param2 <- subset(list.param,rownames(list.param) %in% select.yr & list.param[, "nb.series"] >= 5)
    if(nrow(list.param2) == 0) {
      stop("the year can not be displayed, as the selected index is available for < 5 series")
    }
    nb.series <- NULL
    list.param2 <- subset(list.param2, select = -nb.series)
    list.param2 <- matrix2long(list.param2, new.ids = FALSE)
    colnames(list.param2) <- c("year","tree","value")
  }
  else {
    yrs <- as.character(select.yr)
    list.param2 <- subset(list.param,rownames(list.param) %in% select.yr & list.param[, "nb.series"] >= 5)
    if(FALSE %in% (yrs %in% rownames(list.param2))) {
      warning("years for which the selected index is available for < 5 series are not displayed")
    }
    if(nrow(list.param2) == 0) {
      stop("years can not be displayed, as the selected index is available for < 5 series in all cases")
    }
    yrs <- rownames(list.param2)
    list.param2 <- subset(list.param2, select = -nb.series)
    list.param2 <- matrix2long(list.param2, new.ids = FALSE)
    colnames(list.param2) <- c("year","tree","value") 
  }

  year <- value <- NULL
  index <- ifelse(param2 == "resist", "resistance index",
                  ifelse(param2 == "recov", "recovery index",
                         ifelse(param2 == "resil", "resilience index", 
                                ifelse(param2 == "rel.resil", "relative resilience index",
                                       ifelse(param2 == "rec.period", "recovery period",
                                              ifelse(param2 == "avg.rec.rate", "average recovery rate",
                                                     ifelse(param2 == "tot.abs.grow.red", "total absolute growth reduction",
                                                            ifelse(param2 == "tot.rel.grow.red", "total relative growth reduction",
                                                                   ifelse(param2 == "avg.abs.grow.red", "average absolute growth reduction",
                                                                          ifelse(param2 == "avg.rel.grow.red", "average relative growth reduction"))))))))))
  plot.param <- ggplot(na.omit(list.param2), aes(x = factor(year), y = value)) + 
    geom_boxplot(fill = "#f0f0f0") + theme_bw() + guides(fill = "none") +
    xlab("year") + ylab(index) + scale_x_discrete(labels = yrs)
    
  plot(plot.param)
  message("Please note that the number of series upon which individual boxplots are based can be queried by calling the component 'nb.series' from the res.comp output. Numbers may deviate between years and indices, depending on tree-specific growth behaviour as well as characteristics inherent to the calculation methods")
}




