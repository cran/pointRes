% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/res.plot.R
\name{res.plot}
\alias{res.plot}
\title{Plot resilience indices}
\usage{
res.plot(list.name, select.yr = NULL,
         param = c("resist", "recov", "resil", "rel.resil",
                   "rec.period", "avg.rec.rate",
                   "tot.abs.grow.red", "tot.rel.grow.red",
                   "avg.abs.grow.red", "avg.rel.grow.red"))
}
\arguments{
\item{list.name}{a \code{list} as produced by \code{\link{res.comp}}.}

\item{select.yr}{an \code{integer} or \code{vector} specifying the year(s) to be plotted (e.g., c(1948, 1992)).}

\item{param}{a \code{character} string specifying the resilience index to be plotted. Argument matching is performed.}
}
\value{
Box plot.
}
\description{
The function creates box plots for selected years of the resilience indices as calculated by \code{\link{res.comp}}, and is intended for quick visualization.
}
\details{
The function creates a box plot for a selected resilience index showing the full range of variation for individual trees in specific years. Box plots are only created for years for which indices are available for >= 5 series, as this value represents the number of statistics that a box plot represents in its' simplest form.
}
\examples{
## Plot the recovery period for three selected years
data(s033)
res <- res.comp(s033)
res.plot(res, select.yr = c(1976, 1992, 2003), param = "resist")

}
\author{
Marieke van der Maaten-Theunissen and Ernst van der Maaten.
}
