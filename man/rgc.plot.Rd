% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgc.plot.R
\name{rgc.plot}
\alias{rgc.plot}
\title{Plot mean relative growth changes and pointer years}
\usage{
rgc.plot(list.name, start.yr = NULL, end.yr = NULL, 
         sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)
}
\arguments{
\item{list.name}{a \code{list} as produced by \code{\link{pointer.rgc}}}

\item{start.yr}{an \code{integer} specifying the first year to be plotted. Defaults to the first year with data if \code{\var{start.yr}} is \code{NULL}.}

\item{end.yr}{an \code{integer} specifying the last year to be plotted. Defaults to the last year with data if \code{\var{end.yr}} is \code{NULL}.}

\item{sd.disp}{a \code{logical} specifying whether error bars (stdev) should be displayed. Defaults to FALSE.}

\item{x.tick.major}{an \code{integer} controlling the major x-axis tick labels. Defaults to 10 years.}

\item{x.tick.minor}{an \code{integer} controlling the minor x-axis ticks. Defaults to 5 years.}
}
\value{
Bar plot.
}
\description{
The function creates a bar plot of mean relative growth changes from a \code{list} of the type as produced by \code{\link{pointer.rgc}} and highlights years identified as pointer years.
}
\details{
The function makes a plot showing mean relative growth changes; pointer years are indicated with dark-gray bars. Error bars can be set.
}
\examples{
## Plot mean relative growth changes and pointer years
data(s033)
py <- pointer.rgc(s033, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, 
                  series.thresh = 75)
rgc.plot(py, start.yr = 1950, end.yr = NULL,  
         sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)

}
\author{
Marieke van der Maaten-Theunissen and Ernst van der Maaten.
}

