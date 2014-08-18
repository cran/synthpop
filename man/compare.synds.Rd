\name{compare.synds}
\alias{compare.synds}
\title{
 Compare synthesised and real data
}
\description{
 Compare synthesised data set with the original (real) data set 
 using percent frequency tables and histograms. When more than one 
 synthetic data set has been generated (object$m>1), only the 
 first one is used for comparison.
}
\usage{
 compare.synds(object, data, vars = NULL, nrows = 1, ncols = 1, breaks = 50, ...)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original (real) data set.}
  \item{vars}{variables to be compared. If \code{vars} is \code{NULL} 
    (the default) all synthesised variables are compared.}
  \item{nrows}{the number of rows for the plotting area.}
  \item{ncols}{the number of columns for the plotting area.}
  \item{breaks}{the number of cells for the histogram.}
  \item{\dots}{additional parameters passed to \code{\link{barplot}}.}
}
\value{
 Comparative percent frequency tables and histograms.
}
\examples{
  rds <- SD2011[ ,c("sex","age","edu","marital","ls","income")]
  s1  <- syn(rds)
  compare.synds(s1, rds, vars = "ls")
  compare.synds(s1, rds, vars = "income")
}