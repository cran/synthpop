\name{compare.synds}
\alias{compare.synds}
\alias{print.compare.synds}
\title{
 Compare synthesised and observed data
}
\description{
 Compare synthesised data set with the original (observed) data set 
 using percent frequency tables and histograms. When more than one 
 synthetic data set has been generated (object$m>1), only the 
 first one is used for comparison.
}
\usage{
compare.synds(object, data, vars = NULL, ncol = 4, breaks = 20, 
plot.na = TRUE, rel.size.x = 1, cols = c("#1A3C5A","#4187BF"), ...)

\method{print}{compare.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original (observed) data set.}
  \item{vars}{variables to be compared. If \code{vars} is \code{NULL} 
    (the default) all synthesised variables are compared.}
  \item{ncol}{the number of columns for the plotting area.}
  \item{breaks}{the number of cells for the histogram.}
  \item{plot.na}{a logical variable indicating whether missing data categories 
    for numeric variables should be ploted.}  
  \item{rel.size.x}{a number representing the relative size of x-axis labels.}
  \item{cols}{bar colors.}
  \item{\dots}{additional parameters passed to \code{\link{barplot}}.}

  \item{x}{an object of class \code{compare.synds.}}
}
\value{
 An object of class \code{compare.synds} which is a list including ggplot 
 object with histograms and a list of comparative percent frequency tables.
}
\examples{
ods <- SD2011[ , c("sex","age","edu","marital","ls","income")]
s1  <- syn(ods)
compare.synds(s1, ods, vars = "ls")
compare.synds(s1, ods, vars = "income")
}
