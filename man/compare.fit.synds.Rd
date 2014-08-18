\name{compare.fit.synds}
\alias{compare.fit.synds}
\alias{print.compare.fit.synds}
\title{
 Compare model estimates based on synthesised and real data
}
\description{
 The same model that was used for the synthesised data set is fitted to the 
 real data set. The coefficients with confidence intervals for both model
 estimates are plotted for comparison. When more than one 
 synthetic data set has been generated (object$m>1) average of the means 
 and standard errors are used.
}
\usage{
compare.fit.synds(object, real.data, plot = "Z", 
  return.result = TRUE, plot.intercept = FALSE, col = 1:2, ...)

\method{print}{compare.fit.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of type 'fit.synds' created by fitting a model to 
   synthesised data set using function \code{\link{glm.synds}} or
   \code{\link{lm.synds}}.}
  \item{real.data}{the original (real) data set.}
  \item{plot}{values to be plotted: "Z", "coefficients" or "both".}
  \item{return.result}{a logical value indicating whether a table of estimates 
   should be printed.}
  \item{plot.intercept}{a logical value indicating whether estimates for 
   intercept should be plotted.}
  \item{col}{vector of colours for lines and points for plotting real and 
   synthesised estimates.}
  \item{\dots}{additional parameters passed to \code{\link{plot}}.}
  \item{x}{an object of class \code{compare.fit.synds.}}
}
\value{
 Plot(s) of the coefficients with confidence intervals for both models. 
 
 An object of class \code{compare.fit.synds} which is a list with the 
 following components:
  \item{fit.synds.call}{the original call to fit the model to the 
   synthesised data set.}
  \item{coefficients}{a data frame including coefficients, Z-scores and  
   standard errors for models fitted to real and synthetic data sets
   and synthesising errors for coefficients and Z-scores.}
} 
\seealso{
  \code{\link{summary.fit.synds}}
}
\examples{
  rds <- SD2011[,c("sex","age","edu","smoke")]
  s1 <- syn(rds, m = 5)
  f1 <- glm.synds(smoke ~ sex + age + edu, object = s1, family = "binomial")
  compare.fit.synds(f1, rds) 
  compare.fit.synds(f1, rds, plot = "both")
}