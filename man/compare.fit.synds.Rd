\name{compare.fit.synds}
\alias{compare.fit.synds}
\alias{print.compare.fit.synds}
\title{
 Compare model estimates based on synthesised and observed data
}
\description{
 The same model that was used for the synthesised data set is fitted to the 
 observed data set. The coefficients with confidence intervals for the 
 observed data is plotted together with their estimates from synthetic data. 
 When more than one synthetic data set has been generated (\code{object$m>1}) 
 combining rules are applied. 
}
\usage{
compare.fit.synds(object, obs.data, plot = "Z", return.result = TRUE,  
                  plot.intercept = FALSE, lwd = 1, lty = 1,
                  lcol = c("#1A3C5A","#4187BF"), dodge.height = .5, 
                  point.size = 2.5, ...)

\method{print}{compare.fit.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of type \code{fit.synds} created by fitting a model to 
   synthesised data set using function \code{\link{glm.synds}} or
   \code{\link{lm.synds}}.}
  \item{obs.data}{the original observed data set.}
  \item{plot}{values to be plotted: \code{"Z"} (Z scores) or \code{"coef"}
  (coefficients).}
  \item{return.result}{a logical value indicating whether a table of estimates 
   should be printed.}
  \item{plot.intercept}{a logical value indicating whether estimates for 
   intercept should be plotted.}
  \item{lwd}{the line type.}  
  \item{lty}{the line width.}
  \item{lcol}{line colours.}
  \item{dodge.height}{size of vertical shifts for confidence intervals to 
  prevent overlaping.} 
  \item{point.size}{size of plotting symbols used to plot point estimates of 
   coefficients.} 
  \item{\dots}{additional parameters passed to \code{\link{ggplot}}.}
  
  \item{x}{an object of class \code{compare.fit.synds.}}
}
\details{This function can be used to evaluate whether the model used for 
synthesis is appropriate for the fitted model. If this is the case the 
estimates from the synthetic data (\code{B.syn} and \code{Z.syn}) should 
not differ from the estimates from the observed data (\code{Beta} and 
\code{Z}) by more than would be expected from the standard errors 
(\code{se(B.syn)} and \code{se(Z.syn)}).
}

\value{
 An object of class \code{compare.fit.synds} which is a list with the 
 following components:
  \item{fit.synds.call}{the original call to fit the model to the 
   synthesised data set.}
  \item{coef.obs}{a data frame including estimates based on the observed
   data: coefficients (\code{Beta}), their standard errors (\code{se(Beta)})
   and Z scores (\code{Z}).}
  \item{coef.syn}{a data frame including (combined) estimates based on 
   the synthesised data: point estimates of observed data coefficients 
   (\code{B.syn}), standard errors of those estimates (\code{se(B.syn)}), 
   estimates of the observed standard errors (\code{se(Beta).syn}), Z scores 
   estimates (\code{Z.syn}) and their standard errors (\code{se(Z.syn)}). 
   Note that \code{se(B.syn)} and \code{se(Z.syn)} give the standard errors 
   of the mean of the \code{m} syntheses and can be made very small by 
   increasing \code{m}.}
  \item{ci.plot}{\code{ggplot} of the the coefficients with confidence 
   intervals for models based on observed and synthetic data.}
} 
\seealso{
  \code{\link{summary.fit.synds}}
}
\examples{
ods <- SD2011[,c("sex","age","edu","smoke")]
s1 <- syn(ods, m = 5)
f1 <- glm.synds(smoke ~ sex + age + edu, object = s1, family = "binomial")
compare.fit.synds(f1, ods) 
compare.fit.synds(f1, ods, plot = "coef")
}