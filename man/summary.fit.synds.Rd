\name{summary.fit.synds}
\alias{summary.fit.synds}
\alias{print.summary.fit.synds}
\title{Inference from synthetic data}
\description{Combines the results of models fitted to each of the \code{m} 
  synthetic data sets.}
\usage{
\method{summary}{fit.synds}(object, populationInference = FALSE, \dots)

\method{print}{summary.fit.synds}(x, \dots)
}

\arguments{
  \item{object}{an object of class \code{fit.synds} created by fitting a model to 
   synthesised data set using function \code{\link{glm.synds}} or
   \code{\link{lm.synds}}.}
  \item{populationInference}{a logical value indicating whether inference 
   should be made to population quantities. If \code{FALSE} inference is made
   to original data quantities.}
  \item{\dots}{additional parameters.}
  \item{x}{an object of class \code{summary.fit.synds}.}
}
\details{The mean of the estimates from each of the \code{m} synthetic data sets yields unbaised estimates of the coefficients. The variance is estimated differently depending whether inference is made to the real data quantities or the population parameters and whether synthetic data were produced using simple or proper synthesis (for details see Raab et al. (2014); expressions used to calculate variance for different cases are presented in Table 1).}

\value{
  An object of class \code{summary.fit.synds} which is a list with the following 
  components:
  \item{call}{the original call to \code{glm.synds} or \code{lm.synds}.}
  \item{proper}{a logical value indicating whether synthetic data were 
    generated using proper synthesis.}  
  \item{fitting.function}{function used to fit the model.} 
  \item{m}{the number of synthetic versions of the real data.} 
  \item{coefficients}{a matrix with combined estimates. For inference 
   to original data quantities it includes coefficients (\code{beta syn}) 
   and their standard errors (\code{se beta syn}), Z scores (\code{Z syn}) 
   and synthesising errors for coefficients (\code{syn err beta}) 
   and Z-scores (\code{syn err Z}). For inference to population quantities
   it includes coefficients (\code{est}), synthesising errors (\code{se})
   and Z scores (\code{Z}).} 
  \item{n}{a number of cases in the original data.} 
  \item{k}{a number of cases in the synthesised data.}
}
\references{
Raab, G.M., Nowok, B., Dibben, C. (2014). 
A simplified approach to synthetic data.
\emph{Submitted}. 

}
\seealso{
\code{\link{summary}},\code{\link{print}}
}
\examples{
  rds <- SD2011[1:2000,c("sex","age","edu","ls","smoke")]
  
  ### simple synthesis
  s1 <- syn(rds, m = 5)
  f1 <- glm.synds(smoke ~ sex + age + edu + ls, object = s1, family = "binomial")
  summary(f1)
  summary(f1, populationInference = TRUE)
  
  ### proper synthesis
  s2 <- syn(rds, m = 5, proper = TRUE)
  f2 <- glm.synds(smoke ~ sex + age + edu + ls, object = s2, family = "binomial")
  summary(f2)
  summary(f2, populationInference = TRUE)
}