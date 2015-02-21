\name{syn.sample}
\alias{syn.sample}
\title{Synthesis by simple random sampling}
\description{
 Generates a random sample from the observed data.
}
\usage{
syn.sample(y, xp, proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{xp}{a target length \code{k} of a synthetic data vector.}
  \item{proper}{if \code{proper = TRUE} values are sampled from  
    a bootstrapped sample of the original data.}  
  \item{\dots}{additional parameters passed to \code{sample}.}
}
\details{
 A simple random sample with replacement is taken from the 
 observed values in \code{y} and used as synthetic values.
}
\value{
 A vector of length \code{k} with synthetic values.
}
\seealso{
  \code{\link{syn}}
}
\keyword{datagen}
