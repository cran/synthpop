\name{summary.synds}
\alias{summary.synds}
\alias{print.summary.synds}
\title{Synthetic data object summaries}
\description{
Produce summaries of the synthesised variables for selected copies of synthetic data.
}
\usage{
\method{summary}{synds}(object, msel = 1, \dots)

\method{print}{summary.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}; a result of a call to \code{\link{syn}}.}
  \item{msel}{index or indices of synthetic data copies for which a summary is desired.}
  \item{\dots}{additional arguments passed to \code{\link{summary}}.}
  \item{x}{an object of class \code{summary.synds}.}
}
\details{See \code{\link[base]{summary}} for more details.}
\value{
An object of class \code{summary.synds}, which is a list with the following components:
  \item{m}{the number of synthetic versions of the real data.}
  \item{msel}{index or indices of synthetic data copies for which a summary is produced.}
  \item{method}{a vector of synthesising methods applied to each variable
    in the saved synthesised data.}
  \item{result}{a table or a list of tabels (if more than one synthetic data set is selected)  
    with summaries of synthesised variables.}
}
\seealso{
\code{\link{summary}},\code{\link{print}}
}
\examples{
s1 <- syn(SD2011[,c("sex","age","edu","marital")], m = 3)
summary(s1)
summary(s1, msel = c(1,3))
}
