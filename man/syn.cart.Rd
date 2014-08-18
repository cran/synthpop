\name{syn.ctree, syn.cart}
\alias{syn.cart}
%\alias{syn.cart.bboot}
\alias{syn.ctree}
\title{Synthesis by classification and regression trees (CART)}
\description{
  Generates univariate synthetic data using classification and
  regression trees (without or with bootstrap).
}
\usage{
syn.ctree(y, x, xp, smoothing, proper = FALSE, minbucket, ...)
syn.cart(y, x, xp, smoothing, proper = FALSE, minbucket, cp = 1e-04, ...)  
%syn.cart.bboot(y, x, xp, minbucket, cp = 1e-04, ...) 
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{smoothing}{smoothing method for continuous variables.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) a CART 
    model is fitted to a bootstrapped sample of the original data.}  
  \item{minbucket}{the minimum number of observations in
    any terminal node. See \code{\link[rpart]{rpart.control}} and
    \code{\link[party]{ctree_control}} for details.}
  \item{cp}{complexity parameter. Any split that does not
    decrease the overall lack of fit by a factor of cp is not
    attempted. See \code{\link[rpart]{rpart.control}} for details.}

  \item{\dots}{additional parameters.}
}
\details{
  The procedure for synthesis by a CART model is 
  as follows: \enumerate{ \item Fit a classification 
  or regression tree by binary recursive partitioning. 
  \item For each \code{xp} find the terminal node. \item Randomly 
  draw a donor from the members of the node and take the observed 
  value of \code{y} from that draw as the synthetic value.}
  
  \code{syn.ctree} uses \code{\link{ctree}} function from the
  \pkg{party} package and \code{syn.cart} uses \code{\link{rpart}}
  function from the \pkg{rpart} package.  They differ, among others, 
  in a stopping rule for the splitting process.
}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
%\references{...}
\seealso{
  \code{\link{syn}}, \code{\link{syn.surv.ctree}},
  \code{\link[rpart]{rpart}}, \code{\link[party]{ctree}}
}
\keyword{datagen}
