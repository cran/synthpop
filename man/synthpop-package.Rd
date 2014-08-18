\name{synthpop-package}
\alias{synthpop-package}
\alias{synthpop}
\docType{package}
\title{Generating synthetic versions of sensitive microdata for statistical disclosure control}
\description{
Generate synthetic versions of a data set using parametric or CART methods.}
\details{
\tabular{ll}{
Package: \tab synthpop\cr
Type: \tab Package\cr
Version: \tab 1.0-0\cr
Date: \tab 2014-08-03\cr
License: \tab GPL-2 | GPL-3\cr
}
Synthetic data are generated from the real data by the function \code{\link{syn}}. 
The package includes also tools to compare synthetic data with the real data 
(\code{\link{compare.synds}}) and to fit (generalized) linear model to synthetic data 
(\code{\link{lm.synds}},\code{\link{glm.synds}}) and compare the estimates with those 
for the real data (\code{\link{compare.fit.synds}}). More extensive documentation with 
illustrative examples is provided in the package vignette.}
\author{
Beata Nowok, Gillian M Raab and Chris Dibben (first two authors in alphabetical order) 
based on package \pkg{mice} (2.18) by Stef van Buuren and Karin Groothuis-Oudshoorn

Maintainer: Beata Nowok <beata.nowok@gmail.com>
}
\keyword{package}
