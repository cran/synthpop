\name{tab.utility}
\alias{tab.utility}
\alias{print.tab.utility}
\title{[EXPERIMENTAL] Tabular utility}
\description{Produce tables from observed and synthesized data and calculates 
  chi-square statistics to compare them.}
\usage{
tab.utility(object, data, vars = NULL, ngroups = 5, with.missings = TRUE, ...)

\method{print}{tab.utility}(x, tables = FALSE, digits = 2, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s).}
  \item{data}{the original (observed) data set.}
  \item{vars}{a single string or a vector of strings with the names of 
    variables to be used to form the table.}
  \item{ngroups}{if numerical (non-factor) variables are included they will be 
    classified into this number of groups to form tables. Repeated identical 
    obsevations may produce a smaller number of groups in some cases.}
  \item{with.missings}{a logical value indicating whether the table is to 
    include rows and columns for missing data categories.}
  \item{\dots}{additional parameters.}   
  \item{x}{an object of class \code{tab.utility}.}
  \item{tables}{a logical value that determines if tables of observed, 
    synthesised and Z differences are to be printed.}
  \item{digits}{an integer indicating the number of decimal places 
    for printing \code{tab.Zdiff}.}
}
\details{Forms tables of observed and synthesised values for the variables 
  specified in \code{vars}. A chi-square statistic is calculated from the cells 
  of the tables, as (observed-synthesied)^2/[(observed + synthesised)/2], 
  ignoring those where observed and synthesised are both zero.
}
\value{An object of class \code{tab.utility} which is a list with the following 
   components: 
 \item{Chisq}{a vector with \code{object$m} values for the test.} 
 \item{df}{a vector with corresponding degrees of freedom.} 
 \item{ratio}{a vector with ratios of \code{chisq} to \code{df}.} 
 \item{nempty}{a vector of length \code{object$m} with number of empty cells 
   not contributing to the chi-square statistic.} 
 \item{pval}{a vector of length \code{object$m} with p-values for the 
   chi-square test.} 
 \item{tab.obs}{a table from the observed data.} 
 \item{tab.syn}{a table or a list of \code{object$m} tables from the synthetic 
   data.}  
 \item{tab.Zdiff}{a table or a list of \code{object$m} tables of Z statistics 
   for differences between observed and synthesised cells of the tables. Large 
   absolute values indicate a large contribution to lack-of-fit.}  
}
\seealso{
 \code{\link{utility.synds}}
}
\examples{
  ods <- SD2011[1:1000, c("sex", "age", "edu", "marital")]
  s1 <- syn(ods, m = 3)
  t1 <- tab.utility(s1, ods, vars = c("marital", "sex"))
  print(t1, tables = TRUE)
}

