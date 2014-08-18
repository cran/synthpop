\name{syn}
\alias{syn}
\alias{print.synds}
\title{Generating synthetic data sets}
\description{
  Generates synthetic version(s) of a data set.
}

\usage{
 syn(data, method = vector("character", length = ncol(data)), 
   visitSequence = (1:ncol(data)), predictorMatrix = NULL,  
   m = 1, k = nrow(data), proper = FALSE, nlevelmax = 5,
   maxfaclevels = 60,
   rules = as.list(rep("", ncol(data))), 
   rvalues = as.list(rep(NA, ncol(data))), 
   contNA = as.list(rep(NA, ncol(data))),
   event = rep(0, ncol(data)),
   smoothing = rep("", ncol(data)), 
   denom = rep(0, ncol(data)),
   minbucket = 5,
   drop.not.used = TRUE,
   drop.pred.only = TRUE,
   defaultMethod = c("normrank", "logreg", "polyreg", "polr"), 
   diagnostics = FALSE, printFlag = TRUE, ...)
   
\method{print}{synds}(x, \dots)   
}

\arguments{
  \item{data}{a data frame or a matrix (\code{n} x \code{p}) containing the original data. 
     Observations are in rows and variables are in columns.}
  \item{method}{a single string or a vector of strings of length
     \code{ncol(data)} specifying the synthesising method to be
     used for each variable in the data. Order of variables is exactly the
     same as in \code{data}. If specified as a single string, the same method
     is used for all variables in a visit sequence unless a data type or
     a position in a visit sequence requires a different method. 
     If \code{method} is set to \code{"parametric"} the
     default synthesising method specified by the \code{defaultMethod} argument
     are applied. Variables that are transformations of other variables can
     be synthesised using a passive method that is specified as a string 
     starting with \code{~}. Variables that need not to be synthesised have the empty
     method \code{""}. By default all variables are synthesised using
     \code{ctree} implementation of a CART model. See details for more
     information.}
  \item{visitSequence}{a vector of integers of arbitrary length, specifying 
     the column indices of the synthesising/visiting sequence.
     The default sequence \code{1:ncol(data)} implies that column variables are
     synthesised/visited from left to right. See details for more information.}
  \item{predictorMatrix}{a square matrix of size \code{ncol(data)} specifying
     the set of column predictors to be used for each target variable in the row.
     Each entry has value 0 or 1. A value of 1 means that the column
     variable is used as a predictor for the row variable. Order of
     variables is exactly the same as in \code{data}. By default all
     variables that are earlier in the visit sequence are used as predictors.
     For the default visit sequence (\code{1:ncol(data)}) the default
     \code{predictorMatrix} will have values of 1 in the lower triangle.
     See details for more information.}
  \item{m}{number of synthetic copies of the real data to be generated.
     The default is \code{m = 1}.}
  \item{k}{a size of the synthetic data set (\code{k x p}),
     which can be smaller or greater than the size of the original data 
     set (\code{n x p}). The default is \code{nrow(data)} which means 
     that the number of individuals in the synthesised data is the same 
     as in the real data (\code{k = n}).}
  \item{proper}{a logical value with default set to \code{FALSE}. 
     If \code{TRUE} proper synthesis is conducted.}
  \item{nlevelmax}{a minimum number of values a numeric variable should
     have to be treated as numeric. Numeric variables with fewer levels
     than \code{nlevelmax} are changed into factors.}
  \item{maxfaclevels}{a maximum number of factor levels that can be handled by
     parametric methods.}
  \item{rules}{a vector or a list of rules for restricted values. Restricted
     values are those that are determined explicitly by values of other
     variables.}
  \item{rvalues}{a vector or a list of the values corresponding to the rules
     specified by \code{rules}.}
  \item{contNA}{a vector or a list of codes for missing values for continuous
     variables if different from the \code{R} missing data code \code{NA}.}
  \item{event}{a vector of integers of length \code{ncol(data)} specifying
     for survival data the column indices for corresponding event indicators.
     For non-survival data it has a value of 0. Order of variables is exactly
     the same as in \code{data}.}
  \item{smoothing}{a vector of length \code{ncol(data)} specifying smoothing 
     method (\code{"density"} or \code{""}) to be used for each variable in 
     the data. Smoothing can be applied to continuous variables only and the 
     empty method \code{""} should be set for all other variables.}
  \item{denom}{a vector of integers of length \code{ncol(data)} specifying
     for variables to be modelled using binomial regression the column indices 
     for corresponding denominator variables. For other variables it has 
     a value of 0.}
  \item{minbucket}{the minimum number of observations in any terminal node of 
     a CART model.}     
  \item{drop.not.used}{a logical value. If \code{TRUE} (default) variables not 
     used in synthesis are not saved in the synthesised data and are not 
     included in the corresponding synthesis parameters.}
  \item{drop.pred.only}{a logical value. If \code{TRUE} (default) variables not 
     synthesised and used as predictors only are not saved in the synthesised 
     data.}
  \item{defaultMethod}{a vector of four strings containing the default
     parametric synthesising methods for numerical variables, factors 
     with two levels, unordered factors with more than two levels 
     and ordered factors with more than two levels respectively. 
     They are used when \code{method} is set to \code{"parametric"} or 
     when there is an inconsistency between variable type and provided method.}
  \item{diagnostics}{a logical value. If \code{TRUE} diagnostic information are
     appended to the value of the function. If \code{FALSE} (default) only the
     synthesised data are saved.}
  \item{printFlag}{if \code{TRUE} (default) \code{syn} will print synthesising 
     history at the console. For silent computation use \code{print = FALSE}.}
  \item{\dots}{additional arguments passed to synthesising functions.}
  \item{x}{an object of class \code{synds}; a result of a call to \code{syn}.}
}

\details{
 Only variables that are in \code{visitSequence} with corresponding non-empty 
 \code{method} are synthesised. The only exceptions are event indicators. They 
 are synthesised along with the corresponding time to event variables and should
 not be included in \code{visitSequence}. All other variables (not in 
 \code{visitSequence} or in \code{visitSequence} with a corresponding blank method) 
 can be used as predictors. Including them in \code{visitSequence} generates 
 a default \code{predictorMatrix} reflecting the order of variables in
 the \code{visitSequence} otherwise \code{predictorMatrix} has to be 
 adjusted accordingly. All predictors of the variables that are not in
 \code{visitSequence} or are in \code{visitSequence} but with a blank method
 are removed from \code{predictorMatrix}.

 Variables to be synthesised that are not synthesised yet cannot be used
 as predictors. Also all variables used in passive synthesis or in restricted
 values rules (\code{rules}) have to be synthesised before the variables they
 apply to.

 Mismatch between data type and synthesising method stops execution and
 print an error message but numeric variables with number of levels less
 than \code{nlevelmax} are changed into factors and methods are changed
 automatically, if necessary, to methods for categorical variables. 
 Methods for variables not in a visit sequence will be changed into blank.

 The built-in elementary synthesising methods include:
 \describe{
   \item{ctree, cart}{classification and regression trees (CART)}
   \item{surv.ctree}{classification and regression trees (CART)
      for duration time data (parametric methods for survival data are 
      not implemented yet)}
   \item{norm}{normal linear regression}
   \item{normrank}{normal linear regression preserving the marginal distribution}
   \item{logreg}{logistic regression}
   \item{polyreg}{unordered polytomous regression}
   \item{polr}{ordered polytomous regression}
   \item{pmm}{predictive mean matching}
   \item{sample}{random sample from the observed data}
   \item{passive}{function of other synthesised data}
  }
 }
\value{
  An object of class \code{synds}, which stands for 'synthesised 
  data set'. It is a list with the following components:

  \item{call}{the original call to \code{syn}.}
  \item{m}{the number of synthetic versions of the real data.}
  \item{syn}{a data frame (for \code{m = 1}) or a list of \code{m} data frames 
    (for \code{m > 1}) with synthetic data set(s).}
  \item{method}{a vector of synthesising methods applied to each variable
     in the saved synthesised data.}
  \item{visitSequence}{a vector of column indices of the visiting sequence.
     The indices refer to the columns in the saved synthesised data.}
  \item{predictorMatrix}{a matrix specifying the set of predictors used
     for each variable in the saved synthesised data.}
  \item{event}{a vector of integers specifying for survival data the column 
     indices for corresponding event indicators. The indices refer to the 
     columns in the saved synthesised data.}
  \item{smoothing}{a vector specifying smoothing methods 
     applied to each variable in the saved synthesised data.}
  \item{denom}{a vector of integers specifying for variables modelled 
     using binomial regression the column indices for corresponding 
     denominator variables. The indices refer to the columns in the 
     saved synthesised data.}
  \item{minbucket}{a minimum number of observations in any terminal 
     node of a CART model.}
  \item{proper}{a logical value indicating whether proper synthesis was
     conducted.}
  \item{n}{a number of cases in the original data.}  
  \item{k}{a number of cases in the synthesised data.}
  \item{rules}{a list of rules for restricted values applied to 
     the synthetic data.}
  \item{rvalues}{a list of the values corresponding to the rules
     specified by \code{rules}.}
  \item{contNA}{a list of codes for missing values for continuous variables.}
  \item{drop.not.used}{a logical value indicating whether variables not used
     in synthesis are saved in the synthesised data and corresponding
     synthesis parameters.}
  \item{drop.pred.only}{a logical value indicating whether variables not
     synthesised and used as predictors only are saved in the synthesised data.}
}
%\references{...}
\note{See package vignette for additional information.}
\seealso{\code{\link{compare.synds}}, \code{\link{summary.synds}}}

\examples{
  ### selection of variables
  vars <- c("sex","age","marital","income","ls","smoke")
  rds  <- SD2011[1:2000,vars]
  
  ### default synthesis
  s1 <- syn(rds)
  s1
  
  ### synthesis with default parametric methods
  s2 <- syn(rds, method = "parametric")
  s2$method
  
  ### multiple synthesis of selected variables with customised methods
  s3 <- syn(rds, visitSequence = c(2, 1, 6, 5), m = 2,
            method = c("logreg","sample","","normrank", "ctree",""))
  summary(s3)
  summary(s3, msel = 1:2)
  
  ### adjustment to the default predictor matrix 
  s4.ini <- syn(data = rds, visitSequence = c(1, 2, 5, 3),
                m = 0, drop.not.used = FALSE)
  pM.cor <- s4.ini$predictorMatrix
  pM.cor["marital","ls"] <- 0
  s4 <- syn(data = rds, visitSequence = c(1, 2, 5, 3),
            predictorMatrix = pM.cor)
  
  ### handling missing values in continuous variables
  contNA.income <- as.list(rep(NA, ncol(rds)))
  contNA.income[[4]] <- c(NA,-8)
  s5 <- syn(rds, contNA = contNA.income)
  
  ### rules for restricted values - marital status of males under 18 should be 'single'
  rules.marital <- list("","","age < 18 & sex == 'MALE'","","","")
  rvalues.marital <- list(NA,NA,'SINGLE',NA,NA,NA)
  s6 <- syn(rds, rules = rules.marital, rvalues = rvalues.marital, method = "parametric")
  with(s6$syn, table(marital[age < 18 & sex == 'MALE']))
  ### results for default parametric synthesis without the rule  
  with(s2$syn, table(marital[age < 18 & sex == 'MALE']))
}

\keyword{datagen}
\keyword{tree}
\keyword{regression}
