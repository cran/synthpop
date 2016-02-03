\name{utility.synds}
\alias{utility.synds}
\alias{print.utility.synds}
\title{
 Distributional comparison of synthesised and observed data
}
\description{
 Distributional comparison of synthesised data set with the original (observed) 
 data set using propensity scores.
}
\usage{
utility.synds(object, data, method = "cart", cp = 0.0001, maxorder = 1,
  deviance = FALSE, null.utility = FALSE, syn.only = FALSE, all.comb = FALSE, ...)

\method{print}{utility.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s).}
  \item{data}{the original (observed) data set.}
  \item{method}{a single string or a vector of strings specifying the 
    method(s) for modeling the propensity scores.  Methods can be selected 
    from \code{"cart"}, \code{"logit"} and \code{"poly"}. 
    See details for more.} %\item{cforest_control}{optional control parameters when modeling propensity scores with a random forest.} \item{cart_control}{optional control parameters when modeling propensity scores with a classification tree.}
  \item{cp}{complexity parameter for classification with method \code{"cart"}.}
  \item{maxorder}{maximum order of interactions to be considered in 
    \code{"logit"} method. For model without interactions \code{0} should be
    provided.}  
  \item{deviance}{a logical value indicating whether deviance statistic for 
    propensity score model should be calculated. This statistic gives 
    similar inference as the default utility value based on the variance.}
  \item{null.utility}{a logical value indicating whether the mean squared error 
    (MSE) of the propensity score value for the methods chosen according to the 
    null distribution should be estimated. See details for more.}
  \item{syn.only}{a logical value indicating whether any non-synthesized 
    variables present in the released data should be included in the estimation 
    of the propensity scores. Default is FALSE which means that non-synthesized 
    variables are included.}
  \item{all.comb}{alternative method for estimating propensity scores where all 
    synthetic data from each \code{m} imputation is combined with the original 
    data. Only one utility value is returned rather than one for each synthetic 
    data set.}  
  \item{\dots}{additional parameters passed to \code{\link{glm}}, 
    \code{\link{rpart}}, or \code{\link{polyclass}}.} %or \link{cforest}}
  \item{x}{an object of class \code{utility.synds.}}
}
\details{
  This function follows the method for evaluating the utility of masked data 
  as given in Snoke et al. (forthcoming) and originally proposed by Woo et al. 
  (2009). Propensity scores, as detailed in Rosenbaum and Rubin (1983), are 
  estimated for probability of membership in the synthetic data set. The value 
  returned then is the mean squared difference between these probabilities and 
  the probability of indistinguishability (perctentage of combined records in 
  the synthetic data set).

  Propensity scores can be modeled in a variey of ways. Commonly a simple 
  logistic regression is used with all variables in the data set as predictors, 
  implemented here as method \code{"logit"}. Alternative modeling options 
  available here are classification and regression trees as method \code{"cart"} 
  (default) and multivariate adaptive regression splines using the 
  \code{"polspline"} packing under function \code{"polyclass"} 
  as method \code{"poly"} (in testing). %Modeling by random forest with method "forest" is also possible (currently in testing).

  If missing values exist, indicator varibales are added and included in the 
  model as recommended by Rosenbaum and Rubin (1984). The missing cell is then 
  either imputed as zero or the mean. For cateogrical variables, \code{NA} is 
  treated as a new category. Thanks to \url{https://github.com/markmfredrickson/optmatch/blob/master/R/fill.NAs.R} 
  for useful code chunks on flagging and filling.

  Null propensity score MSEs can be estimated and returned for each selected 
  method. These values are estimated as proposed in Snoke et. al. (forthcoming) 
  and give an estimate of the expected values under a correct synthesis model. 
  These are used to produce a difference or ratio between the observed and null 
  MSEs to improve accuracy of the utility estimate.
}
\value{
  An object of class \code{utility.synds} which is a list include the mean and 
  sd of the propensity score utility values (if \code{m > 1}) and also the raw 
  utility values for each synthetic set.

  The list also contains the chosen method for modeling the propensity scores. 
  If multiple methods were provided, the list will contain a matrix of utility 
  value results for each of the different methods.

  If \code{null.utility} is set to \code{TRUE}, two more elements are given in 
  the output list. \code{nullstats.summary} provides both the differences and 
  ratios between observed and null mean utility statistics (MSEs) for each 
  method. \code{nullstats.raw} gives utility statistics from each pair of 
  synthetic data sets compared.
  
  If \code{"deviance"} is set to \code{TRUE}, the function also returns the 
  deviance statistic. Calculated as -2[log likelihood (estimated model)] / N, 
  where the estimated model uses the predicted propensity scores and the null 
  model uses the true proportion of records in the synthetic data (i.e. no
  distinguishability between original and synthetic data).
}
\references{
  Woo, M-J., Reiter, J.P., Oganian, A., and Karr, A.F. (2009). Global Measures 
  of Data Utility for Microdata Masked for Disclosure Limitation. 
  \emph{Journal of Privacy and Confidentiality}, \bold{1}(1), 111-124.
  
  Rosenbaum, P.R. and Rubin, D.B. (1984). Reducing Bias in Observational Studies 
  Using Subclassification on the Propensity Score.
  \emph{Journal of the American Statistical Association},
  \bold{79}(387), 516-524.
  
  Snoke, J., Nowok, B., Raab, G.M., Dibben, C. and Slavkovic, A. (forthcoming). 
  General and specific utility measures for synthetic data. 
}
\examples{
  ods <- SD2011[1:1000, c("age", "bmi", "depress", "alcabuse", "englang")]
  s1 <- syn(ods, m = 5)
  utility.synds(s1, ods)
}
