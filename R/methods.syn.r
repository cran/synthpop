###-----print.synds--------------------------------------------------------

print.synds <- function(x, ...){
  cat("Call:\n($call) ")
  print(x$call)
  cat("\nNumber of synthesised data sets: \n($m) ",x$m,"\n")  
  if (x$m==1){
    cat("\nFirst rows of synthesised data set: \n($syn)\n")
    print(head(x$syn))
  } else {
    cat("\nFirst rows of first synthesised data set: \n($syn)\n")
    print(head(x$syn[[1]]))
  }    
  cat("...\n")
  cat("\nSynthesising methods: \n($method)\n")
  print(x$method)
  cat("\nOrder of synthesis: \n($visit.sequence)\n")
  print(x$visit.sequence)
  cat("\nMatrix of predictors: \n($predictor.matrix)\n")
  print(x$predictor.matrix)     
  invisible(x)
}


###-----summary.synds------------------------------------------------------

summary.synds <- function(object, msel=1, ...){
  if (!all(msel %in% (1:object$m))) stop("Invalid synthesis number(s)")
  sy <- list(m=object$m,msel=msel,method=object$method)
  if (object$m==1){
    sy$result <- summary(object$syn,...)
  } else if (length(msel)==1){
    sy$result <- summary(object$syn[[msel[1]]],...)
  } else {
    for (i in (1:length(msel))){
      sy$result[[i]] <- summary(object$syn[[msel[i]]],...)
    }
  }
  class(sy) <- "summary.synds"
  sy
}


###-----print.summary.synds------------------------------------------------

print.summary.synds <- function(x, ...){
 if (x$m==1){
   cat("Synthetic object with one synethesis using methods:\n")
   print(x$method)
   cat("\n")
   print(x$result)
 } else if (length(x$msel)==1){
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   cat("\nSummary for synthetic data set ",x$msel[1],":\n",sep="")
   print(x$result)
 } else {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   for (i in (1:length(x$msel))){
     cat("\nSummary for synthetic data set ",x$msel[i],":\n",sep="")
     print(x$result[[i]])
   }
 }
 invisible(x)
}


###-----print.fit.synds----------------------------------------------------

print.fit.synds <- function(x, msel=1, ...)
{
  if (!all(msel %in% (1:x$m))) stop("Invalid synthesis number(s)")
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  for(i in msel) {          
    cat("\nsyn =",i,"\n")
    print(x$analyses[[i]]$coefficients)
  }
  invisible(x)
}


###-----lm.synds-----------------------------------------------------------

lm.synds <- function(formula, object, ...)
{
 if (!class(object)=="synds") stop("Object must have class synds\n")
 call <- match.call()
 fitting.function <- "lm"
 analyses <- as.list(1:object$m)

 # do the repated analysis, store the result without data
 if (object$m==1) {
   analyses[[1]] <- summary(lm(formula,data=object$syn,...))
 } else {
   for(i in 1:object$m) {
     analyses[[i]] <- summary(lm(formula,data=object$syn[[i]],...))
   }
 }
 # return the complete data analyses as a list of length m
 object <- list(call=call, proper=object$proper, m=object$m, 
                analyses=analyses, fitting.function=fitting.function,
                n=object$n, k=object$k)
 class(object) <- "fit.synds"
 return(object)
}


###-----glm.synds----------------------------------------------------------

glm.synds <- function(formula, family="binomial", object, ...)
{
 if (!class(object)=="synds") stop("Object must have class synds\n")
 call <- match.call()
 fitting.function <- "glm"
 analyses <- as.list(1:object$m)
 
 # do the repated analysis, store the result without data
 if (object$m==1) {
   analyses[[1]] <- summary(glm(formula,data=object$syn,family=family,...))
 } else {
   for(i in 1:object$m) {
     analyses[[i]] <- summary(glm(formula,data=object$syn[[i]],family=family,...))
   }
 }
 # return the complete data analyses as a list of length m
 object <- list(call=call, proper=object$proper, m=object$m,
                analyses=analyses, fitting.function=fitting.function,
                n=object$n, k=object$k)
 class(object) <- "fit.synds"
 return(object)
}


###-----summary.fit.synds--------------------------------------------------

summary.fit.synds <- function(object, population.inference = FALSE, ...)
{ # df.residual changed to df[2] because didn't work for lm - check if that's ok
  if (!class(object) == "fit.synds") stop("Object must have class fit.synds\n")
  m <- object$m
  k <- object$k
  n <- object$n
  if (m == 1) {
    coefficients  <- object$analyses[[1]]$coefficients[,1]
    vars          <- object$analyses[[1]]$coefficients[,2]^2
  } else {
    namesbyfit <- lapply(lapply(object$analyses,coefficients),rownames)
    allnames <- Reduce(union,namesbyfit)
    matcoef <- matvar <- matrix(NA,m,length(allnames))
    dimnames(matcoef)[[2]] <- dimnames(matvar)[[2]] <- allnames
    for (i in 1:m){
      pos <- match(namesbyfit[[i]],allnames)
      matcoef[i,pos] <- object$analyses[[i]]$coefficients[,1]
      matvar [i,pos] <- object$analyses[[i]]$coefficients[,2]^2
    }
    coefficients <- apply(matcoef,2,mean,na.rm = TRUE)
    vars <- apply(matvar,2,mean,na.rm = TRUE)
    #bm <- apply(matcoef,2,var) not needed xpt for partial synthesis
  }

  if (population.inference == F){ ## inf to Q hat

    if (object$proper == F){
      ## simple synthesis
      result <- cbind(coefficients,
                      sqrt(vars/m),
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    } else {
      ## proper synthesis
      result <- cbind(coefficients,
                      sqrt(vars*(1+k/n)/m), 
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      sqrt((1 + k/n + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    }

  } else { ## pop inference to Q

    if (object$proper == F){
      ## simple synthesis
      Tf <- vars*(k/n+1/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
    else {
      ## proper synthesis
      Tf <- vars*(k/n+(1+k/n)/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
  }
  res <- list(call = object$call, proper = object$proper,
              population.inference = population.inference,
              fitting.function = object$fitting.function,
              m = m, coefficients = result, n = n, k = k)
  class(res) <- "summary.fit.synds"
  return(res)
}


###-----print.summary.fit.synds--------------------------------------------

print.summary.fit.synds <- function(x, ...) {
  if (x$m==1) {
    cat("\nFit to synthetic data set with a single synthesis.\n")
  } else {
    cat("\nFit to synthetic data set with ",x$m," syntheses.\n",sep="")
  }
  if (x$population.inference) {
    cat("Inference to population coefficients.\n")
  } else {
    cat("Inference to coefficients and standard errors\nthat would be obtained from the observed data.\n")
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nCombined estimates:\n")
  if (x$population.inference){
    print(x$coefficients[,c("B.syn","se(B.syn)")])
  } else {
    print(x$coefficients[,c("B.syn","se(Beta).syn")])
  }      
  invisible(x)
}


###-----print.compare.fit.synds--------------------------------------------

print.compare.fit.synds <- function(x, ...){
  cat("\nCall used to fit models to the synthetised data set(s):\n")
  print(x$fit.synds.call)
  cat("\nEstimates for the observed data set:\n")
  print(x$coef.obs)
  cat("\nCombined estimates for the synthetised data set(s):\n")
  print(x$coef.syn)
  print(x$ci.plot)
  invisible(x)
}


###-----print.compare.synds------------------------------------------------

print.compare.synds <- function(x, ...){
  cat("\nComparing percentages observed with synthetic.\n")
  if(x$plot.na) cat("For numeric variables missing data categories are presented seperately.\n")
  cat("\n")
  print(x$freq.table)
  print(x$p)
  invisible(x)
}
