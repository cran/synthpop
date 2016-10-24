##-----utility.synds-------------------------------------------------------

utility.synds <- function(object, data, method = "cart", cp = 0.0001, 
                          maxorder = 1, deviance = FALSE, null.utility = FALSE, 
                          syn.only = FALSE, all.comb = FALSE, ...){
                 
 if(any(is.na(match(method, c('cart','logit','poly','forest'))))){
   stop('Invalid method type.', call. = FALSE)}
 if(!(maxorder %in% 0:4)) stop('Invalid maximum order of interactions.\n', call. = FALSE)  
 if(is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
 if(!class(object)=="synds") stop("Object must have class 'synds'.\n", call. = FALSE)
 if(null.utility & object$m == 1) stop("Null utility measures cannot be calculated for a single synthetic data set.\n", call. = FALSE)
  
 m <- object$m
 syn.vars <- names(object$method)[!object$method==""] 
 obs.vars <- object$obs.vars
 
 ##Check if not used or used for prediction only vars removed 
 # and if they are needed for utility calculation 
 if((length(setdiff(obs.vars, names(object$method))) > 0) & 
   syn.only == FALSE) stop("If unchanged variables are to be used in computing the utility they cannot be removed during synthesis.\n Set both 'drop.not.used' and 'drop.pred.only' to FALSE.\n", call. = FALSE)

 ##Remove non-synthesized variables  
 if(syn.only == TRUE){
   rm.vars   <- setdiff(obs.vars, syn.vars)
   keep.vars <- names(data)[!(names(data) %in%  rm.vars)]
   realnames <- keep.vars
   df.obs    <- data[ , keep.vars, drop = FALSE] 
   synnames  <- syn.vars
 } else {
   realnames <- names(data)
   if(m > 1) synnames <- names(object$syn[[1]]) else synnames <- names(object$syn)
   df.obs <- data[, realnames, drop = FALSE]
 }  

 ##Check similarity of datasets
 commonnames <- synnames[match(realnames,synnames)]
 uncommonnames <- commonnames[is.na(commonnames)]
 if((length(uncommonnames) != 0) && (typeof(commonnames) == "character")) stop("Original and synthetic data do not correspond.", call. = FALSE)

 if (syn.only == TRUE && length(rm.vars)!=0) cat("Non-synthesized variables not used in computing the utility.\n")
 if (syn.only == FALSE & length(obs.vars) > length(syn.vars)) cat("NOTE: Non-synthesized variables also used in computing the utility.\n")

 ##Convert any character variables into factors
 ischar = sapply(df.obs,is.character)
 if (length(df.obs[ischar]) > 0) {
   cat("Variable(s): ",colnames(df.obs)[ischar],
     " have been changed from character to factor.\n",sep=" ")
   for (j in colnames(df.obs[ischar])){ df.obs[,j] = as.factor(df.obs[,j]) }
 }
 ##Note Factors
 facVar = sapply(df.obs,is.factor) 

 ##Check for missing data
 if(any(is.na(df.obs)) == TRUE) missing = TRUE else missing = FALSE

##-------------------------------------------------------------------------
## Matrices / lists to store the results
 if (all.comb == TRUE) {
   m0 <- 1
   n0 <- object$n + m*object$k
   dimnms0 <- "syn.comb" 
 } else {
   m0 <- m
   n0 <- object$n + object$k 
   dimnms0 <- paste0("syn=",1:m)
 }  
  
 noMeth  <- length(method)
 propStore <- array(NA, dim = c(n0,noMeth,m0), 
   dimnames = list(1:n0,method,dimnms0))
 
 utilVal <- matrix(NA, ncol = noMeth, nrow = m0, 
   dimnames = list(dimnms0, method))
  
 if (deviance == TRUE) dev <- matrix(NA, nrow = m0, ncol = 2*noMeth,
   dimnames = list(dimnms0, paste0(rep(method, each = 2), 
   rep(c(" Stat"," DF"), times = noMeth))))
##-------------------------------------------------------------------------


 for(i in 1:m0){
   if(m0 == 1) synds <- rbind.fill(object$syn) else synds <- object$syn[[i]]
	 if(syn.only == TRUE){
     df.syn <- synds[, syn.vars, drop = FALSE]
   } else {
     df.syn <- synds[, commonnames, drop = FALSE]
	 }
	 temp.comb <- rbind(df.obs, df.syn)

	 ##Flag and fill missing values
   if(missing == TRUE || method == 'forest'){
     df.comb <- fill.NA(temp.comb, 'zero')
   } else {df.comb <- temp.comb}

   ##Original or Synthetic classification variable and 
   ##ratio of total rows in synthetic data
   t <- c(rep(0,nrow(df.obs)),rep(1,nrow(df.syn)))
	 df.prop  <- data.frame(cbind(df.comb, t))
   synRatio <- nrow(df.syn)/length(t)


   ##Fit models, get propensity scores
	 ##---------------------------------
   
   ##Logit
   if (is.element("logit", method)) {
	   if (maxorder >= 1) {
	     logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
	     propLogit <- suppressWarnings(glm(logit.int, data = df.prop, family = "binomial"))
	   } else {  
       propLogit <- suppressWarnings(glm(t ~ ., data = df.prop, family = "binomial"))
	   }
     logitPhat = predict(propLogit, type = 'response')
		 utilVal[i, "logit"] = (sum((logitPhat - synRatio)^2, na.rm = T) / length(t))
     propStore[, "logit", i] <- logitPhat
     
     if(deviance == TRUE){
     	 #nullVal = -2*(length(t) * log(1 - synRatio) + nrow(df.syn) * log(synRatio / (1 - synRatio)))
		   #dev[i, "logit Stat"] = propLogit$null.deviance - propLogit$deviance
		   dev[i, "logit Stat"] = propLogit$deviance / length(t)
		   dev[i, "logit DF"] = propLogit$df.null - propLogit$df.residual
     }
   }	  
  
   ##Polynomial   
   if (is.element("poly", method)) {
     propPoly = tryCatch({polyclass(data = t, cov = df.comb, penalty = 0, ...)}, 
		   error = function(cond){
		     message("Avoided Self-Comparison")
		  	 return(NULL)})
		   if(is.null(propPoly) == FALSE){
		     polyPhat = ppolyclass(fit = propPoly, cov = df.comb)[, 2]
		   } else {polyPhat = rep(synRatio, length(t))}
		 	#### From all.comb 
		  ### propPoly = polyclass(data = t, cov = df.comb, penalty = 0, ...)
      ### polyPhat = ppolyclass(fit = propPoly, cov = df.comb)[, 2]  
       
       #propPoly = polymars(t, df.comb, classify = TRUE, factors = colnames(df.comb[facVar]), gcv = gcv)
		   #polyPhat = predict(propPoly, df.comb, classify = FALSE)[,2]
     utilVal[i, "poly"] = (sum((polyPhat - synRatio)^2, na.rm = T) / length(t))
     propStore[, "poly", i] <- polyPhat
     
     		##propPoly = suppressWarnings(glm(t ~ poly(data.matrix(df.comb), degree = degree, raw = T)))
			  ##propScore = suppressWarnings(predict(propPoly, type = 'response'))
			  ##utilVal[i] = (sum((propScore - synRatio)^2, na.rm = T) / length(t))*10000
	
     
     if(deviance == TRUE){
       #nullVal = -2*(length(t) * log(1 - synRatio) + nrow(df.syn) * log(synRatio / (1 - synRatio)))
       #dev[i, "poly Stat"] = nullVal - (-2*(sum(log(1 - polyPhat)) + 
		   #				sum(log(polyPhat[(nrow(df.obs)+1):length(t)] / (1 - polyPhat[(nrow(df.obs)+1):length(t)]))) ))
		   #dev[i, "poly Stat"] = -2*( sum( log( ( 1 - synRatio ) / ( 1 - polyPhat[ 1 : nrow(df.obs) ] ) ) ) + 
			 #sum( log( synRatio / polyPhat[ ( nrow(df.obs)+1 ) : length( t ) ] ) ) )
		   dev[i, "poly Stat"] = -2*( sum( log( 1 - polyPhat[ 1 : nrow(df.obs) ] ) ) + 
			 sum( log( polyPhat[ ( nrow(df.obs) + 1 ) : length( t ) ] ) ) ) / length(t)
		   if (is.null(propPoly) == FALSE) dev[i, "poly DF"] = propPoly$ndim else dev[i, "poly DF"] = 0
		   #### From all.comb
       ### dev[1, "poly DF"] = propPoly$ndim - 1
       #dev[i, "poly DF"] = propPoly$model.size - 1
     }
   }

   ##CART
   if (is.element("cart", method)) { 		  
     propCART = rpart(t ~ ., data = df.prop, method = 'class', control = rpart.control(cp = cp, ...))
     cartPhat = predict(propCART)[,2]
		 utilVal[i, "cart"] = (sum((cartPhat - synRatio)^2, na.rm = T) / length(t))
     propStore[, "cart", i] <- cartPhat
     
     if(deviance == TRUE){
       #nullCART = rpart(t ~ ., data = df.prop, method = 'class', control = rpart.control(cp = 1))
		   ##nullVal = -2*(length(t) * log(1 - synRatio) + nrow(df.syn) * log(synRatio / (1 - synRatio)))
       #dev[i, "cart Stat"] = nullVal - sum(residuals(propCART, type = "deviance"))
		   #dev[i, "cart Stat"] = -2*( sum( log( ( 1 - synRatio ) / ( 1 - cartPhat[ 1 : nrow(df.obs) ] ) ) ) + 
			 #sum( log( synRatio / cartPhat[ ( nrow(df.obs)+1 ) : length( t ) ] ) ) )
		 	 dev[i, "cart Stat"] = -2*( sum( log( 1 - cartPhat[ 1 : nrow(df.obs) ] ) ) + 
			 sum( log( cartPhat[ ( nrow(df.obs) + 1 ) : length( t ) ] ) ) ) / length(t)
		 	 dev[i, "cart DF"] = nrow(propCART$frame) * 2 - 1
     }
   }

   ##Random forest
   #if (is.element("forest", method)) { 
   #	 if(is.null(cforest_control) == TRUE){
	 #   cforest_control = cforest_unbiased(mtry = sqrt(ncol(df.comb)), ntree = 50)
	 #	 }
   #propForest = cforest(t ~ ., data = df.prop, ...)
   #propScore = predict(propForest)
	 #utilVal[i, "forest Stat"] = (sum((propScore - synRatio)^2, na.rm = T) / length(t))
   #}
 }

 tempMatrix <- matrix(c(apply(utilVal, 2, mean), apply(utilVal, 2, sd)), 
   byrow = FALSE, ncol = 2)
 dimnames(tempMatrix) <- list(method, c('Mean','SD'))
 
 if(null.utility == TRUE){
   null.summary <- null.compute(tempMatrix, object, method, cp, maxorder, ...)
 } else {
   null.summary <- NULL
 }
 
 if(deviance == T){ 
   tempDev <- matrix(c(apply(dev[, paste0(method, " Stat"), drop = FALSE], 2, mean), 
     apply(dev[, paste0(method, " Stat"), drop = FALSE], 2, sd)), byrow = F, ncol = 2)
   dimnames(tempDev) <- list(method, c('Mean', 'SD'))
 } else {
   tempDev <- dev <- NULL
 }
      
 output <- list("method" = method, 
   "utility.summary" = tempMatrix, "utility.raw" = utilVal, 
   "deviance.summary" = tempDev, "deviance.raw" = dev,
   "nullstats.summary" = null.summary$MSEDiffRatio,
   "nullstats.raw" = null.summary$NullUtility, 
   "ind.propensity.scores" = propStore)

 class(output) <- "utility.synds"
 return(output)

}


###-----fill.NA------------------------------------------------------------
#######################
#Date: 04/06/15
#Author: Joshua Snoke
#Purpose: Handles missing data as an internal function for the propensity
# score utility estimation. Code junks borrowed from 
# https://github.com/markmfredrickson/optmatch/blob/master/R/fill.NAs.R
# and methods follow Rosenbaum and Rubin (1984).
###fill.NAs function from package 'optmatch'

fill.NA <- function (data, fill.func) {
    
 ###Checks and setup
 if (is.matrix(data)) {
   dat = as.data.frame(data)
 } else if(is.data.frame(data)){ 
   dat = data
 } else stop("Data must be matrix or data frame")

 ###Split into numeric and factors
 ###Factors
 isfactor = sapply(dat,is.factor)
 facData = dat[isfactor]
   if(length(facData) > 0){
     for (j in (1:ncol(facData))){
       facData[,j] <- addNA(facData[,j], ifany=TRUE)
     } 
   }
 ###Numeric
 tempData = dat[!isfactor]
 if(length(tempData) > 0){
   tempX <- as.formula(~.)
   original.NAs <- sapply(tempData, function(i) any(is.na(i)))
   original.names <- colnames(tempData)[original.NAs]
   modmat <- model.matrix(tempX, model.frame(tempX, tempData, na.action = na.pass))
   modmat <- as.data.frame(modmat)
   modmat["(Intercept)"] <- NULL

   if (!any(original.NAs)) return(modmat)
   
   NA.columns <- sapply(tempData[original.names], function(column) is.na(column))
   colnames(NA.columns) <- paste(colnames(NA.columns), "NA", sep = ".")
   expanded.NAs <- colnames(modmat)[apply(modmat, 2, function(i) any(is.na(i)))]

   ###
   if(fill.func == 'mean'){
     modmat[expanded.NAs] <- sapply(modmat[expanded.NAs], fill.column.numeric, simplify = FALSE)
   } else if (fill.func == 'zero'){
     modmat[expanded.NAs] <- sapply(modmat[expanded.NAs], fill.column.zero, simplify = FALSE)
   } else {stop("Unknown fill function.")}
 }
    
 if(length(facData) > 0 && length(tempData) > 0){
   result <- cbind(facData, modmat, NA.columns)
 } else if(length(facData) == 0){
	 result <- cbind(modmat, NA.columns)
 } else if(length(tempData) == 0){
   result <- facData
 } else stop("There's no data here.")
    
 return(result)
}

###-----fill.column.numeric------------------------------------------------

fill.column.numeric <- function(column) {
  nas <- is.na(column)
  cm <- mean(column[!nas])
  column[nas] <- cm
  return(column)
}


###-----fill.column.zero---------------------------------------------------

fill.column.zero <- function(column) {
  nas <- is.na(column)
  column[nas] <- 0
  return(column)
}


###-----null.compute-------------------------------------------------------

null.compute <- function(util, syn, method, cp, maxorder, ...){

 m0 <- syn$m
  
 ##list to hold each null values
 self <- vector("list", length(method))
 names(self) <- paste(method, "null", sep = '.')
 for(i in 1:length(self)){
   self[[i]] <- matrix(NA, nrow = m0, ncol = m0, 
                       dimnames = list(paste0("syn=",1:m0), paste0("syn=",1:m0)))
 }
    
 for(j in 1:m0){
   if(any(is.na(syn$syn[[j]])) == TRUE) missing = TRUE else missing = FALSE
        
   for(ind in 1:m0){
     df1 <- rbind(syn$syn[[j]], syn$syn[[ind]])
     t   <- rep(0:1, each = nrow(syn$syn[[ind]]))
     if(missing == TRUE){
       df2 <- fill.NA(df1, 'zero')
     } else {df2 = df1}
     df3 <- cbind.data.frame(df2, t)
     synRatio <- 0.5
            
     if(is.element("logit", method)){
       if (maxorder >= 1) {
  	     logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
  	     propLogit <- suppressWarnings(glm(logit.int, data = df3, family = "binomial"))
  	   } else {  
         propLogit <- suppressWarnings(glm(t ~ ., data = df3, family = "binomial"))
  	   }
       logitPhat <- predict(propLogit, type = 'response')
       #utilVal[i, "logit"] = (sum((logitPhat - synRatio)^2, na.rm = T) / length(t))
       self$logit.null[ind, j] = (sum((logitPhat - synRatio)^2, na.rm = T) / length(t))
       #propStore[, "logit", i] = logitPhat ## should we save propensity scores again here?
     }
            
     if(is.element("poly", method)){
       propPoly <- tryCatch({polyclass(data = t, cov = df2, penalty = 0, ...)}, 
         error = function(cond){
           #message("Avoided Self-Comparison")
           return(NULL)})
     if(is.null(propPoly) == FALSE){
       polyPhat = ppolyclass(fit = propPoly, cov = df2)[, 2]
       } else {polyPhat = rep(synRatio, length(t))}
       self$poly.null[ind, j] = (sum((polyPhat - synRatio)^2, na.rm = T) / length(t))
       #propStore[, "poly", i] = polyPhat
     }
       
     if(is.element("cart", method)){
       propCART = rpart(t ~ ., data = df3, method = 'class', control = rpart.control(cp = cp, ...)) ## does this carry over?
       cartPhat = predict(propCART)[,2]
       self$cart.null[ind, j] = (sum((cartPhat - synRatio)^2, na.rm = T) / length(t))
       #propStore[, "cart", i] = cartPhat
     }
   }
 }
 
 mseDiffRat <- matrix(NA, length(method), ncol=2)
 dimnames(mseDiffRat) <- list(method, c('Mean diff','Mean ratio'))

 for(i in 1:length(method)){
   tempSelf = mean(c(self[[i]][lower.tri(self[[i]])], self[[i]][upper.tri(self[[i]])]))
   mseDiffRat[i, 1] = util[i, 1] - tempSelf
   mseDiffRat[i, 2] = util[i, 1] / tempSelf
 }
    
 return(list("MSEDiffRatio" = mseDiffRat, "NullUtility" = self))
}



###-----tab.utility--------------------------------------------------------
# utility measure for tabular data  

tab.utility <- function(object, data, vars = NULL, ngroups = 5, 
                        with.missings = TRUE, ...) 
{
 # CHECKS #---                                
 if (is.null(data)) 
   stop("Requires parameter 'data' to give name of the real data.\n", 
   call. = FALSE)
 if (!class(data) == "data.frame") 
   stop("Data must have class 'data.frame'.\n", call. = FALSE)
 if (!class(object) == "synds") 
   stop("Object must have class 'synds'.\n", call. = FALSE)

 if (is.null(vars)) stop("Need to set variables with vars parameter.\n", call. = FALSE)
 else if (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ", 
   paste(vars[!(vars %in% names(data))],collapse=", "), call. = FALSE)
 #--- 
 
 data  <- data[, vars, drop=FALSE]
 nvars <- ncol(data)  
 data.orig <- data
 
 m <- object$m
 if (m==1) syndata <- list(object$syn) else syndata <- object$syn   
 syndata <- lapply(syndata,'[',vars)
 
 dfs <- ratios  <- chis <- pvals <- nempties <- vector("numeric", m)
 tabs <- Ztab <- vector("list", m) 

 for (i in 1:m) {
   data <- data.orig
      
   # make all variables into factors
   for (j in 1:nvars) {
 	   if (is.numeric(data[,j])){
	     grpd <- group_num(data[,j], syndata[[i]][,j], ngroups)
		   data[,j] <- grpd[[1]]; syndata[[i]][,j] <- grpd[[2]]
	   } else if (is.character(data[,j])) {
	     data[,j] <- factor(data[,j])
	     syndata[[i]][,j] <- factor(syndata[[i]][,j],levels=levels(data[,j]))
	   }
     if (with.missings==TRUE & any(is.na(data[,j]))){
       data[,j] <- addNA(data[,j]) ### makes missings into part of factors 
       syndata[[i]][,j] <- addNA(syndata[[i]][,j]) ### makes missings into part of factors
     }
   }

   tabs[[i]] <- xtabs(~., data = syndata[[i]])
   if (i==1) {
     tabd <- xtabs(~., data = data)  #! can tabd be different for different m???
     totcells <- length(tabd)        #!
   }
   dfs[i] <- totcells - sum(tabd + tabs[[i]] == 0) - 1
   diffs <- tabs[[i]] - tabd
   expect <- (tabs[[i]] + tabd)/2
   tabsq  <- diffs^2/expect
   tabsq[expect==0] <- 0
   chis[i] <- sum(tabsq)
   Ztab[[i]] <- diffs/sqrt(expect)
   # pvals[i] <- 1 - pchisq(chis[i],dfs[i])
   ratios[i] <- chis[i] / dfs[i]  
   nempties[i] <- sum(tabd + tabs[[i]] == 0)
 }  
  
 if (m==1) {
   tabs <- tabs[[1]]
   Ztab <- Ztab[[1]]
 } 
   
 res <- list(Chisq = chis,
   df = dfs,
   ratio = ratios,
   # pval = pvals,
   nempty = nempties,
   tab.obs = tabd,
   tab.syn = tabs,
   tab.Zdiff = Ztab)

 class(res) <- "tab.utility"
 return(res)

}


###-----group_num----------------------------------------------------------

# function to make groups of any continuous variables

group_num <- function(x1, x2, n = 5, cont.na = NULL){
#
# makes 2 cont variable into a factor of n groups with same groupings
# determined by first one
# groups will often be of uneven size because of missing values
# there may even be fewer than n groups if there are big sets of repeats
#
#
 if (length(cont.na)==0) cc <- NULL
 if (!is.numeric(x1)) stop ("x must be numeric\n")
 xn <- x1[!(x1 %in% cont.na | is.na(x1))]
 #print(table(x1 %in% cont.na) )
 #print(length(xn))
 xr <- rank(xn)
 ptiles <- round(c(length(xr)*(1:(n-1))/n,max(xr)))
 #ptiles
 ptiles <- xn[order(xn)][ptiles]
  ptiles <- ptiles[!duplicated(ptiles)]##  to allow for repeated values
 xnew <- rep(1,length(xn))
 for (i in 2:length(ptiles)) xnew[xn>ptiles[i-1] & xn<=ptiles[i]] <- i
 xnew <- factor(xnew)
 res1 <- x1
 res1[!(x1 %in% cont.na | is.na(x1))]<-xnew
# cat("\npercentiles ",ptiles,"\n")############################################################################

# now next one
xn<-x2[!(x2 %in% cont.na | is.na(x2))  ]
xnew=rep(1,length(xn))
for (i in 2:length(ptiles)) xnew[xn>ptiles[i-1] & xn<=ptiles[i]]<-i
xnew<-factor(xnew)
res2<-x2
res2[!(x2 %in% cont.na | is.na(x2))]<-xnew

# now make into factors including missing and cont.na data
nlev<-length(table(c(res1[!(x1 %in% cont.na | is.na(x1))],res2[!(x2 %in% cont.na | is.na(x2))])))
labels=paste("Gp",1:nlev)
#print(labels)
if (any(is.na(c(res1,res2)))) {
    res1[is.na(res1)]<-nlev+1
    res2[is.na(res2)]<-nlev+1
	labels<-c(labels,"missing")
}
if (!is.null(cont.na)) {
    exlevs<-cont.na[!is.na(cont.na)]
    for ( i in 1:length(exlevs)){
      res1[res1==exlevs[i]]<-nlev+i
      res2[res2==exlevs[i]]<-nlev+i
    }
	labels<-c(labels,as.character(exlevs))
	#print(labels)
	#print(table(res1))
}
#print(labels)
res1<-factor(res1,labels=labels)
res2<-factor(res2,labels=labels)
#print(table(res1))
for(i in 1:length(labels)) labels[i]<-paste(range(c(x1[res1==labels[i]],x2[res2==labels[i]])),collapse="-")
#print(labels)
res1<-factor(res1,labels=labels)
res2<-factor(res2,labels=labels)

list(res1,res2)
}          



