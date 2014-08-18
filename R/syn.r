# Source functions for synthpop library to create synthetic data
# for the SYLLS project.

# Structure and some functions based on code from MICE package
# by S. van Buuren and K. Groothuis-Oudshoorn

syn <- function(data,
                method = vector("character",length=ncol(data)),
                visitSequence = (1:ncol(data)),
                predictorMatrix = NULL,  
                m = 1, 
                k = nrow(data),
                proper = FALSE, 
                nlevelmax = 5,
                maxfaclevels = 60,
                rules = as.list(rep("",ncol(data))),                                    
                rvalues = as.list(rep(NA,ncol(data))),                                  
                contNA = as.list(rep(NA,ncol(data))),
                event = rep(0,ncol(data)),
                smoothing = rep("",ncol(data)), 
                denom = rep(0,ncol(data)),
                minbucket = 5,    
                drop.not.used = TRUE, 
                drop.pred.only = TRUE,                                            
                defaultMethod = c("normrank","logreg","polyreg","polr"),
                diagnostics = FALSE,
                printFlag = TRUE,
                ...
                )
{

#----------------------------------------------------------------------
# the code for the following checking functions is included within syn
# so as to allow them to access the local objects
#----------------------------------------------------------------------
# set default method for everything to ctree and to blank (which will get defaults) if method is "parametric"
if (all(method=="")) method="ctree"
#else if (length(method)==1 && method=="parametric") method=rep("",dim(data)[2])
# change vectors to lists
if (!is.list(rules))   rules <- as.list(rules)
if (!is.list(contNA))  contNA <- as.list(contNA)
if (!is.list(rvalues)) rvalues <- as.list(rvalues)

##-----------------------check.visitSequence.syn-----------------------

check.visitSequence.syn <- function(setup){

 vis      <- setup$visitSequence
 nvar     <- setup$nvar
 varnames <- setup$varnames
 method   <- setup$method

 # visitSequence can include column indices only
 # not all columns have to be given - variables
 # not in the visit sequence won't be synthesised
 
 # stop if variable in visitsequence more than once
 if (any(duplicated(vis)))  stop("\n Visit sequence includes repeated variable numbers\n")

 # remove any visitSequnce members outside 1:nvar
 if (any(!(vis %in% 1:nvar))) {
   cat("Element(s): {",paste(vis[!(vis %in% 1:nvar)],
       collapse=", "),"} of the 'visitSequence' removed as not valid. No such column.\n\n",sep="")
   vis <- as.numeric(vis[vis %in% 1:nvar])
 }

 # remove event indicator(s) from visitSequence, if present
 event.in.vis <- !is.na(match(vis,event))
 if (!is.null(event) & any(event.in.vis)) {
   cat("Element(s) {",paste(vis[event.in.vis],collapse=", "),
       "} of the 'visitSequence' with method(s) {",
       paste(method[vis[event.in.vis]],collapse=", "),
       "} removed because used as event indicator(s).\nAny event indicators will be synthesised along with the corresponding survival time(s). \n\n")
   vis <- vis[!event.in.vis]
   if (length(vis)<2) stop("Visit sequence now of length ",
       length(vis),". No synthesis done.")
 } 
                                                            #GRdenom new code
 #!BN adjusted to allow visit sequence with selected vars only 
 #! have to add a condition when denominator is not in visit seq at all;
 #! sampler has to be changed still
 #!---                                                                   
 #  check that denominator comes before the count for a binomial with denom
 #if (any(denom>0)) {
 #   denom.in.vis<-(1:nvar)[denom>0]
 #       for (j in denom.in.vis){
 #          posj<-(1:length(vis))[match(j,vis)]
 #          kj <-denom[j]
 #          posk<-(1:length(vis))[match(kj,vis)]
 #      if (posj<=posk) 
 #         stop("\n Denominator ",varnames[j]," for ",varnames[kj]," must be synthesisied before it\n")
 #   }
 #
 # }                                                               
 
 # check that denominator comes before the count for a binomial with denom 
 for (j in which(denom[vis]>0)){ 
   denom.pos <- match(denom[vis][j],vis)
   if (j < denom.pos) stop("Denominator ",varnames[denom[vis][j]]," for ",
                           varnames[vis[j]],"must be synthesisied before it.\n")
 }
 #!---
                                                                        
 # stop if visitSequence is empty
 if (length(vis)==0) stop(paste("Seems no variables being synthesised.\nCheck parameter 'visitSequence'."))

 # All variables used in passive synthesis have to be synthesised BEFORE
 # the variables they apply to
 for (j in which(is.passive(method[vis]))){  #  
   var.present <- match(all.vars(as.formula(method[vis][j])),varnames) 
   var.in.vis  <- match(var.present,vis)
   if(j < max(var.in.vis) | any(is.na(var.in.vis))) stop("Variable(s) used in passive synthesis for ",
     varnames[vis][j]," has/have to be synthesised BEFORE the variables they apply to.")
 }

 setup$visitSequence <- vis
 return(setup)
}
##-----------------end of--check.visitSequence.syn---------------------


##-----------------------check.predictorMatrix.syn---------------------

check.predictorMatrix.syn <- function(setup){
 ## checks the predictorMatrix
 ## makes consistency edits of the predictorMatrix

 pred     <- setup$predictorMatrix
 nvar     <- setup$nvar
 varnames <- setup$varnames
 vis      <- setup$visitSequence
 method   <- setup$method
 denom     <-setup$denom                     #GRdenom new
    
 # set up default predictor matrix (if not provided by a user)
 # to lower diagonal in order of visitSequnce but with
 # elements for variables not to be synthesised set to 0
 
 pred.dt          <- matrix(0,nvar,nvar)
 pred.dt[vis,vis] <- outer(1:length(vis),1:length(vis),">")
 if(is.null(pred)) pred <- pred.dt

 # basic corrections for a default matrix or the one provided by a user
 dimnames(pred)   <- list(varnames, varnames)
 diag(pred)       <- 0

 # select from visitSequence variables that are synthesised
 # (=method different than "")
 vis.syn <- vis
 if (!all(method=="") & length(method)>1) vis.syn <- intersect(vis,which(method!=""))
 # removing predictors for variables with "" method
 if (length(vis.syn) < length(vis)){
   vis.blank        <- setdiff(vis,vis.syn)
   pred[vis.blank,] <- 0
 }
 # removing predictors for variables not in visitSequence
 pred[setdiff(1:nvar,vis),] <- 0
 
 # removing predictors for variables with "sample" method
 for (j in which(method=="sample")) pred[j,] <- 0
 
 # removing survival time from predictors
 for (j in which(method=="surv.ctree")) pred[,j] <- 0

 # removing event indicator from predictors
 for (j in which(method=="surv.ctree" & event>0)) pred[,event[j]] <- 0
                                                                     #GRdenom new lines
 #  remove denom from prediction of its numerator
    for  (j in which(method=="logreg"))  {
       if (denom[j]>0) pred[j,denom[j]]<-0
    }
                                                                    # to here
 # checking consistency between visitSequence and predictor matrix
 # provided by a user: dropping from predictors variables that are
 # not already synthesised

 preddel <- which((pred[,vis.syn]==1 & pred.dt[,vis.syn]==0),arr.ind=TRUE)

if (length(vis)>1) {
 	pred[,vis.syn]  <- ifelse((pred[,vis.syn]==1 & pred.dt[,vis.syn]==0),
                           0,pred[,vis.syn])
 	if(nrow(preddel)>0) cat(paste("Not synthesised predictor ",
                         varnames[vis.syn][preddel[,2]],
                         " removed from predictorMatrix for variable ",
                         varnames[preddel[,1]],".\n",sep=""))

}
setup$predictorMatrix <- pred
 setup$visitSequence   <- vis
 return(setup)
}
##-----------------end of--check.predictorMatrix.syn-------------------


##-----------------------check.method.syn------------------------------

check.method.syn <- function(setup, data, proper) {
 ## check method, set defaults if appropriate

 method        <- setup$method
 defaultMethod <- setup$defaultMethod
 vis           <- setup$visitSequence
 nvar          <- setup$nvar
 varnames	     <- setup$varnames
 pred		       <- setup$predictorMatrix
 event         <- setup$event
 denom         <- setup$denom                    #GRdenom new

 # check if var has predictors
 if (sum(pred)>0) has.pred <- apply(pred!=0,1,any)   # GR condition added
 else has.pred<-rep(0,nvar)

 if(any(method=="parametric")){ # the default argument
   # set method for first in visitsequence to "sample"
   # method[vis[1]] <- "sample"  
   # change to default methods for variables with predictors

   if (length(vis)>1){
	   for(j in vis[-1]){
	     if (has.pred[j]){
	       y <- data[,j]
	       if(is.numeric(y))        method[j] <- defaultMethod[1]
	       else if(nlevels(y) == 2) method[j] <- defaultMethod[2]
	       else if(is.ordered(y) & nlevels(y) > 2) method[j] <- defaultMethod[4]
	       else if(nlevels(y) > 2)  method[j] <- defaultMethod[3]
	       else if(is.logical(y))   method[j] <- defaultMethod[2]
	       else if(nlevels(y)!=1) stop("Variable ",j," ",varnames[j],
		 " type not numeric or factor.") # to prevent a constant values failing
	     }
	   }
	}
 }

 # check method for variables without predictors
 # set to "sample" if method is not valid
 for (j in vis) {
   if (!has.pred[j] & is.na(any(match(method[j],
      c("","sample","sample.proper",
      "norm","norm.proper","logreg","logreg.proper")))))
   {
     cat('Method "',method[j],
     '" is not valid for a variable without predictors (',
     names(data)[j],').\nMethod has been changed to "sample".\n\n',sep="")
     method[j] <- "sample"
   }
 }
 
 # check whether the elementary syhthesising methods are available
 # on the search path
 active    <- !is.passive(method) & !(method=="")
 fullNames <- paste("syn", method[active], sep=".")
 notFound  <- !sapply(fullNames, exists, mode="function", inherit=TRUE)
 if (any(notFound)) stop(paste("The following functions were not found:",
                         paste(fullNames[notFound],collapse=", ")))

 # type checks on built-in  methods 

for(j in vis) {
 	 y     <- data[,j]
   vname <- colnames(data)[j]
   mj    <- method[j]
 mlist <- list(m1 = c("logreg","polyreg","polr"),
               m2 = c("norm","normrank","surv.ctree"),
               m3 = c("norm","normrank","surv.ctree","logreg"))
   # In case of type mismatch stop execution

#                                                       #GRdenom lines changed
# check for logistic with denominator
# 
   if (denom[j]>0) {
   if (!(mj %in% c("logreg"))) stop("Variable ", vname," has denominator ",colnames(data[denom[j]]),
                               " but method ",mj," should be logreg\n")
#
#  check all integers
#
  if (!( (is.integer(y) | all((y-round(y))==0)) & 
            (is.integer(data[denom[j]]) | all((data[denom[j]]-round(data[denom[j]])==0))) )) 
     stop("Variable ", vname," and denominator ",colnames(data[denom[j]]),
         " must be integers\n")
  if (any((data[denom[j]]-y)< 0)) stop("Variable ", vname," must be less than or equal denominator ",colnames(data[denom[j]]),
                            "\n")
   }
   else{
      if (is.numeric(y) & (mj %in% mlist$m1)){
   	 stop('Type mismatch for variable ', vname,
      		'.\n  Syhthesis method "', mj, '" is for categorical data.',sep="")
      }
      else if (is.factor(y) & nlevels(y) == 2 & (mj %in% mlist$m2)){
        stop('Type mismatch for variable ', vname,
        			'.\n  Syhthesis method "', mj, '" is not for factors.')
      }
      else if (is.factor(y) & nlevels(y) > 2 & (mj %in% mlist$m3)){
        stop('Type mismatch for variable ', vname,
        			'.\n  Syhthesis method "', mj,
        			'" is not for factors with three or more levels.')
      }                                               # to here
   }

 }

 # remove constant variables but leave passive variables untouched
 # factors with missing data won't count as NA is made into an additional level
 for(j in 1:nvar) {
 #
 #  GR change    #NA replaced with length(table(data[,j]))-1
 #
   if (!is.passive(method[j])){
     v <- ifelse(is.character(data[,j]),length(table(data[,j]))-1,var(data[,j],na.rm=TRUE))
     if (!is.na(v)) constant <- (v < 1000 * .Machine$double.eps) else
     constant <- is.na(v) | v < 1000 * .Machine$double.eps
     if (constant) {
       if (any(vis==j) & any(pred[,j]!=0)) cat("Variable ",varnames[j],
         " removed from visit sequence and as predictor because only one value.\n\n",sep="")
       else if (any(vis==j)) cat("Variable ",varnames[j],
         " removed from visit sequence  because only one value.\n\n",sep="")
       else if (any(pred[,j]!=0)) cat("Variable ",varnames[j],
         " removed as predictor because only one value.\n\n",sep="")
     pred[,j]  <- 0
     pred[j,]  <- 0
     method[j] <- ""
     }
   }
 }
 ############   this bit moved after constants out
 # check method for variables without predictors
 # set to "sample" if method is not valid
 
 # check if var has predictors  re compute it
 if (sum(pred)>0) has.pred <- apply(pred!=0,1,any)  #  GR condition added
 else has.pred<-rep( 0,sqrt( length(pred) ) )           # this needed in case pred now has dimension 1
 
 for (j in vis) {
    if (!has.pred[j] & is.na(any(match(method[j],
                                       c("","sample","sample.proper",
                                         "norm","norm.proper","logreg","logreg.proper")))))
    {
       cat('Method "',method[j],
           '" is not valid for a variable without predictors (',
           names(data)[j],').\nMethod has been changed to "sample".\n\n',sep="")
       method[j] <- "sample"
    }
 }
 # chck survival method	and events are consistent
 if (any(method=="surv.ctree")) {
   for(j in vis){   # checks for survival variables
     y     <- data[,j]
     vname <- colnames(data)[j]
     mj    <- method[j]
     if (mj=="surv.ctree") {
   	   if (!is.numeric(y)) stop("Variable ",vname,
         " should be a numeric survival time.")
    	 if (any(y < 0)) stop("Variable ",vname,
         " should be a non-negative survival time.")

       if (is.na((match(event[j],1:nvar)))) {
         cat("Variable ",vname," is a survival time. Corresponding event not in data, assuming no censoring.\n\n",sep="")
         event[j] <- -1      # used to indicate no censoring
       }
    	 else {
         tabEI <- table(data[,event[j]])
         if (length(tabEI)!=2) {
           if (length(tabEI)==1 & all(tabEI==1)) cat("Variable ",vname,
             " is a survival time with all cases having events.\n",sep="")
           else if (length(tabEI)==1 & all(tabEI==0)) stop("Variable ",
             vname," is a survival time with no cases having events.\n",
             "Estimation not possible.",sep="")
           else stop("Event must be binary 0/1 but has values {", paste(
                names(tabEI),collapse=", "),"}.\nNo data synthesised.")
         }
         if (!all(as.character(names(tabEI))==c("0","1")) &&
             !is.logical(data[,event[j]])){
           stop("Event must be binary 0/1 but it is a ",class(data[,event[j]]),
                " variable with values {", paste(names(tabEI)[1:2],
                collapse=", "),"}.",sep="")
         }
       }
     }
     else { #checks for non-survival variables
       if (event[j]!=0){
         cat("Event for variable ",vname," set to ",event[j],
 	           ' although method is "',mj,'". Event reset to 0.\n',sep="")
         event[j]<-0
       }
     }# end checking events for non  survival
   }# end of j loop
 }# end of checking when some surv variables
 else if (!all(event==0)){
   cat("No variables have a survival method so event vector which was \n{",
   paste(event,collapse=","),"} set to 0s.\n\n",sep="")
   event <- rep(0,nvar)
 }

 # change names for proper imputations and check
 #for(j in unique(vis)){
 #  if(proper==T & method[j]!="") method[j] <- paste(method[j],
 #                                                   ".proper",sep="")
 #}

 # check collinearity of variables
 if (sum(pred>0)) {                                     # GR added condition
	 inpred <- apply(pred!=0,1,any) | apply(pred!=0,2,any)
	 if (any(inpred)) collmx <- find.collinear(data[,inpred,drop=FALSE]) else
	   collmx <- NULL
	 if (length(collmx)>0) stop("Variables ", paste(collmx,collapse=", "),
	   " are pairwise collinear and cannot both be used in synthesis.")
 }
 setup$event           <- event
 setup$method          <- method
 setup$predictorMatrix <- pred
 setup$visitSequence   <- vis
 setup$denom           <-denom                  #GRdenom new
 
 return(setup)
}
##--------------------end of--check.method.syn-------------------------
 

##------------------check.rules.syn------------------------------------

check.rules.syn <- function(setup, data) {

 rules      <- setup$rules
 rvalues    <- setup$rvalues
 pred       <- setup$predictorMatrix
 nvar       <- setup$nvar
 varnames   <- setup$varnames
 method     <- setup$method
 vis        <- setup$visitSequence
  
 # Check the syntax
 #------------------
 # check the length of the rules and corresponding values
 if (length(rules)!=nvar | length(rvalues)!=nvar)
   stop(paste('Data rules (',length(rules),') and corresponding values (',
     length(rvalues),') have to be provided for each column in the data (',
     nvar,').\n  Set a rule to "" and a corresponding value to NA if no rules are to be applied for a variable.',sep=""))
 if (any(sapply(rules,length)!=sapply(rvalues,length)))
   stop("The number of data rules for each variable should equal the number of corresponding values.\n  Check variable(s): ",
     paste(varnames[sapply(rules,length)!=sapply(rvalues,length)],collapse=", "),".")

 # special characters 
 char.allowed <- c("|","||","&","&&","==",">=","<=","<",">","!=","","==-",">=-","<=-","<-",">-","!=-","'","=='",".",")","(",";","-") #### . ( and ) added
 char.present <- paste(gsub("\\w"," ",unlist(rules)),collapse=" ") # remove word characters and concatenate
 char.present <- strsplit(char.present,"[[:space:]]+")[[1]]    # split into seperate characters
 char.wrong   <- !(char.present %in% char.allowed)             # identify unxepected characters
 #if (any(char.wrong)) stop("Unexpected character(s) in rules: ",paste(char.present[char.wrong],collapse=" "),".")

 # variables names (=a string before a special character) must be in varnames 
 rule.sep <- lapply(sapply(rules,strsplit,"[|&]"),unlist)       # split into seperate conditions
 get.vars <- lapply(rule.sep,function(x) gsub("[<>=!].*","",x)) # remove evrything after a special character
 get.vars <- lapply(get.vars,function(x) gsub(" ","",x))        # remove spaces
 get.vars <- lapply(get.vars,function(x) gsub("[\\(\\)]","",x)) # remove brackets      #################  line added
 get.vars <- lapply(get.vars,function(x) gsub("is.na","",x))    # remove function name #####################  line added
 get.vars <- lapply(get.vars,function(x) x[x!=""])              # remove empty strings  ?? why this
 
 vars.in.rules <- unique(unlist(get.vars))
 vars.wrong <- !(vars.in.rules %in% varnames)                   # identify unxepected variables
 #if (any(vars.wrong)) stop("Unexpected variable(s) in rules: ",paste(vars.in.rules[vars.wrong],collapse=" "),".")
 
 if (any(char.wrong) | any(vars.wrong)) {
   cat("One of rules may not be correct. Please check and correct if necessary.") 
   rs <- unlist(rules); names(rs) <- varnames
   rs <- cbind(rs[rs!=""]); colnames(rs)<-""
   cat("\nYour rules are:")
   print(rs); cat("\n")
 }


 # Check that missingness in the data obeys the rules in rules
 nonmissing    <- vector("list",nvar)
 isfactor      <- sapply(data,is.factor)
 yes.rules <- sapply(rules,function(x) any(x!=""))
 lth.rules <- sapply(rules,length)
 for (i in 1:nvar){
   if (yes.rules[i]){
     for (r in 1:lth.rules[i]){
       if (is.na(rvalues[[i]][r]) & !isfactor[i]){
         nonmissing[[i]][r] <- with(data,sum(!is.na(data[eval(parse(text=rules[[i]][r])),i])))
       } else if (is.na(rvalues[[i]][r]) & isfactor[i]){    # different for factors because <NA> is treated as a level
         nonmissing[[i]][r] <- with(data,sum(!is.na(as.character(data[eval(parse(text=rules[[i]][r])),i]))))
       } else {
         nonmissing[[i]][r] <- with(data,sum(data[eval(parse(text=rules[[i]][r])),i]!=rvalues[[i]][r] |
                               is.na(data[eval(parse(text=rules[[i]][r])),i])))
       }
     }
   }
 }
 any.nonmissing <- sapply(nonmissing, function(x) any(x>0))
 if (any(any.nonmissing)>0) cat("Unexpected values (not obeying the rules) found for variable(s): ",
     paste(varnames[any.nonmissing>0],collapse=" "),".")

 # Check visit sequence 
 # all variables used in missing data rules have to be synthesised BEFORE 
 # the variables they apply to
 var.position <- lapply(get.vars, function(x) match(unique(x),varnames))
 var.in.vis   <- lapply(var.position, function(x) if (length(x)==0){
                                        x <- 0
                                        } else if(any(is.na(match(x,vis)))) {
                                        x[!is.na(match(x,vis))] <- match(x,vis)
                                        x[is.na(match(x,vis))]  <- nvar
                                        } else {
                                        x <- match(x,vis)})
 max.seq      <- sapply(var.in.vis,max,na.rm=T)
 not.synth    <- match(1:nvar,vis)[!is.na( match(1:nvar,vis))] <= max.seq[!is.na( match(1:nvar,vis))]
 if (any(not.synth,na.rm=TRUE)) stop("Variable(s) used in missing data rules for ",
       paste(varnames[!is.na( match(1:nvar,vis))][not.synth & !is.na(not.synth)],collapse=" "),
       " have to be synthesised BEFORE the variables they apply to.")

 # Check if a variable with missing values predicts other variables only if its
 # missing values are a subset of the missing values of the predicted variables
 # and remove from a prediction matrix if not. 
 # It refers to missing values coded as NA, otherwise variable can be used as 
 # a predictor without restrictions.
  
 #for (i in 1:nvar){
 #  if (!is.na(rvalues[i])) data[with(data,eval(parse(text=rules[i]))),i] <- NA
 #}
 patternRules <- matrix(0,nrow=nrow(data),ncol=ncol(data))
 for (i in 1:nvar){
   if (yes.rules[i]){
     for (r in 1:lth.rules[i]){
       if(is.na(rvalues[[i]][r])) patternRules[with(data,eval(parse(text=rules[[i]][r]))),i] <- 1
     }
   }
 }
 patternNA <- is.na(data)+0
 patternNA <- ifelse(patternRules==patternNA,patternNA,0)
 diffNAij  <- function(i,j,dataNA) sum(dataNA[,i]-dataNA[,j]<0)
 diffNA    <- Vectorize(diffNAij,vectorize.args=list("i","j"))
 predNA    <- outer(1:nvar,1:nvar,diffNA,dataNA=patternNA)
   
 # predNAwrong <- which ((pred==1 & predNA>0),arr.ind=TRUE)
 # pred        <- ifelse((pred==1 & predNA>0),0,pred)
 # if(nrow(predNAwrong)>0) cat(paste("Missing values of variable ",
 # varnames[predNAwrong[,2]]," are not a subset of missing values of variable ",
 # varnames[predNAwrong[,1]]," and cannot be used as its predictor (removed).\n",sep=""),
 # "\n",sep="")

 # check length of missing data codes for continuous variable
  if (length(contNA)!=nvar)
    stop(paste("The length of 'contNA' (",length(contNA),
     ") does not match the number of columns in the data (",nvar,").",sep=""))

 setup$predictorMatrix <- pred

 return(setup)
}
##-----------------end of--check.rules.syn----------------------------



#----------------------- now syn continues here ----------------------
# Basic checks of provided parameters:
# dimensions, consistency, replication, ...

 call <- match.call()
 nvar <- ncol(data)

 if(!(is.matrix(data) | is.data.frame(data)))
    stop("Data should be a matrix or data frame.")
 if(nvar < 2) stop ("Data should contain at least two columns.")
 if(length(event)!=nvar) stop("'event' is of length ",length(event),
                              " and should be of length ",nvar,".",sep="")
 if(length(denom)!=nvar) stop("'denom' is of length ",length(denom),                  
                              " and should be of length ",nvar,".",sep="")
 if(length(smoothing)!=nvar) stop("'smoothing' is of length ",length(smoothing),                  
                                  " and should be of length ",nvar,".",sep="")
                              
 # S U B S A M P L E   S I Z E
 if(k!=nrow(data)) {
  # if (k > nrow(data)) {
  #   cat("Warning: Subpopulation size (k=",k,") cannot be greater than the population size (",
  #       nrow(data),").\n","Synthetic data sets of same size as data will be produced.\n\n",sep="")
  #       k <- nrow(data)
  # } else
   cat("Sample(s) of size ",k," will be generated from original data of size ",
         nrow(data),".\n\n",sep="")
 }

 # M E T H O D S
 method <- gsub(" ","",method) # remove any spaces in or around method
 # expand user's syhthesising method (single string) to all variables
 if (length(method)==1){
   if(is.passive(method)) stop("Cannot have a passive syhthesising method for every column.")
   method <- rep(method, nvar)
   method[visitSequence[1]] <- "sample" # set 1st in visit seq to "sample"
   # set method to "" for vars not in visitSequence
   method[setdiff(1:nvar,visitSequence)] <- ""
 }
 
 # if user specifies multiple methods, check the length of the argument
 # methods must be given for all columns in the data
 if (length(method)!=nvar) stop(paste("The length of method (",length(method),
    ") does not match the number of columns in the data (",nvar,").",sep=""))

 # P R E D I C T O R   M A T R I X
 if(!is.null(predictorMatrix)){
   if (!is.matrix(predictorMatrix)) {
     stop("Argument 'predictorMatrix' is not a matrix.")
   } else if (nvar!=nrow(predictorMatrix)| nvar!=ncol(predictorMatrix))
     stop(paste("The 'predictorMatrix' has ",nrow(predictorMatrix)," row(s) and ",ncol(predictorMatrix),
          " column(s). \nBoth should match the number of columns in the data (",nvar,").",sep=""))
 }

 data     <- as.data.frame(data)
 varnames <- dimnames(data)[[2]]

 # Perform various validity checks on the specified arguments
 setup <- list(visitSequence = visitSequence,
               method = method,
               defaultMethod = defaultMethod,
               predictorMatrix = predictorMatrix,
               nvar = nvar,
               varnames = varnames, 
               rules = rules,
               rvalues = rvalues,
               contNA = contNA,                                         
               event = event,
               denom=denom)                   #GRdenom new
              
 setup <- check.visitSequence.syn(setup)
 setup <- check.predictorMatrix.syn(setup)


 # C H A N G E  D A T A  T Y P E  &  M O D I F Y  F A C T O R  L E V E L S
#---
 # apply only if in predictor matrix 
                                                                               # GR added condition and else
 if (!is.null(setup$predictorMatrix) & sum(setup$predictorMatrix>0)) {
             inpred   <- apply(setup$predictorMatrix!=0,1,any)*(!(method %in% c("","sample") )) |                # GR added to allow null methods not affected
             apply(setup$predictorMatrix!=0,2,any)  # if anywhere in predictorMatrix
        }
  else {
    inpred<-rep(FALSE,sqrt(length(setup$predictorMatrix)))
  }
 notevent <- is.na(match(1:nvar,setup$event))       # if not in event list

 # Convert any character variables into factors for variables in pred
 ischar    <- sapply(data,is.character)
 chartofac <- (ischar * inpred)>0
 if (sum(chartofac)>0) {
   cat("Variable(s): ",varnames[chartofac],
       " have been changed from character to factor.\n",sep=" ")
   for (j in (1:nvar)[chartofac]) data[,j] <- as.factor(data[,j]) 
 }
 # Changing numeric variables with fewer than 'nlevelmax' into factors
 #  Default for this now set to 5, 20 too many as picked up months
 #  Also only need to do this if variable in predictionMatrix
 #  and any inappropriate methods are changed to the default for factors
 nlevel      <- sapply(data, function(x) length(table(x)))
 ifnum       <- sapply(data, is.numeric)
 vartofactor <- which(nlevel<=nlevelmax & ifnum & inpred & notevent)
 for (j in vartofactor) data[,j] <- as.factor(data[,j])
 if (length(vartofactor)>0) {
   cat("Variable(s): ",paste0(varnames[vartofactor],collapse=", "),
       "numeric but with fewer than", nlevelmax,"levels turned into factor(s).\n\n",sep=" ")
   for (j in vartofactor) {
     if (setup$method[j] %in% c("norm","norm.proper",
                                "normrank","normrank.proper")) {
       if (nlevels[j]==2) setup$method[j] <- defaultMethod[2]
       else setup$method[j] <- defaultMethod[3]
  	   cat("Method for ",varnames[j]," changed to ",setup$method[j],"\n\n")
     }
     }
 }

 # Modifies a factor by turning NA into an extra level
 isfactor  <- sapply(data,is.factor)
 for (j in (1:nvar)[isfactor & inpred & notevent]){
   data[,j] <- addNA(data[,j],ifany=TRUE)
 } 

 # Identify any factors with > maxfaclevels levels that are in visitSequence
 too.many.levels <- sapply(data, function(x) length(levels(x))) > maxfaclevels
 if (any(inpred & too.many.levels)) {
   cat("Factor(s) with more than",maxfaclevels,"levels:",
     paste0(varnames[inpred & too.many.levels],collapse=", "),
     "\nIn case of computational problems, consider removing them from prediction matrix or combining categories.\n\n")
 }
 
#---

 setup                     <- check.method.syn(setup, data, proper)
 if (any(rules!="")) setup <- check.rules.syn (setup, data)

 method          <- setup$method
 predictorMatrix <- setup$predictorMatrix
 visitSequence   <- setup$visitSequence
 event           <- setup$event
 rules           <- setup$rules
 rvalues         <- setup$rvalues
 contNA          <- setup$contNA
 defaultMethod   <- setup$defaultMethod
 denom           <- setup$denom                       #GRdenom new

 # Not used variables are identified and dropped if drop.not.used==T
 # reclculate inpred & notevent in case they have changed after
 # check.method and check.data
 if (sum(predictorMatrix)>0){                            # GR condition added
   inpred      <- apply(predictorMatrix!=0,1,any) | apply(predictorMatrix!=0,2,any) # if anywhere in predictorMatrix
	 ispredictor <- apply(predictorMatrix!=0,2,any)    # if used as predictor
 }
 else inpred<-ispredictor<-rep( 0,sqrt( length(predictorMatrix) ) ) 

 notinvs     <- is.na(match(1:nvar,visitSequence)) # if not in visitSequence
 notsynth    <- notinvs | (!notinvs & method=="")  # if not synthesised
 notevent    <- is.na(match(1:nvar,event))         # if not in event list

 # identify columns not used as events or predictors or in visitSequnce
 out <- !inpred & notevent & notsynth

 if (any(out)) {
   cat("Variable(s):",paste0(varnames[out],collapse=", "),
       "not synthesised or used in prediction.\n",sep=" ")
   if (drop.not.used==T) cat("The variable(s) will be removed from data and not saved in synthesised data.\n\n")
   else cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n\n")
 }
   
# remove columns not used from data and replace predictor matrix, visit sequence, nvar and others
 if (any(out) & drop.not.used==T){
   if(sum(!out)==0) stop("No variables left to be synthesised") ######to stop if all data excluded 
   newnumbers       <- rep(0,nvar)
   newnumbers[!out] <- 1:sum(!out)
   visitSequence    <- newnumbers[visitSequence]
   visitSequence    <- visitSequence[!visitSequence==0]
   predictorMatrix  <- predictorMatrix[!out,!out]
   event[event!=0]  <- newnumbers[event[event!=0]]
   event            <- event[!out]
   data             <- data[,!out]
   nvar             <- sum(!out)
   method           <- method[!out]
   varnames         <- varnames[!out]

   if (nvar==1) {                             #  GR added  note having to reassign character vector
   	 cl<-class(data)
  	 data<-data.frame(data)
  	 if (cl=="character") data[,1]<-as.character(data[,1])
  	 names(data)<-varnames
   }
   contNA           <- contNA[!out]
   rules            <- rules[!out]
   rvalues          <- rvalues[!out]
   # recalculate these
   if (sum(predictorMatrix>0)) {                                 # GR condition added
	   inpred <- apply(predictorMatrix!=0,1,any) |
	             apply(predictorMatrix!=0,2,any)      # if anywhere in predictorMatrix
	   ispredictor <- apply(predictorMatrix!=0,2,any) # if used as predictor
   }
   else inpred <- ispredictor<-rep(0,sqrt(length(predictorMatrix))) 

   notinvs  <- is.na(match(1:nvar,visitSequence))  # if not in visitSequence
   notsynth <- notinvs | (!notinvs & method=="")   # if not synthesised
   notevent <- is.na(match(1:nvar,event))          # if not in event list
 }

 # Print out info on variables not synthesised but used in prediction
 pred_not_syn <- (ispredictor & notsynth)
 if (sum(pred_not_syn )>0) {
   cat("Variable(s):", paste0(varnames[ispredictor & notsynth],collapse=", "),
       "used as predictors but not synthesised.\n",sep=" ")
   if (drop.pred.only==T) cat("The variable(s) will not be saved with the synthesised data.\n\n")
   else  {
     cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n\n")
     pred_not_syn[pred_not_syn==TRUE]<-FALSE
   }
 } 
                                                      #  GR condition added
 if (sum(predictorMatrix)>0) {
	 pm <- padMis.syn(data, method, predictorMatrix, visitSequence,
			   nvar, rules, rvalues, defaultMethod, contNA, smoothing, denom)

	 # Pad the Syhthesis model with dummy variables for the factors
	 # p <- padModel.syn(data, method, predictorMatrix, visitSequence,
	 #                   nvar, rules, rvalues)
	 p  <- padModel.syn(pm$data, pm$method, pm$predictorMatrix, pm$visitSequence,
			   pm$nvar, pm$rules, pm$rvalues, pm$factorNA, pm$smoothing, pm$denom)
   
   if (k != dim(data)[1]){
   	 p$syn <- p$syn[sample(1:dim(data)[1],k,replace=TRUE),]
   }
 	 if(sum(duplicated(names(p$data))) > 0)
     stop("Column names of padded data should be unique.")
 }
 else {p <- list(data = data,                             #  GR this added
              syn = data,
              predictorMatrix = predictorMatrix, 
              method = method, 
              visitSequence = visitSequence, 
              rules = rules,
              rvalues = rvalues, 
              event=event,                  #GRdenom new
              denom=denom,                  #GRdenom new
              categories = NULL,
              smoothing=smoothing)         
 }

 if (m > 0) {
   syn <- list(m)
	 for(i in 1:m){
     syn[[i]] <- data
     if (k != dim(data)[1]) syn[[i]] <-
       syn[[i]][sample(1:dim(data)[1],k,replace=TRUE),]
   }    
 }
 else syn <- NULL

 syn <- sampler.syn(p, data, m, syn, visitSequence, rules, rvalues, 
                    event, proper, minbucket, printFlag, k, pred_not_syn)
 if (m==1) syn <- syn[[1]]

#-----------------------------------------------------------------------
# restore the original NA's in the data
# for(j in p$visitSequence) p$data[(!r[,j]),j] <- NA
  
 names(method)        <- varnames
 names(visitSequence) <- varnames[visitSequence]

 # save, and return, but don't return data
 syndsobj <- list(call = call,
                 m = m,
                 syn = syn,
                 method = method,
                 visitSequence = visitSequence,
                 predictorMatrix = predictorMatrix,
                 event = event,
                 smoothing=smoothing,
                 denom=denom,                  #GRdenom new
                 minbucket = minbucket,
                 proper = proper,
                 n = nrow(data),
                 k = k,
                 rules = rules,
                 rvalues = rvalues,
                 contNA = contNA,
                 drop.not.used = drop.not.used,
                 drop.pred.only = drop.pred.only)
 if (diagnostics) syndsobj <- c(syndsobj, list(pad = p))
 class(syndsobj) <- "synds"
 return(syndsobj)
}

