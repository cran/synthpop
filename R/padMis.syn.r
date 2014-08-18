padMis.syn <- function(data, method, predictorMatrix, visitSequence,
                       nvar, rules, rvalues, defaultMethod, contNA, smoothing, denom) {

 # Function called by syn to make dummy/factor variable for missing values
 # in continuous variables. Data is augmented by columns for dummy/factor 
 # variables when they are used in sythesis. 

 # Check presence of missing values not covered by missing rules
 # (missing values for non-numeric variables are not counted)   
 No.NA <- vector("list",nvar)
 yes.rules <- sapply(rules, function(x) any(x!=""))
 com.rules <- lapply(rules, paste, collapse=" | ")
 for (j in 1:nvar){
   if (yes.rules[j]){
     No.NA[j] <- with(data,sum(data[!eval(parse(text=com.rules[[j]])),j] %in% contNA[[j]]))      
   } else {
     No.NA[j] <- sum(data[,j] %in% contNA[[j]]) 
   }
 }
 No.NA    <- sapply(No.NA,function(x) x>0)
 inpred   <- apply(predictorMatrix!=0,1,any)|apply(predictorMatrix!=0,2,any)
 factorNA <- rep(FALSE,nvar)

 for(j in 1:nvar){
    if (No.NA[j] & is.numeric(data[,j]) & inpred[j]==TRUE){
 
    # augment the data with a column for the original continuous variable with 
    # missing values replaced by zeros and a column for a new factor for 
    # missing values 
      y.0  <- ifelse(data[,j] %in% c(contNA[[j]],rvalues[[j]]),0,data[,j])
      y.NA <- ifelse(data[,j] %in% c(contNA[[j]],rvalues[[j]]),data[,j],0)
      y.NA <- addNA(y.NA,ifany=TRUE) 
      data <- cbind(data,y.0,y.NA)           
      name.0  <- paste(attr(data,"names")[j],0,sep=".")
      name.NA <- paste(attr(data,"names")[j],NA,sep=".")
      names(data)[(ncol(data)-1):ncol(data)] <- c(name.0,name.NA)
      factorNA[(ncol(data)-1):ncol(data)] <- c(FALSE,TRUE) 

    # predictorMatrix is given two extra rows and columns for the new variables
    # rows and columns are copied from an original variable j in predictorMatrix
      predictorMatrix <- rbind(predictorMatrix,  matrix(rep(predictorMatrix[j,], 
                               times=2),byrow=TRUE,nrow=2))
      predictorMatrix <- cbind(predictorMatrix, matrix(rep(predictorMatrix[,j], 
                               times=2),ncol=2))
    # the original variable is removed from predictors (=insert zeros) 
      predictorMatrix[,j] <- 0
    # the original variable is imputed passively so its predictors can be removed as well
      predictorMatrix[j,] <- 0

    # add methods for new variables
      method[ncol(data)-1] <- method[j]
      if (method[j] %in% c("ctree","ctree.proper","cart","cart.proper")) {
        method[ncol(data)] <- method[j]  
      } else {
        method[ncol(data)] <- ifelse(nlevels(data[,ncol(data)])==2,
                                     defaultMethod[2],defaultMethod[3])
      }   
    
    # pass smoothing to new variables
      smoothing[ncol(data)-1] <- smoothing[j]
      smoothing[ncol(data)]   <- ""
    
    # add denom for new variables
      denom[(ncol(data)-1):ncol(data)] <- 0
    
    # insert the column numbers for the new variables into the visit sequence 
    # before the jth column
      if (any(visitSequence==j)){                    
        newcols <- c(ncol(data),ncol(data)-1)
        idx <- (1:length(visitSequence))[visitSequence==j]-1
        visitSequence <- append(visitSequence,newcols,idx)
      # modify method for the original variable
        method[j] <- paste0("~(ifelse(",name.0,"==0 | is.na(",name.0,
            "),as.numeric(levels(",name.NA,"))[",name.NA,"],",name.0,"))")
      }

    # update missing rules and values for the new variables
      if (any(rules[[j]]!="")) rules[[ncol(data)-1]] <-
        c(rules[[j]][rules[[j]]!=""],paste(name.NA,"!=0",sep=""))
      else rules[[ncol(data)-1]] <- paste(name.NA,"!=0",sep="") 
      rules[[ncol(data)]]        <- rules[[j]]                      
      rules[[j]]                 <- ""
      #!BN1513 rule "year_death.NA!=0" should have only one correspnding 
      #! rvalue equal to 0; before a vector c(NA,0) was assigned instead of 0 
      if (length(rules[[j]])==1) rvalues[[ncol(data)-1]] <- 0           
      else rvalues[[ncol(data)-1]] <- c(rvalues[[j]],0)                 
      rvalues[[ncol(data)]]        <- rvalues[[j]]
    }
  }
   
  varnames <- dimnames(data)[[2]]  # now includes new names
  dimnames(predictorMatrix) <- list(varnames,varnames)
  names(method) <- varnames
  names(visitSequence) <- varnames[visitSequence]
  return(list(data = as.data.frame(data), 
              nvar = ncol(data),
              predictorMatrix = predictorMatrix, 
              method = method, 
              visitSequence = visitSequence, 
              rules = rules,
              rvalues = rvalues,
              factorNA = factorNA,
              smoothing = smoothing,
              denom = denom))
}
