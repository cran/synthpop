padModel.syn <- function(data, method, predictorMatrix, visitSequence,
                         nvar, rules, rvalues, factorNA, smoothing, denom) {

 # Function called by syn to make dummy variables data frame data is 
 # augmented by columns for dummy variables when they are used as predictors. 
 # This is returned as part of a list which also contains structures 
 # to tell sampler.syn how to select columns

  categories <- data.frame(yes.no.categorical = factor(rep(FALSE,nvar), 
                                                levels=c("TRUE","FALSE")), 
                           number.of.dummies  = rep(0,nvar), 
                           yes.no.dummy       = factor(rep(FALSE,nvar), 
                                                levels=c("TRUE","FALSE")), 
                           corresponding.column.dummy = rep(0,nvar))
 
 # This is a data frame with a row for each variable and extra rows for 
 # each of the dummy variables added at the end of the j loop with
 # col1 = T/F if category and predictor in parametric model                     #!BN1605
 # col2 = number of dummy variables for factors with col1==TRUE, else 0         #!BN1605
 # col3 = yes.no.dummy TRUE for dummy variables
 # col4 = corresponding.column.dummy original column number for dummies, else 0

  pred.with.cart <- method %in% c("ctree","ctree.proper","cart","cart.proper")         #!BN1605

  for(j in 1:nvar){
    if ((is.factor(data[,j]) & any(predictorMatrix[1:nvar,j]!=0 & !pred.with.cart)) |  #!BN1605
        (factorNA[j]==TRUE & !pred.with.cart[j])){                                     #!BN1605
      categories[j, 1] <- TRUE

      # all factors defined to have treatment contrasts
      data[, j] <- C(data[, j], contr.treatment)
      n.dummy   <- length(levels(data[, j])) - 1
      categories[j, 2] <- n.dummy

      # predictorMatrix is given extra rows and columns for the dummy variables
      # rows are set to zero initially
      predictorMatrix <- rbind(predictorMatrix, matrix(0,
                               ncol=ncol(predictorMatrix), nrow=n.dummy))
      
      # columns are set to zero and then for vars with non-CART method          #!BN1605
      # copied from an original variable j in predictorMatrix for               #!BN1605
      # -> 1 for all the rows for which this variable is being used as          #!BN1605
      # a predictor in a non-CART model, 0 otherwise                            #!BN1605
      predictorMatrix <- cbind(predictorMatrix, matrix(0, ncol=n.dummy,         #!BN1605
                               nrow=nrow(predictorMatrix)))                     #!BN1605
      predictorMatrix[!pred.with.cart,(ncol(predictorMatrix)-n.dummy+1):        #!BN1605
        ncol(predictorMatrix)] <- matrix(rep(predictorMatrix[!pred.with.cart,j],times=n.dummy)) #!BN1605
                               
      # the original categorical variable is removed from predictors (=insert zeros)
      # for variables with non-CART method 
      predictorMatrix[!pred.with.cart,j] <- 0                                   #!BN1605


 # insert the column number for first of this set of dummies into
 # the visit sequence immediately after the jth column is predicted
      if (any(visitSequence == j)){
        # set an original categorical variable as predictor for its dummies  
        predictorMatrix[(ncol(predictorMatrix) - n.dummy + 1):
                         ncol(predictorMatrix), j] <- rep(1, times = n.dummy)
        # insert dummies into visit sequence 
        newcol <- ncol(predictorMatrix) - n.dummy + 1
        nloops <- sum(visitSequence == j)
          for (ii in 1:nloops){
            idx <- (1:length(visitSequence))[visitSequence==j][ii]
            visitSequence <- append(visitSequence, newcol, idx)
          }
      }

 # augment the data with columns for the new dummies
      data <- (cbind(data, matrix(0, ncol = n.dummy, nrow = nrow(data))))

 # set dummies to missing when variable is missing
      data[is.na(data[, j]), (ncol(predictorMatrix) - n.dummy + 1):
                              ncol(predictorMatrix)] <- NA
      cat.column <- data[!is.na(data[, j]), j]  # these are the non missing values of this factor

 # next bit sets the colums for the dummies to the dummy variables 
 # when data are not missing and labels columns
      data[!is.na(data[, j]),(ncol(predictorMatrix) - n.dummy + 1):
           ncol(predictorMatrix)] <- model.matrix(~cat.column - 1)[,-1]
      names(data)[(ncol(predictorMatrix) - n.dummy + 1):
                   ncol(predictorMatrix)] <- paste(attr(data,"names")[j],
                                                   (1:n.dummy),sep=".")
      method     <- c(method, rep("dummy", n.dummy))
      rules      <- c(rules,rep(rules[j],n.dummy))
      rvalues    <- c(rvalues,rep(rvalues[j],n.dummy))
      categories <- rbind(categories, 
                          data.frame(yes.no.categorical = rep(FALSE,n.dummy),
                          number.of.dummies = rep(0, n.dummy), 
                          yes.no.dummy      = rep(TRUE, n.dummy), 
                          corresponding.column.dummy = rep(j,n.dummy)))
    }
  }                                                       
   
  varnames <- dimnames(data)[[2]]  # now includes dummy names
  dimnames(predictorMatrix) <- list(varnames, varnames)
  names(method) <- varnames
  names(visitSequence) <- varnames[visitSequence]
  dimnames(categories)[[1]] <- dimnames(data)[[2]]
  
  #print(predictorMatrix)
  #print(visitSequence)
  #print(categories)
  
  return(list(data = as.data.frame(data), 
              syn = as.data.frame(data),
              predictorMatrix = predictorMatrix, 
              method = method, 
              visitSequence = visitSequence, 
              rules = rules,
              rvalues = rvalues, 
              categories = categories,
              smoothing = smoothing,
              denom=denom))
}

