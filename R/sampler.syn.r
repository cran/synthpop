#-----------------------------sampler.syn-------------------------------

sampler.syn <- function(p, data, m, syn, visit.sequence,
                        rules, rvalues, event, proper,
                        minbucket, print.flag, 
                        k, pred.not.syn, ...){
 # The sampler controls the generation of conditional distributions
 # This function is called by syn
 #
 # Authors: G Raab & B Nowok 2013-14
 # 
  if (m > 0){
	  if (print.flag) cat("syn  variables")
    for (i in 1:m){  # begin i loop : repeated synthesising loop
      if (print.flag) cat("\n",i,"   ",sep="")
## augment the data with the actual dummy variables  maybe not needed now??  GR
 # This next code replaces the dummy variables for a synthesised variable with 
 # the new values 
 #     for (j in setdiff(p$visit.sequence,visit.sequence)){
 #       cat.columns <- p$syn[, p$categories[j, 4]]
 #       p$syn[,(j:(j+p$categories[p$categories[j,4],2]-1))] <- 
 #                        matrix((model.matrix(~cat.columns-1)[,-1]),
 #                        ncol=p$categories[p$categories[j,4],2],nrow=nrow(p$data)) 
 #       p$syn[,j]<-as.numeric(p$syn[,j])
 #     }

      for(j in p$visit.sequence) {
        theMethod <- p$method[j]
        vname     <- dimnames(p$data)[[2]][j]
        
        if(print.flag & theMethod!="dummy"  & j<=ncol(data)) cat(" ",vname,sep="")
        if(j%%10==0 & j<=ncol(data)) cat("\n    ")
        
        ya <-  1:nrow(p$data) 
        ypa <- 1:k                  ############  new
        # ya=yavailable, ym=ymissing                                            
        if(any(p$rules[[j]]!="")) {
          com.rules  <- paste(p$rules[[j]],collapse=" | ")
          evalrul.y  <- with(p$data,eval(parse(text=com.rules)))
          ym         <- which(evalrul.y==TRUE & !is.na(evalrul.y))
          ya         <- setdiff(1:nrow(p$data),ym)                                  
          evalrul.yp <- with(p$syn,eval(parse(text=com.rules)))         
          ypm        <- which(evalrul.yp==TRUE & !is.na(evalrul.yp))        
          ypa        <- setdiff(1:nrow(p$syn),ypm)       
                    
        }                                                                       
           
        if(theMethod!="" & (!is.passive(theMethod)) & theMethod!="dummy" ){
          if (theMethod %in% c("sample","sample.proper")) {
            ##### new code for method sample
            y   <- p$data[ya, j]
            if (is.factor(y)) y <- y[,drop=TRUE]
            xp  <- length(ypa)
            x   <- length(ya)
            nam <- vname
            f   <- paste("syn", theMethod, sep = ".")
            p$syn[ypa, j]  <- do.call(f, args = list(y,xp,proper=proper, ...)) 
          }
          else
          {
            x    <- p$data[ya, p$predictor.matrix[j, ] == 1, drop = FALSE]
            xp   <- p$syn [ypa, p$predictor.matrix[j, ] == 1, drop = FALSE]
            y    <- p$data[ya, j]
            if (is.factor(y)) y <- y[,drop=TRUE]
            nam  <- vname
            f    <- paste("syn", theMethod, sep = ".")
            keep <- remove.lindep.syn(x, y, ...)
            x    <- x[, keep, drop = FALSE]
            xp   <- xp[, keep, drop = FALSE]
         
            if (theMethod=="surv.ctree") {
              if (event[j]==-1) yevent <- rep(1,length(y))
              else yevent  <- p$data[ya,event[j]]
              survres      <- do.call(f, args = list(y,yevent,x,xp,proper=proper,minbucket=minbucket,...))
              p$syn[ypa,j] <- survres[[1]]# synthetic data survival goes to p$syn
              if (event[j]!=-1) p$syn[ypa,event[j]] <- survres[[2]] # synthetic data event goes to p$syn
            }
            else if (theMethod=="logreg" & j<=ncol(data) & p$denom[j]!=0) {                   #GRdenom new
              p$syn[ypa, j] <- do.call(f, args = list(y,x,xp,denom=p$data[,p$denom[j]],denomp=p$syn[,p$denom[j]],proper=proper, ...))
            }   
            else {
              p$syn[ypa, j] <- do.call(f, args = list(y,x,xp,smoothing=p$smoothing[j],proper=proper,minbucket=minbucket,...))
            }

          }

          if(any(p$rules[[j]]!="")){
            if(length(p$rules[[j]])==1 & length(ypm)>0){
              p$syn[ypm,j] <- p$rvalues[[j]] 
            } else {
              for (r in 1:length(p$rules[[j]])){
                revalrul.yp  <- with(p$syn,eval(parse(text=p$rules[[j]][r])))  
                rypm <- which(revalrul.yp==TRUE & !is.na(revalrul.yp))
                if (length(rypm)>0) p$syn[rypm,j] <- p$rvalues[[j]][r]
              }
            }                 
          }  
        }

        else if (is.passive(theMethod)) {
          class0 <- class(p$syn[,j])
          p$syn[,j] <- model.frame(as.formula(theMethod), p$syn, na.action=na.pass)	#RJ - FIXED passive synthesising: as.formula()
          class(p$syn[,j]) <- class0
        }

        else if (theMethod=="dummy") {    # replace dummy variables in p$syn
          # getting dummy values from a synthesised categorical variable
          cat.columns <- p$syn[,p$categories[j,4]]  # this is the single column with the data for which this is the dummy
          model.frame(~cat.columns-1,data=p$syn) 
          p$syn[,(j:(j+p$categories[p$categories[j,4],2]-1))] <-   # replaces all the dummies for this variable with
          matrix((model.matrix(~cat.columns-1)[,-1]),              # dummies calculated from the synthesised data
                  ncol=p$categories[p$categories[j,4],2],
                  nrow=nrow(p$syn))
          p$syn[,j] <- as.numeric(p$syn[,j])
          remove("cat.columns")
        }

      } # end j loop 
    


      #if (k==dim(data)[1]) syn[[i]] <- p$syn[,1:dim(data)[2]]
      #else syn[[i]] <- p$syn[sample(1:dim(data)[1],k),1:dim(data)[2]]
      syn[[i]] <- p$syn[,1:dim(data)[2]]
      nms<-names(data)
      # exclude unsynthesised if drop.pred.only set to true
      if (sum(pred.not.syn )>0) {
        syn[[i]] <- syn[[i]][,!pred.not.syn]
        nms<-nms[!pred.not.syn]                               #GR save names to use below if data just one column
      }
                                          # GR changes extra lines needed # to prevent a single character column being changed to a factor
      chgetochar<- (sum(!pred.not.syn)==1 & class(syn[[i]][,1])=="character")       
  
      syn[[i]] <- as.data.frame(syn[[i]])
      if (chgetochar) {
         syn[[i]][,1]<-as.character(syn[[i]][,1])
         names(syn[[i]])<-nms
      }
   

      #turn NA level in factors to missing NA's    # to delete - replaced by the code below
  #    for (j in (1:ncol(syn[[i]]))){
  #      if(is.factor(syn[[i]][,j])) {
  #        syn[[i]][is.na(as.character(syn[[i]][,j])),j] <- NA
  #        syn[[i]][,j] <- factor(syn[[i]][,j])
  #      }
  #    }

      #turn NA level in factors to missing NA's
      for (j in (1:ncol(syn[[i]]))){
        if(is.factor(syn[[i]][,j])) {
          syn[[i]][,j] <- factor(syn[[i]][,j],exclude=NA,levels=levels(syn[[i]][,j]))
        }
      }

    } # end i loop
 } # end synthesising
  
 if (print.flag) cat("\n")
 return(syn)
}


remove.lindep.syn <- function(x, y, eps=0.00001, maxcor=0.99999, 
                              allow.na=FALSE, ...) {
  if (ncol(x)==0) return(NULL) 
  if (eps <= 0) stop("\n Argument 'eps' must be positive.")
  xobs <- sapply(x,as.numeric)                                       #!BN1605 was xobs <- x
  yobs <- as.numeric(y)
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor)) # if y includes NA -> NAs error
  if (all(!keep)) warning("All predictors are constant or have too high correlation.")
  ksum <- sum(keep)
  cx   <- cor(xobs[, keep, drop=FALSE], use = "all.obs")
  eig  <- eigen(cx, symmetric = TRUE)
  ncx  <- cx
  while(eig$values[ksum]/eig$values[1] < eps) {
    j   <- (1:ksum)[order(abs(eig$vectors[, ksum]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx  <- cx[keep[keep], keep[keep], drop = FALSE]
    ksum <- ksum - 1
    eig  <- eigen(ncx)
  }
  # if (!all(keep)) cat("\tVariable(s): ", paste(dimnames(x)[[2]][!keep], collapse = ", "),
  #   " removed due to linear dependency",sep="")
  return(keep)
}

# make list of collinear variables
find.collinear <- function(x, threshold=0.99999, ...) {
  nvar      <- ncol(x)
  x         <- data.matrix(x)
  varnames  <- dimnames(x)[[2]]
  z         <- suppressWarnings(cor(x, use="pairwise.complete.obs"))
  hit       <- outer(1:nvar,1:nvar,"<") & (abs(z)>=threshold)
  collvar   <- which (hit==1 ,arr.ind=TRUE)
  collvar[] <- varnames[collvar]
  collvarpairs <- paste(collvar[,1],collvar[,2],sep=" & ")
  return(collvarpairs)
}
