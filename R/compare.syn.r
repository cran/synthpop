###-----compare.synds------------------------------------------------------

compare.synds <- function(object,data,vars=NULL,nrows=1,ncols=1,breaks=50,...){

# vars selects just variables wanted

 #print(object$call)
 if (is.null(data)) stop("Requires parameter data to give name of the real data.\n")
 if (!class(object)=="synds") stop("Object must have class synds.\n")
 if (object$m==1) synds <- object$syn else synds  <- object$syn[[1]]
 synnames    <- names(synds)
 realnames   <- names(data)
 commonnames <- synnames[match(realnames,synnames)]
 if (!is.null(vars)) commonnames <- commonnames[match(vars,commonnames)]
 commonnames <- commonnames[!is.na(commonnames)]
 
 posinreal   <- match(commonnames,realnames)
 posinsyn    <- match(commonnames,synnames)
   
 cat("Data types of variables selected for comparison:\n\n")
 if (length(vars)==1) {
   classes <- class(synds[,posinsyn])
   tc <- cbind(classes); rownames(tc)<-""; colnames(tc)<-vars 
   print(tc)
 } else {
   classes <- sapply(synds[,posinsyn],class)
   print(classes)
 }
 
 par(mfrow=c(nrows,ncols))
  #if (object$m==1){
   for (i in 1:length(commonnames)){
     if (classes[i] %in% c("numeric","integer") & length(table(data[,posinreal[i]]))>breaks){  #GR changed
       rr <- data[,posinreal[i]]
       ss <- synds[,posinsyn[i]]
       if (any(rr %in% object$contNA[[posinreal[i]]])){
         missrr <- rr %in% object$contNA[[posinreal[i]]]
         missss <- ss %in% object$contNA[[posinreal[i]]]
         t1 <- table(rr[missrr],exclude=NULL)/length(rr)*100
         t2 <- table(ss[missss],exclude=NULL)/length(ss)*100
         labelsNA <- unique(c(names(t1),names(t2)))
         tt <- matrix(0,2,length(labelsNA))
	       tt[1,match(names(t1),labelsNA)] <- t1
         tt[2,match(names(t2),labelsNA)] <- t2
         dimnames(tt)[[1]] <- c("real","synthetic")
         labelsNA[is.na(labelsNA)] <- "NA"
         dimnames(tt)[[2]] <- labelsNA
         cat("\nComparing actual (black) with synthetic (grey) for % missing data in ",
             commonnames[i],"\n\n")
		     print(tt)
		     if (nrows==1 & ncols==1) layout(matrix(c(2,1),nrow=1), widths=c(2,1))
         barplot(tt,beside=T,main=paste("% missing\nfor",commonnames[i]),space=c(0,0.3),ylab="%",las=1,...)
         #cat("\nPress return to view comparison of non-missing values\n\n")
         #bringToTop(-1)
         #readline()
        
         rr <- rr[!(rr %in% object$contNA[[posinreal[i]]])]
         ss <- ss[!(ss %in% object$contNA[[posinreal[i]]])]
       }

	     hh  <- hist(c(rr,ss),plot=F,breaks=breaks)
	     #print(hh$breaks)
	     hrr <- hist(rr,breaks=hh$breaks,plot=F)
	     hss <- hist(ss,breaks=hh$breaks,plot=F)
       tt  <- rbind(hrr$counts/length(rr)*100,
                    hss$counts/length(ss)*100)
       cat("\nComparing actual (black) with synthetic (grey) for ",commonnames[i],"\n\n")
       dimnames(tt)[[1]] <- c("real","synthetic")
       dimnames(tt)[[2]] <- hh$breaks[-length(hh$breaks)]
		   print(tt)
       barplot(tt,beside=T,main=commonnames[i],space=c(0,0.3),ylab="%",las=1,...)
		   if (nrows==1 & ncols==1) par(mfrow=c(nrows,ncols))
       if (i < length(commonnames)) {
         cat("\nPress return to view next variable\n\n")
         #bringToTop(-1)
         readline()
       }
     }

     else if (classes[i]=="factor" || length(table(data[,posinreal[i]]))<=breaks){ #GR changed
       if (any(is.na(synds[,posinsyn[i]]))) {        # GR data changed to synds and next bit changed
         t1 <- table(data[,posinreal[i]],exclude=NULL)/dim(data)[1]*100
         t2 <- table(synds[,posinsyn[i]],exclude=NULL)/dim(synds)[1]*100
       } else {
       	 t1 <- table(data[,posinreal[i]])/dim(data)[1]*100
         t2 <- table(synds[,posinsyn[i]])/dim(synds)[1]*100
       }
       #labels <- sort(unique(c(names(t1),names(t2))),na.last=TRUE) 
       labels <- unique(c(names(t1),names(t2)))      # this needed if some categories don't exists in either real or syn
       tt <- matrix(0,2,length(labels))
       tt[1,match(names(t1),labels)] <- t1
       tt[2,match(names(t2),labels)] <- t2
       dimnames(tt)[[1]] <- c("real","synthetic")
       labels[is.na(labels)] <- "NA"
       dimnames(tt)[[2]] <- labels
       cat("\nComparing percentages actual (black) with synthetic (grey) for ",
           commonnames[i],"\n\n")
       print(tt)
       par(mar = c(7.5, 7, 4, 2) + 0.1)
       barplot(tt,beside=T,main=commonnames[i],space=c(0,0.3),ylab="%",las=2,xaxt="n",...)
       text(seq(1.5,2.3*dim(tt)[[2]],2.3), par("usr")[3] - 0.25, srt=30, adj=1,
            labels=dimnames(tt)[[2]], xpd=TRUE)
       if (i < length(commonnames)) {
         cat("\nPress return to view next variable\n\n")
         #bringToTop(-1)
         readline()
       }
     }
     else cat("\nClass of type ",classes[i]," doesn't have a method.\n\n")
	 }
 #}
}


###-----compare.fit.synds----------------------------------------------------
compare.fit.synds <- function(object, real.data, plot="Z", return.result=TRUE,
                              plot.intercept=FALSE, col=1:2, ...) {

 # compares and plots fits to synthetic and real data
 # first parameter must be a fit to synthetic data from glm.synds()

 call <- match.call()
 if (!class(object)=="fit.synds") stop("Object must have class fit.synds\n")
 if (!is.data.frame(real.data)) stop("real.data must be a data frame\n")  # theoretically can be a matrix
 #m <- object$m
 sum.fit.object   <- summary.fit.synds(object) #synthpop:::
 fitting.function <- object$fitting.function
 
 # get fit to real data
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=object$call$formula,
                     family=object$call$family,data=call$real.data)))

 # combine with fit to synthetic
 res <- cbind(real.fit$coefficients[,-4],sum.fit.object$coefficients[match(rownames(real.fit$coefficients),rownames(sum.fit.object$coefficients)),])
 dimnames(res)[[2]][1:3] <- c("beta real","se beta real","Z real")
 result <- res

 # save result to return
 # if (return.result==TRUE) result=res else result=invisible()

 # plotting
 #----------
 if (plot %in% c("both","coef","Z")) {

   # add baselines for factors for showing on plots
   
   yvar <- as.character(object$call$formula[2])
   # get vector of explanatory variable numbers in formula
   vars_in_form <- get.names(object$call$formula,names(real.data)) #synthpop:::

   # first coefficient is intercept
   newres <- res[1,]
   bases  <- ""
   j <- 2; l <- 2
   for (i in vars_in_form){
	   nl <- levels(real.data[,i])
	   if (!is.null(nl)) {
	     newres <- rbind(newres,rep(0,8))   # line of zeros for baseline
	     dimnames(newres)[[1]][j] <- nl[1]
	     bases  <- c(bases,names(real.data)[i])
	     j <- j+1
	     for (kl in 2:length(nl)) {
	       newres <- rbind(newres,res[l,])
		     dimnames(newres)[[1]][j] <- nl[kl]
		     bases <- c(bases,"")
		     j <- j+1; l <- l+1
	     }
     } else {
	     newres <- rbind(newres,res[l,])
		   dimnames(newres)[[1]][j] <- dimnames(res)[[1]][l]
		   bases <- c(bases,"")
		   j <- j+1; l <- l+1
	   }
   }
   dimnames(newres)[[1]][1] <- dimnames(res)[[1]][1]
   res <- newres
   
   if (plot.intercept==FALSE) {res <- res[-1,]; bases <- bases[-1]}
   #print(res) 
   bases[bases!=""] <- paste0("baseline for ", bases[bases!=""])

   if (plot=="both") par(mfrow=c(2,1))

   # coefficient plot
   #-----------------
   if (plot %in% c("both","coef")) {
	   yvals <- dim(res)[[1]]
	   xlims <- c(min(res[,c(1,4)]-2*res[,c(2,5)]),
                max(res[,c(1,4)]+2*res[,c(2,5)]))
	   par(mar=c(5.1,2.1,4.1,18))
	   plot(xlims,c(0,1+yvals),type="n",yaxt="n",
          ylab="",xlab="coefficient",bty="n")
	   title(main=paste("Compare coefficients for fit to",yvar),xpd=T)
	   #lines(xlims,c(-0.4,-0.4),xpd=T,lwd=2)      
	   
     points(res[bases=="","beta real"],(1:yvals)[bases==""],pch=15,col=col[1],...)
	   segments(res[bases=="","beta real"]-1.96*res[bases=="","se beta real"],
              (1:yvals)[bases==""],
              res[bases=="","beta real"]+1.96*res[bases=="","se beta real"],
              (1:yvals)[bases==""],col=col[1],...)
	   
     points(res[bases=="","beta syn"],(1:yvals+0.3)[bases==""],pch=16,col=col[2],...)
	   segments(res[bases=="","beta syn"]-1.96*res[bases=="","se beta syn"],
              (1:yvals)[bases==""]+0.3,
              res[bases=="","beta syn"]+1.96*res[bases=="","se beta syn"],
              (1:yvals)[bases==""]+0.3,col=col[2],...)
	   
     points(res[bases!="",4],(1:yvals)[bases!=""],pch=17,...)
	   lines(c(0,0),c(0,yvals+2),lty=3)
     
     text(rep(xlims[2],yvals)+.2,(1:yvals)+.1,dimnames(res)[[1]],adj=0,xpd=T)
	   text(rep(xlims[1],yvals),(1:yvals)+.1,bases,xpd=T,adj=0)
     
     legend(xlims[2],0.2,c("synthesised","real"),
       lty=1,col=col[2:1],pch=16:15,xpd=T,bty="n",...)
   }
   
   # Z plot
   #-------
   if (plot %in% c("both","Z")){
	   yvals <- dim(res)[[1]]
	   xlims <- c(min(res[,c(3,6)])-2,
                max(res[,c(3,6)])+2)
	   par(mar=c(5.1,2.1,4.1,18))
	   plot(xlims,c(0,1+yvals),type="n",yaxt="n",ylab="",xlab="Z value",bty="n")
	   title(main=paste("Compare Z values for fit to",yvar),xpd=T)
	   #lines(xlims,c(-0.4,-0.4),xpd=T,lwd=2)
	   
     points(res[bases=="","Z real"],(1:yvals)[bases==""],pch=15,col=col[1],...)
	   segments(res[bases=="","Z real"]-1.96,(1:yvals)[bases==""],
              res[bases=="","Z real"]+1.96,(1:yvals)[bases==""],col=col[1],...)
	   
     points(res[bases=="","Z syn"],(1:yvals+0.3)[bases==""],pch=16,col=col[2],...)
  	 segments(res[bases=="","Z syn"]-1.96,(1:yvals)[bases==""]+0.3,
              res[bases=="","Z syn"]+1.96,(1:yvals)[bases==""]+0.3,col=col[2],...)
     
     points(res[bases!="","Z syn"],(1:yvals)[bases!=""],pch=17,...)
	   lines(c(0,0),c(0,yvals+2),lty=3)
     
     text(rep(xlims[2],yvals)+.2,(1:yvals)+.1,dimnames(res)[[1]],adj=0,xpd=T)
	   text(rep(xlims[1],yvals),(1:yvals)+.1,bases,xpd=T,adj=0)
	   
     legend(xlims[2],0.2,c("synthesised","real"),
       lty=1,col=col[2:1],pch=16:15,xpd=T,bty="n",...)
   }
 }
 par(mfrow=c(1,1))
 res <- list(fit.synds.call=object$call,coefficients=result)
 class(res) <- "compare.fit.synds"
 if (return.result==TRUE) return(res) else invisible()
}  

