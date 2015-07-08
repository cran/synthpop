###-----compare------------------------------------------------------------
compare <- function(object, data, ...) UseMethod("compare")

###-----compare.default----------------------------------------------------
compare.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----compare.synds------------------------------------------------------
compare.synds <- function(object, data, vars = NULL, 
  msel = NULL, breaks = 20, nrow = 2, ncol = 2, rel.size.x = 1,
  cols = c("#1A3C5A","#4187BF"), ...){ 
                                                                         ## GR6 & GR7 extra parameter nrow and drop plot.na
 if (is.null(data)) stop("Requires parameter data to give name of the real data\n", call.=FALSE)
 if (!is.data.frame(data)) stop("Argument data must be a data frame\n", call.=FALSE)          ## GR2 extra check
 #if (class(object)!="synds") stop("Object must have class synds\n", call.=FALSE )                                                                    
 if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s)", call.=FALSE)                                                                        

 # browser()
 # single / pooled synthetic data sets                                  
 if (object$m==1) {
   syndsall <- object$syn 
 } else if (length(msel) == 1) {
   syndsall <- object$syn[[msel]]
 } else if (length(msel) > 1 | is.null(msel)) {                               
   syndsall <- do.call(rbind,object$syn)
 }
 # list of synthetic data sets for non-pooled results 
 if (length(msel) > 1) {
   synds <- vector("list",length(msel))
   for (i in 1:length(msel)) synds[[i]] <- object$syn[[msel[i]]]
 }
 synnames    <- names(syndsall)
 realnames   <- names(data)
 commonnames <- synnames[match(realnames,synnames)]
 commonnames <- commonnames[!is.na(commonnames)]

 if (!is.null(vars)){
   if (!(all(vars %in% synnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,synnames))], collapse=", "),
     " not in synthetic data \n", call.=FALSE)
   if (!(all(vars %in% realnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,realnames))], collapse=", "),
     " not in observed data \n", call.=FALSE)
   commonnames <- commonnames[match(vars,commonnames)]
 }
 
 if (!(all(synnames %in% realnames))) cat("Warning: Variable(s)", 
   paste0(synnames[is.na(match(synnames,realnames))], collapse=", "),
   "in synthetic object but not in observed data\n",
   " Looks like you might not have the correct data for comparison\n")
 
 if ((length(commonnames) == 0) && (typeof(commonnames) == "character"))        #! when would it apply?
   stop("None of variables selected for comparison in data", call. = FALSE)

 df.obs    <- data[, commonnames, drop = FALSE]
 df.synall <- syndsall[, commonnames, drop = FALSE]
 if (length(msel) > 1) {
   df.syn <- vector("list",length(msel))
   for (i in 1:length(msel)) df.syn[[i]] <- synds[[i]][, commonnames, drop = FALSE]
 }
 num <- sapply(df.synall, is.numeric) | sapply(df.synall, is.integer)  
 fac <- sapply(df.synall, is.factor)  

 # to exclude from summaries if no missing in data
 if (sum(num) > 0) any.num.na <- unlist(apply(df.obs[num,,drop=FALSE],2,function(x) any(is.na(x))))  

 # frequency tables for factors
 if (sum(fac) > 0) {
   per.obs.fac <- ggfac(df.obs[, fac, drop = FALSE])
   if (length(msel) <= 1) per.syn.facall <- ggfac(df.synall[, fac, drop = FALSE], name = "synthetic")
   if (length(msel) > 1) {
     per.syn.fac <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.fac[[i]] <- ggfac(df.syn[[i]][, fac, drop = FALSE], 
       name = paste0("syn=", msel[i]))
   }
 } else {
   per.obs.fac    <- NULL
   per.syn.facall <- NULL
 }

 # frequency tables for numeric variables
 if (sum(num) > 0) {
   na.index <- match(colnames(df.obs[, num, drop = FALSE]), colnames(syndsall))    
   any.num.na.index <- match(colnames(df.obs[, num, drop = FALSE]), names(any.num.na))
   na <- object$cont.na[na.index]
   any.na <- any.num.na[any.num.na.index]                               
   lbreaks  <- as.list(rep(breaks, length(na)))                   
   df.both  <- rbind.data.frame(df.obs, df.synall)                ## GR3 ## gets limits from both
   per.both <- ggnum(df.both[, num, drop = FALSE], na = na, 
     breaks = lbreaks, anyna = any.na) 
   per.obs.num <- ggnum(df.obs[, num, drop = FALSE], na = na, 
     breaks = per.both$hbreaks, anyna = any.na)
   
   if (length(msel) <= 1) {
     per.syn.numall <- ggnum(df.synall[, num, drop = FALSE], name = "synthetic", 
       na = na, breaks = per.both$hbreaks, anyna = any.na) 
   } 
   if (length(msel) > 1) {
     per.syn.num <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.num[[i]] <- ggnum(df.syn[[i]][, num, drop = FALSE], 
       name =  paste0("syn=",msel[i]), 
       na = na, breaks = per.both$hbreaks, anyna = any.na)
   } 
 } else {
   per.obs.num<-NULL
   per.syn.numall <- NULL
 }

 # data frame for ploting 
 if (length(msel) <= 1) per.fac <- rbind.data.frame(per.obs.fac$perdf, 
   per.obs.num$perdf, per.syn.facall$perdf, per.syn.numall$perdf )

 if (length(msel) > 1) {
   per.fac <- rbind.data.frame(per.obs.fac$perdf, per.obs.num$perdf)
   for (i in 1:length(msel)) {
     if (sum(fac)>0) temp.fac <- per.syn.fac[[i]]$perdf else temp.fac <- NULL
     if (sum(num)>0) temp.num <- per.syn.num[[i]]$perdf else temp.num <- NULL  
     per.fac <- rbind.data.frame(per.fac, temp.fac, temp.num)
   }
 }
  per.fac$Variable <- factor(per.fac$Variable, levels = commonnames, 
   ordered = T, exclude = NULL)
  per.fac$Value <- factor(per.fac$Value, levels = unique(per.fac$Value), 
   ordered = T, exclude = NULL)
 
 # list of result tables
 if (length(msel) <= 1) {
   os.table.fac <- mapply(rbind, obs = per.obs.fac$perlist, 
     syn = per.syn.facall$perlist, SIMPLIFY = FALSE)
   os.table.num <- mapply(rbind, obs = per.obs.num$perlist, 
     syn = per.syn.numall$perlist, SIMPLIFY = FALSE)
 }
 if (length(msel) > 1) {
   os.table.fac <- per.obs.fac$perlist 
   os.table.num <- per.obs.num$perlist 
   for (i in 1:length(msel)) {
     if (sum(fac)>0) temp.fac <- per.syn.fac[[i]]$perlist else temp.fac <- NULL
     if (sum(num)>0) temp.num <- per.syn.num[[i]]$perlist else temp.num <- NULL
     os.table.fac <- mapply(rbind, os.table.fac, temp.fac, SIMPLIFY = FALSE)
     os.table.num <- mapply(rbind, os.table.num, temp.num, SIMPLIFY = FALSE)     
   }
  }
 os.table <- c(os.table.fac, os.table.num)
 if (is.null(msel)){
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed","synthetic")
 } else {
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed",paste0("syn=", msel))
 }
 Value <- Percent <- Data <- NULL                               ## otherwise 'no visible binding for global variables'
 # sorts the factor labels in the right order for numeric vars
 per.fac$Value <- as.character(per.fac$Value)
 vals <- unique(per.fac$Value)
 valsnum <- unique(per.fac$Value[per.fac$Variable %in% names(num[num==TRUE])])
 valsnum.nonmiss <- sort(as.numeric(vals[vals %in% valsnum & substr(vals,1,4)!="miss"]))
 valsnum.miss <- sort(vals[vals %in% valsnum & substr(vals,1,4)=="miss"])
 vals[vals %in% valsnum] <- c(valsnum.nonmiss,valsnum.miss)
 per.fac$Value <- factor(as.character(per.fac$Value),levels=vals)

 # get different plots in order of data
 nperplot <- nrow*ncol
 nplots   <- ceiling(length(commonnames)/nperplot)
 plots    <- vector("list",nplots)
 tables   <- vector("list",nplots)
 
 for (i in 1:nplots) {
   min <-(i-1)*nperplot+1
   max <- min(length(commonnames),(i-1)*nperplot+nperplot)
 # tables
   ttables <- vector("list",(max-min+1))
   names(ttables) <- commonnames[min:max]
   for (j in commonnames[min:max]) {
     ttables[[j]] <- os.table[[j]]
   } 
   tables[[i]] <- ttables

  # plots
   per.fact <- per.fac[per.fac$Variable %in% commonnames[min:max],]
   # per.fact <- per.fact[order(match(per.fact$Variable,commonnames[min:max])),]
   p <- ggplot(data=per.fact, aes(x=Value,y=Percent,fill=Data))
   p <- p + geom_bar(position="dodge",colour=cols[1], stat="identity") + 
      facet_wrap(~ Variable, scales="free", ncol=ncol)  
   p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
        theme(axis.text.x=element_text(angle=-30, hjust=0, vjust=1,size=rel(rel.size.x)), 
              legend.position="top", 
              legend.key = element_rect(colour = cols[1])) 
   p <- p + theme(legend.title=element_blank())
   if (length(msel) > 1) p <- p + scale_fill_manual(values = c(cols[1],rep(cols[2],length(msel))))
   if (length(msel) <= 1) p <- p + scale_fill_manual(values = cols)
   plots[[i]] <- p
  }
 
 if (length(tables)==1) {
   tables <- tables[[1]]
   plots  <- plots[[1]]
 }  
 #browser()
 res <- list(tables = tables, plots = plots)
 class(res) <- "compare.synds"
 return(res)
}


###-----pertable-----------------------------------------------------------
# calculate percentages
pertable <- function(x) {
 if (any(is.na(x))) res <- table(x, exclude=NULL)*100/sum(table(x, exclude=NULL)) else
 res <- table(x)*100/sum(table(x, exclude=NULL))
 return(res)
}


###-----ggfac--------------------------------------------------------------
# calculate percentages for factors and store in a data frame (a long format)
ggfac <- function(data, name = "observed"){ 
  data <- as.data.frame(data)
  perlist  <- lapply(data, pertable)
  Percent  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),lapply(perlist,length))
  perdf    <- cbind.data.frame(Percent, Value, Variable)
  perdf$Data <- name
  return(list(perdf = perdf, perlist = perlist))
} 


###-----ggnum--------------------------------------------------------------
# calculate percentages for numeric variables (non-missing values and 
# missing data categories seperately) and store in a data frame (a long format)
ggnum <- function(data, name = "observed", na = as.list(rep(NA,ncol(data))), 
                  breaks = as.list(rep(30,ncol(data))), anyna){ 
  data <- as.data.frame(data)

# non-missing values  
  nvar <- ncol(data)
  perlist <- vector("list", nvar)
  hbreaks <- vector("list", nvar)
  
  for (i in 1:nvar){
    vardata <- data[!(data[,i] %in% na[[i]]),i]                  
    hh <- hist(vardata, breaks = breaks[[i]], plot=F)
    counts <- hh$counts
    names(counts) <- hh$breaks[-length(hh$breaks)]
    hbreaks[[i]] <- hh$breaks
    dataNA <- data[(data[,i] %in% na[[i]]),i] ### new bits in this
     NAcounts<-table(dataNA,exclude=NULL)
     names(NAcounts)<-paste("miss",names(NAcounts),sep=".") 
     counts<-c(counts,NAcounts)
     if (!anyna[i]) counts<-counts[!(names(counts)=="miss.NA")] ## drop NA category if not in data
     perlist[[i]]<-counts*100/length(data[,i])
  }
  names(perlist) <- colnames(data)

# create data frame in a long format  
  Percent  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),lapply(perlist,length))
  perdf    <- cbind.data.frame(Percent, Value, Variable)
  perdf$Data <- name
  
return(list(perdf = perdf, perlist = perlist, hbreaks = hbreaks))
}


###-----dfNA---------------------------------------------------------------
dfNA <- function(data, na){
  all  <- length(data)
  data <- data[data %in% unlist(na)]
  NA.counts <- table(data, exclude = NULL)*100/all
  return(NA.counts)
}


###-----compare.fit.synds--------------------------------------------------
compare.fit.synds <- function(object, data, plot = "Z", 
  return.result = TRUE, return.plot = TRUE, plot.intercept = FALSE, 
  lwd = 1, lty = 1, lcol = c("#1A3C5A","#4187BF"), 
  dodge.height = .5, point.size = 2.5, ...) {   # c("#132B43", "#56B1F7")
           
 # compares and plots fits to synthetic and real data
 # first parameter must be a fit to synthetic data from glm.synds() 
 # or lm.synds()
 value <- "Value"
 coefficient <- c("Coefficient", "Model")
 
 call <- match.call()
 # if (!class(object)=="fit.synds") stop("Object must have class fit.synds\n")
 if (!is.data.frame(data)) stop("Data must be a data frame\n")  # theoretically can be a matrix (?)

 syn.fit          <- summary.fit.synds(object) 
 fitting.function <- object$fitting.function
 
 # get fit to real data
 if (fitting.function=="lm") {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=object$call$formula,data=call$data)))
 } else {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=object$call$formula,
                     family=object$call$family,data=call$data)))
 }                     

 if (return.plot == TRUE){
   yvar <- as.character(object$call$formula[2])
  
   # prepare data for plotting confidence intervals (one data frame)
   if (plot=="Z"){
     BetaCI <- dfCI(real.fit, Z = TRUE)
     BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn", 
                 model.name = "synthetic")
     xlab = "Z value"            
     title = paste0("Z values for fit to ",yvar)
   } else {
     BetaCI <- dfCI(real.fit)
     BsynCI <- dfCI(syn.fit, names.est.se = c("B.syn","se(Beta).syn"),
                 model.name = "synthetic")
     xlab = "Value"            
     title = paste0("Coefficients for fit to ",yvar)
   }
   
   modelCI <- rbind.data.frame(BetaCI, BsynCI)
   rownames(modelCI) <- 1:nrow(modelCI)
  
   if(!plot.intercept) modelCI <- modelCI[modelCI$Coefficient!="(Intercept)",]
  
   CI.geom <- geom_errorbarh(aes_string(xmin = "LowCI", xmax = "HighCI",
     color = "Model", linetype = "Model"),lwd = lwd, lty = lty, height = 0,
     position = coefplot:::position_dodgev(height = dodge.height))   #! :::  
  
   point.geom <- geom_point(aes_string(xmin = value, xmax = value, 
     color = "Model", shape = "Model"), size = point.size, 
     position = coefplot:::position_dodgev(height = dodge.height))   #! :::   
  
   col.geom <- scale_colour_manual(values = lcol, breaks = c("synthetic","observed"))
   # col.geom <- scale_colour_manual(values = rev(brewer.pal(3,"Blues")))
   # col.geom <- scale_colour_grey(start = 0, end = .6)
  
   p <- ggplot(data=modelCI, aes_string(x = value, y = coefficient))  
   p <- p + geom_vline(xintercept=0, colour="grey", linetype=2, lwd=1)
   p <- p + CI.geom + point.geom + labs(title = title, x = xlab)
   p <- p + col.geom
   p <- p + scale_shape_manual(values=c(16:17),breaks = c("synthetic","observed"))
   #p <- p + theme_bw()
   p
 } else p <- NULL
 
 if (return.result == TRUE){ 
   res.obs <- real.fit$coefficients[,-4]
   colnames(res.obs) <- c("Beta","se(Beta)","Z")
   res.syn <- syn.fit$coefficients[,c("B.syn","se(Beta).syn","se(B.syn)","Z.syn","se(Z.syn)")]
   ##JS Addition, re-order synthetic results according to real results
   res.syn = res.syn[order(match(rownames(res.syn), rownames(res.obs) ) ), ]
 } else res.obs <- res.syn <- NULL
 
 res <- list(fit.synds.call = object$call, coef.obs = res.obs,
   coef.syn = res.syn, ci.plot = p)
 class(res) <- "compare.fit.synds"
 return(res)
}


###-----dfCI---------------------------------------------------------------
# extract info for plotting confidence intervals
dfCI <- function(modelsummary, names.est.se = c("Estimate","Std. Error"),
          model.name = "observed", CI = 1.96, Z = FALSE, 
          name.Z = colnames(modelsummary$coefficients)[3]){
  
  if (!Z) {
    msCI <- as.data.frame(modelsummary$coefficients[,names.est.se])
    names(msCI) <- c("Value", "SE")
  } else {
    msCI <- as.data.frame(modelsummary$coefficients[,name.Z])
    names(msCI) <- c("Value")
    msCI$SE <- 1
  }  
  msCI$Coefficient <- rownames(msCI)

  msCI$HighCI <- msCI$Value + CI*msCI$SE
  msCI$LowCI  <- msCI$Value - CI*msCI$SE
  msCI$SE     <- NULL
  msCI$Model  <- model.name
  msCI$Coefficient <- factor(msCI$Coefficient, levels = msCI$Coefficient)

  return(msCI)
}
