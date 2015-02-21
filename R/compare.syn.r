###-----compare.synds------------------------------------------------------
compare.synds <- function(object, data, vars = NULL, ncol = 4, breaks = 20, 
                   plot.na = TRUE, rel.size.x = 1, 
                   cols = c("#1A3C5A","#4187BF"), ...){ 

 if (is.null(data)) stop("Requires parameter data to give name of the real data.\n")
 if (!class(object)=="synds") stop("Object must have class synds.\n")
 if (object$m==1) synds <- object$syn else synds  <- object$syn[[1]]
 synnames    <- names(synds)
 realnames   <- names(data)
 commonnames <- synnames[match(realnames,synnames)]
 if (!is.null(vars)) commonnames <- commonnames[match(vars,commonnames)]
 commonnames <- commonnames[!is.na(commonnames)]
 if ((length(commonnames) == 0) && (typeof(commonnames) == "character")) stop("Variables selected for comparison do not exist.", call. = FALSE)
 
 df.obs <- data[, commonnames, drop=FALSE]
 df.syn <- synds[, commonnames, drop=FALSE]
 num <- sapply(df.obs, is.numeric) | sapply(df.obs, is.integer)
 fac <- sapply(df.obs, is.factor) 

 # frequency tables for factors
 if (sum(fac)>0) {
   per.obs.fac <- ggfac(df.obs[,fac,drop=FALSE])
   per.syn.fac <- ggfac(df.syn[,fac,drop=FALSE], name = "synthetic")
 } else {
   per.obs.fac <- per.syn.fac <- NULL
 }
 
 # frequency tables for numeric variables
 if (sum(num)>0) {
   na.index <- match(colnames(df.obs[,num,drop=FALSE]),colnames(data))
   na <- object$cont.na[na.index] 
   lbreaks  <- as.list(rep(breaks,length(na)))
   per.obs.num <- ggnum(df.obs[,num,drop=FALSE], 
                    na = na, breaks = lbreaks, plot.na = plot.na)
   per.syn.num <- ggnum(df.syn[,num,drop=FALSE], name = "synthetic", 
                    na = na, breaks = per.obs.num$hbreaks, plot.na = plot.na) 
 } else { 
   per.obs.num <- per.syn.num <- NULL
   plot.na = FALSE
 }
                              
 # data frame for ploting 
 per.fac <- rbind.data.frame(per.obs.fac$perdf, per.syn.fac$perdf, 
                             per.obs.num$perdf, per.syn.num$perdf)
 per.fac$value <- factor(per.fac$value, levels = unique(per.fac$value), 
                         ordered = T, exclude = NULL)
 per.fac <- per.fac[per.fac$percent!=0, ]

 # list of result tables
 os.table.fac <- mapply(rbind, observed = per.obs.fac$perlist, 
                   synthetic = per.syn.fac$perlist, SIMPLIFY = FALSE)
 os.table.num <- mapply(rbind, observed = per.obs.num$perlist, 
                   synthetic = per.syn.num$perlist, SIMPLIFY = FALSE)
 os.table <- c(os.table.fac,os.table.num)

 value <- percent <- NULL

 p <- ggplot(na.omit(per.fac), aes(x=value,y=percent,fill=data))
 p <- p + geom_bar(position="dodge", stat="identity") + 
      facet_wrap(~ variable, scales="free", ncol=ncol)  
 p <- p + theme(axis.text.x=element_text(angle=-30, hjust=0, 
      vjust=1,size=rel(rel.size.x)), legend.position="top") 
 p <- p + theme(legend.title=element_blank())
 p <- p + scale_fill_manual(values = cols)
 p
 res <- list(freq.table = os.table, p = p, plot.na = plot.na)
 class(res) <- "compare.synds"
 return(res)
}


###-----pertable-----------------------------------------------------------
# calculate percentages
pertable <- function(x) table(x, exclude=NULL)*100/sum(table(x, exclude=NULL))


###-----ggfac--------------------------------------------------------------
# calculate percentages for factors and store in a data frame (a long format)
ggfac <- function(data, name = "observed"){ 
  data <- as.data.frame(data)
  perlist  <- lapply(data, pertable)
  percent  <- unlist(perlist, use.names = F)
  value    <- unlist(lapply(perlist,names), use.names = F)
  variable <- rep(names(perlist),lapply(perlist,length))
  perdf    <- cbind.data.frame(percent, value, variable)
  perdf$data <- name
  return(list(perdf = perdf, perlist = perlist))
} 


###-----ggnum--------------------------------------------------------------
# calculate percentages for numeric variables (non-missing values and 
# missing data categories seperately) and store in a data frame (a long format)
ggnum <- function(data, name = "observed", na = as.list(rep(NA,ncol(data))), 
                  breaks = as.list(rep(30,ncol(data))), plot.na = TRUE){ 
  data <- as.data.frame(data)

# non-missing values  
  nvar <- ncol(data)
  datanonNA <- vector("list", nvar)
  hbreaks <- vector("list", nvar)
  names(datanonNA) <- colnames(data)
  for (i in 1:nvar){
    vardata <- data[!(data[,i] %in% na[[i]]),i]                  
    hh <- hist(vardata, breaks = breaks[[i]], plot=F)
    datanonNA[[i]] <- hh$counts*100/length(vardata)
    names(datanonNA[[i]]) <- hh$breaks[-length(hh$breaks)]
    hbreaks[[i]] <- hh$breaks
  }

# calculate % for missing values if they are to be ploted  
  if(plot.na){
    dataNA  <- lapply(data,dfNA,na=na)     
    names(dataNA) <- paste0(names(dataNA),".missing")
    mth <- c(t(matrix(1:(2*nvar),nrow=nvar)))
    perlist  <- c(datanonNA,dataNA)[mth]
  } else {
    perlist  <- c(datanonNA)
  }

# create data frame in a long format  
  percent  <- unlist(perlist, use.names = F)
  value    <- unlist(lapply(perlist,names), use.names = F)
  variable <- rep(names(perlist),lapply(perlist,length))
  perdf    <- cbind.data.frame(percent, value, variable)
  perdf$data <- name
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
compare.fit.synds <- function(object, obs.data, plot = "Z", 
  return.result = TRUE, plot.intercept = FALSE, lwd = 1, lty = 1, 
  lcol = c("#1A3C5A","#4187BF"), dodge.height = .5, point.size = 2.5, ...) {   # c("#132B43", "#56B1F7")
           
 # compares and plots fits to synthetic and real data
 # first parameter must be a fit to synthetic data from glm.synds() 
 # or lm.synds()
 value <- "Value"
 coefficient <- c("Coefficient", "Model")
 
 call <- match.call()
 if (!class(object)=="fit.synds") stop("Object must have class fit.synds\n")
 if (!is.data.frame(obs.data)) stop("obs.data must be a data frame\n")  # theoretically can be a matrix

 syn.fit          <- summary.fit.synds(object) 
 fitting.function <- object$fitting.function
 
 # get fit to real data
 if (fitting.function=="lm") {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=object$call$formula,data=call$obs.data)))
 } else {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=object$call$formula,
                     family=object$call$family,data=call$obs.data)))
 }                     

 yvar <- as.character(object$call$formula[2])

 # prepare data for plotting confidence intervals (one data frame)
 if (plot=="Z"){
   BetaCI <- dfCI(real.fit, Z = TRUE)
   BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn", 
               model.name = "Synthetic")
   xlab = "Z value"            
   title = paste0("Z values for fit to ",yvar)
 } else {
   BetaCI <- dfCI(real.fit)
   BsynCI <- dfCI(syn.fit, names.est.se = c("B.syn","se(Beta).syn"),
               model.name = "Synthetic")
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

 col.geom <- scale_colour_manual(values = lcol)
 # col.geom <- scale_colour_manual(values = rev(brewer.pal(3,"Blues")))
 # col.geom <- scale_colour_grey(start = 0, end = .6)

 p <- ggplot(data=modelCI, aes_string(x = value, y = coefficient))  
 p <- p + geom_vline(xintercept=0, colour="grey", linetype=2, lwd=1)
 p <- p + CI.geom + point.geom + labs(title = title, x = xlab)
 p <- p + col.geom
 #p <- p + theme_bw()
 p
 
 res.obs <- real.fit$coefficients[,-4]
 colnames(res.obs) <- c("Beta","se(Beta)","Z")
 res.syn <- syn.fit$coefficients[,c("B.syn","se(Beta).syn","se(B.syn)","Z.syn","se(Z.syn)")]

 res <- list(fit.synds.call = object$call, coef.obs = res.obs,
   coef.syn = res.syn, ci.plot = p)
 class(res) <- "compare.fit.synds"
 if (return.result==TRUE) return(res) else return(p)

}


###-----dfCI---------------------------------------------------------------
# extract info for plotting confidence intervals
dfCI <- function(modelsummary, names.est.se = c("Estimate","Std. Error"),
          model.name = "Observed", CI = 1.96, Z = FALSE, 
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