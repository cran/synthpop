### R code from vignette source 'utility.Rnw'

###################################################
### code chunk number 1: utility.Rnw:61-62
###################################################
options(prompt="R> ", width=77, digits=4, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: utility.Rnw:81-95 (eval = FALSE)
###################################################
## library("synthpop")
## ods <- SD2011[, c("sex", "income", "age", "edu" , "socprof", "trust",
##   "height", "weight", "smoke", "region")]
## syn_para <- syn(ods, method = "parametric", cont.na = list(income = -8),
##   seed = 34567, print.flag = FALSE)
## syn_cart <- syn(ods, method = "cart", cont.na = list(income = -8),
##   seed = 34567, print.flag = FALSE)
## u_ods_para <- utility.gen(syn_para, ods, method = "cart",
##   resamp.method = "none", print.flag = FALSE)
## u_ods_cart <- utility.gen(syn_cart, ods, method = "cart",
##   resamp.method = "none", print.flag = FALSE)
## cat("\nParametric pMSE ", u_ods_para$pMSE,
##   "\nCART pMSE ", u_ods_cart$pMSE,
##   "\nUtility ratio parametric to CART: ", u_ods_para$pMSE/u_ods_cart$pMSE)


###################################################
### code chunk number 3: utility.Rnw:111-112 (eval = FALSE)
###################################################
## utility.tables(syn_para, ods, tables = "twoway", nworst = 4)


###################################################
### code chunk number 4: utility.Rnw:145-154 (eval = FALSE)
###################################################
## u_3way <- utility.tables(syn_para, ods, tab.stats = "all",
##   tables = "threeway")
## cors <- cor(u_3way$tabs)
## cat("Correlations:\nVW with pMSE = ", cors["VW", "pMSE"],
##   ", SPECKS with MabsDD = ", cors["SPECKS", "MabsDD"],
##   ", and SPECKS with PO50 = ", cors["SPECKS", "MabsDD"], "." , sep = "")     
## toplot <- u_3way$tabs[, c(1:4, 6, 5, 7)]
##   dimnames(toplot)[[2]][c(1, 4)] <- c("VW\npMSE", "SPECKS\nMabsDD\nPO50")
## pairs(toplot)   


###################################################
### code chunk number 5: utility.Rnw:232-233 (eval = FALSE)
###################################################
## compare(syn_para, ods, utility.stats = c("S_pMSE", "df"))


###################################################
### code chunk number 6: utility.Rnw:270-274 (eval = FALSE)
###################################################
## syn_para2 <- syn(ods, method = "parametric", cont.na = list(income = -8),       
##   visit.sequence = c(1, 3, 7:9, 2, 4:6, 10), seed = 34567,
##   print.flag = FALSE)
## compare(syn_para2, ods, utility.stats = c("S_pMSE", "df"), plot = FALSE)


###################################################
### code chunk number 7: utility.Rnw:303-321 (eval = FALSE)
###################################################
## syn_para3 <- syn.strata(ods, method = "parametric", 
##   strata = ods$age > 55 & !is.na(ods$age), cont.na = list(income = -8), 
##   visit.sequence = c(1, 3, 7:9, 2, 4:6, 10), seed = 34567, 
##   print.flag = FALSE) 
## 
## u.para <- utility.tables(syn_para, ods, max.scale = 31, 
##   plot.title = "(a) parametric synthesis\n")
## u.para2 <- utility.tables(syn_para2, ods, max.scale = 31, 
##   plot.title = "(b) reordered parametric synthesis\n")
## u.para3 <- utility.tables(syn_para3, ods, max.scale = 31, 
##   plot.title = "(c) reordered and age startified\nparametric synthesis")
## u.cart <- utility.tables(syn_cart, ods, max.scale = 31, 
##   plot.title = "(d) CART synthesis\n")
## 
## list.plots <- list(u.para$utility.plot, u.para3$utility.plot,  
##   u.para2$utility.plot, u.cart$utility.plot)
## gridExtra::marrangeGrob(list.plots, nrow = 2, ncol = 2, 
##   top = "Two-way pMSE ratios")


###################################################
### code chunk number 8: utility.Rnw:432-462 (eval = FALSE)
###################################################
## m <- 1000
## ods_cat <- numtocat.syn(ods, cont.na = list(income = -8), 
##   style.groups = "quantile", catgroups = 5)$data
## syn_bad <- syn(ods_cat, method = "sample",  print.flag = FALSE, m = m, 
##   seed = 12345)
## 
## calc_power <- function(nvars, m = m) {
##   syn_good <- syn(ods_cat[, 1:nvars], method = "catall", 
##     print.flag = FALSE, m = m, seed = 12345)
##   u_good <- utility.tab(syn_good, ods_cat,
##     vars = names(syn_good$syn[[1]]), print.flag = FALSE) 
##   u_bad <- utility.tab(syn_bad, ods_cat, 
##     vars = names(syn_bad$syn[[1]]), print.flag = FALSE) 
##   results_bad <- with(u_bad, data.frame(pMSE, FT, JSD, SPECKS, U, 
##     WMabsDD, G, df, dfG))
##   results_good <- with(u_good, data.frame(pMSE, FT, JSD, SPECKS, U,
##     WMabsDD, G, df, dfG))
##   power <- (apply(results_bad[, 1:7], 2, "mean") - 
##     apply(results_good[, 1:7], 2, "mean"))/sqrt(apply(results_good[,1:7], 2, "var"))
##   list(power = c(power, df_good = median(results_good[, 8]), 
##          dfG_good = median(results_good[, 9]), 
##          df_bad = median(results_bad[, 8]), 
##          dfG_bad = median(results_bad[, 9])),
##        expect = apply(with(u_good, data.frame(S_pMSE, S_FT, 
##          S_JSD, S_WMabsDD, S_G)), 2 ,"mean"))
##   }
## pe <- as.list(2:6)
## for (i in 2:6) {pe[[i - 1]] <- calc_power(i, m)}
## table1 <- t(sapply(pe, '[[', "power"))
## table5 <- t(sapply(pe, '[[', "expect"))


###################################################
### code chunk number 9: utility.Rnw:511-529 (eval = FALSE)
###################################################
## calc_power2 <- function(nvars) {
##   syn_good1 <- syn(ods_cat[, 1:nvars], method = "catall", m = 1000, 
##     print.flag = FALSE, seed = 12345)
##   u_good1 <- utility.gen(syn_good1, ods_cat[, 1:nvars],  
##     vars = names(syn_good1$syn[[1]]), resamp.method = "perm", 
##     utility.stats = "S_pMSE", print.flag = FALSE) 
##   syn_good2 <- syn(ods_cat[, 1:nvars], method = "catall", m = 46, 
##     print.flag = FALSE, seed = 12345)
##   u_good2 <- utility.gen(syn_good2, ods_cat, 
##     vars = names(syn_good2$syn[[1]]), resamp.method = "pairs",
##     utility.stats = c("S_pMSE", "S_SPECKS", "S_U"), print.flag = FALSE)
##   table6 <- c(S_pMSE1 = table5[nvars - 1, 1], 
##     S_pMSE2 = mean(u_good1$S_pMSE), S_pMSE3 = mean(u_good2$S_pMSE), 
##     S_SPECKS = mean(u_good2$S_SPECKS), S_U = mean(u_good2$S_U)) 
##   }
## pe <- as.list(2:6)
## for (i in 2:6) {pe[[i - 1]] <- calc_power2(i)}
## table6 <- do.call(rbind, pe)


