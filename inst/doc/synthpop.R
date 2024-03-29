### R code from vignette source 'synthpop.Rnw'

###################################################
### code chunk number 1: synthpop.Rnw:64-65
###################################################
options(prompt="R> ", width=77, digits=4, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: synthpop.Rnw:168-169
###################################################
library("synthpop")


###################################################
### code chunk number 3: ods
###################################################
vars <- c("sex", "age", "edu", "marital", "income", "ls", "wkabint")
ods <- SD2011[, vars]
head(ods)


###################################################
### code chunk number 4: synthpop.Rnw:200-202 (eval = FALSE)
###################################################
## my.seed <- 17914709
## sds.default <- syn(ods, seed = my.seed)


###################################################
### code chunk number 5: synthpop.Rnw:207-208 (eval = FALSE)
###################################################
## sds.default


###################################################
### code chunk number 6: synthpop.Rnw:259-260 (eval = FALSE)
###################################################
## names(sds.default)


###################################################
### code chunk number 7: synthpop.Rnw:280-281 (eval = FALSE)
###################################################
## sds.parametric <- syn(ods, method = "parametric", seed = my.seed)


###################################################
### code chunk number 8: synthpop.Rnw:283-284 (eval = FALSE)
###################################################
## sds.parametric$method


###################################################
### code chunk number 9: synthpop.Rnw:300-302 (eval = FALSE)
###################################################
## sds.selection <- syn(ods, visit.sequence = c(1, 2, 6, 4, 3), 
##   seed = my.seed, drop.not.used = TRUE)


###################################################
### code chunk number 10: synthpop.Rnw:307-308 (eval = FALSE)
###################################################
## sds.selection


###################################################
### code chunk number 11: synthpop.Rnw:361-365
###################################################
visit.sequence.ini <- c(1, 2, 5, 6, 4, 3)
method.ini <- c("sample", "ctree", "ctree", "polyreg", "", "ctree", "")
sds.ini <- syn(data = ods, visit.sequence = visit.sequence.ini,
  method = method.ini, m = 0, drop.not.used = FALSE)


###################################################
### code chunk number 12: synthpop.Rnw:367-371
###################################################
sds.ini$predictor.matrix
predictor.matrix.corrected <- sds.ini$predictor.matrix
predictor.matrix.corrected["marital", "ls"] <- 0
predictor.matrix.corrected


###################################################
### code chunk number 13: synthpop.Rnw:373-376 (eval = FALSE)
###################################################
## sds.corrected <- syn(data = ods, visit.sequence = visit.sequence.ini,
##   method = method.ini, predictor.matrix = predictor.matrix.corrected,
##   seed = my.seed)


###################################################
### code chunk number 14: synthpop.Rnw:382-384 (eval = FALSE)
###################################################
## sds.income <- syn(ods, cont.na = list(income = c(NA, -8)), 
##   smoothing = list(income = "density"), seed = NA)


###################################################
### code chunk number 15: synthpop.Rnw:390-398 (eval = FALSE)
###################################################
## M18.ods <- table(subset(ods,
##   age < 18 & sex == "MALE", marital))
## M18.default <- table(subset(sds.default$syn,
##   age < 18 & sex == "MALE", marital))
## M18.parametric <- table(subset(sds.parametric$syn,
##   age < 18 & sex == "MALE", marital))
## cbind("Observed data" = M18.ods, CART = M18.default,
##   Parametric = M18.parametric)


###################################################
### code chunk number 16: synthpop.Rnw:411-417 (eval = FALSE)
###################################################
## rules.marital <- list(marital = "age < 18 & sex == 'MALE'")
## rvalues.marital <- list(marital = "SINGLE")
## sds.rmarital <- syn(ods, rules = rules.marital,
##   rvalues = rvalues.marital, seed = my.seed)
## sds.rmarital.param <- syn(ods, rules = rules.marital,
##   rvalues = rvalues.marital, method = "parametric", seed = my.seed)


###################################################
### code chunk number 17: synthpop.Rnw:420-426 (eval = FALSE)
###################################################
## rM18.default <- table(subset(sds.rmarital$syn,
##   age < 18 & sex == "MALE", marital))
## rM18.parametric <- table(subset(sds.rmarital.param$syn,
##   age < 18 & sex == "MALE", marital))
## cbind("Observed data" = M18.ods, CART = rM18.default,
##   Parametric = rM18.parametric)


###################################################
### code chunk number 18: synthpop.Rnw:443-448
###################################################
ods$wkabint <- as.character(ods$wkabint)
ods$wkabint[ods$wkabint == "YES, TO EU COUNTRY" |
  ods$wkabint == "YES, TO NON-EU COUNTRY"] <- "YES"
ods$wkabint <- factor(ods$wkabint) 
ods$income[ods$income == -8] <- NA


###################################################
### code chunk number 19: synthpop.Rnw:453-454 (eval = FALSE)
###################################################
## sds <- syn(ods, method = "ctree", m = 5, seed = my.seed)


###################################################
### code chunk number 20: synthpop.Rnw:459-460
###################################################
summary(ods)


###################################################
### code chunk number 21: synthpop.Rnw:465-466 (eval = FALSE)
###################################################
## summary(sds)


###################################################
### code chunk number 22: synthpop.Rnw:498-500 (eval = FALSE)
###################################################
## summary(sds, msel = 2)
## summary(sds, msel = 1:5)


###################################################
### code chunk number 23: synthpop.Rnw:505-506 (eval = FALSE)
###################################################
## compare(sds, ods, vars = "income")  


###################################################
### code chunk number 24: synthpop.Rnw:525-526 (eval = FALSE)
###################################################
## compare(sds, ods, vars = "ls", msel = 1:3)


###################################################
### code chunk number 25: synthpop.Rnw:545-548
###################################################
model.ods <- glm(wkabint ~ sex + age + edu + log(income), 
  family = "binomial", data = ods) 
model.ods


###################################################
### code chunk number 26: synthpop.Rnw:550-553 (eval = FALSE)
###################################################
## model.sds <- glm.synds(wkabint ~ sex + age + edu + log(income), 
##   family = "binomial", data = sds) 
## model.sds


###################################################
### code chunk number 27: synthpop.Rnw:578-579 (eval = FALSE)
###################################################
## summary(model.sds)


###################################################
### code chunk number 28: synthpop.Rnw:612-613 (eval = FALSE)
###################################################
## compare(model.sds, ods)


