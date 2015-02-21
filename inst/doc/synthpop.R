### R code from vignette source 'synthpop.Rnw'

###################################################
### code chunk number 1: synthpop.Rnw:61-62
###################################################
options(prompt="R> ", width=77, digits=4, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: synthpop.Rnw:154-155
###################################################
library("synthpop")


###################################################
### code chunk number 3: ods
###################################################
vars <- c("sex", "age", "edu", "marital", "income", "ls", "wkabint")
ods <- SD2011[, vars]
head(ods)


###################################################
### code chunk number 4: synthpop.Rnw:186-188
###################################################
my.seed <- 17914709
sds.default <- syn(ods, seed = my.seed)


###################################################
### code chunk number 5: synthpop.Rnw:193-194
###################################################
sds.default


###################################################
### code chunk number 6: synthpop.Rnw:199-200
###################################################
names(sds.default)


###################################################
### code chunk number 7: synthpop.Rnw:209-210
###################################################
sds.parametric <- syn(ods, method = "parametric", seed = my.seed)


###################################################
### code chunk number 8: synthpop.Rnw:212-213
###################################################
sds.parametric$method


###################################################
### code chunk number 9: synthpop.Rnw:223-224
###################################################
sds.selection <- syn(ods, visit.sequence = c(1, 2, 6, 4, 3), seed = my.seed)


###################################################
### code chunk number 10: synthpop.Rnw:229-230
###################################################
sds.selection


###################################################
### code chunk number 11: synthpop.Rnw:245-249
###################################################
visit.sequence.ini <- c(1, 2, 5, 6, 4, 3)
method.ini <- c("sample", "ctree", "ctree", "polyreg", "", "ctree", "")
sds.ini <- syn(data = ods, visit.sequence = visit.sequence.ini,
  method = method.ini, m = 0, drop.not.used = FALSE)


###################################################
### code chunk number 12: synthpop.Rnw:251-255
###################################################
sds.ini$predictor.matrix
predictor.matrix.corrected <- sds.ini$predictor.matrix
predictor.matrix.corrected["marital", "ls"] <- 0
predictor.matrix.corrected


###################################################
### code chunk number 13: synthpop.Rnw:257-260
###################################################
sds.corrected <- syn(data = ods, visit.sequence = visit.sequence.ini,
  method = method.ini, predictor.matrix = predictor.matrix.corrected,
  seed = my.seed)


###################################################
### code chunk number 14: synthpop.Rnw:267-269
###################################################
cont.na.income <- as.list(rep(NA, ncol(ods)))
cont.na.income[[5]] <- c(NA, -8)


###################################################
### code chunk number 15: synthpop.Rnw:275-282
###################################################
maritalM18.ods <- table(ods[ods$age < 18 & ods$sex == 'MALE', "marital"])
maritalM18.default <- table(sds.default$syn[sds.default$syn$age < 18 &
  sds.default$syn$sex == 'MALE', "marital"])
maritalM18.parametric <- table(sds.parametric$syn[sds.default$syn$age < 18 &
  sds.parametric$syn$sex == 'MALE', "marital"])
cbind("Observed data" = maritalM18.ods, CART = maritalM18.default,
  Parametric = maritalM18.parametric)


###################################################
### code chunk number 16: synthpop.Rnw:286-298
###################################################
rules.marital <- list("", "", "", "age < 18 & sex == 'MALE'", "", "", "")
rvalues.marital <- list(NA, NA, NA, 'SINGLE', NA, NA, NA)
sds.rmarital <- syn(ods, rules = rules.marital, 
  rvalues = rvalues.marital, seed = my.seed)
sds.rmarital.param <- syn(ods, rules = rules.marital, 
  rvalues = rvalues.marital, method = "parametric", seed = my.seed)

rmaritalM18.default <- table(sds.rmarital$syn[sds.rmarital$syn$age < 18
  & sds.rmarital$syn$sex == 'MALE', "marital"])
rmaritalM18.parametric <- table(sds.rmarital.param$syn[
  sds.rmarital.param$syn$age < 18
  & sds.rmarital.param$syn$sex == 'MALE', "marital"])


###################################################
### code chunk number 17: synthpop.Rnw:300-302
###################################################
cbind("Observed data" = maritalM18.ods, CART = rmaritalM18.default,
  Parametric = rmaritalM18.parametric) 


###################################################
### code chunk number 18: synthpop.Rnw:310-315
###################################################
ods$wkabint <- as.character(ods$wkabint)
ods$wkabint[ods$wkabint == 'YES, TO EU COUNTRY' |
  ods$wkabint == 'YES, TO NON-EU COUNTRY'] <- 'YES'
ods$wkabint <- factor(ods$wkabint) 
ods$income[ods$income == -8] <- NA


###################################################
### code chunk number 19: synthpop.Rnw:320-321
###################################################
sds <- syn(ods, m = 5, seed = my.seed)


###################################################
### code chunk number 20: synthpop.Rnw:326-327
###################################################
summary(ods)


###################################################
### code chunk number 21: synthpop.Rnw:332-333
###################################################
summary(sds)


###################################################
### code chunk number 22: synthpop.Rnw:335-337
###################################################
summary(sds, msel = 2)
summary(sds, msel = 1:5)


###################################################
### code chunk number 23: synthpop.Rnw:342-343 (eval = FALSE)
###################################################
## compare.synds(sds, ods)


###################################################
### code chunk number 24: synthpop.Rnw:345-346
###################################################
compare.synds(sds, ods, vars = "ls")


###################################################
### code chunk number 25: synthpop.Rnw:355-356
###################################################
compare.synds(sds, ods, vars = "income")  


###################################################
### code chunk number 26: synthpop.Rnw:366-373
###################################################
model.ods <- glm(wkabint ~ sex + age + edu + log(income), 
  family = "binomial", data = ods) 
summary(model.ods)
                             
model.sds <- glm.synds(wkabint ~ sex + age + edu + log(income), 
  family = "binomial", object = sds) 
model.sds


###################################################
### code chunk number 27: synthpop.Rnw:376-377
###################################################
print(model.sds, msel = 3)


###################################################
### code chunk number 28: synthpop.Rnw:384-385
###################################################
summary(model.sds)


###################################################
### code chunk number 29: synthpop.Rnw:390-391
###################################################
compare.fit.synds(model.sds, ods)


