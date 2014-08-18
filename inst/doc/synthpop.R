### R code from vignette source 'synthpop.Rnw'

###################################################
### code chunk number 1: synthpop.Rnw:62-63
###################################################
options(width=77,digits=4,useFancyQuotes=FALSE)


###################################################
### code chunk number 2: synthpop.Rnw:164-165
###################################################
library(synthpop)


###################################################
### code chunk number 3: rds
###################################################
vars <- c("sex","age","edu","marital","income","ls","wkabint")
rds <- SD2011[,vars]
head(rds)


###################################################
### code chunk number 4: synthpop.Rnw:195-197
###################################################
set.seed(17914709)
sds.default <- syn(rds)


###################################################
### code chunk number 5: synthpop.Rnw:202-203
###################################################
sds.default


###################################################
### code chunk number 6: synthpop.Rnw:208-209
###################################################
names(sds.default)


###################################################
### code chunk number 7: synthpop.Rnw:218-220
###################################################
set.seed(17914709)
sds.parametric <- syn(rds, method = "parametric")


###################################################
### code chunk number 8: synthpop.Rnw:222-223
###################################################
sds.parametric$method


###################################################
### code chunk number 9: synthpop.Rnw:233-235
###################################################
set.seed(17914709)
sds.selection <- syn(rds, visitSequence = c(1, 2, 6, 4, 3))


###################################################
### code chunk number 10: synthpop.Rnw:240-241
###################################################
sds.selection


###################################################
### code chunk number 11: synthpop.Rnw:256-261
###################################################
visitSequence.ini <- c(1, 2, 5, 6, 4, 3)
method.ini <- c("sample", "ctree", "ctree", "polyreg", "", "ctree", "")
set.seed(17914709)
sds.ini <- syn(data = rds, visitSequence = visitSequence.ini,
  method = method.ini, m = 0, drop.not.used = FALSE)


###################################################
### code chunk number 12: synthpop.Rnw:263-267
###################################################
sds.ini$predictorMatrix
predictorMatrix.corrected <- sds.ini$predictorMatrix
predictorMatrix.corrected["marital","ls"] <- 0
predictorMatrix.corrected


###################################################
### code chunk number 13: synthpop.Rnw:269-272
###################################################
set.seed(17914709)
sds.corrected <- syn(data = rds, visitSequence = visitSequence.ini,
  method = method.ini, predictorMatrix = predictorMatrix.corrected)


###################################################
### code chunk number 14: synthpop.Rnw:279-281
###################################################
contNA.income <- as.list(rep(NA, ncol(rds)))
contNA.income[[5]] <- c(NA,-8)


###################################################
### code chunk number 15: synthpop.Rnw:287-294
###################################################
maritalM18.rds <- table(rds[rds$age < 18 & rds$sex == 'MALE' ,"marital"])
maritalM18.default <- table(sds.default$syn[sds.default$syn$age < 18 &
  sds.default$syn$sex == 'MALE',"marital"])
maritalM18.parametric <- table(sds.parametric$syn[sds.default$syn$age < 18 &
  sds.parametric$syn$sex == 'MALE',"marital"])
cbind("Real data" = maritalM18.rds, CART = maritalM18.default,
  Parametric = maritalM18.parametric)


###################################################
### code chunk number 16: synthpop.Rnw:298-311
###################################################
rules.marital <- list("","","","age < 18 & sex == 'MALE'","","","")
rvalues.marital <- list(NA,NA,NA,'SINGLE',NA,NA,NA)
set.seed(17914709)
sds.rmarital <- syn(rds, rules = rules.marital, rvalues = rvalues.marital)
set.seed(17914709)
sds.rmarital.param <- syn(rds, rules = rules.marital, 
  rvalues = rvalues.marital, method = "parametric")

rmaritalM18.default <- table(sds.rmarital$syn[sds.rmarital$syn$age < 18
  & sds.rmarital$syn$sex == 'MALE',"marital"])
rmaritalM18.parametric <- table(sds.rmarital.param$syn[
  sds.rmarital.param$syn$age < 18
  & sds.rmarital.param$syn$sex == 'MALE',"marital"])


###################################################
### code chunk number 17: synthpop.Rnw:313-315
###################################################
cbind("Real data" = maritalM18.rds, CART = rmaritalM18.default,
  Parametric = rmaritalM18.parametric) 


###################################################
### code chunk number 18: synthpop.Rnw:323-328
###################################################
rds$wkabint <- as.character(rds$wkabint)
rds$wkabint[rds$wkabint=='YES, TO EU COUNTRY' |
  rds$wkabint=='YES, TO NON-EU COUNTRY'] <- 'YES'
rds$wkabint <- factor(rds$wkabint) 
rds$income[rds$income==-8] <- NA


###################################################
### code chunk number 19: synthpop.Rnw:333-335
###################################################
set.seed(17914709)
sds <- syn(rds, m = 5)


###################################################
### code chunk number 20: synthpop.Rnw:340-341
###################################################
summary(rds)


###################################################
### code chunk number 21: synthpop.Rnw:346-347
###################################################
summary(sds)


###################################################
### code chunk number 22: synthpop.Rnw:349-351
###################################################
summary(sds, msel = 2)
summary(sds, msel = 1:5)


###################################################
### code chunk number 23: synthpop.Rnw:356-357 (eval = FALSE)
###################################################
## compare.synds(sds,rds)


###################################################
### code chunk number 24: synthpop.Rnw:359-360
###################################################
compare.synds(sds,rds,vars="ls")


###################################################
### code chunk number 25: synthpop.Rnw:368-369
###################################################
compare.synds(sds,rds,vars="income")  


###################################################
### code chunk number 26: synthpop.Rnw:385-392
###################################################
model.rds <- glm(wkabint ~ sex + age + edu + log(income), 
  family = "binomial", data = rds) 
summary(model.rds)
                             
model.sds <- glm.synds(wkabint ~ sex + age + edu + log(income), 
  family = "binomial", object = sds) 
model.sds


###################################################
### code chunk number 27: synthpop.Rnw:395-396
###################################################
print(model.sds, msel = 3)


###################################################
### code chunk number 28: synthpop.Rnw:401-402
###################################################
summary(model.sds)


###################################################
### code chunk number 29: synthpop.Rnw:407-408
###################################################
compare.fit.synds(model.sds,rds)


