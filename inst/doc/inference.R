### R code from vignette source 'inference.Rnw'

###################################################
### code chunk number 1: inference.Rnw:48-49
###################################################
options(prompt="R> ", width=77, digits=4, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: inference.Rnw:109-117
###################################################
library(synthpop)
ods <- SD2011[, c("smoke", "sex", "age", "edu")]
levels(ods$edu) <- c("NONE", "VOC", "SEC", "HIGH")
s1 <- syn(ods, seed = 1234)
summary(glm(smoke ~ sex + age + edu + sex * edu, 
  data = s1$syn, family = "binomial"))
summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
  data = s1, family = "binomial"))


###################################################
### code chunk number 3: inference.Rnw:122-125
###################################################
s2 <- syn(ods, seed = 1234,visit.sequence=c("smoke","edu","age"))
summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
data = s2, family = "binomial"))


###################################################
### code chunk number 4: inference.Rnw:134-142
###################################################
s3 <- syn(ods, seed = 1234, k = 500)

## analysing synthetic data with just glm()
summary(glm(smoke ~ sex + age + edu + sex * edu, 
            data = s3$syn, family = "binomial"))
## using glm.synds()
summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
                  data = s3, family = "binomial"))


###################################################
### code chunk number 5: inference.Rnw:153-155
###################################################
summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
  data = s1, family = "binomial"), population.inference = TRUE)


###################################################
### code chunk number 6: inference.Rnw:163-166
###################################################
s4 <- syn(ods, seed = 5678, proper = TRUE)
summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
  data = s4, family = "binomial"), population.inference = TRUE)


###################################################
### code chunk number 7: inference.Rnw:173-176 (eval = FALSE)
###################################################
## s5 <- syn(ods, seed = 5678, m = 10, proper = TRUE, print.flag = FALSE)
## summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
##   data = s5, family = "binomial"), population.inference = TRUE)


###################################################
### code chunk number 8: inference.Rnw:209-214 (eval = FALSE)
###################################################
## s6 <- syn(ods, seed = 910011, m = 12, 
##   method = c("", "", "", "cart"), print.flag = FALSE) 
## 
## summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
##   data = s6, family = "binomial"), population.inference = TRUE)


###################################################
### code chunk number 9: inference.Rnw:243-249
###################################################
s7 <- syn(ods, seed = 910011, m = 1, 
  method = c("", "", "", "cart"), print.flag = FALSE) 

summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
  data = s7, family = "binomial"), population.inference = TRUE)



###################################################
### code chunk number 10: inference.Rnw:257-260 (eval = FALSE)
###################################################
## summary(glm.synds(smoke ~ sex + age + edu + sex * edu, 
##   data = s6, family = "binomial"), population.inference = TRUE, 
##   incomplete = TRUE, msel = 1:2)


###################################################
### code chunk number 11: inference.Rnw:340-347 (eval = FALSE)
###################################################
## ods <- ods[!is.na(ods$smoke), ] # remove 10 observations with missing "smoke"
## s8 <- syn.strata(ods, m = 5, method = "parametric", strata = "smoke", 
##   seed = 5678, visit.sequence = c("smoke", "sex", "edu", "age"), 
##   print.flag = FALSE, tab.strataobs = FALSE)
## f8 <- glm.synds(smoke ~ sex  + edu + sex * edu, data = s8,
##   family = "binomial")
## compare(f8, ods, plot.intercept = TRUE, plot = "coef") 


###################################################
### code chunk number 12: inference.Rnw:380-386 (eval = FALSE)
###################################################
## ods <- ods[!is.na(ods$smoke), ] # remove 10 observationswith missing "smoke"
## s9 <- syn(ods, m = 5, seed = 5678, 
##   visit.sequence = c("sex", "edu", "age", "smoke"), print.flag = FALSE)
## f9 <- glm.synds(smoke ~ sex  + age + edu + sex * age, data = s9, 
##   family = "binomial")
## compare(f9, ods)  


###################################################
### code chunk number 13: inference.Rnw:423-428 (eval = FALSE)
###################################################
## s10 <- syn(ods, m = 5, seed = 5678, method = "sample", 
##   visit.sequence = c("sex", "edu", "age", "smoke"), print.flag = FALSE)
## f10 <- glm.synds(smoke ~ sex  + age + edu + sex * edu, data = s10, 
##   family = "binomial")
## compare(f10, ods)


###################################################
### code chunk number 14: inference.Rnw:461-466 (eval = FALSE)
###################################################
## s11 <- syn(ods, seed = 910011, m = 20, method = c("", "", "cart", "cart"), 
##   print.flag = FALSE) 
## f11 <- glm.synds(smoke ~ sex + age + edu + sex * edu, data = s11, 
##   family = "binomial")
## compare(f11, ods)


###################################################
### code chunk number 15: inference.Rnw:513-519 (eval = FALSE)
###################################################
## ods <- ods[!is.na(ods$smoke), ]
## s12 <- syn.strata(ods, m = 5, visit.sequence = c(4, 1, 2), 
##   method = "parametric", strata = "smoke", seed = 5678, print.flag = FALSE)
## s13 <- syn.strata(ods, m = 5, visit.sequence = c(4, 1, 2), 
##   method = "parametric", strata = "smoke", seed = 1234, proper = TRUE, 
##   print.flag = FALSE, tab.strataobs = FALSE)


###################################################
### code chunk number 16: inference.Rnw:521-524 (eval = FALSE)
###################################################
## f12 <- glm.synds(smoke ~ sex  + edu + sex * edu, data = s12, 
##   family = "binomial")
## compare(f12, ods, plot.intercept = TRUE, plot = "coef") 


###################################################
### code chunk number 17: inference.Rnw:551-554 (eval = FALSE)
###################################################
## f13 <- glm.synds(smoke ~ sex  + edu + sex * edu, data = s13, 
##   family = "binomial")
## compare(f13, ods, plot.intercept = TRUE)


###################################################
### code chunk number 18: inference.Rnw:588-592 (eval = FALSE)
###################################################
## s14 <- syn(ods, m = 5, seed = 9101112, method = "parametric",
##   print.flag = FALSE)
## s15 <- syn(ods, m = 5, seed = 1415, method = "parametric", 
##   proper = TRUE, print.flag = FALSE)


###################################################
### code chunk number 19: inference.Rnw:594-597 (eval = FALSE)
###################################################
## f14 <- glm.synds(smoke ~ sex  + edu + age + sex * edu, data = s14, 
##   family = "binomial")
## compare(f14, ods, plot.intercept = TRUE, plot = "coef") 


###################################################
### code chunk number 20: inference.Rnw:625-628 (eval = FALSE)
###################################################
## f15 <- glm.synds(smoke ~ sex  + edu + age + sex * edu, data = s15, 
##   family = "binomial")
## compare(f15, ods, plot.intercept = TRUE) 


###################################################
### code chunk number 21: inference.Rnw:672-675 (eval = FALSE)
###################################################
## s16 <- syn(ods, proper = TRUE, print.flag = FALSE)
## f16 <- glm.synds(smoke ~ sex  + edu  + sex * edu, data = s16, 
##   family = "binomial")


###################################################
### code chunk number 22: inference.Rnw:677-678 (eval = FALSE)
###################################################
## compare(f16, ods, plot.intercept = TRUE, plot = "coef") 


###################################################
### code chunk number 23: inference.Rnw:705-707 (eval = FALSE)
###################################################
## compare(f16, ods, plot.intercept = TRUE, population.inference = TRUE, 
##   plot = "coef")


