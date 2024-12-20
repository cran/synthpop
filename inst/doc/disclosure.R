### R code from vignette source 'disclosure.Rnw'

###################################################
### code chunk number 1: disclosure.Rnw:44-45
###################################################
options(prompt="R> ", width=77, digits=4, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: synth
###################################################
library("synthpop")
ods <- SD2011[, c("sex", "age", "region","placesize","depress",
   "income","ls","marital" , "workab")]
s1 <- syn(ods, seed = 8564, print.flag = FALSE, cont.na = list(income = -8))
t1 <- multi.disclosure(s1, ods, print.flag = FALSE, plot = FALSE, 
 keys = c("sex", "age", "region", "placesize"),   ngroups_targets = c(0,20,0,0,0))
print(t1,to.print = "ident")


###################################################
### code chunk number 3: mof5
###################################################
s5 <- syn(ods, seed = 8564, m = 5, print.flag = FALSE)
t5 <- disclosure( s5, ods, keys = c("sex", "age", "region",
   "placesize"), target = "depress", print.flag = FALSE)
print(t5, to.print = c("ident"))


###################################################
### code chunk number 4: mof5a
###################################################
print(t5, to.print = c("attrib"))


###################################################
### code chunk number 5: mof5
###################################################
multi.disclosure(s1, ods, print.flag = FALSE, plot = FALSE,
   keys = c("sex", "age", "region", "placesize"),
   denom_lim =1, exclude_ov_denom_lim = TRUE)


###################################################
### code chunk number 6: workab
###################################################
d1_workab <- disclosure(s1, ods, print.flag = FALSE, target = "workab",
    keys = c("sex", "age", "region", "placesize"),plot = FALSE)
print(d1_workab, to.print = c("check_1way"))


###################################################
### code chunk number 7: marital
###################################################
d1_marital <- disclosure(s1, ods, print.flag = FALSE, target = "marital",
    keys = c("sex", "age", "region", "placesize"),plot = FALSE)
print(d1_marital, to.print = c("check_2way"))


###################################################
### code chunk number 8: disclosure.Rnw:446-447
###################################################
print(t5, to.print = "allCAPs")


