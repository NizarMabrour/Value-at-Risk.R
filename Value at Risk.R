


#################################################################################
#########   ETUDE DE LA VALUE AT RISK SUR ACTION EUROFINS SCIENTIFICS    #######
#################################################################################


rm(list=ls(all=TRUE))

library(QuantTools)
library(xts)
library(forecast)
library(moments)
library(BatchGetSymbols)
library(rugarch)
library(ghyp)
first.date <- "2010-01-04" #date debut
last.date <- "2019-12-31"#date fin
freq.data <- 'daily'#frequence
type.return<-'log'#type de rendement
tickers <- 'ERF.PA' #symbole de l'action sur yahoo/finance
tab <- BatchGetSymbols(tickers = tickers,
                          first.date = first.date,
                         last.date = last.date,
                         freq.data = freq.data,
                         type.return=type.return,
                         cache.folder = file.path(tempdir(),
                                                    'BGS_Cache') ) # cache in tempdir()
pt<-tab$df.tickers$price.adjusted

dates<-tab$df.tickers$ref.date[-1]
rendement=tab$df.tickers$ret.adjusted.prices[-1]
N<-length(rendement)
rt<-xts(x=rendement,order.by=dates)
length(rt)
rte=rendement[1:1510]#rt sur ensemble estimation de 2010 a 2015 inclus
datesrte<-dates[1:1510]
datesrtt = dates[1511:N]
rtt=rendement[1511:N]#rt sur l'ensemble de test de 2016 a 2019 inclus
alpha=0.95
Ne=length(rte)
Nt=length(rtt)

#############################################
#######ESTIMATION DE LA DISTRIBUTION
############################################
qqnorm(rte)
qqline(rte, col = 2)

#estimation d'une distribution normale
fitnorm<-fit.gaussuv(data=rte)
summary(fitnorm)

#estimation d'une distribution student asymetrique.
fitstudentasy<-fit.tuv(data=rte, silent =T)
summary(fitstudentasy)

#estimation d'une distribution gaussienne inverse asymetrique.
nig<-fit.NIGuv(data=rte,silent=T)
summary(nig)

#estimation d'une distribution asymetrique hyperbolique.

fithyp<-fit.hypuv(rte,silent=T)
summary(fithyp)

#Hyperbolique generalise asymetrique 
fitghypuv<-fit.ghypuv(rte,silent=T)
summary(fitghypuv)

#Skewness Hyperbolique
library(SkewHyperbolic)
skewhyper<-skewhypFit(rte, print = FALSE, plot =FALSE, hessian = TRUE)
summary(skewhyper)

plot(density(rte),col = 1)
lines(fitstudentasy,col=2)#student asymetrique
lines(nig,col=3)# nig
lines(fithyp,col=4)#hyperbolique 
lines(fitghypuv,col=5)#hyperbolique generalisee
legend("topleft",legend =c("rte","student asy","nig","hyper","hyper gen"), col =1:5,lty=rep(1,5))

library(rugarch)


#############################################
####### BONNE DISTRIBUTION
############################################


# student asymetrique.
fitstudentasy<-fit.tuv(data=rte, silent =T)
summary(fitstudentasy)

plot(density(rte),col = 1)
lines(fitstudentasy,col=2)
legend("topleft",legend =c("rte","student asymetrique"), col =1:5,lty=rep(1,5))


#############################################
####### VAR normal, Historique, Cornish Fisher 
############################################


backTestVaR <- function(x, p = alpha) {
  normal.VaR = as.numeric(VaR(x, p=p, method="gaussian"))
  historical.VaR = as.numeric(VaR(x, p=p, method="historical"))
  modified.VaR = as.numeric(VaR(x, p=p, method="modified"))
  ans = c(normal.VaR, historical.VaR, modified.VaR)
  names(ans) = c("Normal", "HS", "Modified")
  return(ans)
}

# rolling 1-step ahead estimates of VaR
VaR.results = rollapply(as.zoo(rt), width=Ne, 
                        FUN = backTestVaR, p=alpha, by.column = FALSE,
                        align = "right")
#VaR.results = lag(VaR.results, k=-1)


violations.mat = matrix(0, 3, 5)
rownames(violations.mat) = c("Normal", "HS", "Modified")
colnames(violations.mat) = c("En1", "n1", "1-alpha", "Percent", "VR")
violations.mat[, "En1"] = (1-alpha)*Nt
violations.mat[, "1-alpha"] = 1 - alpha

# Show Normal VaR violations
normalVaR.violations = as.numeric(as.zoo(rt[index(VaR.results)])) < VaR.results[, "Normal"]
violation.dates = index(normalVaR.violations[which(normalVaR.violations)])

for(i in colnames(VaR.results)) {
  VaR.violations = as.numeric(as.zoo(rt[index(VaR.results)])) < VaR.results[, i]
  violations.mat[i, "n1"] = sum(VaR.violations)
  violations.mat[i, "Percent"] = sum(VaR.violations)/Nt
  violations.mat[i, "VR"] = violations.mat[i, "n1"]/violations.mat[i, "En1"]
}
violations.mat

resultats<-data.frame(matrix(NA,ncol=5,nrow=3))
colnames(resultats)<-c("expected.exceed","actual.exceed","Kupiecpv","Christoffersenpv","Violation_rate")
rownames(resultats)<-c("Normale","HS","CF")

# normale
VaR.test1 = VaRTest(1-alpha,actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"Normal"]))
resultats[1,1]=VaR.test1$expected.exceed
resultats[1,2]=VaR.test1$actual.exceed
resultats[1,3]=VaR.test1$uc.LRp
resultats[1,4]=VaR.test1$cc.LRp
resultats[1,5]=(VaR.test1$actual.exceed/length(VaR.results))

# historique
VaR.test2 = VaRTest(1-alpha,actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"HS"]))
resultats[2,1]=VaR.test2$expected.exceed
resultats[2,2]=VaR.test2$actual.exceed
resultats[2,3]=VaR.test2$uc.LRp
resultats[2,4]=VaR.test2$cc.LRp
resultats[2,5]=VaR.test2$actual.exceed/length(VaR.results)


# modifie
VaR.test3 = VaRTest(1-alpha, actual=coredata(rt[index(VaR.results)]), VaR=coredata(VaR.results[,"Modified"]))

resultats[3,1]=VaR.test3$expected.exceed
resultats[3,2]=VaR.test3$actual.exceed
resultats[3,3]=VaR.test3$uc.LRp
resultats[3,4]=VaR.test3$cc.LRp
resultats[3,5]=VaR.test3$actual.exceed/length(VaR.results)



###################################################
#### VAR PARAMETRIQUE ROLLING WINDOW APARCH STUDENT ASYMETRIQUE AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_ap = ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_ap, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)

roll_ap=ugarchroll(spec_ap, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                refit.window="moving",solver = "hybrid",solver.control=list(tol=1e-6, trace=1),
                cluster=cl,fit.control=list(scale=1),calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_ap <- zoo(roll_ap@forecast$VaR[, 1])
tail(roll_VaR_ap)

report(roll_ap,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_ap = ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(gamma1 = 0 ))
fit = ugarchfit(spec =spec_ap, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_ap
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("sstd",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_ap=fitted(filt)+sigma(filt)*qdist("sstd",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])
index(VaR_ap) = dates[1511:N]
# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_ap, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_ap, boot = TRUE))

stopCluster(cl)


###################################################
####### VAR PARAMETRIQUE ROLLING WINDOW  gjrGARCH SSTD AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_gjr = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(alpha1 = 0))
fit = ugarchfit(spec =spec_gjr, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)
roll_gjr=ugarchroll(spec_gjr, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                refit.window="moving",solver = "hybrid",solver.control=list(tol=1e-6, trace=1),
                cluster=cl,fit.control=list(scale=1),calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_gjr <- zoo(roll_gjr@forecast$VaR[, 1])

tail(roll_VaR_gjr)

report(roll_gjr,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_gjr = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(alpha1 = 0))
fit = ugarchfit(spec =spec_gjr, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_gjr
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("sstd",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_gjr=fitted(filt)+sigma(filt)*qdist("sstd",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])

# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_gjr, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_gjr, boot = TRUE))
stopCluster(cl)


###################################################
### VAR PARAMETRIQUE ROLLING WINDOW eGARCH SSTD AND BACKTESTING
###################################################


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
spec_e = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(gamma1 = 0))
fit = ugarchfit(spec =spec_e, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit)
roll_e=ugarchroll(spec_e, data=rt,n.ahead=1,forecast.length=length(rtt),refit.every=10,
                refit.window="moving",solver = "hybrid",solver.control=list(tol=1e-6, trace=1),
                cluster=cl,fit.control=list(scale=1),calculate.VaR=TRUE,VaR.alpha=0.05,keep.coef = TRUE)




roll_VaR_e <- zoo(roll_e@forecast$VaR[,1])
tail(roll_VaR_e)

report(roll_e,type="VaR",VaR.alpha=0.05,conf.level=0.95)

# calcul ES
spec_e = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(2,2)),distribution.model="sstd",fixed.pars = list(gamma1 = 0))
fit = ugarchfit(spec =spec_e, data = rt,out.sample=length(rtt),solver="hybrid")
spec_a = spec_e
setfixed(spec_a)<-as.list(coef(fit))
filt = ugarchfilter(spec_a, data=rtt)
f = function(x) qdist("sstd",p=x,mu=0,sigma=1,skew=coef(fit)["skew"],shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
VaR_e=fitted(filt)+sigma(filt)*qdist("sstd",p=0.05,mu=0,sigma=1,skew = coef(fit)["skew"],shape=coef(fit)["shape"])


# EXPECTED SHORTFALL FOR ROLLING WINDOW VAR METHOD
print(ESTest(0.05, rtt, ES, roll_VaR_e, boot = TRUE))
# EXPECTED SHORTFALL FOR ONE ESTIMATION VAR METHOD
print(ESTest(0.05, rtt, ES, VaR_e, boot = TRUE))
stopCluster(cl)


#############################
###########Ploting of All Var
############################

plot(dates[1511:N],rtt,type='c',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting model GJR GARCH sstd")
lines(dates[1511:N],roll_VaR_gjr,type='l',col="red")
legend("topright", inset=.05, c("ERF.PA return","VaR (gjrGARCH - sstd)"),col = c("black","red"), lty = c(1,1))

plot(dates[1511:N],rtt,type='b',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting Model APARCH sstd")
lines(dates[1511:N],roll_VaR_ap,type='l',col="brown")
legend("topright", inset=.05, c("ERF.PA return","VaR (APGARCH - sstd)"),col = c("black","red"), lty = c(1,1))

plot(dates[1511:N],rtt,type='b',xlab="Dates",ylab="Return/VaR",main = "95% 1 day VaR Backtesting Model EGARCH sstd")
lines(dates[1511:N],roll_VaR_e,type='l',col="blue")
legend("topright", inset=.05, c("ERF.PA return","VaR (EGARCH - sstd)"),col = c("black","red"), lty = c(1,1))

############# END
