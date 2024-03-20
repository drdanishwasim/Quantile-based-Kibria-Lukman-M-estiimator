#R-Code for Robust Ridge Regression in presence of outliers (For real data)
#################Muhammad Suhail########################
#############Roll no P2-14########################
rm(list=ls())
set.seed(1990)
library(MASS)
p=5                                      #No of explanatory variables
I=diag(p)
setwd("F:/M.Phil & PhD/Ph.D Course work/PhD Research Work/PhD work/R-Codes practice/Real data/Sports Data")      #To change working directory
data=read.csv("sports data.csv",header=TRUE)
data=data[1:30,]
x1=data$x1                # use $ to add variables
x2=data$x2
x3=data$x3
x4=data$x4
x5=data$x5
n=length(x1)
x=cbind(x1,x2,x3,x4,x5)
as.matrix(x)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y

y=data$y
y=(y-mean(y))/sd(y)

#y=scale(x,center = TRUE,scale = TRUE)
eigen(cor(x))
ev=eigen(cor(x))$values
vec=eigen(cor(x))$vectors
which.max(ev)
CN=max(ev)/min(ev)        #Computing condition number
CI=sqrt(max(ev)/min(ev))
beta=vec[,which.max(ev)]
C=t(x)%*%(x)
D=eigen(C)$vectors
Z=x%*%D
lam=t(Z)%*%Z
lam=round(lam,4)
lamda=diag(lam)
alpha=t(D)%*%beta

betahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=x%*%betahat
sigmahat=(sum((y-yhat)^2))/(n-p)
alphahat=solve(lam)%*%t(Z)%*%y
y[5]=y[5]+200*sigmahat
y[10]=y[10]+300*sigmahat
y[15]=y[15]+400*sigmahat
#............................................
#y[35]=y[35]-700*sigmahat
#y[40]=y[40]+800*sigmahat

#Estimation of ridge parameter
#OLS estimator
alphahat.ols=alphahat
alphahat.ols=c(alphahat.ols)


#Ridge Regression estimator of HK method (1970)

HK=sigmahat/(max(alphahat^2))
K1=HK
T1=solve(lam+I*K1)%*%lam
alphahat.RR=T1%*%alphahat.ols
alphahat.RR = c(alphahat.RR)

# M-estimator denoted by alphahat.m of Huber 1981
r.fit=rlm(y~x-1,psi=psi.huber,k=1.345,maxit = 1000)
alphahat.m=as.vector(r.fit$coefficients)
alphahat.m=c(alphahat.m)
a.huber=r.fit$k

# To obtain A2
e=r.fit$residuals
s=r.fit$s
u=e/s
w.psi=psi.huber(u)
u.psi=u*w.psi
deriv.psi=psi.huber(u,deriv = 1)
num.a2=((s^2)*sum(u.psi^2))/(n-p)
denum.a2=(sum(deriv.psi)/n)^2
A2=num.a2/denum.a2

# Silvapulle (1991) Proposed M-estimator of HKB(1975) & LW(1976)

MHKB=(p*A2)/(sum(alphahat.m^2))
K2=MHKB
T2=solve(lam+I*K2)%*%lam
alphahat.MHKB=T2%*%alphahat.m
alphahat.MHKB=c(alphahat.MHKB)

MLW=(p*A2)/(sum(ev*(alphahat.m^2)))
K3=MLW
T3=solve(lam+I*K3)%*%lam
alphahat.MLW=T3%*%alphahat.m
alphahat.MLW=c(alphahat.MLW)

# Majid Kibria-Lukman M-estimator (2022) M-estimator  

KLME=A2/((2*alphahat.m^2)+(A2/ev))
KLME1=sum(KLME)/p
K4=KLME1
T4=solve(lam+I*K4)%*%lam
alphahat.KLME1=T4%*%alphahat.m
alphahat.KLME1=c(alphahat.KLME1)

KLME2=(prod(KLME))^(1/p)
K5=KLME2
T5=solve(lam+I*K5)%*%lam
alphahat.KLME2=T5%*%alphahat.m
alphahat.KLME2=c(alphahat.KLME2)

KLME3=p/(sum(1/KLME))
K6=KLME3
T6=solve(lam+I*K6)%*%lam
alphahat.KLME3=T6%*%alphahat.m
alphahat.KLME3=c(alphahat.KLME3)

# Danish Modified M-estimator of Majid Kibria Lukman M-estimators(2022)

DSM.01=quantile(KLME,probs = 0.01)
DSM.05=quantile(KLME,probs = 0.05)
DSM.25=quantile(KLME,probs = 0.25)
DSM.50=quantile(KLME,probs = 0.50)
DSM.75=quantile(KLME,probs = 0.75)
DSM.95=quantile(KLME,probs = 0.95)
DSM.99=quantile(KLME,probs = 0.99)

K7=DSM.01
T7=solve(lam+I*K7)%*%lam
alphahat.DSM.01=T7%*%alphahat.m
alphahat.DSM.01=c(alphahat.DSM.01)

K8=DSM.05
T8=solve(lam+I*K8)%*%lam
alphahat.DSM.05=T8%*%alphahat.m
alphahat.DSM.05=c(alphahat.DSM.05)

K9=DSM.25
T9=solve(lam+I*K9)%*%lam
alphahat.DSM.25=T9%*%alphahat.m
alphahat.DSM.25=c(alphahat.DSM.25)

K10=DSM.50
T10=solve(lam+I*K10)%*%lam
alphahat.DSM.50=T10%*%alphahat.m
alphahat.DSM.50=c(alphahat.DSM.50)

K11=DSM.75
T11=solve(lam+I*K11)%*%lam
alphahat.DSM.75=T11%*%alphahat.m
alphahat.DSM.75=c(alphahat.DSM.75)

K12=DSM.95
T12=solve(lam+I*K12)%*%lam
alphahat.DSM.95=T12%*%alphahat.m
alphahat.DSM.95=c(alphahat.DSM.95)

K13=DSM.99
T13=solve(lam+I*K13)%*%lam
alphahat.DSM.99=T13%*%alphahat.m
alphahat.DSM.99=c(alphahat.DSM.99)



comb.coeff=rbind(alphahat.ols,alphahat.RR,alphahat.m,alphahat.MHKB,alphahat.MLW,alphahat.KLME1,
                 alphahat.KLME2,alphahat.KLME3,alphahat.DSM.01,alphahat.DSM.05,alphahat.DSM.25,
                 alphahat.DSM.50,alphahat.DSM.75,alphahat.DSM.95,alphahat.DSM.99)
comb.coeff=round(comb.coeff,6)

###########################################################################################################
#MSE of all estimators
#MSE of OLS and RR

K<-c(0,K1) 
a=rep(0,length(K)); b=rep(0,length(K))                            #to obtain MSE of ridge estimators
for(i in 1:length(K)){
  
  a[i]=sum(ev/(ev+K[i])^2)
  b[i]=sum(((K[i]^2)*(alphahat.ols^2))/(ev+K[i])^2)
  
}
MSE1=sigmahat*a+b
#MSE of M and Ridge M-estimators
Km=c(0,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11)
a.m=rep(0,length(Km)); b.m=rep(0,length(Km))
for (i in 1:length(Km)) {
  a.m[i]=sum(((ev^2)*A2)/(ev*(ev+Km[i])^2))
  b.m[i]=sum(((Km[i]^2)*(alphahat.m^2))/(ev+Km[i])^2)
}
MSE2=a.m+b.m

#a.m1=sum(((ev^2)*A2)/(ev*(ev+KLME)^2))
#b.m1=sum(((KLME^2)*(alphahat.m^2))/(ev+KLME)^2)
#MSE3=a.m1+b.m1

#a.m2=sum(((ev^2)*A2)/(ev*(ev+DSS11)^2))
#b.m2=sum(((DSS11^2)*(alphahat.m^2))/(ev+DSS11)^2)

# MSE4=a.m2+b.m2

MSE=c(MSE1,MSE2)
MSE=as.matrix(MSE)
MSE=round(MSE,4)
K.hat=c(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11)
write.csv(K.hat,file = "F:/M.Phil & PhD/Ph.D Course work/PhD Research Work/PhD work/R-Codes practice/M-estimators/K.hat.PCRWR Report 2015-16.csv")
row.names(MSE)=c("OLS","RR","HUB.M","MHKB","MLW","KLME1","KLME2","KLME3","DSM.25","DSM.50","DSM.75","DSM.95","DSM.99")
colnames(MSE)=c("MSE")
write.csv(MSE,file = "F:/M.Phil & PhD/Ph.D Course work/PhD Research Work/PhD work/R-Codes practice/M-estimators/MSE.PCRWR Report 2015-16.csv")
write.csv(comb.coeff,file = "F:/M.Phil & PhD/Ph.D Course work/PhD Research Work/PhD work/R-Codes practice/M-estimators/coefficents.PCRWR Report 2015-16.csv")

cor.x=round(cor(x),4)
cor.y=round(cor(x,y),4)
cor.xy=cbind(cor.x,cor.y)
write.csv(cor.xy,file = "F:/M.Phil & PhD/Ph.D Course work/PhD Research Work/PhD work/R-Codes practice/M-estimators/correlation.PCRWR Report 2015-16.csv")
CN
CI
ev
xnam=paste0("x", 1:p) 
fmla=as.formula(paste("y ~ ", paste(xnam, collapse= "+"))) 
model=lm(fmla)

#Checking outliers
#Influence measures
inflm.model= influence.measures(model)
which(apply(inflm.model$is.inf, 1, any))      #outliers
summary(inflm.model)                          #Summary of outlying observations
#--------------------------------------------------------------------------------
#Prediction Summ of Square 
yhat.ols.cv=rep(0,n)
yhat.ols.cv1=rep(0,n)
yhat.RR.cv=rep(0,n)
yhat.m.cv=rep(0,n)
yhat.MHKB.cv=rep(0,n)
yhat.MLW.cv=rep(0,n)
yhat.KLME1.cv=rep(0,n)
yhat.KLME2.cv=rep(0,n)
yhat.KLME3.cv=rep(0,n)
yhat.DSM.01.cv=rep(0,n)
yhat.DSM.05.cv=rep(0,n)
yhat.DSM.25.cv=rep(0,n)
yhat.DSM.50.cv=rep(0,n)
yhat.DSM.75.cv=rep(0,n)
yhat.DSM.95.cv=rep(0,n)
yhat.DSM.99.cv=rep(0,n)

for (i in 1:n) {
  betahat.cv=solve(t(x[-i,])%*%x[-i,])%*%t(x[-i,])%*%y[-i]
  yhat.ols.cv[i]=x[i,]%*%betahat.cv
  
  alphahat.cv=solve(t(Z[-i,])%*%Z[-i,])%*%t(Z[-i,])%*%y[-i]
  alphahat.ols.cv=c(alphahat.cv)
  yhat.ols.cv1[i]=Z[i,]%*%alphahat.cv
  
  #Ridge Regression estimator of HK method (1970)
  K1.cv=sigmahat/(max(alphahat.cv^2))
  alphahat.RR.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K1.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.ols.cv)
  yhat.RR.cv[i]=Z[i,]%*%alphahat.RR.cv
  
  # M-estimator denoted by alphahat.m of Huber 1981
  r.fit.cv=rlm(y[-i]~x[-i,]-1,psi=psi.huber,k=1.345,maxit = 1000)
  alphahat.m.cv=as.vector(r.fit.cv$coefficients)
  alphahat.m.cv=c(alphahat.m.cv)
  a.huber=r.fit.cv$k
  yhat.m.cv[i]=Z[i,]%*%alphahat.m.cv
  
  # To obtain A2
  e=r.fit.cv$residuals
  s=r.fit.cv$s
  u=e/s
  w.psi=psi.huber(u)
  u.psi=u*w.psi
  deriv.psi=psi.huber(u,deriv = 1)
  num.a2=((s^2)*sum(u.psi^2))/(n-p)
  denum.a2=(sum(deriv.psi)/n)^2
  A2.cv=num.a2/denum.a2
  # Silvapulle (1991) Proposed M-estimator of HKB(1975) & LW(1976)
  K2.cv=(p*A2.cv)/(sum(alphahat.m.cv^2))
  alphahat.MHKB.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K2.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.MHKB.cv[i]=Z[i,]%*%alphahat.MHKB.cv
  
  ev=eigen(cor(x[-i,]))$values
  
  K3.cv=(p*A2.cv)/(sum(ev*(alphahat.m.cv^2)))
  alphahat.MLW.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K3.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.MLW.cv[i]=Z[i,]%*%alphahat.MLW.cv
  
  # Majid Kibria-Lukman M-estimator (2022) M-estimator  
  
  KLME.cv=A2.cv/((2*alphahat.m.cv^2)+(A2.cv/ev))
  K4.cv=sum(KLME.cv)/p
  K5.cv=(prod(KLME.cv))^(1/p)
  K6.cv=p/(sum(1/KLME.cv))
  
  alphahat.KLME1.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K4.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.KLME1.cv[i]=Z[i,]%*%alphahat.KLME1.cv
  
  alphahat.KLME2.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K5.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.KLME2.cv[i]=Z[i,]%*%alphahat.KLME2.cv
  
  alphahat.KLME3.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K6.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.KLME3.cv[i]=Z[i,]%*%alphahat.KLME3.cv
 
  # Danish Modified M-estimator of Majid Kibria Lukman M-estimators(2022)
  
  K7.cv=quantile(KLME.cv,probs = 0.01)
  K8.cv=quantile(KLME.cv,probs = 0.05)
  K9.cv=quantile(KLME.cv,probs = 0.25)
  K10.cv=quantile(KLME.cv,probs = 0.50)
  K11.cv=quantile(KLME.cv,probs = 0.75)
  K12.cv=quantile(KLME.cv,probs = 0.95)
  K13.cv=quantile(KLME.cv,probs = 0.99)
  
  alphahat.DSM.01.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K7.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.01.cv[i]=Z[i,]%*%alphahat.DSM.01.cv
  
  alphahat.DSM.05.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K8.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.05.cv[i]=Z[i,]%*%alphahat.DSM.05.cv
  
  alphahat.DSM.25.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K9.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.25.cv[i]=Z[i,]%*%alphahat.DSM.25.cv
  
  alphahat.DSM.50.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K10.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.50.cv[i]=Z[i,]%*%alphahat.DSM.50.cv
  
  alphahat.DSM.75.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K11.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.75.cv[i]=Z[i,]%*%alphahat.DSM.75.cv
  
  alphahat.DSM.95.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K12.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.95.cv[i]=Z[i,]%*%alphahat.DSM.95.cv
  
  alphahat.DSM.99.cv=c(solve(t(Z[-i,])%*%Z[-i,]+I*K13.cv)%*%t(Z[-i,])%*%Z[-i,]%*%alphahat.m.cv)
  yhat.DSM.99.cv[i]=Z[i,]%*%alphahat.DSM.99.cv
  
}
psse.ols=sum((y-yhat.ols.cv)^2)
psse.ols1=sum((y-yhat.ols.cv1)^2)
psse.RR=sum((y-yhat.RR.cv)^2)
psse.m=sum((y-yhat.m.cv)^2)
psse.MHKB=sum((y-yhat.MHKB.cv)^2)
psse.MLW=sum((y-yhat.MLW.cv)^2)
psse.KLME1=sum((y-yhat.KLME1.cv)^2)
psse.KLME2=sum((y-yhat.KLME2.cv)^2)
psse.KLME3=sum((y-yhat.KLME3.cv)^2)
psse.DSM.01=sum((y-yhat.DSM.01.cv)^2)
psse.DSM.05=sum((y-yhat.DSM.05.cv)^2)
psse.DSM.25=sum((y-yhat.DSM.25.cv)^2)
psse.DSM.50=sum((y-yhat.DSM.50.cv)^2)
psse.DSM.75=sum((y-yhat.DSM.75.cv)^2)
psse.DSM.95=sum((y-yhat.DSM.95.cv)^2)
psse.DSM.99=sum((y-yhat.DSM.99.cv)^2)

psse=c(psse.ols,psse.RR,psse.m,psse.MHKB,psse.MLW,psse.KLME1,psse.KLME2,psse.KLME3,
       psse.DSM.01,psse.DSM.05,psse.DSM.25,psse.DSM.50,psse.DSM.75,psse.DSM.95,
       psse.DSM.99)

psse
############################VIF##############################################

