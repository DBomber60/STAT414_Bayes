rm(list=ls())
setwd("~/Documents/STAT414_Bayes/HW 10")
data = read.table('mathstandard.dta',header=T)
y = data[,2]
n = length(y)
X = cbind(rep(1,n),data[,3])

# compute log posterior distribution
logpi = function(beta,y,X) {
  c = 100
  a = X%*%beta
  val=sum(y*a-log(1+exp(a)))-(0.5/(c^2))*crossprod(beta)
  return(val)
}

RWM_logRegr = function(Niter, y, X) {
  p=length(X[1,]); n = length(y)
  Res=matrix(NA, ncol=p, nrow=Niter)
  lsig=0 #initialize lsig
  beta=double(p)
  lpi=logpi(beta,y,X)
  for(jj in 1:Niter){
    beta_prop=beta+exp(lsig)*rnorm(p)
    lpi_prop=logpi(beta_prop,y,X)
    Acc=min(1,exp(lpi_prop-lpi))
    if(runif(1)<=Acc){
      beta=beta_prop
      lpi=lpi_prop
    }
    lsig = lsig + (1/jj^0.7)*(Acc-0.4)
    Res[jj,]=beta
  }
  print(lsig)
  return(list(lsig,Res))
}

Niter=100000
Res= as.matrix(RWM_logRegr(Niter,y,X)[[2]])
lsig = 
Output=Res[10001:Niter,2]
par(mfrow=c(1,3))
plot(Output,type="l",col='blue')
acf(Output,lag.max=100,col='blue')
hist(Output,nclass=50,prob=T,col='blue')
# generate confidence intervals
CI = apply(Res,2, function(x) {quantile(x, c(0.025, 0.975))})
colnames(CI)=c('Intercept','y')
CI

a = glm(y~X-1, family="binomial")
summary(a)

##### PROBLEM 2

# create function to calculate T given a vector of y values
model_mat=model.matrix(~data$county-1)
sums = apply(model_mat,2,sum)

T_fun = function(y_vec) {
  ones = t(y_vec)%*%model_mat
  prop = t(ones/sums)
  return(1/ncol(model_mat)*sum((prop-mean(prop))^2))
}

# T statistic for the data set = 0.08317
t_data = T_fun(y)

post_beta_m0 = Res[10001:Niter,]
n_sample = nrow(post_beta_m0)
T_stat = c()
for(i in 1:n_sample){
  a = X%*%post_beta_m0[i,]
  y_sample = rbinom(n,1,exp(a)/(1+exp(a)))
  T_stat[i] = T_fun(y_sample)
}

hist(T_stat)
abline(v=t_data, col='blue')
mean(T_stat) # mean of T for first model = 0.0625

##### PROBLEM 3 

# Create new model with a different intercept for each county

X = cbind(model_mat, data[,3])

ARWM_logRegr = function(Niter, y, X) {
  p=length(X[1,]); n = length(y)
  Res=matrix(NA, ncol=p, nrow=Niter)
  C = diag(p); K=t(chol(C))
  mu=numeric(p)
  cptUpdate = 0; LUpdate = 1000
  lsig=-1 #initialize lsig
  beta=numeric(p)
  lpi=logpi(beta,y,X)
  for(jj in 1:Niter){
    beta_prop=beta+exp(lsig)*(rnorm(p))*(K%*%rnorm(p))
    lpi_prop=logpi(beta_prop,y,X)
    Acc=min(1,exp(lpi_prop-lpi))
    if(runif(1)<=Acc){
      beta=beta_prop
      lpi=lpi_prop
    }
    lsig = lsig + (1/jj^0.7)*(Acc-0.4)
    mu=mu+(1/jj)*(beta-mu)
    bmu=as.vector(beta-mu)
    C=C+(1/jj)*(outer(bmu,bmu)-C)
    if(cptUpdate==LUpdate){
      K=t(chol(C))
      cptUpdate=0
    }else{
      cptUpdate=cptUpdate+1
    }
    Res[jj,]=beta
  }
  output=list(Res,lsig,C)
}

Niter=100000
Res= ARWM_logRegr(Niter,y,X)
Output=Res[[1]][10001:Niter,3]
par(mfrow=c(1,3))
plot(Output,type="l",col='blue')
acf(Output,lag.max=100,col='blue')
hist(Output,nclass=50,prob=T,col='blue')
# generate confidence intervals
CI = apply(Res[[1]],2, function(x) {quantile(x, c(0.025, 0.975))})
CI

post_beta_m1 = Res[[1]][10001:Niter,]
n_sample = nrow(post_beta_m1)
T_stat = c()
for(i in 1:n_sample){
  a = X%*%post_beta_m1[i,]
  y_sample = rbinom(n,1,exp(a)/(1+exp(a)))
  T_stat[i] = T_fun(y_sample)
}

hist(T_stat)
abline(v=t_data, col='blue')
mean(T_stat) # mean of T for second model = 0.122



##### PROBLEM 4
ma_pct = data[,3]
X = cbind(rep(1,n),ma_pct*model_mat)

ARWM_logRegr = function(Niter, y, X) {
  p=length(X[1,]); n = length(y)
  Res=matrix(NA, ncol=p, nrow=Niter)
  C = diag(p); K=t(chol(C))
  mu=numeric(p)
  cptUpdate = 0; LUpdate = 1000
  lsig=-1 #initialize lsig
  beta=numeric(p)
  lpi=logpi(beta,y,X)
  for(jj in 1:Niter){
    beta_prop=beta+exp(lsig)*(rnorm(p))*(K%*%rnorm(p))
    lpi_prop=logpi(beta_prop,y,X)
    Acc=min(1,exp(lpi_prop-lpi))
    if(runif(1)<=Acc){
      beta=beta_prop
      lpi=lpi_prop
    }
    lsig = lsig + (1/jj^0.7)*(Acc-0.4)
    mu=mu+(1/jj)*(beta-mu)
    bmu=as.vector(beta-mu)
    C=C+(1/jj)*(outer(bmu,bmu)-C)
    if(cptUpdate==LUpdate){
      K=t(chol(C))
      cptUpdate=0
    }else{
      cptUpdate=cptUpdate+1
    }
    Res[jj,]=beta
  }
  output=list(Res,lsig,C)
}

Niter=100000
Res= ARWM_logRegr(Niter,y,X)
Output=Res[[1]][10001:Niter,3]
par(mfrow=c(1,3))
plot(Output,type="l",col='blue')
acf(Output,lag.max=100,col='blue')
hist(Output,nclass=50,prob=T,col='blue')
# generate confidence intervals
CI = apply(Res[[1]],2, function(x) {quantile(x, c(0.025, 0.975))})
CI

post_beta_m1 = Res[[1]][10001:Niter,]
n_sample = nrow(post_beta_m1)
T_stat = c()
for(i in 1:n_sample){
  a = X%*%post_beta_m1[i,]
  y_sample = rbinom(n,1,exp(a)/(1+exp(a)))
  T_stat[i] = T_fun(y_sample)
}

hist(T_stat)
abline(v=t_data, col='blue')
mean(T_stat) # mean of T for second model = 0.120






