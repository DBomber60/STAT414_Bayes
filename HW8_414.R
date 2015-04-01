##### PROBLEM 1

# generate data
n=200
p=5
X=cbind(rep(1,n), matrix(rnorm((p-1)*n),ncol=p-1))
beta_true=c(-2,0,1,2,1)
mp=X%*%beta_true
prob=exp(mp)/(1+exp(mp))
y=rbinom(n,1,prob)
c=100
sigma= 1/4

# compute log-posterior
logpi = function(beta,y,X){
  a=X%*%beta
  val=sum(y*a-log(1+exp(a)))-(0.5/(c^2))*crossprod(beta)
  return(val)
}

# gradient function
grad = function(beta,y,X){
  a=X%*%beta
  return(-(1/c^2)*beta+t(X)%*%(y-(exp(a)/(1+exp(a)))))
}

# propose a new beta using gradient
beta_prop = function(beta,y,X){
  Sigma= sigma^2*diag(p)
  L = t(chol(Sigma))
  Z = rnorm(p,0,1)
  mu = beta+0.5*(sigma^2)*grad(beta,y,X)
  beta_prop = mu+L%*%Z
  return(beta_prop)
}

# q function
log_q = function(rv, mean) {
  mu = mean+0.5*sigma^2*grad(mean,y,X)
  val = -0.5*sigma^2*crossprod(rv-mu)
  return(val)
}

MetroAdj_logRegr = function(Niter, y, X) {
  p=length(X[1,]); n = length(y)
  Res=matrix(NA, ncol=p, nrow=Niter)
  beta=double(p)
  lpi=logpi(beta,y,X)
  for(jj in 1:Niter){
    beta_prop=beta_prop(beta,y,X)
    lpi_prop=logpi(beta_prop,y,X)
    Acc=min(1,exp(lpi_prop+log_q(rv=beta_prop,beta)-lpi-log_q(rv=beta,beta_prop)))
    if(runif(1)<=Acc){
      beta=beta_prop
      lpi=lpi_prop
    }
    Res[jj,]=beta
  }
  return(Res)
}

Niter=20000
Res= MetroAdj_logRegr(Niter,y,X)
Output=Res[2001:Niter,1]
par(mfrow=c(1,3))
plot(Output,type="l",col='blue')
acf(Output,lag.max=100,col='blue')
hist(Output,nclass=50,prob=T,col='blue')
# generate confidence intervals
CI = apply(Res,2, function(x) {quantile(x, c(0.025, 0.975))})
colnames(CI)=c('beta0','beta1','beta2','beta3','beta4')
CI

##### PROBLEM 2

# Read data
setwd("~/Documents/STAT414_Bayes")
polls = read.table('polls.dt',header=T)
X = cbind(rep(1,nrow(polls)),polls[,2:5])
names(X) = c('Intercept',names(polls)[2:5])
X = as.matrix(X)
y = polls[,1]
p=ncol(X)
c=100


# compute log-posterior
logpi = function(beta,y,X){
  a=X%*%beta
  val=sum(log((pnorm(a)^y)*((1-pnorm(a))^(1-y))))-(0.5/(c^2))*crossprod(beta)
  return(val)
}

RWM_probRegr = function(Niter, y, X) {
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
  return(list(lsig,Res))
}

Niter=10000
Res= as.matrix(RWM_probRegr(Niter,y,X)[[2]])
Output=Res[1001:Niter,1]
par(mfrow=c(1,3))
plot(Output,type="l",col='blue')
acf(Output,lag.max=100,col='blue')
hist(Output,nclass=50,prob=T,col='blue')
# generate confidence intervals
CI = apply(Res,2, function(x) {quantile(x, c(0.025, 0.975))})
colnames(CI)=c('Intercept',names(polls)[2:5])
CI

# check using R's probit model
probit_test <- glm(y~X-1, family=binomial(link="probit"))
summary(probit_test)