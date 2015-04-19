###SET-UP###
rm(list=ls())
setwd("~/Documents/STAT503_MvAnalysis")

#Load necessary packages
library(lubridate)

#Load data
data = read.csv('train-1.csv', header=T)

#Extract hour, day, and year from datetime to create separate variables for each
data$datetime = format(as.POSIXct((data$datetime), format = "%m/%d/%y %H:%M"))
data$hour = hour(data$datetime)             #hour of the day
data$day = wday(data$datetime)              #day of the week
data$month = month(data$datetime)           #month of the year
data$year = year(data$datetime)             #year

#Change the only weather=4 observation to weather=3
data$weather[data$weather==4] = 3

#Redefine categorical predictors as factors
data$hour = as.factor(data$hour)             #(24-hour clock, 0=12:00AM,..., 23=11:00PM)
data$day = as.factor(data$day)               #(1=Sunday, 2=Monday,...,7=Saturday)
data$month = as.factor(data$month)           #(1=January,...,12=December)
data$year = as.factor(data$year)           
data$season = as.factor(data$season)         #(1=spring, 2=summer, 3=fall, 4=winter)
data$holiday = as.factor(data$holiday)       #(1=holiday, 0=non-holiday)
data$workingday = as.factor(data$workingday) #(1=neither weekend nor holiday, 0=weekend or holiday)
data$weather = as.factor(data$weather)       #(1=clear, 2=mist/cloudy, 3=light precip.)

# processed dataset
X = model.matrix(~., data[,c(2:9, 13:16)])
y = data[,12]
n=length(y)
p=dim(X)[2]

# testing with reduced model
# X=X0[,c(1,9)]
# p=dim(X)[2]

###Implement Gibbs Sampler

# prior parameters
kappa=0.01
a=b=1.1
Niter=100
Res=matrix(NA, ncol=p+1, nrow=Niter)


# initialize Gibbs sampler
beta = rep(0,p)
M = kappa*diag(p)+t(X)%*%X
invM=solve(M)
R=chol(M)
for(jj in 1:Niter){
  # update sigmasq
  A=a+0.5*(n+p)
  B=b+0.5*kappa*crossprod(beta)+0.5*crossprod(y-X%*%beta)
  sigmasq=1/rgamma(1,shape=A,rate=B)
  # update beta
  beta=invM%*%(t(X)%*%y)+sqrt(sigmasq)*solve(R,rnorm(p))
  # record result
  Res[jj,]=c(sigmasq, as.vector(beta))
}

Output=Res[,3]
par(mfrow=c(1,3))
plot(Output, type='l',col='blue')
hist(Output, breaks=50, prob=T, col='blue')
acf(Output, lag=50, col='red')
CI = apply(Res,2, function(x) {quantile(x, c(0.025, 0.975))})
CI

## tester ##
m0 = lm(y~X-1)
summary(m0)

## test statistic ##
T_fun = function(y_vec, X_mat) {
  hour_average = c()
  for(i in 13:35) {
    hour_average = c(hour_average, mean(y_vec[which(X_mat[,i]==1)]))
  }
  return (var(hour_average))
}

t_data = T_fun(y,X) # 17689.3

# sample from the posterior
post_beta = Res[,-1]
post_sigsq = Res[,1]
n_sample = nrow(post_beta)
T_stat = c()
for(i in 1:n_sample){
  mu = X%*%post_beta[i,]
  sig = sqrt(post_sigsq[i])
  y_sample = rnorm(n,mean=mu,sd=sig)
  T_stat[i] = T_fun(y_sample, X)
}

qplot(T_stat, geom="histogram", binwidth=10)
abline(v=t_data, col='blue')
mean(T_stat)

Rsq=numeric(Niter)
ss_tot=sum((y-mean(y))^2)
for(i in 1:Niter){
  Rsq[i]=1-crossprod(y-X%*%Res[i,-1],y-X%*%Res[i,-1])/ss_tot
}
hist(Rsq,nclass=30,prob=T)
CI_Rsq=quantile(Rsq[1:Niter], c(0.025,0.975))
CI_Rsq