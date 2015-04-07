rm(list=ls())
setwd("~/Documents/STAT414_Bayes/HW9")
y = read.table('stop.dta',header=F)
# of stop and frisk
y = c(y[,1],y[,2],y[,3])
n=length(y)
mat = read.table('arrest.dta',header=F)
# log number of arrests
X = cbind(rep(1,n),c(mat[,1],mat[,2],mat[,3]))
mod=glm(y~X-1, family=poisson)
summary(mod)

logpi = function(beta,y,X) {
  c = 100
  a = X%*%beta
  val = sum(-exp(a)+y*a)-(0.5/(c^2))*crossprod(beta)
  return(val)
}

RWM_pois = function(Niter, y, X, lsig) {
  p = dim(X)[2]
  Res = matrix(NA, ncol=p, nrow=Niter)
  beta = double(p)
  lpi = logpi(beta, y, X)
  for (jj in 1:Niter) {
    beta_prop = beta + exp(lsig) * rnorm(p)
    lpi_prop = logpi(beta_prop, y, X)
    acc = min(1, exp(lpi_prop - lpi))
    u = runif(1)
    if (u <= acc) {
      beta = beta_prop
      lpi = lpi_prop
    }
    Res[jj,] = beta
  }
  return (Res)
}

#Call sampler
Niter = 100000
lsig = -2
Res = RWM_pois(Niter, y, X, lsig)

#Evaluate results
Output = Res[10001:Niter,1]
par(mfrow=c(1,3))
plot(Output, type='l')
acf(Output, lag.max=100)
hist(Output, nclass=50)
CI = apply(Res,2, function(x) {quantile(x, c(0.025, 0.975))})
CI

###### PROBLEM 2

# reformat data for clarity
y_race1 = y[1:(n/3)]
y_race2 = y[((n/3)+1):((n/3)*2)]
y_race3 = y[((n/3)*2+1):n]

x_race1 = X[1:(n/3),2]
x_race2 = X[((n/3)+1):((n/3)*2),2]
x_race3 = X[((n/3)*2+1):n,2]

# compute 'true' statistics based on data
T_1 = sum(y_race1*exp(-x_race1)-y_race3*exp(-x_race3)) # 884 
T_2 = sum(y_race2*exp(-x_race2)-y_race3*exp(-x_race3)) # 527

# simulate T1 and T2 based on log linear model
B_0 = mean(Res[10001:Niter,1])
B_1 = mean(Res[10001:Niter,2])

Iter=10000
# simulate y's
y_race1_sim = c()
y_race2_sim = c()
y_race3_sim = c()
for (i in 1:Iter) {
  x1 = sample(x_race1,1)
  y_race1_sim[i] = rpois(1, exp(B_0+B_1*x1))
  x2 = sample(x_race2,1)
  y_race2_sim[i] = rpois(1, exp(B_0+B_1*x2))
  x3 = sample(x_race3,1)
  y_race3_sim[i] = rpois(1, exp(B_0+B_1*x3))
}

# simulate T1, T2
T_1 = c()
T_2 = c()
for (i in 1:Iter) {
  y1 = sample(y_race1_sim,27)
  y3 = sample(y_race3_sim,27)
  x1 = sample(x_race1,27)
  x3 = sample(x_race3,27)
  T_1[i] = sum(y1*exp(-x1)-y3*exp(-x3))
  y2 = sample(y_race2_sim,27)
  y3 = sample(y_race3_sim,27)
  x2 = sample(x_race2,27)
  x3 = sample(x_race3,27)
  T_2[i] = sum(y2*exp(-x2)-y3*exp(-x3))
}

hist(T_1, prob=TRUE)
mean(T_1)
hist(T_2, prob=TRUE)
mean(T_2)
