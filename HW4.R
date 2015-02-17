library(MASS)

## PROBLEM 1a
# sample function below samples randomly from f(x) specified in part a
sample = function(n) tan(pi/2*runif(n))

## PROBLEM 1b
r_sample = function(n) {
  Vy=numeric(n); Vcpt=integer(n); # initialize vectors
  for(j in 1:n){
    y = rexp(1,rate=1); u=runif(1); cpt=1;
      while(u>exp((-.5)*(y-1)^2)) {
        u=runif(1); y = rexp(1,rate=1);
        cpt=cpt+1
      }
    Vy[j]=y; Vcpt[j]=cpt
  }
  return(list(Vy,Vcpt))
}
X = r_sample(50000)
# visualize 5000 draws from this distribution
hist(X[[1]], breaks=100)
# visually make sure 'g' function completely envelopes 'f'
fun1 = function(x) 0.5*exp(-(x^2)/2)
fun2 = function(x) exp(-x)
fun3 = function(x) fun1(x)/fun2(x)
curve(fun1, from=0, to=5, ylim=c(0,1))
curve(fun2, from=0, to=5, add=TRUE)
curve(fun3, from=0, to=5, add=TRUE)


## PROBLEM 2
m=c(2.16, 0.74, 1.87, 3.03, 3.11, 2.74, 1.23, 3.64, 1.57, 2.12)
#exponential function when x = m_bar
fun = function(x) exp(-.5*sum((rep(x,length(m))-m)^2))
c = fun(mean(m))

gamma_sample = function(n) {
  Vy=numeric(n); Vcpt=integer(n); # initialize vectors
  for(j in 1:n){
    y = rgamma(1,shape=2, rate=0.5); u=runif(1); cpt=1;
    while(u>((1/c)*y*exp(-0.5*y)*fun(y))) {
      u=runif(1); y = rgamma(1,shape=2, rate=0.5);
      cpt=cpt+1
    }
    Vy[j]=y; Vcpt[j]=cpt
  }
  return(list(Vy,Vcpt))
}

X = gamma_sample(1000)
hist(X[[1]], breaks=100)

# monte carlo integral and confidence interval
ihat = mean(X[[1]])
se = sqrt(var(X[[1]])/1000)
z= qnorm(0.975)
CI= c(ihat-z*se,ihat+z*se) # 95% conf int

## PROBLEM 3
mu=c(0,0)
Sigma = matrix(c(1,0.8,0.8,1.0),nrow=2)
d=length(mu)
bivn = matrix(0,nrow=100,ncol=2)
for (i in 1:100) {
  L=t(chol(Sigma)); Z=rnorm(d,0,1)
  X=mu+L%*%Z
  bivn[i,1]=X[1,1]
  bivn[i,2]=X[2,1]
}

# plot bivariate normal
bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 100)
contour(bivn.kde)

# importance sampling using bivariate normal
x = bivn[,1]
y = bivn[,2]
sxx = sum(x^2)
syy = sum(y^2)
sxy = sum(x*y)

n = 5000
# correlation density function assuming n = 100
corr_fun = function(p) { 
  ((1/(1-p^2))^(100/2))*exp(-1/(2*(1-p^2))*(sxx+syy-2*p*sxy))
}
sple = runif(n, min=-1, max=1)
om = corr_fun(sple)
est=sum(sple*om)/sum(om)
v1=var((om/mean(om))*sple)
se=sqrt(v1/n)
z=qnorm(0.975)
ci = c(est-z*se, est+z*se)