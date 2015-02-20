# Read Data
setwd("~/Documents/STAT414_Bayes")
data=read.table('studentdata.txt',header=T)
dt=cbind(data$ToSleep, data$WakeUp-data$ToSleep)
dt=dt[complete.cases(dt),]
x = dt[,1]
y = dt[,2]
n=length(y)

# Set up prior
sigmasq_0=5
tau_0 = tau_1 = 1e-3

# log posterior B0, B1
mod = lm(y~x)
sigmasq = sum(mod$residuals^2)/(n-2)
sigma_n_inv = matrix(c(tau_0+n/sigmasq, sum(x)/sigmasq, sum(x)/sigmasq, tau_1+sum(x^2)/sigmasq), ncol=2)
sigma_n = solve(sigma_n_inv)
mu_n=(1/sigmasq)*sigma_n%*%c(sum(y), sum(y*x))
mu_n = as.vector(mu_n)

lpi = function(y) -t(y-mu_n)%*%sigma_n_inv%*%(y-mu_n)

# Metropolis algorithm for B0, B1
Niter = 15000
k = diag(1,2); C=chol(k)
sig = 0.0001
Res = matrix(NA, nrow=Niter, ncol=2)
beta = numeric(2)
log_pi=lpi(beta)
for (jj in 1:Niter){
  beta_prop = beta+sqrt(sig)*t(C)%*%rnorm(2,0,1)
  log_pi_prop = lpi(beta_prop)
  alpha=min(1, exp(log_pi_prop-log_pi))
  if (runif(1)<=alpha){
    beta=beta_prop
    log_pi=log_pi_prop
  }
  Res[jj,]=beta
}

# Intercept
output_b0 = Res[5000:Niter,1]
par(mfrow=c(1,3))
plot(output_b0, type="l",col="blue")
hist(output_b0, breaks=50, prob=T, col="blue")
acf(output_b0, lag=50, col="red")
mean(output_b0)
CI = quantile(output_b0, c(0.025, 0.975))

# Slope
output_b1 = Res[5000:Niter,2]
par(mfrow=c(1,3))
plot(output_b1, type="l",col="blue")
hist(output_b1, breaks=50, prob=T, col="blue")
acf(output_b1, lag=50, col="red")
mean(output_b1)
CI = quantile(output_b1, c(0.025, 0.975))

