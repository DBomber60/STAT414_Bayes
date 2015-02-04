library(MASS)
setwd("~/Documents/STAT414_Bayes")

### PROBLEM ONE
student_data = read.table("studentdata.txt",header=T,sep='')
x_0 = student_data$ToSleep
z_0 = student_data$WakeUp
incompletes = unique(c((which(is.na(student_data$ToSleep))), (which(is.na(student_data$WakeUp)))))
x = x_0[-incompletes]
z = z_0[-incompletes]
# define a variable y that holds number of hours of sleep for each student
y = z-x

# part b (set up indices for group 1 and group 2)
tmp = 1:length(y)
group_a = sample(1:length(y), size=500)
group_b = tmp[-group_a]
group_1 = y[group_a] # total sleep
group_2 = y[group_b]
group_1_bedTime = x[group_a]
group_2_bedTime = x[group_b]

# part c/d: set up the prior
n = length(group_1)
mu_0 = 8
kappa_0 = 1
nu_0 = 100
sigmasq_0 = 8

# posterior parameters
kappa_n = n + kappa_0
nu_n = n + nu_0
mu_n = (kappa_0*mu_0+sum(group_1))/(kappa_n)
sigmasq_n = (sigmasq_0*nu_0+kappa_0*mu_0^2+sum(group_1^2)-mu_n^2)/nu_n

# part e: posterior inference
nMC = 153
post_sigsq_sample = 1/rgamma(nMC,0.5*nu_n, 0.5*nu_n*sigmasq_n)
post_mu_sample = rnorm(nMC, mu_n, sqrt(post_sigsq_sample/kappa_n))
predictions = rnorm(nMC, mean=post_mu_sample, sd= sqrt(post_sigsq_sample))
sum(abs(predictions-group_2))/length(group_2) # average prediction error: ~6 hours

### PROBLEM TWO
n = length(group_1)
m = length(group_2)
tau_0 = 1e-3
tau_1 = 1e-3
sigmasq = mean(post_sigsq_sample)
sigma_n_inv = matrix(c(tau_0+n/sigmasq, 1/sigmasq*sum(group_1_bedTime), 
                        1/sigmasq*sum(group_1_bedTime), 1/sigmasq*sum(group_1_bedTime^2)+tau_1), nrow=2)
sigma_n = solve(sigma_n_inv)
mu_n =(1/sigmasq)*sigma_n%*%c(sum(group_1),sum(group_1_bedTime*group_1))
mu_n = as.vector(mu_n)
mean_b0 = mu_n[1]; mean_b1 = mu_n[2]
# sample from the multivariate normal to get intercept and slope coefficient
coefficients = mvrnorm(n = m, mu_n, sigma_n)
predictions = coefficients[,1]+ coefficients[,2]*group_2_bedTime
sum(abs(predictions-group_2))/length(group_2) # average prediction error: ~1 hour