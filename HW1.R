list.of.packages = c("ggplot2")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)

### PROBLEM 2
prior = function(theta) {
  ifelse(theta>=0 & theta <=0.5, 4*theta, ifelse(theta > 0.5 & theta <= 1, 4-4*theta,0))
}

likelihood = function(theta, n, success){theta^success*(1-theta)^(n-success)}

marginal_y = function(theta) {
  ifelse(theta>=0 & theta <=0.5, choose(82,45)*4*theta*likelihood(theta,82,45), 
  ifelse(theta > 0.5 & theta <= 1, choose(82,45)*(4-4*theta)*likelihood(theta,82,45),0))
}

# evaluate marginal of y at y=45
integrate(marginal_y,0,1)

posterior_R1 = function(theta) {dbeta(theta,45+2,82-45+1)}
a = integrate(posterior_R1,0,.5)$value
posterior_R2 = function(theta) {dbeta(theta,45+1,82-45+2)}
b = integrate(posterior_R2,0.5,1)$value
c1 = a/(a+b)
c2 = 1-c1

posterior = function(theta) {
  ifelse(theta>=0 & theta <=0.5, (c1/a)*dbeta(theta,45+2,82-45+1),
  ifelse(theta > 0.5 & theta <= 1, (c2/b)*dbeta(theta,45+1,82-45+2),0))
}

posterior_mean = function(theta) {
  ifelse(theta>=0 & theta <=0.5, theta*(c1/a)*dbeta(theta,45+2,82-45+1),
  ifelse(theta > 0.5 & theta <= 1, theta*(c2/b)*dbeta(theta,45+1,82-45+2),0))
}

# posterior mean of theta
mean = integrate(posterior_mean, 0,1)$value

posterior_second_moment = function(theta) {
  ifelse(theta>=0 & theta <=0.5, theta^2*(c1/a)*dbeta(theta,45+2,82-45+1),
  ifelse(theta > 0.5 & theta <= 1, theta^2*(c2/b)*dbeta(theta,45+1,82-45+2),0))
}

variance = integrate(posterior_second_moment,0,1)$value - mean^2

ggplot(data.frame(x=c(0, 1)), aes(x)) + stat_function(fun=posterior, aes(colour="Posterior")) + stat_function(fun=prior, aes(colour="Prior"))+
  scale_colour_manual("Distributions", values = c("red", "blue"))

### PROBLEM 3
# prior distribution of the mean - exponential with rate parameter = 1/5
x = seq(0, 20, length=1000)
density = dexp(x, rate = 1/5)
m = ((120*14.8)/4 - .2)/(120/4)
v = 4/120
posterior = dnorm(x, mean=m, sd = sqrt(v))
df = data.frame(cbind(x,density, posterior))
ggplot(df, aes(x)) + geom_line(aes(y=density, colour="prior")) + geom_line(aes(y=posterior, color="posterior"))+
  scale_colour_manual("Distributions", values = c("red", "blue"))

