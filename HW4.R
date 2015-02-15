## PROBLEM 1


sample = function(n) tan(pi/2*runif(n))
# visualize resulting density, excluding values over 10
hist(a[a<10], breaks=100)

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
mean(X[[2]])
hist(X[[1]], breaks=100)

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
mean(X[[1]])
mean(X[[2]])


x = seq(from=0, to=10, by=1/1000)
y = dgamma(x, shape=2, rate = 0.5)
plot(x,y)

