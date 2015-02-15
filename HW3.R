### PROBLEM 1
log_series = function(p, n) {
  X = vector('numeric', n)
  for (i in 1:n){
    U = runif(1)
    k = 1; S = -(1-p)/log(p)
    while(S<U){
      k = k+1
      S = S+((-(1-p)^k)/(k*log(p)))
    }
    X[i] = k
  }
  return(X)
}

inversion_sample = log_series(0.3, 5000)
y_inversion = as.data.frame(table(inversion_sample)/5000)

#  sample directly from the probability mass function to compare distribution
test_fun = function(x) {(-.7^x)/(x*log(0.3))}
x = as.vector(as.numeric(y_inversion[,1])) # use x values from our inversion sample
y_direct = test_fun(x)

# compare graphically
barplot(rbind(y_direct, as.vector(as.numeric(y_inversion[,2]))), beside=T, col=c('red','blue'),
        main = 'Comparison of Inversion Sampling (Blue) \n versus Direct Sampling from PMF (red)')

