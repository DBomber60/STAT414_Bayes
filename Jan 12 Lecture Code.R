# data frame on student survey
wd = "~/Documents/STATS414/Jan 12"
setwd(wd)
data = read.table("studentdata.txt", header=T)
boxplot(data$Shoes ~ data$Gender)
boxplot(data$Job ~ data$Gender)

# Can we predict the cost of haircut from gender, shoes, number of hours on job?
model1 = lm(data$Haircut ~ data$Gender + data$Shoes + data$Job)
summary(model1)

# Can we explain how many shoes are purchased based on job, gender?
nb.Shoes = floor(data$Shoes)
model2 = glm(nb.Shoes ~ data$Job + data$Gender, family=poisson)
summary(model2)

# congress data - each column is a bill; each row is a senator; value=1 is the sponsor; value=2 is cosponsor
# let's measure how popular certain members of congress are
senmat = read.table("senmatrix.txt", header=F, sep=',')
ptm <- proc.time()
n.sen = length(senmat[,1])  # number of senators - 102
n.bill = length(senmat[1,]) # number of bills - 10327
bill_sp = numeric(n.sen)    # number of bills sponsored by each senator
V = double(n.sen)           # each entry of V is a senator
for (i in 1:n.sen) {
  bill_sp[i] = length(which(senmat[i,]==1))
  if (bill_sp[i]==0) {V[i]=0}
  else {
    cosponsors = length(which(senmat[,which(senmat[i,]==1)]==2))
    V[i] = cosponsors/bill_sp[i]
  }
}
proc.time() - ptm

# Slightly more vectorized approach to the same problem (use apply to find the vector of total bills sponsored for each senator)
# Runs roughly twice as fast as for loop approach above
ptm <- proc.time()
sponsor = function(x) length(which(x==1))
bills_sponsored = apply(senmat,1,sponsor)
b = rep(0,102)
for(i in 1:102) {
  if (bills_sponsored[i]==0) {b[i]=0}
  else {b[i] = length(which(senmat[,which(senmat[i,]==1)]==2))}
}
print (b/bills_sponsored)
proc.time() - ptm