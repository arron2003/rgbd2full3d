library(fitdistrplus)

dtnorm<- function(x, mean, sd, a, b) {
dnorm(x, mean, sd)/(pnorm(b, mean, sd)-pnorm(a, mean, sd))
}

ptnorm <- function(x, mean, sd, a, b) {
(pnorm(x,mean,sd) - pnorm(a,mean,sd)) / 
  (pnorm(b,mean,sd) - pnorm(a,mean,sd))
}

mydf = read.table("tmp.txt")
mydata = as.vector(as.matrix(mydf));
para = fitdist( mydata, dtnorm, method="mle", start=list(mean=0, sd=1), fix.arg=list(a=0,b=10));

write.table(para$estimate, 'tmp.txt',  sep = " ", row.names = F, col.names = F);