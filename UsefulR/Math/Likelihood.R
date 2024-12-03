#Plot likelihood
likelihood<-function(n,z,theta){return(theta^z*(1-theta)^(n-z))}
theta<-seq(from=.01,to=.99,by=.01)
plot(theta,likelihood(400,12,theta))
abline(v=12/400)
plot(theta,likelihood(400,12,theta),type='l')
abline(v=12/400)

#Plot log-likelihood
loglikelihood<-function(n,z,theta){return(z*log(theta)+(n-z)*log(1-theta))}
plot(theta,loglikelihood(400,12,theta),type='l')
abline(v=12/400)
