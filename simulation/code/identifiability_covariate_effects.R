b<-0.003
p<-0.0004
duration.max<-100
duration.infectious.period<-2

x<-rbinom(duration.max, 1, 0.5)  # time-dependent covariates
beta0 <- log(2)  # exp(beta0) is the odds ratio for the covariate  
                 # effect on environment-to-human transmission
beta1 <- log(1.5)  # exp(beta1) is the odds ratio for the covariate  
                   # effect on human-to-human transmission
pop.size <- 1000
n.sim <- 500

logit <- function(p) {log(p/(1-p))}
inv.logit <- function(x) {exp(x)/(1+exp(x))}
log.likelihood <- function(para, covariate, n_susceptible, 
                           n_infectious, n_new_infection)
{
   b <- inv.logit(para[1])
   p <- inv.logit(para[2])
   beta0 <- para[3]
   beta1 <- para[4]
   b.effective <- inv.logit(logit(b) + beta0 * covariate)
   p.effective <- inv.logit(logit(p) + beta1 * covariate)
   log.prob.escape <- log(1 - b.effective) + 
                      n_infectious * log(1 - p.effective)
   log.prob.inf <- log(1-exp(log.prob.escape))
   sum( n_new_infection * log.prob.inf + 
       (n_susceptible - n_new_infection) * log.prob.escape)
}
est <- NULL 
size.epidemic <- rep(NA, n.sim)  
for(iter in 1:n.sim)
{
   print(iter)
   set.seed(iter * 10)
   n.susceptible <- c(pop.size, rep(0, duration.max-1))
   n.infectious <- n.new.infection <- n.recovered <- rep(0, duration.max)
   duration.epidemic <- duration.max
   for(t in 1:duration.max)
   {
      b.effective <- inv.logit(logit(b) + beta0 * x[t])
      p.effective <- inv.logit(logit(p) + beta1 * x[t])
      prob.escape <- (1-b.effective) * ((1-p.effective)^n.infectious[t]) 
      n.new.infection[t] <- rbinom(1, n.susceptible[t], 1-prob.escape) 
      if(t < duration.max)
      {      
         n.susceptible[t+1] <- n.susceptible[t] - n.new.infection[t]
         n.infectious[t+1] <- n.infectious[t] + n.new.infection[t]
         if(t > duration.infectious.period)  
         {
            n.infectious[t+1]<-n.infectious[t+1] - 
                               n.new.infection[t-duration.infectious.period]
            n.recovered[t+1]<-n.recovered[t] + 
                              n.new.infection[t-duration.infectious.period]
         }   
         if(n.susceptible[t+1] == 0)  {duration.epidemic <- t; break}
      }
   }
   size.epidemic[iter] <- sum(n.new.infection)
   est.ini <- c(logit(b), logit(p), beta0, beta1)
   fit <- optim(est.ini, log.likelihood, covariate=x[1:duration.epidemic],  
                n_susceptible=n.susceptible[1:duration.epidemic], 
                n_infectious=n.infectious[1:duration.epidemic], 
                n_new_infection=n.new.infection[1:duration.epidemic], 
                method="BFGS", 
                control=list(fnscale=-1,reltol=1e-8,maxit=200,trace=0))
   est <- rbind(est, c(inv.logit(fit$par[1]), inv.logit(fit$par[2]), 
                exp(fit$par[3]), exp(fit$par[4])))
}
summary(size.epidemic)
apply(est, 2, mean)
apply(est, 2, sd)
