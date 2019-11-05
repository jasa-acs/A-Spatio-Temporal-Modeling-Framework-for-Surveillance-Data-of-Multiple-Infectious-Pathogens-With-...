
# This file generates one replicate of simulation data with two infectious weeks in the balanced setting, two imbalanced setting, and real data setting
# This file also prepares the data in .txt format as inputs of C code

# data generator for balanced setting with 2 infectious weeks
# gamma0: pathogen-specific baseline environment-to-human transmission rate (gamma_E in paper)
# gamma1: pathogen-specific baseline human-to-human transmission rate within the same prefecture (gamma_{H1} in paper)
# gamma2: pathogen-specific baseline human-to-human transmission rate between adjacent prefectures (gamma_{H2} in paper)
# eta_E: covariate coefficients for environment-to-human transmission (eta_E in paper)
# eta_S: covariate coefficients for human-to-human transmission (eta_H in paper)
# eta_P0: g(p_0) in paper
# eta_P: covariate coefficients for pathogenicity (eta_P in paper)
# c_2: rho in paper

library(spdep)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) 
{
	stop("Please provide an index for the simulation\n", call.=FALSE)
} else {
	index.rep <- args[1]
}

load("data/2009_south5.RData") # need neighborhood structure from real data
load("simulation/parameters/south5_simdata_2infweek_epi_parameters.RData")

cloglog <- function(x) log(-log(1-x))
inv.cloglog <- function(x) 1 - exp(-exp(x))


# covariate effects
cov.E <- array(0, dim=c(n.region, n.week, n.virus))
cov.S <- array(0, dim=c(n.region, n.week, n.virus))
cov.P <- array(0, dim=c(n.region, n.week, n.virus))

for (r in 1:n.region)
{
	for (t in 1:n.week)
	{
		cov.E[r,t,] <- covariates[r,t,]%*%eta_E
		cov.S[r,t,] <- covariates[r,t,]%*%eta_S
		cov.P[r,t,] <- covariates[r,t,]%*%eta_P
		cov.P[r,t,] <- cov.P[r,t,] + eta_P0
	}
}

# spatial effects
alpha <- array(0, dim=c(n.region, n.region, n.virus))   
nu <- 0
for(v in 1:n.virus)
{
	for(j in 1:n.region)
	{
		for(r in 1:n.region)
		{
			alpha[r,j,v] <- alpha_element[r,v] + nu * alpha_element[j,v]
		}
	}
}

# calculate temporal effects #

beta <- array(0, dim=c(n.region, n.week, n.virus))

for (r in 1:n.region)
{
	for (v in 1:n.virus)
	{
		beta[r,,v] <- basis %*% beta_fix[, v]
	}
}

# transmission rates 
W <- nb2mat(nb.prefecture.south5, style = "B", zero.policy=TRUE)
gammaP2P <- array(0, dim=c(n.region,n.region,n.week,n.virus)) # the effective transmission rate 
for(v in 1:n.virus)
{
  for(r in 1:n.region)
  {
    for(j in 1:n.region)
    {
      for(t in 1:n.week)
      {
           if (r==j) gammaP2P[r, j, t, v] <- exp(log(gamma1[v]) + alpha[r,j,v] + beta[r,t,v] + cov.S[r,t,v])
           else gammaP2P[r, j, t, v] <- exp(log(gamma2[v]) + alpha[r,j,v] + beta[r,t,v] + cov.S[r,t,v])*W[r,j]
      } 
    }
  }
}

gammaRT <- array(0, dim=c(n.region, n.week, n.virus))
for (r in 1:n.region)
{
	for (t in 1:n.week)
	{
		for (v in 1:n.virus)
		{
			gammaRT[r,t,v] <- gamma0[v]*exp(cov.E[r,t,v])
		}
	}
}


# number of cases #
Y.init <- array(rbinom(n.region*n.virus,1,0.5), dim=c(n.region, n.virus))   # initial number of infected people
Lambda <- array(0, dim=c(n.region, n.week, n.virus))
Y.total <- array(0, dim=c(n.region, n.week, n.virus)) # cases

Y.total[,1,] <- Y.init

for(v in 1:n.virus)
{
    for(r in 1:n.region)
    { 
    	Br <- nb.prefecture.south5[[r]]
		Lambda[r,2,v] <- gammaRT[r,2,v] + Y.total[r,1,v] * gammaP2P[r,r,2,v] + sum(gammaP2P[r,Br,2,v]*Y.total[Br,1,v])
		Y.total[r,2,v] <- rpois(1, Lambda[r,2,v])  
    }
}

for(v in 1:n.virus)
{
  for(t in 3:n.week)
  {
    for(r in 1:n.region)
    { 
    	Br <- nb.prefecture.south5[[r]]
		Lambda[r,t,v] <- gammaRT[r,t,v] + (c_2 * Y.total[r,(t-2),v] + Y.total[r,(t-1),v]) * gammaP2P[r,r,t,v] + sum(gammaP2P[r,Br,t,v]*(Y.total[Br,(t-1),v] + c_2 * Y.total[Br,(t-2),v]))
		Y.total[r,t,v] <- rpois(1, Lambda[r,t,v])  
    }
  }
}

Y <- array(0, dim=c(n.region, n.severity, n.week, n.virus)) # cases (virus specific)

ps <- 1 - inv.cloglog(cov.P) # probability of getting a severe case
for(v in 1:n.virus){
  for(r in 1:n.region){
    for(t in 1:n.week){
      Y[r,2,t,v] <- rbinom(1, round(Y.total[r,t,v]), ps[r,t,v])
      Y[r,1,t,v] <-  Y.total[r,t,v] - Y[r,2,t,v]
    }
  }
}

Yv <- apply(Y,c(1,2,3),sum)

for (t in 1:n.week)
{
	for (s in 1:n.severity)
	{
		cat(Yv[,s,t], "\n", file=paste("simulation/sim_data/south5_Yv_simdata_2infweek_rep", index.rep, ".txt", sep=""), append=TRUE)
	}	
}

## sample lab data

# balanced design

load("simulation/parameters/south5_simdata_lab_parameters_balanced.RData")

p_m <- lab_para_m_cand[1]
Z <- array(0, dim=c(n.region,n.severity,n.week,n.virus)) # virus specific lab tests
for(r in 1:n.region)
{
	for(t in 1:n.week)
    {
    	for(v in 1:n.virus)
    	{		
    		Z[r,1,t,v] <- rbinom(1,Y[r,1,t,v], p_m)
    		Z[r,2,t,v] <- rbinom(1,Y[r,2,t,v], p_s)
		}
  	}
}

for (v in 1:n.virus)
{
	for (t in 1:n.week)
	{
		for (s in 1:n.severity)
		{
			cat(Z[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_simdata_2infweek_balanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}
save(c_2, gamma0, gamma1, gamma2, sigma2_alpha, covariates, eta_E, eta_S, eta_P, eta_P0, alpha_element, beta_fix, Y, Yv, Z, 
     file=paste("simulation/sim_data/south5_simdata_2infweek_balanced_combo1_rep", index.rep, ".RData", sep=""))


# imbalanced design
load("simulation/parameters/south5_simdata_lab_parameters_imbalanced.RData")

p_m <- lab_para_m_cand[1, 2]
cutoff <- lab_para_m_cand[1, 1]

Z <- array(0, dim=c(n.region,n.severity,n.week,n.virus)) # virus specific lab tests
for(r in 1:n.region)
{
	for(t in 1:n.week)
    {
    	for(v in 1:n.virus)
    	{		
    		Z[r,1,t,v] <- ifelse(lab_zero_prob_m[r,t] < cutoff, 1, 0) * rbinom(1,Y[r,1,t,v], p_m)
		    Z[r,2,t,v] <- rbinom(1,Y[r,2,t,v], p_s)
		}
  	}
}

for (v in 1:n.virus)
{
	for (t in 1:n.week)
	{
		for (s in 1:n.severity)
		{
			cat(Z[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_simdata_2infweek_imbalanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}
save(c_2, gamma0, gamma1, gamma2, sigma2_alpha, covariates, eta_E, eta_S, eta_P, eta_P0, alpha_element, beta_fix, Y, Yv, Z, 
     file=paste("simulation/sim_data/south5_simdata_2infweek_imbalanced_combo1_rep", index.rep, ".RData", sep=""))


