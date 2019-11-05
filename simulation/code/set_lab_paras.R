
# This file set the parameters for generating HFMD epidemics in the simulation study

# gamma0: pathogen-specific baseline environment-to-human transmission rate (gamma_E)
# gamma1: pathogen-specific baseline human-to-human transmission rate within the same prefecture (gamma_{H1})
# gamma2: pathogen-specific baseline human-to-human transmission rate between adjacent prefectures (gamma_{H2})
# eta_E: covariate coefficients for environment-to-human transmission (eta_E)
# eta_S: covariate coefficients for human-to-human transmission (eta_H)
# eta_P0: g(p_0)
# eta_P: covariate coefficients for pathogenicity (eta_P)
# c_2: coefficient for the second infectious week (rho)


## balanced sampling proportion; same across pathogens
p_s <- 0.7 # sampling proportion for severe cases
lab_para_m_cand <- c(0.02, 0.05, 0.10) # sampling proportion for mild cases

save(p_s, lab_para_m_cand, file="simulation/parameters/south5_simdata_lab_parameters_balanced.RData")

# imbalanced sampling proportion
p_s <- 0.7 # sampling proportion for severe cases
# parameters for sampling mild cases

# fit zero-lab-case probability curve for mild cases from real data
library(mgcv)

load("data/2009_south5.RData")

n.region <- nregion
n.week <- nweek
n.severity <- ns
n.virus <- nv

Zv <- apply(Z, 1:3, sum)

mdata2009 <- data.frame(Yv=as.vector(Yv[,1,]), Zv=as.vector(Zv[,1,]), Week=rep(1:n.week, each=n.region), Region=rep(1:n.region, time=n.week))
mdata2009_nonzeroY <- mdata2009[mdata2009$Yv != 0, ]
mdata2009_nonzeroY$Zv_ind <- as.numeric(mdata2009_nonzeroY$Zv == 0)
fit_m <- gam(Zv_ind ~ as.factor(Region) + s(Week, bs="cr", k=13), family="binomial", data=mdata2009_nonzeroY)
lab_zero_prob_m <- matrix(predict(fit_m, 
newdata=data.frame(Week=rep(1:n.week, each=n.region), Region=as.factor(rep(1:n.region, time=n.week))), 
type="response"), nrow=n.region, ncol=n.week)

# set parameters
lab_para_m_cand <- matrix(c(0.4, 0.4, 0.68, 0.4, 0.905,
					        0.133, 0.325, 0.133, 0.65, 0.133), nrow=5, ncol=2)
# column 1: threshold for sampling lab cases; if zero-lab-case probility > threshold, no lab cases will be sampled
# column 2: sampling proportion
save(p_s, lab_para_m_cand, lab_zero_prob_m, file="simulation/parameters/south5_simdata_lab_parameters_imbalanced.RData")

# sampling proportion resembling the real data
load("data/2009_south5.RData")
n.region <- nregion
n.week <- nweek
n.severity <- ns
n.virus <- nv

Zv <- apply(Z, 1:3, sum)

lab_prop <- array(0, dim=c(n.region, n.severity, n.week))

for (r in 1:n.region)
{
	for (s in 1:n.severity)
	{
		for (t in 1:n.week)
		{
			if (Yv[r,s,t] != 0) lab_prop[r,s,t] <- Zv[r,s,t]/Yv[r,s,t]
		}
	}
}
save(lab_prop, file="simulation/parameters/south5_simdata_lab_parameters_real.RData")


