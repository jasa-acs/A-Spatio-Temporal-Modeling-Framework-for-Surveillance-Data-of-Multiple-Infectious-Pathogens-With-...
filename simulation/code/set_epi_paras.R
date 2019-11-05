
# This file set the parameters for generating HFMD epidemics in the simulation study

# gamma0: pathogen-specific baseline environment-to-human transmission rate (gamma_E)
# gamma1: pathogen-specific baseline human-to-human transmission rate within the same prefecture (gamma_{H1})
# gamma2: pathogen-specific baseline human-to-human transmission rate between adjacent prefectures (gamma_{H2})
# eta_E: covariate coefficients for environment-to-human transmission (eta_E)
# eta_S: covariate coefficients for human-to-human transmission (eta_H)
# eta_P0: g(p_0)
# eta_P: covariate coefficients for pathogenicity (eta_P)
# c_2: coefficient for the second infectious week (rho)

library(caTools)
library(splines)
library(spdep)

load("data/2009_south5.RData")

n.region <- nregion # 69 # number of regions (R)
n.week <- nweek # 53 # number of weeks (T)
n.virus <- nv # 3 # number of pathogens (V)
n.severity <- ns # 2 # number of severity categories (S)

n.cov <- 1 # number of covariates

##### parameters in one-infectious-week setting #####

# baseline transmission rates
gamma0 <- rep(1.0, n.virus)
gamma1 <- rep(0.5, n.virus)
gamma2 <- rep(0.05, n.virus)

# covariates
covariates <- array(0, dim=c(n.region, n.week, n.cov))
covariates[,,1] <- array(scan("data/covariates/2009_south5_temp_std.txt"), dim=c(n.region, n.week))

# covariate coefficients
eta_P0 <- c(1.61, 1.79, 1.71)
eta_P <- c(-0.15, -0.22, -0.06)
eta_E <- c(-0.05, 0.15, 0.12)
eta_S <- c(0.18, 0.35, -0.18)

# fixed temporal effects
knots <- c(20,30,40)
n.knots <- length(knots)

beta_fix <- matrix(0, n.knots + 3, n.virus)
beta_fix[, 1] <- c(0.0, 1.9, -1.0, 1.0, 0.3, -0.5)
beta_fix[, 2] <- c(0.2, 1.7, -0.7, 1.0, 0.5, -0.2)
beta_fix[, 3] <- c(-0.1, 1.8, -0.5, 0.5, 0.3, 0.0)

basis <- bs(1:n.week, knots=knots, degree=3)
# center the basis w.r.t integral
t.grid <- seq(from=1, to=n.week, by=0.001)
basis.integrand <- bs(t.grid, knots=knots, degree=3)
basis.integral <- rep(0, n.knots + 3)
for (k in 1:(n.knots + 3))
{
	basis.integral[k] <- trapz(t.grid,basis.integrand[, k])
	basis[, k] <- basis[, k] - basis.integral[k] / n.week
}

# spatial effects
sigma2_alpha <- rep(0.01, n.virus)

source("simulation/multinorm_lc.R") # function for sampling from multivariate normal with linear constraint
set.seed(8915)
rho <- 1
W <- nb2mat(nb.prefecture.south5, style = "B", zero.policy=TRUE)
D <- diag(rowSums(W))
alpha_element <- array(0, dim = c(n.region, n.virus)) 
for(v in 1:n.virus) # generate spatial effects from CAR model
{
	Prec <- (D - rho * W) / sigma2_alpha[v]
	alpha_element[, v] <- multiNorm(1, rep(0, n.region), Prec)
}

save(n.region, n.week, n.severity, n.virus, n.cov,
     gamma0, gamma1, gamma2, 
     covariates, eta_P0, eta_P, eta_E, eta_S, 
     n.knots, knots, beta_fix, basis, 
     sigma2_alpha, alpha_element, 
     file="simulation/parameters/south5_simdata_epi_parameters.RData")



##### parameters in two-infectious-week setting #####
c_2 <- 0.2

# baseline transmission rates
gamma0 <- rep(0.92, n.virus)
gamma1 <- rep(0.4, n.virus)
gamma2 <- rep(0.05, n.virus)

# covariates
covariates <- array(0, dim=c(n.region, n.week, n.cov))
covariates[,,1] <- array(scan("data/covariates/2009_south5_temp_std.txt"), dim=c(n.region, n.week))

# covariate coefficients
eta_P0 <- c(1.61, 1.79, 1.71)
eta_P <- c(-0.15, -0.22, -0.06)
eta_E <- c(-0.05, 0.15, 0.12)
eta_S <- c(0.18, 0.35, -0.18)

# fixed temporal effects
knots <- c(20, 30, 40)
n.knots <- length(knots)

beta_fix <- matrix(0, n.knots + 3, n.virus)
beta_fix[, 1] <- c(0.0, 1.9, -1.0, 1.0, 0.3, -0.5)
beta_fix[, 2] <- c(0.2, 1.7, -0.7, 1.0, 0.5, -0.2)
beta_fix[, 3] <- c(-0.1, 1.8, -0.5, 0.5, 0.3, 0.0)

basis <- bs(1:n.week, knots=knots, degree=3)
# center the basis w.r.t integral
t.grid <- seq(from=1, to=53, by=0.001)
basis.integrand <- bs(t.grid, knots=knots, degree=3)
basis.integral <- rep(0, n.knots + 3)
for (k in 1:(n.knots + 3))
{
	basis.integral[k] <- trapz(t.grid,basis.integrand[, k])
	basis[, k] <- basis[, k] - basis.integral[k] / n.week
}

# spatial effects
sigma2_alpha <- rep(0.01, n.virus)

source("simulation/multinorm_lc.R") # function for sampling from multivariate normal with linear constraint
set.seed(8915)
rho <- 1
W <- nb2mat(nb.prefecture.south5, style = "B", zero.policy=TRUE)
D <- diag(rowSums(W))
alpha_element <- array(0, dim = c(n.region, n.virus)) 
for(v in 1:n.virus) # generate spatial effects from CAR model
{
	Prec <- (D - rho * W) / sigma2_alpha[v]
	alpha_element[, v] <- multiNorm(1, rep(0, n.region), Prec)
}

save(c_2, n.region, n.week, n.severity, n.virus, n.cov,
     gamma0, gamma1, gamma2, 
     covariates, eta_P0, eta_P, eta_E, eta_S, 
     n.knots, knots, beta_fix, basis, 
     sigma2_alpha, alpha_element, 
     file="simulation/parameters/south5_simdata_2infweek_epi_parameters.RData")

