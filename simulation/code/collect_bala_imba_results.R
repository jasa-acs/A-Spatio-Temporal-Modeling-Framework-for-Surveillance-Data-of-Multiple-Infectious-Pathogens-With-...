rm(list = ls())

nregion <- 69
nweek <- 53
ns <- 2
nv <- 3
nknots <- 3
ncov <- 1

nepi <- 100
nchain <- 5
nsize <- 10000
npara <- 8

load("simulation/parameters/south5_simdata_epi_parameters.RData")

eta_E <- array(eta_E, dim=c(1, nv))
eta_S <- array(eta_S, dim=c(1, nv))
eta_P <- array(eta_P, dim=c(1, nv))

bmse <- function(mylist)
{
	x0 <- mylist$true
	x_chain <- mylist$mean
	x0_ex <- array(as.vector(x0), dim=dim(x_chain))
	sd_chain <- mylist$sd
	nd <- length(dim(x0))
	if (nd==0) nd <- 1
	x_overall <- apply(x_chain, 1:nd, mean)
	
	part1 <- (x_overall - x0)^2
	part2 <- (1-1/nchain)*apply(x_chain, 1:nd, var)
	part3 <- apply(sd_chain^2*(1-1/nsize), 1:nd, mean)
	
	part1+part2+part3
}


alpha_element.bmse <- array(0, dim=c(nregion, nv, npara, nepi))
beta_fix.bmse <- array(0, dim=c(nknots+3, nv, npara, nepi))
eta_P0.bmse <- array(0, dim=c(nv, npara, nepi))
eta_P.bmse <- array(0, dim=c(ncov, nv, npara, nepi))
eta_E.bmse <- array(0, dim=c(ncov, nv, npara, nepi))
eta_S.bmse <- array(0, dim=c(ncov, nv, npara, nepi))
gamma0.bmse <- array(0, dim=c(nv, npara, nepi))
gamma1.bmse <- array(0, dim=c(nv, npara, nepi))
gamma2.bmse <- array(0, dim=c(nv, npara, nepi))
sigma2_alpha.bmse <- array(0, dim=c(nv, npara, nepi))

alpha_element.stat <- array(0, dim=c(nregion, nv, nepi, npara))
beta_fix.stat <- array(0, dim=c(nknots+3, nv, nepi, npara))
eta_P.stat <- array(0, dim=c(ncov, nv, nepi, npara))
eta_P0.stat <- array(0, dim=c(nv, nepi, npara))
eta_S.stat <- array(0, dim=c(ncov, nv, nepi, npara))
eta_E.stat <- array(0, dim=c(ncov, nv, nepi, npara))
gamma0.stat <- array(0, dim=c(nv, nepi, npara))
gamma1.stat <- array(0, dim=c(nv, nepi, npara))
gamma2.stat <- array(0, dim=c(nv, nepi, npara))
sigma2_alpha.stat <- array(0, dim=c(nv, nepi, npara))

alpha_element.means <- array(0, dim=c(nregion, nv, npara, nepi))
beta_fix.means <- array(0, dim=c(nknots+3, nv, npara, nepi))
eta_P0.means <- array(0, dim=c(nv, npara, nepi))
eta_P.means <- array(0, dim=c(ncov, nv, npara, nepi))
eta_E.means <- array(0, dim=c(ncov, nv, npara, nepi))
eta_S.means <- array(0, dim=c(ncov, nv, npara, nepi))
gamma0.means <- array(0, dim=c(nv, npara, nepi))
gamma1.means <- array(0, dim=c(nv, npara, nepi))
gamma2.means <- array(0, dim=c(nv, npara, nepi))
sigma2_alpha.means <- array(0, dim=c(nv, npara, nepi))


for (index.epi in 1:nepi)
{
	for (index.para in 1:npara)
	{
		alpha_element.mean <- array(0, dim=c(nregion, nv, nchain))
		beta_fix.mean <- array(0, dim=c(nknots+3, nv, nchain))
		eta_P.mean <- array(0, dim=c(ncov, nv, nchain))
		eta_P0.mean <- array(0, dim=c(nv, nchain))
		eta_S.mean <- array(0, dim=c(ncov, nv, nchain))
		eta_E.mean <- array(0, dim=c(ncov, nv, nchain))
		gamma0.mean <- array(0, dim=c(nv, nchain))
		gamma1.mean <- array(0, dim=c(nv, nchain))
		gamma2.mean <- array(0, dim=c(nv, nchain))
		sigma2_alpha.mean <- array(0, dim=c(nv, nchain))
		
		
		alpha_element.sd <- array(0, dim=c(nregion, nv, nchain))
		beta_fix.sd <- array(0, dim=c(nknots+3, nv, nchain))
		eta_P.sd <- array(0, dim=c(ncov, nv, nchain))
		eta_P0.sd <- array(0, dim=c(nv, nchain))
		eta_S.sd <- array(0, dim=c(ncov, nv, nchain))
		eta_E.sd <- array(0, dim=c(ncov, nv, nchain))
		gamma0.sd <- array(0, dim=c(nv, nchain))
		gamma1.sd <- array(0, dim=c(nv, nchain))
		gamma2.sd <- array(0, dim=c(nv, nchain))
		sigma2_alpha.sd <- array(0, dim=c(nv, nchain))

		
		for (index_chain in 1:nchain)
		{
			if (index.para <= 3) load(paste("simulation/output_chains/south5_simdata_balanced_combo", index.para, "_block_", index.epi, "_", index.chain, "_results.RData", sep=""))
			else load(paste("simulation/output_chains/south5_simdata_imbalanced_combo", index.para-3, "_block_", index.epi, "_", index.chain, "_results.RData", sep=""))
				
			alpha_element.mean[ , , index.chain] <- apply(alpha_element.chain, 1:2, mean)
			beta_fix.mean[ , , index.chain] <- apply(beta_fix.chain, 1:2, mean)
			eta_P.mean[ , , index.chain] <- apply(eta_P.chain, 1:2, mean)
			eta_P0.mean[ , index.chain] <- apply(eta_P0.chain, 1, mean)
			eta_S.mean[ , , index.chain] <- apply(eta_S.chain, 1:2, mean)
			eta_E.mean[ , , index.chain] <- apply(eta_E.chain, 1:2, mean)
			gamma0.mean[ , index.chain] <- apply(gamma0.chain, 1, mean)
			gamma1.mean[ , index.chain] <- apply(gamma1.chain, 1, mean)
			gamma2.mean[ , index.chain] <- apply(gamma2.chain, 1, mean)
			sigma2_alpha.mean[ , index.chain] <- apply(sigma2_alpha.chain, 1, mean)
			
		
			alpha_element.sd[ , , index.chain] <- apply(alpha_element.chain, 1:2, sd)
			beta_fix.sd[ , , index.chain] <- apply(beta_fix.chain, 1:2, sd)
			eta_P.sd[ , , index.chain] <- apply(eta_P.chain, 1:2, sd)
			eta_P0.sd[ , index.chain] <- apply(eta_P0.chain, 1, sd)
			eta_S.sd[ , , index.chain] <- apply(eta_S.chain, 1:2, sd)
			eta_E.sd[ , , index.chain] <- apply(eta_E.chain, 1:2, sd)
			gamma0.sd[ , index.chain] <- apply(gamma0.chain, 1, sd)
			gamma1.sd[ , index.chain] <- apply(gamma1.chain, 1, sd)
			gamma2.sd[ , index.chain] <- apply(gamma2.chain, 1, sd)
			sigma2_alpha.sd[ , index.chain] <- apply(sigma2_alpha.chain, 1, sd)

		}
				
		alpha_element.bmse[,,index.para, index.epi] <- bmse(list(true=alpha_element, mean=alpha_element.mean, sd=alpha_element.sd))
		beta_fix.bmse[,,index.para, index.epi] <- bmse(list(true=beta_fix, mean=beta_fix.mean, sd=beta_fix.sd))
		eta_P0.bmse[,index.para, index.epi] <- bmse(list(true=eta_P0, mean=eta_P0.mean, sd=eta_P0.sd))
		eta_P.bmse[,,index.para, index.epi] <- bmse(list(true=eta_P, mean=eta_P.mean, sd=eta_P.sd))
		eta_E.bmse[,,index.para, index.epi] <- bmse(list(true=eta_E, mean=eta_E.mean, sd=eta_E.sd))
		eta_S.bmse[,,index.para, index.epi] <- bmse(list(true=eta_S, mean=eta_S.mean, sd=eta_S.sd))
		gamma0.bmse[,index.para, index.epi] <- bmse(list(true=gamma0, mean=gamma0.mean, sd=gamma0.sd))
		gamma1.bmse[,index.para, index.epi] <- bmse(list(true=gamma1, mean=gamma1.mean, sd=gamma1.sd))
		gamma2.bmse[,index.para, index.epi] <- bmse(list(true=gamma2, mean=gamma2.mean, sd=gamma2.sd))
		sigma2_alpha.bmse[,index.para, index.epi] <- bmse(list(true=sigma2_alpha, mean=sigma2_alpha.mean, sd=sigma2_alpha.sd))
				
		alpha_element.stat[ , , index.epi, index.para] <- apply(alpha_element.mean, 1:2, var)/apply(alpha_element.sd^2, 1:2, mean)
		beta_fix.stat[ , , index.epi, index.para] <- apply(beta_fix.mean, 1:2, var)/apply(beta_fix.sd^2, 1:2, mean)
		eta_P.stat[ , , index.epi, index.para] <- apply(eta_P.mean, 1:2, var)/apply(eta_P.sd^2, 1:2, mean)
		eta_P0.stat[ , index.epi, index.para] <- apply(eta_P0.mean, 1, var)/apply(eta_P0.sd^2, 1, mean)
		eta_S.stat[ , , index.epi, index.para] <- apply(eta_S.mean, 1:2, var)/apply(eta_S.sd^2, 1:2, mean)
		eta_E.stat[ , , index.epi, index.para] <- apply(eta_E.mean, 1:2, var)/apply(eta_E.sd^2, 1:2, mean)
		gamma0.stat[ , index.epi, index.para] <- apply(gamma0.mean, 1, var)/apply(gamma0.sd^2, 1, mean)
		gamma1.stat[ , index.epi, index.para] <- apply(gamma1.mean, 1, var)/apply(gamma1.sd^2, 1, mean)
		gamma2.stat[ , index.epi, index.para] <- apply(gamma2.mean, 1, var)/apply(gamma2.sd^2, 1, mean)
		sigma2_alpha.stat[ , index.epi, index.para] <- apply(sigma2_alpha.mean, 1, var)/apply(sigma2_alpha.sd^2, 1, mean)

		alpha_element.means[,,index.paras, index.epi] <- apply(alpha_element.mean, 1:2, mean)
		beta_fix.means[,,index.paras, index.epi] <- apply(beta_fix.mean, 1:2, mean)
		eta_P0.means[,index.paras, index.epi] <- apply(eta_P0.mean, 1, mean)
		eta_P.means[,,index.paras, index.epi] <- apply(eta_P.mean, 1:2, mean)
		eta_E.means[,,index.paras, index.epi] <- apply(eta_E.mean, 1:2, mean)
		eta_S.means[,,index.paras, index.epi] <- apply(eta_S.mean, 1:2, mean)
		gamma0.means[,index.paras, index.epi] <- apply(gamma0.mean, 1, mean)
		gamma1.means[,index.paras, index.epi] <- apply(gamma1.mean, 1, mean)
		gamma2.means[,index.paras, index.epi] <- apply(gamma2.mean, 1, mean)
		sigma2_alpha.means[,index.paras, index.epi] <- apply(sigma2_alpha.mean, 1, mean)
	}
}


bmse_res <- list(alpha_element=alpha_element.bmse, 
				 beta_fix=beta_fix.bmse,
				 eta_P0=eta_P0.bmse,
				 eta_P=eta_P.bmse,
				 eta_E=eta_E.bmse,
				 eta_S=eta_S.bmse,
				 gamma0=gamma0.bmse,
				 gamma1=gamma1.bmse,
				 gamma2=gamma2.bmse,
				 sigma2_alpha=sigma2_alpha.bmse)


stat_res <- list(alpha_element=alpha_element.stat, 
				 beta_fix=beta_fix.stat,
				 eta_P0=eta_P0.stat,
				 eta_P=eta_P.stat,
				 eta_E=eta_E.stat,
				 eta_S=eta_S.stat,
				 gamma0=gamma0.stat,
				 gamma1=gamma1.stat,
				 gamma2=gamma2.stat,
				 sigma2_alpha=sigma2_alpha.stat) 	


save(bmse_res, stat_res, file="simulation/results_summary/summary_simulation_bala_imba_results.RData")

save(alpha_element.means, beta_fix.means, eta_P0.means, eta_P.means, eta_E.means, eta_S.means, gamma0.means, gamma1.means, gamma2.means, sigma2_alpha.means, file="simulation/results_summary/summary_means_simulation_bala_imba_results.RData")



##### Y
impute_Y <- function(Yv, Z)
{
	Y.imputed <- array(0, dim=c(nregion, ns, nweek, nv))
	for (r in 1:nregion)
	{
		for (s in 1:ns)
		{
			for (t in 1:nweek)
			{
				if (Yv[r,s,t] != 0)
				{
					Zsum <- sum(Z[r,s,t,])
					if (Zsum != 0) {
						Y.imputed[r,s,t,] <- floor(Yv[r,s,t] * Z[r,s,t,]/Zsum)
					} else {
						Y.imputed[r,s,t,] <- floor(Yv[r,s,t] / 3)
					}
					Y.imputed[r,s,t,nv] <- Yv[r,s,t] - sum(Y.imputed[r,s,t,-nv])
				}
			}
		}
	}
	
	Y.imputed
}

Y.true.all <- array(0, dim=c(nweek, nv, nepi))
Y.est.all <- array(0, dim=c(npara, nweek, nv, nepi))
Y.imputed.all <- array(0, dim=c(npara, nweek, nv, nepi))

r <- 1
s <- 1

for (index.epi in 1:nepi)
{
	load(paste("simulation/sim_data/south5_simdata_balanced_combo1_rep", index.epi, ".RData", sep=""))
	Y.true.all[,,index.epi] <- Y[r,s,,]
}

for (index.epi in 1:nepi)
{
	for (index.para in 1:npara)
	{
		if (index.para <=3) {
			load(paste("simulation/output_chains/south5_simdata_balanced_combo", index.para, "_rep", index.epi, ".RData", sep=""))
		} else {
			load(paste("simulation/output_chains/south5_simdata_imbalanced_combo", index.para-3, "_rep", index.epi, ".RData", sep=""))
		}
		
		Y.imputed <- impute_Y(Yv, Z)
		Y.imputed.all[index.para, , , index.epi] <- Y.imputed[r,s,,]
	}
}

for (index.epi in 1:nepi)
{
	for (index.para in 1:npara)
	{
		Y.est <- array(0, dim=c(nregion, ns, nweek, nv, nchain))
		for (index.chain in 1:nchain)
		{
			if (index.para <= 3) load(paste("simulation/output_chains/south5_simdata_balanced_combo", index.paras, "_block_", index.epi, "_", index.chain, "_results.RData", sep=""))
			else load(paste("simulation/output_chains/south5_simdata_imbalanced_combo", index.paras-3, "_block_", index.epi, "_", index.chain, "_results.RData", sep=""))
			Y.est[,,,,index.chain] <- apply(Y.chain, 1:4, mean)
		}
		Y.est.all[index.para, , , index.e] <- apply(Y.est, 1:4, mean)[r,s,,]
	}
}

plotvar_est <- array(0, dim=c(npara, nepi, nv))
plotvar_imp <- array(0, dim=c(npara, nepi, nv))

for (v in 1:nv)
{
	for (index.paras in 1:npara)
	{
		plotvar_est[index.para,,v] <- apply(Y.est.all[index.para,,v,] - Y.true.all[,v,], 2, mean)
		plotvar_imp[index.para,,v] <- apply(Y.imputed.all[index.para,,v,] - Y.true.all[,v,], 2, mean)
	}
}
save(plotvar_est, plotvar_imp, file="simulation/results_summary/summary_Y_mse_simulation_bala_imba_results.RData")
