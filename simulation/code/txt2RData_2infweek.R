###################################
##### Xueying Tang 12/01/2018 #####
###################################

# This file read txt outputs from C into RData files

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) 
{
	stop("Please provide data label, index_epi, index_chain as extra arguments!\n", call.=FALSE)
} else {
	data_label <- args[1]
	index_epi <- as.integer(args[2])
	index_chain <- as.integer(args[3])
}

nsim <- 25

nregion <- 69
nweek <- 53
nv <- 3
ns <- 2

nknots <- 3
ncov <- 1

eta_S.chain <- array(scan(paste("simulation/output_chains/eta_S_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")),dim=c(ncov,nv,nsim))
eta_E.chain <- array(scan(paste("simulation/output_chains/eta_E_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")),dim=c(ncov,nv,nsim))
eta_P.chain <- array(scan(paste("simulation/output_chains/eta_P_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")),dim=c(ncov,nv,nsim))
eta_P0.chain <- matrix(scan(paste("simulation/output_chains/eta_P0_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
Y.chain <- array(scan(paste("simulation/output_chains/Y_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), dim=c(nregion, ns, nweek, nv, nsim))
gamma0.chain <- matrix(scan(paste("simulation/output_chains/gamma0_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
gamma1.chain <- matrix(scan(paste("simulation/output_chains/gamma1_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
gamma2.chain <- matrix(scan(paste("simulation/output_chains/gamma2_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
sigma2_alpha.chain <- matrix(scan(paste("simulation/output_chains/sigma2_alpha_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
alpha_element.chain <- array(scan(paste("simulation/output_chains/alpha_element_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), dim=c(nregion, nv, nsim))
beta_fix.chain <- array(scan(paste("simulation/output_chains/beta_fix_output_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, ".txt", sep="")), dim=c(nknots+3, nv, nsim))

save(alpha_element.chain, sigma2_alpha.chain, Y.chain, eta_P0.chain, eta_P.chain, eta_S.chain, eta_E.chain, gamma0.chain, gamma1.chain, gamma2.chain, beta_fix.chain, 
     file=paste("simulation/output_chains/south5_simdata_2infweek_", data_label, "_", index_epi, "_", index_chain, "_results.RData", sep=""))

