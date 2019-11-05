###################################
##### Xueying Tang 12/01/2018 #####
###################################

# This file read txt outputs from C into RData files

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) 
{
	stop("Please provide index_chain as an extra argument!\n", call.=FALSE)
} else {
	index_chain <- as.integer(args[1])
}

nsim <- 25

nregion <- 69
nweek <- 53
nv <- 3
ns <- 2

nknots <- 3
ncov <- 12

eta_S.chain <- array(scan(paste("real_data_analysis/results/eta_S_output_realdata_", index_chain, ".txt", sep="")), dim=c(ncov,nv,nsim))
eta_E.chain <- array(scan(paste("real_data_analysis/results/eta_E_output_realdata_", index_chain, ".txt", sep="")), dim=c(ncov,nv,nsim))
eta_P0.chain <- matrix(scan(paste("real_data_analysis/results/eta_P0_output_realdata_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
Y.chain <- array(scan(paste("real_data_analysis/results/Y_output_realdata_", index_chain, ".txt", sep="")), dim=c(nregion, ns, nweek, nv, nsim))
gamma0.chain <- matrix(scan(paste("real_data_analysis/results/gamma0_output_realdata_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
gamma1.chain <- matrix(scan(paste("real_data_analysis/results/gamma1_output_realdata_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
gamma2.chain <- matrix(scan(paste("real_data_analysis/results/gamma2_output_realdata_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
sigma2_alpha.chain <- matrix(scan(paste("real_data_analysis/results/sigma2_alpha_output_realdata_", index_chain, ".txt", sep="")), nrow=nv, ncol=nsim)
alpha_element.chain <- array(scan(paste("real_data_analysis/results/alpha_element_output_realdata_", index_chain, ".txt", sep="")), dim=c(nregion, nv, nsim))
beta_fix.chain <- array(scan(paste("real_data_analysis/results/beta_fix_output_realdata_", index_chain, ".txt", sep="")), dim=c(nknots+3, nv, nsim))

save(alpha_element.chain, sigma2_alpha.chain, Y.chain, eta_P0.chain, eta_S.chain, eta_E.chain, gamma0.chain, gamma1.chain, gamma2.chain, beta_fix.chain, 
     file=paste("real_data_analysis/results/south5_realdata_", index_chain, "_results.RData", sep=""))

