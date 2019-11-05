################################################################
# This file contains R programs for plotting all the figures for
# simulation study
################################################################
rm(list = ls())

nregion <- 69
nweek <- 53
ns <- 2
nv <- 3
nknots <- 3
ncov <- 1

nepi <- 100

##### Compare balanced and imbalanced lab sampling

npara <- 8

bala_index <- c(1, 2, 3)
imba_c_index <- c(4, 5, 7) # fixed c
imba_p_index <- c(4, 6, 8) # fixed p

adj_dist <- 0.2
mycol <- c("white", "gray80", "gray40")

name.measure <- c("Mixing Statistic", "PMSE")

load("simulation/results_summary/provided/summary_simulation_bala_imba_results.RData")

############################################
##### Figure 1: transmission rates + eta_P0
############################################
pdf(file="simulation/figures/plot_gamma+P0_south5_simulation_cloglog.pdf", width=9, height=6)
par(mfcol=c(2, 4), mar=c(2,2,1,1), oma=c(4,2,2,1))
## gamma0
# mixing statistic
plotvar <- apply(stat_res$gamma0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$gamma0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
# plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

## gamma1
# mixing statistic
plotvar <- apply(stat_res$gamma1, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 9*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$gamma1, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

## gamma2
# mixing statistic
plotvar <- apply(stat_res$gamma2, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 10*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$gamma2, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

# eta_P0
# mixing statistic
plotvar <- apply(stat_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 8.5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}


mtext(c(expression(gamma[0]), expression(gamma[1]), expression(gamma[2]), expression(cloglog(p[0]))), side=3, line=0, outer=TRUE, at=seq(from=0.12, to=0.88, length=4))
mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("balanced", "imbalanced I", "imbalanced II"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")

dev.off()

########################################
##### Figure S3 : covariate coefficents
########################################
pdf(file="simulation/figures/plot_covcoef_south5_simulation_cloglog.pdf", width=8, height=6)
par(mfcol=c(2, 3), mar=c(2,2,1,1), oma=c(4,2,2,1))
## E2P
# mixing statistic
plotvar <- apply(stat_res$eta_E, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_E, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

## P2P
# mixing statistic
plotvar <- apply(stat_res$eta_S, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 6*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$eta_S, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

## Patho
# mixing statistic
plotvar <- apply(stat_res$eta_P, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 10*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$eta_P, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

mtext(c(expression(eta[E]), expression(eta[S]), expression(eta[P])), side=3, line=0, outer=TRUE, at=seq(from=0.18, to=0.85, length=3))
mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("balanced", "imbalanced I", "imbalanced II"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")

dev.off()

############################
##### Figure S4: alpha+beta
############################
pdf(file="simulation/figures/plot_alpha+beta_south5_simulation_cloglog.pdf", width=7, height=6)
par(mfcol=c(2, 2), mar=c(2,2,1,1), oma=c(4,2,2,1))
## alpha
# mixing statistic
plotvar <- apply(stat_res$alpha_element, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(1.0*IQR(plotvar), 7*IQR(plotvar[,-bala_index])), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}
# rmse
plotvar <- apply(bmse_res$alpha_element, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

## beta
# mixing statistic
plotvar <- apply(stat_res$beta_fix, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 10*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[,bala_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_c_index[index.p]], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[,imba_p_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$beta_fix, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 3.5), ylim=c(0, 5.5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
for (index.p in 1:3)
{
	boxplot(plotvar[bala_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_c_index[index.p],], at=index.p, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
	boxplot(plotvar[imba_p_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[3], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(paste("Spatial Effects ", alpha)), expression(paste("Temporal Effects ", eta[B]))), side=3, line=0, outer=TRUE, at=seq(from=0.25, to=0.75, length=2))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("balanced", "imbalanced I", "imbalanced II"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")

dev.off()

####################
##### Figure S5-S7
####################
load("simulation/results_summary/provided/summary_means_simulation_bala_imba_results.RData")
load("simulation/parameters/provided/south5_simdata_epi_parameters.RData")

temporal_trend <- basis%*%beta_fix
temporal_trend.means <- array(0, dim=c(nweek, nv, npara, nepi))
for (i in 1:nepi)
{
	for (j in 1:npara)
	{
		temporal_trend.means[,,j,i] <- basis%*%beta_fix.means[,,j,i]
	}
}
temporal_trend.mean <- apply(temporal_trend.means, 1:3, mean)

titles_bala <- c("2% balanced", "5% balanced", "10% balanced")
titles_imba_p <- c("2% imbalanced II", "5% imbalanced II", "10% imbalanced II")
titles_imba_c <- c("2% imbalanced I", "5% imbalanced I", "10% imbalanced I")

for (v in 1:nv)
{
	pdf(paste("simulation/figures/plot_beta_south5_simulation_cloglog_v=",v,".pdf", sep=""), width=9, height=10)
	
	par(mfcol=c(3,3), mar=c(2,2,4,1), oma=c(2,2,0,0))
	
	for (index_para in 1:3)
	{
		plot(0, xlim=c(1, nweek), ylim=range(temporal_trend.means[,v,,], temporal_trend), type="n", xlab="", ylab="")
		for (i in 1:nepi)
		{
			lines(1:nweek, temporal_trend.means[,v,bala_index[index_para],i], col="grey")
		}
		lines(1:nweek, temporal_trend.mean[,v,bala_index[index_para]], col="black", lwd=1.5)
		lines(1:nweek, temporal_trend[,v], col="red", lwd=1.5)
		title(main=titles_bala[index_para])
	}
	
	for (index_para in 1:3)
	{
		plot(0, xlim=c(1, nweek), ylim=range(temporal_trend.means[,v,,], temporal_trend), type="n", xlab="", ylab="")
		for (i in 1:nepi)
		{
			lines(1:nweek, temporal_trend.means[,v,imba_c_index[index_para],i], col="grey")
		}
		lines(1:nweek, temporal_trend.mean[,v,imba_c_index[index_para]], col="black", lwd=1.5)
		lines(1:nweek, temporal_trend[,v], col="red", lwd=1.5)
		title(main=titles_imba_c[index_para])
	}
	
	for (index_para in 1:3)
	{
		plot(0, xlim=c(1, nweek), ylim=range(temporal_trend.means[,v,,], temporal_trend), type="n", xlab="", ylab="")
		for (i in 1:nepi)
		{
			lines(1:nweek, temporal_trend.means[,v,imba_p_index[index_para],i], col="grey")
		}
		lines(1:nweek, temporal_trend.mean[,v,imba_p_index[index_para]], col="black", lwd=1.5)
		lines(1:nweek, temporal_trend[,v], col="red", lwd=1.5)
		title(main=titles_imba_p[index_para])
	}
	
	
	mtext("Week", side=1, outer=TRUE, at=0.5, line=0.5)
	mtext(expression(beta), side=2, outer=TRUE, at=0.45)
	dev.off()
}

###############################################
# Figure S8: compare imputed Y and estimated Y
###############################################

load("simulation/results_summary/provided/summary_Y_mse_simulation_bala_imba_results.RData")

mycol <- c("red1", "goldenrod1")

pdf("simulation/figures/plot_Y.pdf", width=8, height=6)
par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(4,2,2,1))

for (v in 1:nv)
{	
	plot(0, type="n", ylim=range(plotvar_est[,,v], plotvar_imp[,,v]), xlim=c(0.5, 3.5), ylab="", xlab="", xaxt="n")
	axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
	for (i in 1:3) 
	{
		boxplot(plotvar_est[bala_index[i],,v], at=i-0.2, add=T, col=mycol[1], pars=list(boxwex=0.5))
		boxplot(plotvar_imp[bala_index[i],,v], at=i+0.2, add=T, col=mycol[2], pars=list(boxwex=0.5))
	}
	abline(h=0, lty=2)
	
	plot(0, type="n", ylim=range(plotvar_est[,,v], plotvar_imp[,,v]), xlim=c(0.5, 3.5), ylab="", xlab="", xaxt="n")
	axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
	for (i in 1:3) 
	{
		boxplot(plotvar_est[imba_c_index[i],,v], at=i-0.2, add=T, col=mycol[1], pars=list(boxwex=0.5))
		boxplot(plotvar_imp[imba_c_index[i],,v], at=i+0.2, add=T, col=mycol[2], pars=list(boxwex=0.5))
	}
	abline(h=0, lty=2)
	
	plot(0, type="n", ylim=range(plotvar_est[,,v], plotvar_imp[,,v]), xlim=c(0.5, 3.5), ylab="", xlab="", xaxt="n")
	axis(side=1, at=1:3, labels=c("2%", "5%", "10%"))
	for (i in 1:3) 
	{
		boxplot(plotvar_est[imba_p_index[i],,v], at=i-0.2, add=T, col=mycol[1], pars=list(boxwex=0.5))
		boxplot(plotvar_imp[imba_p_index[i],,v], at=i+0.2, add=T, col=mycol[2], pars=list(boxwex=0.5))
	}
	abline(h=0, lty=2)
}
mtext(c("balanced", "imbalanced I", "imbalanced II"), side=3, line=0, outer=TRUE, at=seq(from=0.17, to=0.84, length=3))
mtext(c("v=1", "v=2", "v=3"), side=2, line=0.7, outer=TRUE, at=seq(from=0.87, to=0.2, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("Posterior Mean", "Imputed"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")
dev.off()


##### Compare no aggregation, aggregation by neighborhood and aggregation by province
load("simulation/results_summary/provided/summary_simulation_agg_results.RData")

orig_index <- c(3, 4)
agg_index <- c(1, 2)
agg_prov_index <- c(5, 6)

adj_dist <- 0.2
orig_col <- c("white", "white")
agg_col <- c("gray80", "gray80")
agg_prov_col <- c("gray40", "gray40")
name.measure <- c("Mixing Statistic", "PMSE")

###############
##### Figure 2
###############
pdf("simulation/figures/plot_gamma_south5_simulation_agg_prov_cloglog.pdf", width=9, height=6)
par(mfcol=c(2, 4), mar=c(2,2,1,1), oma=c(4,2,2,1))
## gamma0
# mixing statistic
plotvar <- apply(stat_res$gamma0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma0, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 4.3*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

## gamma1
# mixing statistic
plotvar <- apply(stat_res$gamma1, c(2,3), mean)

plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 8*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma1, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 6*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

## gamma2
# mixing statistic
plotvar <- apply(stat_res$gamma2, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma2, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 6.2*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# eta_P0
# mixing
plotvar <- apply(stat_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 10*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P0, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 2.3*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(gamma[E]), expression(gamma[H1]), expression(gamma[H2]), expression(cloglog(p[0]))), side=3, line=0, outer=TRUE, at=seq(from=0.12, to=0.88, length=4))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("no aggregation", "agg. by neighborhood", "agg. by province"), xpd=TRUE, horiz=TRUE, fill=c(orig_col[1], agg_col[1], agg_prov_col[1]), inset=c(0,0), cex=1.5, bty="n")

dev.off()

#################
##### Figure S17
#################
pdf("simulation/figures/plot_covcoef_south5_simulation_agg_prov_cloglog.pdf", width=8, height=6)
par(mfcol=c(2, 3), mar=c(2,2,1,1), oma=c(4,2,2,1))
## E2P
# mixing statistic
plotvar <- apply(stat_res$eta_E, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 6*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_E, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 4.5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

## P2P
# mixing statistic
plotvar <- apply(stat_res$eta_S, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_S, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

## Patho
# mixing statistic
plotvar <- apply(stat_res$eta_P, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 2.2*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(eta[E]), expression(eta[H]), expression(eta[P])), side=3, line=0, outer=TRUE, at=seq(from=0.15, to=0.85, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("no aggregation", "agg. by neighborhood", "agg. by province"), xpd=TRUE, horiz=TRUE, fill=c(orig_col[1], agg_col[1], agg_prov_col[1]), inset=c(0,0), cex=1.5, bty="n")
dev.off()


##### Figure S18
pdf("simulation/figures/plot_alpha+beta_south5_simulation_agg_prov_cloglog.pdf", width=7, height=6)
par(mfcol=c(2, 2), mar=c(2,2,1,1), oma=c(4,2,2,1))
## alpha
# mixing statistic
plotvar <- apply(stat_res$alpha_element, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(IQR(plotvar), 6*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$alpha_element, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 7*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

## beta
# mixing statistic
plotvar <- apply(stat_res$beta_fix, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 8.5*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,orig_index[index.p]], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_index[index.p]], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[,agg_prov_index[index.p]], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$beta_fix, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 5.1*IQR(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
	boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
}


mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(paste("Spatial Effects ", alpha)), expression(paste("Temporal Effects ", eta[B]))), side=3, line=0, outer=TRUE, at=seq(from=0.25, to=0.75, length=2))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("no aggregation", "agg. by neighborhood", "agg. by province"), xpd=TRUE, horiz=TRUE, fill=c(orig_col[1], agg_col[1], agg_prov_col[1]), inset=c(0,0), cex=1.2, bty="n")
dev.off()

#################
##### Figure S19
#################

pdf("simulation/figures/plot_pathogenicity_pmse_south5_simulation_agg_prov_cloglog.pdf", width=8, height=6)
par(mfrow=c(2, 3), mar=c(2,2,1,1), oma=c(4,4,2,1))

## eta_P0
for (v in 1:3)
{
	plotvar <- bmse_res$eta_P0[v,,]
	plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
	axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
	for (index.p in 1:2)
	{
		boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
		boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
		boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
	}
}

## eta_P

for (v in 1:3)
{
	plotvar <- bmse_res$eta_P[,v,,]
	plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
	axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
	for (index.p in 1:2)
	{
		boxplot(plotvar[orig_index[index.p],], at=index.p - adj_dist, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
		boxplot(plotvar[agg_index[index.p],], at=index.p, add=TRUE, col=agg_col[index.p], pars=list(boxwex=0.2))
		boxplot(plotvar[agg_prov_index[index.p],], at=index.p + adj_dist, add=TRUE, col=agg_prov_col[index.p], pars=list(boxwex=0.2))
	}
}


mtext(c(expression(cloglog(p[0]^(v))), expression(eta[P]^(v))), side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(paste("v=", 1:3, sep=""), side=3, line=0, outer=TRUE, at=seq(from=0.17, to=0.85, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("no aggregation", "agg. by neighborhood", "agg. by province"), xpd=TRUE, horiz=TRUE, fill=c(orig_col[1], agg_col[1], agg_prov_col[1]), inset=c(0,0), cex=1.5, bty="n")

dev.off()


##### Compare block sampling and element-wise sampling 

load("simulation/results_summary/provided/summary_simulation_block_element_results.RData")
adj_dist <- 0.2
orig_col <- c("white", "gray50")
name.measure <- c("Mixing Statistic", "PMSE")

##################
##### Figure S20
##################
pdf("simulation/figures/plot_gamma+P0_south5_real_simulation_block.pdf", width=9, height=5)
par(mfcol=c(2, 4), mar=c(2,2,1,1), oma=c(2,2,2,1))
## gamma0
# mixing statistic
plotvar <- apply(stat_res$gamma0, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma0, c(2,3), mean)

plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}


## gamma1
# mixing statistic
plotvar <- apply(stat_res$gamma1, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma1, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## gamma2
# mixing statistic
plotvar <- apply(stat_res$gamma2, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma2, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}
# eta_P0
# mixing statistic
plotvar <- apply(stat_res$eta_P0, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(gamma[E]), expression(gamma[H1]), expression(gamma[H2]), expression(cloglog(p[0]))), side=3, line=0, outer=TRUE, at=seq(from=0.12, to=0.88, length=4))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("MBS approach", "MIS approach"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), cex=1.5, bty="n")

dev.off()

#################
##### Figure S21
#################
pdf("simulation/figures/plot_covcoef_south5_real_simulation_block.pdf", width=8, height=5)
par(mfcol=c(2, 3), mar=c(2,2,1,1), oma=c(2,2,2,1))
## E2P
# mixing statistic
plotvar <- apply(stat_res$eta_E, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_E, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## P2P
# mixing statistic
plotvar <- apply(stat_res$eta_S, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_S, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## Patho
# mixing statistic
plotvar <- apply(stat_res$eta_P, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(eta[E]), expression(eta[H]), expression(eta[P])), side=3, line=0, outer=TRUE, at=seq(from=0.15, to=0.85, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("MBS approach", "MIS approach"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), cex=1.5, bty="n")

dev.off()

#################
##### Figure S22
#################
pdf("simulation/figures/plot_alpha+beta_south5_real_simulation_block.pdf", width=7, height=5)
par(mfcol=c(2, 2), mar=c(2,2,1,1), oma=c(2,2,2,1))
## alpha
# mixing statistic
plotvar <- apply(stat_res$alpha_element, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$alpha_element, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## beta
# mixing statistic
plotvar <- apply(stat_res$beta_fix, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$beta_fix, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(paste("Spatial Effects ", alpha)), expression(paste("Temporal Effects ", eta[B]))), side=3, line=0, outer=TRUE, at=seq(from=0.25, to=0.75, length=2))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("MBS approach", "MIS approach"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), bty="n")

dev.off()


##### Compare results for one infectious week and two infectious weeks

load("simulation/results_summary/provided/summary_simulation_2infweek_results.RData")

two_week_index <- c(1,2)
one_week_index <- c(3,4)
mycol <- c("white", "gray50")
name.measure <- c("Mixing Statistic", "PMSE")
adj_dist <- 0.2

#################
##### Figure S29
#################
pdf("simulation/figures/plot_gamma_south5_simulation_cloglog_2infweek.pdf", width=9, height=6)
par(mfcol=c(2, 4), mar=c(2,2,1,1), oma=c(4,2,2,1))
## gamma0
# mixing statistic
plotvar <- apply(stat_res$gamma0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma0, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(range(plotvar)), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

## gamma1
# mixing statistic
plotvar <- apply(stat_res$gamma1, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma1, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

## gamma2
# mixing statistic
plotvar <- apply(stat_res$gamma2, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma2, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# eta_P0
# mixing
plotvar <- apply(stat_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P0, 2:3, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(gamma[E]), expression(gamma[H1]), expression(gamma[H2]), expression(cloglog(p[0]))), side=3, line=0, outer=TRUE, at=seq(from=0.12, to=0.88, length=4))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("one infectious week", "two infectious weeks"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")
dev.off()

#################
##### Figure S30
#################
pdf("simulation/figures/plot_covcoef_south5_simulation_cloglog_2infweek.pdf", width=8, height=6)
par(mfcol=c(2, 3), mar=c(2,2,1,1), oma=c(4,2,2,1))
## E2P
# mixing statistic
plotvar <- apply(stat_res$eta_E, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_E, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

## P2P
# mixing statistic
plotvar <- apply(stat_res$eta_S, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_S, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

## Patho
# mixing statistic
plotvar <- apply(stat_res$eta_P, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P, 3:4, mean)

plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(eta[E]), expression(eta[H]), expression(eta[P])), side=3, line=0, outer=TRUE, at=seq(from=0.15, to=0.85, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("one infectious week", "two infectious weeks"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")

dev.off()

#################
##### Figure S31
#################
pdf("simulation/figures/plot_alpha+beta_south5_simulation_cloglog_2infweek.pdf", width=7, height=6)
par(mfcol=c(2, 2), mar=c(2,2,1,1), oma=c(4,2,2,1))
## alpha
# mixing statistic
plotvar <- apply(stat_res$alpha_element, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$alpha_element, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

## beta
# mixing statistic
plotvar <- apply(stat_res$beta_fix, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[,one_week_index[index.p]], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[,two_week_index[index.p]], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$beta_fix, 3:4, mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
axis(side=1, at=1:2, labels=c("balanced", "imbalanced"))
for (index.p in 1:2)
{
	boxplot(plotvar[one_week_index[index.p],], at=index.p - adj_dist, add=TRUE, col=mycol[1], pars=list(boxwex=0.2))
	boxplot(plotvar[two_week_index[index.p],], at=index.p + adj_dist, add=TRUE, col=mycol[2], pars=list(boxwex=0.2))
}


mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(paste("Spatial Effects ", alpha)), expression(paste("Temporal Effects ", eta[B]))), side=3, line=0, outer=TRUE, at=seq(from=0.25, to=0.75, length=2))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("one infectious week", "two infectious weeks"), xpd=TRUE, horiz=TRUE, fill=mycol, inset=c(0,0), cex=1.5, bty="n")

dev.off()

##### Compare different sampling proportion for severe cases

load("simulation/results_summary/provided/summary_simulation_severe_results.RData")
adj_dist <- 0.2
orig_col <- c("white", "gray50")
name.measure <- c("Mixing Statistic", "PMSE")

##################
##### Figure S32
##################

pdf("simulation/figures/plot_gamma+P0_south5_severe_simulation.pdf", width=9, height=5)
par(mfcol=c(2, 4), mar=c(2,2,1,1), oma=c(2,2,2,1))
## gamma0
# mixing statistic
plotvar <- apply(stat_res$gamma0, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0.0, 0.10), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## gamma1
# mixing statistic
plotvar <- apply(stat_res$gamma1, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.85), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma1, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## gamma2
# mixing statistic
plotvar <- apply(stat_res$gamma2, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.4), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$gamma2, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}
# eta_P0
# mixing statistic
plotvar <- apply(stat_res$eta_P0, c(2, 3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P0, c(2,3), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(gamma[E]), expression(gamma[H1]), expression(gamma[H2]), expression(cloglog(p[0]))), side=3, line=0, outer=TRUE, at=seq(from=0.12, to=0.88, length=4))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("70% Severe Sampling", "20% Severe Sampling"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), cex=1.5, bty="n")

dev.off()

#################
##### Figure S33
#################
pdf("simulation/figures/plot_covcoef_south5_severe_simulation.pdf", width=8, height=5)
par(mfcol=c(2, 3), mar=c(2,2,1,1), oma=c(2,2,2,1))
## E2P
# mixing statistic
plotvar <- apply(stat_res$eta_E, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.025), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_E, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## P2P
# mixing statistic
plotvar <- apply(stat_res$eta_S, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.015), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_S, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## Patho
# mixing statistic
plotvar <- apply(stat_res$eta_P, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$eta_P, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(eta[E]), expression(eta[H]), expression(eta[P])), side=3, line=0, outer=TRUE, at=seq(from=0.15, to=0.85, length=3))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("70% Severe Sampling", "20% Severe Sampling"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), cex=1.5, bty="n")

dev.off()

####################
##### Figure S34
####################
pdf("simulation/figures/plot_alpha+beta_south5_severe_simulation.pdf", width=7, height=5)
par(mfcol=c(2, 2), mar=c(2,2,1,1), oma=c(2,2,2,1))
## alpha
# mixing statistic
plotvar <- apply(stat_res$alpha_element, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.0075), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$alpha_element, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

## beta
# mixing statistic
plotvar <- apply(stat_res$beta_fix, c(3, 4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=c(0, 0.6), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[,index.p], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

# rmse
plotvar <- apply(bmse_res$beta_fix, c(3,4), mean)
plot(0, type="n", xlim=c(0.5, 2.5), ylim=range(plotvar), xlab="", ylab="", xaxt="n")
for (index.p in 1:2)
{
	boxplot(plotvar[index.p,], at=index.p, add=TRUE, col=orig_col[index.p], pars=list(boxwex=0.2))
}

mtext(name.measure, side=2, line=0.7, outer=TRUE, at=seq(from=0.77, to=0.27, length=length(name.measure)))
mtext(c(expression(paste("Spatial Effects ", alpha)), expression(paste("Temporal Effects ", eta[B]))), side=3, line=0, outer=TRUE, at=seq(from=0.25, to=0.75, length=2))

par(fig=c(0,1,0,1), oma = c(0, 2, 2, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c("70% Severe Sampling", "20% Severe Sampling"), xpd=TRUE, horiz=TRUE, fill=orig_col, inset=c(0,0), bty="n")

dev.off()

##################################
##### Figure S9-S16: density plot
##################################
load("simulation/results_summary/provided/summary_simulation_density_results.RData")
load("simulation/parameters/provided/south5_simdata_epi_parameters.RData")

dinvgamma <- function(x, alpha, beta)
{
	x^(-alpha-1) * exp(-beta/x)
}

n_plot <- 20
title_names <- c("Balanced 2%", "Balanced 5%", "Balanced 10%", "Imbalanced I/II 2%", "Imbalanced I 5%", "Imbalanced I 10%", "Severe Balanced 2%", "Imbalanced II 5%", "Imbalanced II 10%")

# gamma0
x <- exp(seq(from=-4, to=1, by=0.01))
y <- dgamma(x, 0.1, 0.1)


pdf(paste("simulation/figures/gamma0_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(0, gamma0_all, gamma0_severe)
yrange <- range(gamma0_den_all, gamma0_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[E]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma0_severe[,index_epi], gamma0_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma0[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[E]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma0_all[,index_epi,index_para], gamma0_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma0[1], col="red", lty=2, lwd=2)
	}
}

dev.off()
	
# gamma1
x <- exp(seq(from=-6, to=1, by=0.01))
y <- dgamma(x, 0.1, 0.1)

pdf(paste("simulation/figures/gamma1_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(0, gamma1_all, gamma1_severe)
yrange <- range(gamma1_den_all, gamma1_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[H1]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma1_severe[,index_epi], gamma1_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma1[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[H1]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma1_all[,index_epi,index_para], gamma1_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma1[1], col="red", lty=2, lwd=2)
	}
}

dev.off()

# gamma2
x <- exp(seq(from=-7, to=1, by=0.01))
y <- dgamma(x, 0.1, 0.1)

pdf(paste("simulation/figures/gamma2_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(0, gamma2_all, gamma2_severe)
yrange <- range(gamma2_den_all, gamma2_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[H2]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma2_severe[,index_epi], gamma2_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma2[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(gamma[H2]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(gamma2_all[,index_epi,index_para], gamma2_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=gamma2[1], col="red", lty=2, lwd=2)
	}
}

dev.off()
		
# eta_P0
x <- seq(from=-20, to=20, by=0.01)
y <- dnorm(x, 0, sqrt(1000))
pdf(paste("simulation/figures/eta_P0_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(eta_P0_all, eta_P0_severe)
yrange <- range(eta_P0_den_all, eta_P0_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(cloglog(p[0])), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_P0_severe[,index_epi], eta_P0_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_P0[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(cloglog(p[0])), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_P0_all[,index_epi,index_para], eta_P0_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_P0[1], col="red", lty=2, lwd=2)
	}
}

dev.off()
	
	
# eta_P
pdf(paste("simulation/figures/eta_P_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(eta_P_all, eta_P_severe)
yrange <- range(eta_P_den_all, eta_P_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[P]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_P_severe[,index_epi], eta_P_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_P[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[P]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_P_all[,index_epi,index_para], eta_P_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_P[1], col="red", lty=2, lwd=2)
	}
}

dev.off()

# eta_S
pdf(paste("simulation/figures/eta_S_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(eta_S_all, eta_S_severe)
yrange <- range(eta_S_den_all, eta_S_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[S]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_S_severe[,index_epi], eta_S_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_S[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[S]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_S_all[,index_epi,index_para], eta_S_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_S[1], col="red", lty=2, lwd=2)
	}
}

dev.off()

# eta_E
pdf(paste("simulation/figures/eta_E_density_compare_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(eta_E_all, eta_E_severe)
yrange <- range(eta_E_den_all, eta_E_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta_E[E]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_E_severe[,index_epi], eta_E_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_E[1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[E]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(eta_E_all[,index_epi,index_para], eta_E_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=eta_E[1], col="red", lty=2, lwd=2)
	}
}

dev.off()
	
# beta_fix
index_k <- 1
pdf(paste("simulation/figures/beta_fix_density_compare_k=", index_k, "_v=",v,".pdf", sep=""), width=8, height=8)
par(mfcol=c(3,3))
xrange <- range(beta_fix_all, beta_fix_severe)
yrange <- range(beta_fix_den_all, beta_fix_den_severe, y)
for (index_para1 in 1:9)
{
	if (index_para1 == 7)
	{
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[B]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(beta_fix_severe[,index_epi], beta_fix_den_severe[,index_epi], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=beta_fix[index_k,1], col="red", lty=2, lwd=2)
	}
	else
	{
		if (index_para1 <= 6) index_para <- index_para1
		else index_para <- index_para1 - 1
		plot(0, xlim=xrange, ylim=yrange, type="n", xlab=expression(eta[B]), ylab="Density", main=title_names[index_para1])
		for (index_epi in 1:n_plot) lines(beta_fix_all[,index_epi,index_para], beta_fix_den_all[,index_epi,index_para], lwd=0.5)
		lines(x, y, col="blue", lty=2, lwd=2)
		abline(v=beta_fix[index_k, 1], col="red", lty=2, lwd=2)
	}
}

dev.off()
