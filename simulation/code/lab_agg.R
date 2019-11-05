
# Lab data aggregation by neighborhood or province
# It loads data generated from data_gen.R with 2% overall sampling proportion (combo1 in balanced and imbalanced design)
# It outputs txt files containing aggregated lab data, which will be used by C code

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) 
{
	stop("Please provide an index for the simulation\n", call.=FALSE)
} else {
	index.rep <- args[1]
}

load("data/2009_south5.RData")

province_code <- scan("data/south5_provincecode.txt")

## imbalanced data
load(paste("simulation/sim_data/south5_simdata_imbalanced_combo1_rep", index.rep, ".RData", sep=""))

Zv <- apply(Z, 1:3, sum)	

# aggregation by neighborhood
Z_agg <- Z
for (r in 1:nregion)
{
	for (s in 1:ns)
	{
		for (t in 1:nweek)
		{
			r_prov <- province_code[r]
			
			# time window for aggregation
			if (t<=2) t_window <- 1:(t+2)
			else if (t >= (nweek-1)) t_window <- (t-1):nweek
			else t_window <- (t-2):(t+2)
			
			# spatial range for aggregation
			r_prov_range <- province_code==r_prov
			r_nb_range <- 1:nregion %in% c(r, nb.prefecture.south5[[r]])
			
			Zv_prov_large <- Zv[r_prov_range,s,t_window]
			Zv_nb_large <- Zv[r_nb_range,s,t_window]
			Z_prov_large_sum <- apply(Z[r_prov_range,s,t_window,], 3, sum)
			Z_nb_large_sum <- apply(Z[r_nb_range,s,t_window,], 3, sum)
			
			if (Yv[r,s,t] !=0 & Zv[r,s,t] != Yv[r,s,t])
			{
				orig_prop <- Zv[r,s,t]/Yv[r,s,t] # original lab test sampling proportion 
				if (sum(Zv_nb_large)!=0) # aggregate by neighborhood if lab sample exists in neighborhood and sampling proportion is higher than the original one. 
				{
					agg_prop <- sum(Zv_nb_large)/sum(Yv[r_nb_range, s, t_window][Zv_nb_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_nb_large_sum / sum(Zv_nb_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
				else if (sum(Zv_prov_large)!=0) # aggregate by province if no lab cases in the neighborhood
				{
					agg_prop <- sum(Zv_prov_large)/sum(Yv[r_prov_range, s, t_window][Zv_prov_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_prov_large_sum / sum(Zv_prov_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
		
				if (Yv[r,s,t] - sum(Z_agg[r,s,t,]) <= 3 & sum(Z_agg[r,s,t,]) != 0) # if aggregated lab cases is no fewer than observed cases - 3, randomly assign the non-lab-tested cases to virus group. This avoids some difficulty in initialization of Y. 
				{
					Z_agg[r,s,t,] <- floor(Z_agg[r,s,t,] * Yv[r,s,t]/sum(Z_agg[r,s,t,]))
					#pick_v <- ceil(runif(1)*nv)
					pick_v <- nv
					Z_agg[r,s,t,pick_v] <- Yv[r,s,t] - sum(Z_agg[r,s,t,-pick_v])
				}

			}	
		}
	}
}

for (v in 1:nv)
{
	for (t in 1:nweek)
	{
		for (s in 1:ns)
		{
			cat(Z_agg[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_agg_simdata_imbalanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}

# aggregation by province
Z_agg <- Z
for (r in 1:nregion)
{
	for (s in 1:ns)
	{
		for (t in 1:nweek)
		{
			r_prov <- province_code[r]
		
			if (t<=2) t_window <- 1:(t+2)
			else if (t >= (nweek-1)) t_window <- (t-1):nweek
			else t_window <- (t-2):(t+2)
			
			r_prov_range <- province_code==r_prov
			
			Zv_prov_large <- Zv[r_prov_range,s,t_window]

			Z_prov_large_sum <- apply(Z[r_prov_range,s,t_window,], 3, sum)
			
			if (Yv[r,s,t] !=0 & Zv[r,s,t] != Yv[r,s,t])
			{
				orig_prop <- Zv[r,s,t]/Yv[r,s,t]
				if (sum(Zv_prov_large)!=0) # aggregate by province
				{
					agg_prop <- sum(Zv_prov_large)/sum(Yv[r_prov_range, s, t_window][Zv_prov_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_prov_large_sum / sum(Zv_prov_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
		
				if (Yv[r,s,t] - sum(Z_agg[r,s,t,]) <= 3 & sum(Z_agg[r,s,t,]) != 0) # if aggregated lab cases is no fewer than observed cases - 3, randomly assign the non-lab-tested cases to virus group. This avoids some difficulty in initialization of Y.
				{
					Z_agg[r,s,t,] <- floor(Z_agg[r,s,t,] * Yv[r,s,t]/sum(Z_agg[r,s,t,]))
					#pick_v <- ceil(runif(1)*nv)
					pick_v <- nv
					Z_agg[r,s,t,pick_v] <- Yv[r,s,t] - sum(Z_agg[r,s,t,-pick_v])
				}

			}	
		}
	}
}

for (v in 1:nv)
{
	for (t in 1:nweek)
	{
		for (s in 1:ns)
		{
			cat(Z_agg[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_agg_prov_simdata_imbalanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}
	
###### balanced data
load(paste("simulation/sim_data/south5_simdata_balanced_combo1_rep", index.rep, ".RData", sep=""))

Zv <- apply(Z, 1:3, sum)

# aggregation by neighborhood
Z_agg <- Z
for (r in 1:nregion)
{
	for (s in 1:ns)
	{
		for (t in 1:nweek)
		{
			r_prov <- province_code[r]
		
			if (t<=2) t_window <- 1:(t+2)
			else if (t >= (nweek-1)) t_window <- (t-1):nweek
			else t_window <- (t-2):(t+2)
			
			r_prov_range <- province_code==r_prov
			r_nb_range <- 1:nregion %in% c(r, nb.prefecture.south5[[r]])
			
			Zv_prov_large <- Zv[r_prov_range,s,t_window]
			Zv_nb_large <- Zv[r_nb_range,s,t_window]
			Z_prov_large_sum <- apply(Z[r_prov_range,s,t_window,], 3, sum)
			Z_nb_large_sum <- apply(Z[r_nb_range,s,t_window,], 3, sum)
			
			if (Yv[r,s,t] !=0 & Zv[r,s,t] != Yv[r,s,t]) 
			{
				orig_prop <- Zv[r,s,t]/Yv[r,s,t]
				if (sum(Zv_nb_large)!=0) # aggregate by neighborhood if lab sample exists in neighborhood and sampling proportion is higher than the original one. 
				{
					agg_prop <- sum(Zv_nb_large)/sum(Yv[r_nb_range, s, t_window][Zv_nb_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_nb_large_sum / sum(Zv_nb_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
				else if (sum(Zv_prov_large)!=0) # aggregate by province if no lab cases in the neighborhood
				{
					agg_prop <- sum(Zv_prov_large)/sum(Yv[r_prov_range, s, t_window][Zv_prov_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_prov_large_sum / sum(Zv_prov_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
		
				if (Yv[r,s,t] - sum(Z_agg[r,s,t,]) <= 3 & sum(Z_agg[r,s,t,]) != 0) # if aggregated lab cases is no fewer than observed cases - 3, randomly assign the non-lab-tested cases to virus group. This avoids some difficulty in initialization of Y.
				{
					Z_agg[r,s,t,] <- floor(Z_agg[r,s,t,] * Yv[r,s,t]/sum(Z_agg[r,s,t,]))
					#pick_v <- ceil(runif(1)*nv)
					pick_v <- nv
					Z_agg[r,s,t,pick_v] <- Yv[r,s,t] - sum(Z_agg[r,s,t,-pick_v])
				}

			}	
		}
	}
}

for (v in 1:nv)
{
	for (t in 1:nweek)
	{
		for (s in 1:ns)
		{
			cat(Z_agg[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_agg_simdata_balanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}

# aggregation by province
Z_agg <- Z
for (r in 1:nregion)
{
	for (s in 1:ns)
	{
		for (t in 1:nweek)
		{
			r_prov <- province_code[r]
		
			if (t<=2) t_window <- 1:(t+2)
			else if (t >= (nweek-1)) t_window <- (t-1):nweek
			else t_window <- (t-2):(t+2)
			
			r_prov_range <- province_code==r_prov
			
			Zv_prov_large <- Zv[r_prov_range,s,t_window]
			Z_prov_large_sum <- apply(Z[r_prov_range,s,t_window,], 3, sum)
			
			if (Yv[r,s,t] !=0 & Zv[r,s,t] != Yv[r,s,t])
			{
				orig_prop <- Zv[r,s,t]/Yv[r,s,t]
				if (sum(Zv_prov_large)!=0) # aggregate by province
				{
					agg_prop <- sum(Zv_prov_large)/sum(Yv[r_prov_range, s, t_window][Zv_prov_large!=0])
					if (agg_prop > orig_prop)
					{
						Z_impute <- agg_prop * Yv[r,s,t] * Z_prov_large_sum / sum(Zv_prov_large)
						Z_agg[r,s,t,] <- Z_impute
					}
				}
		
				if (Yv[r,s,t] - sum(Z_agg[r,s,t,]) <= 3 & sum(Z_agg[r,s,t,]) != 0) # if aggregated lab cases is no fewer than observed cases - 3, randomly assign the non-lab-tested cases to virus group. This avoids some difficulty in initialization of Y.
				{
					Z_agg[r,s,t,] <- floor(Z_agg[r,s,t,] * Yv[r,s,t]/sum(Z_agg[r,s,t,]))
					#pick_v <- ceil(runif(1)*nv)
					pick_v <- nv
					Z_agg[r,s,t,pick_v] <- Yv[r,s,t] - sum(Z_agg[r,s,t,-pick_v])
				}

			}	
		}
	}
}

for (v in 1:nv)
{
	for (t in 1:nweek)
	{
		for (s in 1:ns)
		{
			cat(Z_agg[,s,t,v], "\n", file=paste("simulation/sim_data/south5_Z_agg_prov_simdata_balanced_combo1_rep", index.rep, ".txt", sep=""), append=TRUE)
		}
	}
}

