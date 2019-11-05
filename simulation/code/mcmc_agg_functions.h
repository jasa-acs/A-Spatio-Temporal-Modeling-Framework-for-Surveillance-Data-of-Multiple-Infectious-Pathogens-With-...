
double maxY_fun (int r, int s, int t, int v, POLYGON *poly)
{
    double part1=0;
    part1 = poly[r].Y[s][t][v] + poly[r].Y[s][t][n_v-1] - poly[r].Z_agg[s][t][n_v-1];
    
    return part1;
}

double lambda_fun (POLYGON *poly, int r, int t, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S) // no index s in lambda_fun
{
    double part1, part2, part3,cov_part1,cov_part2,cov_part3;
    int index_nb, j;
    
    cov_part1 = 0.0;
    cov_part2 = 0.0;
    for (j=0; j < n_cov; j++)
    {
        cov_part1 = cov_part1 + poly[r].X[t][j] * eta_E[j][v];
        cov_part2 = cov_part2 + poly[r].X[t][j] * eta_S[j][v];
    }
    
    part1 = gamma0[v] * exp(cov_part1);
    
    part2 = gamma1[v] * exp(poly[r].alpha[r][v] + poly[r].beta[t][v] + cov_part2) * poly[r].Ys[t-1][v];
    part3 = 0.0;
    for (index_nb = 0; index_nb < poly[r].n_nb; index_nb++)
    {
        cov_part3 = 0.0;
        for (j=0; j < n_cov;j++)
        {
            cov_part3 = cov_part3 + poly[r].X[t][j] * eta_S[j][v];
        }
        part3 = part3 + gamma2[v]*exp(poly[r].alpha[poly[r].nb[index_nb]][v] + poly[r].beta[t][v] + cov_part3) * poly[poly[r].nb[index_nb]].Ys[t-1][v];
    }
    
    return part1 + part2 + part3;

}

double eta_E_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S, double **eta_E_mean, double **eta_E_var)
{
    int index_r, index_t, index_j;
    double part1=0.0, lambda=0.0, part2=0.0;
    
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=1; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly, index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S);
            part1 = part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    for (index_j = 0; index_j < n_cov; index_j++)
    {
        part2 = part2 - pow(eta_E[index_j][v]-eta_E_mean[index_j][v],2)/eta_E_var[index_j][v]/2.0;
    }
    return part1 + part2;
}

double logpmfY_fun (int y, POLYGON *poly, int r, int s, int t, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S)
{
    double part1, part2, part3, part4;
    double current_Lambda, *future_Lambdas, future_Lambda;
    int index_nb, old_Ys;
    
    old_Ys = poly[r].Ys[t][v];
    poly[r].Ys[t][v] = old_Ys - poly[r].Y[s][t][v] + y;
    
    part1 = -1.0 * lgamma(y - poly[r].Z_agg[s][t][v] + 1.0);
    
    current_Lambda = lambda_fun(poly, r, t, v, gamma0, gamma1, gamma2, eta_E, eta_S);
    future_Lambda = lambda_fun(poly, r,t+1,v,gamma0,gamma1,gamma2, eta_E, eta_S);
    make_1d_array_double(&future_Lambdas, poly[r].n_nb, 0.0);
    for (index_nb=0; index_nb < poly[r].n_nb; index_nb++)
    {
        future_Lambdas[index_nb] = lambda_fun(poly, poly[r].nb[index_nb], t+1, v, gamma0, gamma1, gamma2, eta_E, eta_S);
    }
    
    part2 = y * log(current_Lambda*poly[r].ps[s][t][v]);
    
    part3 = poly[r].Ys[t+1][v] * log(future_Lambda) - future_Lambda;
    for (index_nb=0; index_nb < poly[r].n_nb; index_nb++)
    {
        part3 = part3 + poly[poly[r].nb[index_nb]].Ys[t+1][v] * log(future_Lambdas[index_nb]) - future_Lambdas[index_nb];
    }
    
    free(future_Lambdas);
    poly[r].Ys[t][v] = old_Ys;
    
    return part1 + part2 + part3;
}

double gamma_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape0, double *shape1, double *shape2, double *rate0, double *rate1, double *rate2, double **eta_E, double **eta_S)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=1; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape0[v]-1) * log(gamma0[v]) - rate0[v]*gamma0[v] + (shape1[v]-1) * log(gamma1[v]) - rate1[v]*gamma1[v] + (shape2[v]-1) * log(gamma2[v]) - rate2[v]*gamma2[v];
    
    return(part1 + part2);
}


double gamma0_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape0, double * rate0, double **eta_E, double **eta_S)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=1; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape0[v]-1.0) * log(gamma0[v]) - rate0[v]*gamma0[v];
    
    return part1 + part2;
}

double gamma1_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape1, double *rate1, double **eta_E, double **eta_S)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=1; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape1[v]-1.0) * log(gamma1[v]) - rate1[v]*gamma1[v];
    
    return part1 + part2;
}

double gamma2_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape2, double *rate2, double **eta_E, double **eta_S)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=1; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape2[v]-1.0) * log(gamma2[v]) - rate2[v]*gamma2[v];
    
    return part1 + part2;
}

int check_nb(int region1, int region2, POLYGON *poly) // check whether region1 and region2 are neighbors
{
    int indicator = 0, index_nb;
    for (index_nb = 0; index_nb < poly[region1].n_nb; index_nb++)
    {
        if (region2 == poly[region1].nb[index_nb])
        {
            indicator = 1;
            break;
        }
    }
    
    return indicator;
}

double alpha_logdensity_fun(POLYGON *poly, int r, int v, double *gamma0, double *gamma1, double *gamma2, double *sigma2_alpha, double **eta_E, double **eta_S)
{
    int index_t,index_nb;
    double part1, part2, lambda_r, lambda_R;
    part1 = 0.0;
    for (index_t=1; index_t < n_week; index_t++)
    {
        lambda_r = lambda_fun(poly,r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S);
        lambda_R = lambda_fun(poly,n_region-1,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S);
        part1 = part1 - lambda_r + log(lambda_r) * poly[r].Ys[index_t][v];
        part1 = part1 - lambda_R + log(lambda_R) * poly[n_region-1].Ys[index_t][v];
    }
    
    part2 = 0.0;
    for (index_nb=0; index_nb < poly[r].n_nb; index_nb++)
    {
        part2 = part2 - pow(poly[r].alpha_element[v]-poly[poly[r].nb[index_nb]].alpha_element[v],2);
    }
    for (index_nb=0; index_nb < poly[n_region-1].n_nb; index_nb++)
    {
        part2 = part2 - pow(poly[n_region-1].alpha_element[v]-poly[poly[n_region-1].nb[index_nb]].alpha_element[v],2);
    }
    
    if (check_nb(r, n_region-1, poly)) part2 = part2 + pow(poly[n_region-1].alpha_element[v] - poly[r].alpha_element[v],2);
    part2 = part2/(sigma2_alpha[v]*2.0);
    
    return part1+part2;
    
}

double alphas_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *sigma2_alpha, double **eta_E, double **eta_S)
{
    int index_r,index_t,index_nb;
    double part1, part2, lambda_r;
    part1 = 0.0;
    for (index_t=1; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda_r = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S);
            part1 = part1 - lambda_r + log(lambda_r) * poly[index_r].Ys[index_t][v];
        }
    }
    
    part2 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_nb=0; index_nb < poly[index_r].n_nb; index_nb++)
        {
            part2 = part2 - pow(poly[index_r].alpha_element[v]-poly[poly[index_r].nb[index_nb]].alpha_element[v],2);
        }
    }
    part2 = part2/(sigma2_alpha[v]*4.0);
    
    return part1+part2;
    
}

double beta_fix_logdensity_fun(POLYGON *poly, int k, int v, double *gamma0, double *gamma1, double *gamma2, double **beta_fix, double **beta_fix_mean, double **beta_fix_var, double **eta_E, double **eta_S)
{
    int index_t, index_r, index_nb;
    double part1, part2, lambda;
    
    part1 = 0.0;
    for (index_t=1; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S);
            part1 = part1 - lambda + log(lambda) * poly[index_r].Ys[index_t][v];
        }
    }
    
    part2 = log(dnorm(beta_fix[k][v],beta_fix_mean[k][v],sqrt(beta_fix_var[k][v])));

    return part1+part2;
}

double beta_fixs_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double **beta_fix, double **beta_fix_mean, double **beta_fix_var, double **eta_E, double **eta_S)
{
    int index_t, index_r, index_nb, index_k;
    double part1, part2, lambda;
    
    part1 = 0.0;
    for (index_t=1; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S);
            part1 = part1 - lambda + log(lambda) * poly[index_r].Ys[index_t][v];
        }
    }
    
    part2 = 0.0;
    for (index_k=0; index_k < n_knots+3; index_k++)
    {
        part2 = part2 + log(dnorm(beta_fix[index_k][v],beta_fix_mean[index_k][v],sqrt(beta_fix_var[index_k][v])));
    }
    
    return part1+part2;
}

void rgmrf (double *y, int n, int n_lc, double* mean, double *eval, double **evec, long *seed) // columns of evec are the eigenvectors of the neighborhood matrix
{
    int i,j;
    double *x;
    
    make_1d_array_double(&x, n-n_lc, 0.0);
    
    for (j=0; j<n-n_lc; j++)
    {
        x[j] = rnorm(0.0, sqrt(1.0/eval[j]), seed);
    }
    
    for (i=0; i<n; i++)
    {
        y[i] = 0.0;
        for (j=0; j<n-n_lc; j++)
        {
            y[i] = y[i] + x[j]*evec[i][j];
        }
        y[i] = y[i] + mean[i];
    }
    free(x);
}




