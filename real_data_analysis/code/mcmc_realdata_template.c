#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define logit(x) ( log( (x) / ( 1.0 - (x) ) ) )
#define inv_logit(x) ( 1.0 / ( 1.0 + exp(-(x)) ) )

#define n_week 53
#define n_region 69
#define n_prov 5
#define n_s 2
#define n_v 3
#define n_cov 12
#define n_comp 1
#define n_sim 50
#define n_burn 50
#define n_knots 3 // number of inner knots
#define thin 2
#define n_batch 10
#define max_batch 5

#include "matrix.h"
#include "mathfunc.h"
#include "distribution.h"
#include "datastructure_agg.h" /*specifying data structure*/

double maxY_fun (int r, int s, int t, int v, POLYGON *poly)
{
    double part1=0.0;
    part1 = poly[r].Y[s][t][v] + poly[r].Y[s][t][n_v-1] - poly[r].Z_agg[s][t][n_v-1];
    
    return part1;
}

double lambda_fun (POLYGON *poly, int r, int t, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S, double c) // no index s in lambda_fun
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
    
    part2 = gamma1[v] * exp(poly[r].alpha[r][v] + poly[r].beta[t][v] + cov_part2) * (poly[r].Ys[t-1][v] + c * poly[r].Ys[t-2][v]);
    
    part3 = 0.0;
    for (index_nb = 0; index_nb < poly[r].n_nb; index_nb++)
    {
        part3 = part3 + gamma2[v] * exp(poly[r].alpha[poly[r].nb[index_nb]][v] + poly[r].beta[t][v] + cov_part2) * (poly[poly[r].nb[index_nb]].Ys[t-1][v] + c * poly[poly[r].nb[index_nb]].Ys[t-2][v]);
    }
    
    return part1 + part2 + part3;

}

double eta_E_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S, double **eta_E_mean, double **eta_E_var, double c)
{
    int index_r, index_t, index_j;
    double part1=0.0, lambda=0.0, part2=0.0;
    
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=2; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly, index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
            part1 = part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    for (index_j = 0; index_j < n_cov; index_j++)
    {
        part2 = part2 - pow(eta_E[index_j][v]-eta_E_mean[index_j][v],2)/eta_E_var[index_j][v]/2.0;
    }
    return part1 + part2;
}

double logpmfY_fun (int y, POLYGON *poly, int r, int s, int t, int v, double *gamma0, double *gamma1, double *gamma2, double **eta_E, double **eta_S, double c)
{
    double part1, part2, part3, part4;
    double current_Lambda, *future_Lambdas, *future2_Lambdas, future_Lambda, future2_Lambda;
    int index_nb, old_Ys;
    
    old_Ys = poly[r].Ys[t][v];
    poly[r].Ys[t][v] = old_Ys - poly[r].Y[s][t][v] + y;
    
    part1 = -1.0 * lgamma(y - poly[r].Z_agg[s][t][v] + 1.0);
    
    current_Lambda = lambda_fun(poly, r, t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
    future_Lambda = lambda_fun(poly, r,t+1,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
    future2_Lambda = lambda_fun(poly, r,t+2,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
    
    make_1d_array_double(&future_Lambdas, poly[r].n_nb, 0.0);
    make_1d_array_double(&future2_Lambdas, poly[r].n_nb, 0.0);
    
    for (index_nb=0; index_nb < poly[r].n_nb; index_nb++)
    {
        future_Lambdas[index_nb] = lambda_fun(poly, poly[r].nb[index_nb], t+1, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
        future2_Lambdas[index_nb] = lambda_fun(poly, poly[r].nb[index_nb], t+2, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
    }
    
    part2 = y * log(current_Lambda*poly[r].ps[s][t][v]);
    
    part3 = poly[r].Ys[t+1][v] * log(future_Lambda) - future_Lambda + poly[r].Ys[t+2][v] * log(future2_Lambda) - future2_Lambda;
    for (index_nb=0; index_nb < poly[r].n_nb; index_nb++)
    {
        part3 = part3 + poly[poly[r].nb[index_nb]].Ys[t+1][v] * log(future_Lambdas[index_nb]) - future_Lambdas[index_nb] + poly[poly[r].nb[index_nb]].Ys[t+2][v] * log(future2_Lambdas[index_nb]) - future2_Lambdas[index_nb];
    }
    
    free(future_Lambdas);
    free(future2_Lambdas);
    poly[r].Ys[t][v] = old_Ys;
    
    return part1 + part2 + part3;
}

double gamma_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape0, double *shape1, double *shape2, double *rate0, double *rate1, double *rate2, double **eta_E, double **eta_S, double c)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=2; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S, c);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape0[v]-1) * log(gamma0[v]) - rate0[v]*gamma0[v] + (shape1[v]-1) * log(gamma1[v]) - rate1[v]*gamma1[v] + (shape2[v]-1) * log(gamma2[v]) - rate2[v]*gamma2[v];
    
    return(part1 + part2);
}


double gamma0_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape0, double * rate0, double **eta_E, double **eta_S, double c)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=2; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S, c);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape0[v]-1.0) * log(gamma0[v]) - rate0[v]*gamma0[v];
    
    return part1 + part2;
}

double gamma1_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape1, double *rate1, double **eta_E, double **eta_S, double c)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=2; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
            part1= part1 - lambda + poly[index_r].Ys[index_t][v] * log(lambda);
        }
    }
    
    part2 = (shape1[v]-1.0) * log(gamma1[v]) - rate1[v]*gamma1[v];
    
    return part1 + part2;
}

double gamma2_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *shape2, double *rate2, double **eta_E, double **eta_S, double c)
{
    double part1, part2, lambda;
    int index_r, index_t;
    
    part1 = 0.0;
    for (index_r=0; index_r < n_region; index_r++)
    {
        for (index_t=2; index_t < n_week; index_t++)
        {
            lambda = lambda_fun(poly,index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
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

double alpha_logdensity_fun(POLYGON *poly, int r, int v, double *gamma0, double *gamma1, double *gamma2, double *sigma2_alpha, double **eta_E, double **eta_S, double c)
{
    int index_t,index_nb;
    double part1, part2, lambda_r, lambda_R;
    part1 = 0.0;
    for (index_t=2; index_t < n_week; index_t++)
    {
        lambda_r = lambda_fun(poly,r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
        lambda_R = lambda_fun(poly,n_region-1,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
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

double alphas_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double *sigma2_alpha, double **eta_E, double **eta_S, double c)
{
    int index_r,index_t,index_nb;
    double part1, part2, lambda_r;
    part1 = 0.0;
    for (index_t=2; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda_r = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
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

double beta_fix_logdensity_fun(POLYGON *poly, int k, int v, double *gamma0, double *gamma1, double *gamma2, double **beta_fix, double **beta_fix_mean, double **beta_fix_var, double **eta_E, double **eta_S, double c)
{
    int index_t, index_r, index_nb;
    double part1, part2, lambda;
    
    part1 = 0.0;
    for (index_t=2; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
            part1 = part1 - lambda + log(lambda) * poly[index_r].Ys[index_t][v];
        }
    }
    
    part2 = log(dnorm(beta_fix[k][v],beta_fix_mean[k][v],sqrt(beta_fix_var[k][v])));

    return part1+part2;
}

double beta_fixs_logdensity_fun(POLYGON *poly, int v, double *gamma0, double *gamma1, double *gamma2, double **beta_fix, double **beta_fix_mean, double **beta_fix_var, double **eta_E, double **eta_S, double c)
{
    int index_t, index_r, index_nb, index_k;
    double part1, part2, lambda;
    
    part1 = 0.0;
    for (index_t=2; index_t < n_week; index_t++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            lambda = lambda_fun(poly,index_r,index_t,v,gamma0,gamma1,gamma2, eta_E, eta_S, c);
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

int main (int argc, char* argv[])
{
    if (argc == 1) {
        printf("\n Please provide index of mcmc chains as an extra argument\n!");
        return 0;
    }
    
    int index_chain;
    index_chain = atoi(argv[1]);
    printf("Chain %d\n", index_chain);
    
    
	FILE *file;
    char filename[100];
    
	int r,s,t,v,i,k,j, index_t, index_k, index_r, index_j, index_nb, index_v;
    
    double c = 0.0;
    
    POLYGON *polygons;
    double *gamma0, *gamma1, *gamma2, *sigma2_alpha;
    int knots[n_knots]={20, 30, 40};
    double **basis;
    long seed, all_seeds[10] = {777429, 351049, 451175, 147253, 372297, 980743, 255947, 572525, 687361, 842243};
    
    double *nb_eval_prov, **nb_evec_prov;
    
    double **eta_S, **eta_E, *eta_P0;
    
    // variables for Y
    double sum_Z=0.0, pmf_Y[n_v]={0.0, 0.0, 0.0}, sum_pmf_Y=0.0;
    int sum_Z_int, new_Y_block[n_v]={0, 0, 0}, old_Y_block[n_v]={0, 0, 0};

    // variables for gamma
    double shape0[n_v]={1.5, 1.5, 1.5},shape1[n_v]={1.5, 1.5, 1.5},shape2[n_v]={1.5, 1.5, 1.5},rate0[n_v]={0.5, 0.5, 0.5},rate1[n_v]={0.5, 0.5, 0.5}, rate2[n_v]={0.5, 0.5, 0.5};// initialize shape and rate parameters of gamma priors
    double logacprate;
    double new_gamma0, new_gamma1, new_gamma2, old_gamma0, old_gamma1, old_gamma2, gamma_lb=1e-6;
    
    // variables for sigma2_alpha
    double shape_sigma2_alpha, rate_sigma2_alpha, prior_shape_alpha[n_v]={2.1, 2.1, 2.1}, prior_rate_alpha[n_v]={1.0, 1.0, 1.0};
    
    // variables for alpha_element
    double nu=0.0, **alpha_element_prov;
    
    // variables for beta
    double old_beta[n_week], diff_beta[n_week]={0.0};
    double **beta_fix, **beta_fix_mean, **beta_fix_var;
    MATRIX design_basis, inv_design_basis, inv_design_beta[n_v];
    double old_beta_fixs[n_knots+3] = {0.0}, new_beta_fixs[n_knots+3]={0.0}, diff_beta_fixs[n_knots+3]={0.0};
    
    // variables for block
    double *new_values_prov, *old_values_prov, *tuned_eval_alpha_prov;
    
    // variables for pathogenicity
    MATRIX design_1, inv_design_1, inv_design_E[n_v], inv_design_S[n_v];
    
    // variables for covariates
    MATRIX design_0, inv_design_0;
    double old_eta_S[n_cov]={0.0}, new_eta_S[n_cov]={0.0}, lambda_new=0.0, lambda_old=0.0;
    double old_eta_E[n_cov]={0.0}, new_eta_E[n_cov]={0.0};
    double **eta_E_mean, **eta_E_var;
    
    double para1_beta = 0.0, para2_beta = 0.0, new_ps0 = 0.0;
    
    // variables for adaption
    int adapt_status = 0, count_batch = 0;
    int n_accept_gammas[n_v] = {0}, n_accept_alpha[n_region][n_v] = {0}, n_accept_beta_fix[n_knots+3][n_v] = {0}, n_accept_eta_S[n_v] = {0}, n_accept_eta_E[n_v] ={0}, n_accept_gamma0[n_v] = {0}, n_accept_gamma1[n_v] = {0}, n_accept_gamma2[n_v] = {0}, n_accept_beta_fixs[n_v] = {0}, n_accept_alphas[n_v] = {0};
    double tuning_gammas[n_v] = {0.0}, tuning_alpha[n_region][n_v] = {0.0}, tuning_beta_fix[n_knots+3][n_v] = {0.0}, tuning_eta_S[n_v] = {0.0}, tuning_eta_E[n_v] = {0.0}, tuning_gamma0[n_v] = {0.0}, tuning_gamma1[n_v] = {0.0}, tuning_gamma2[n_v] = {0.0}, tuning_beta_fixs[n_v]={0.0}, tuning_alphas[n_v] = {0};
    double single_lower = 0.40, single_upper = 0.50, group_lower = 0.20, group_upper = 0.30, accept_rate = 0.0, tuning_adjust = 0.1;
    
    time_t t1, t2;
    
    seed = all_seeds[index_chain-1];
    t1=time(NULL);
    
    mt_init(seed);
    
    make_1d_array_double(&gamma0, n_v, 1.0);
    make_1d_array_double(&gamma1, n_v, 0.1);
    make_1d_array_double(&gamma2, n_v, 0.01);
    make_2d_array_double(&alpha_element_prov, n_prov, n_v, 0.0);
    make_1d_array_double(&sigma2_alpha, n_v, 0.01);
    make_2d_array_double(&basis, n_knots+3, n_week, 0.0);
    make_2d_array_double(&beta_fix, n_knots+3, n_v, 0.0);
    make_2d_array_double(&beta_fix_mean, n_knots+3, n_v, 0.0);
    make_2d_array_double(&beta_fix_var, n_knots+3, n_v, 1000.0);
    
    make_2d_array_double(&eta_E, n_cov, n_v, 0.0);
    make_2d_array_double(&eta_S, n_cov, n_v, 0.0);
    make_1d_array_double(&eta_P0, n_v, 0.0);
    make_2d_array_double(&eta_E_mean, n_cov, n_v, 0.0);
    make_2d_array_double(&eta_E_var, n_cov, n_v, 1000.0);
    
    make_1d_array_double(&old_values_prov, n_prov, 0.0);
    make_1d_array_double(&new_values_prov, n_prov, 0.0);
    make_1d_array_double(&nb_eval_prov, n_prov, 0.0);
    make_1d_array_double(&tuned_eval_alpha_prov, n_prov, 0.0);
    make_2d_array_double(&nb_evec_prov, n_prov, n_prov, 0.0);
    
    polygons = (POLYGON *)malloc((size_t) n_region * sizeof(POLYGON)); // Create one polygons for each region
    
    /* Read basis function */
    file = fopen("data/basis_modified.txt", "r");
    if (file == NULL)
    {
        printf("Can't open the basis function file!\n");
    }
    else
    {
        rewind(file);
        for (k=0; k < n_knots+3; k++)
        {
            for (t=0; t < n_week; t++)
            {
                fscanf(file, "%lf", &basis[k][t]);
            }
        }
    }
    fclose(file);
    
    /* Read region information and initialize the polygons*/
    file = fopen("data/south5_region2prefecture.txt", "r");
    if (file==NULL)
    {
        printf("Can't open the prefecture code file!\n");
    }
    else
    {   
        rewind(file);
        for (r=0; r < n_region;r++)
        {
            fscanf(file, "%d", &polygons[r].id);
            fscanf(file, "%d", &polygons[r].code);
            polygons[r].id--;
            
            make_2d_array_int(&polygons[r].Yv, n_s, n_week, 0);
            make_3d_array_int(&polygons[r].Z, n_s, n_week, n_v, 0);
            make_3d_array_double(&polygons[r].Z_agg, n_s, n_week, n_v, 0.0);
            make_3d_array_int(&polygons[r].Y, n_s, n_week, n_v, 0);
            make_2d_array_int(&polygons[r].Ys, n_week, n_v, 0);
            
            make_1d_array_double(&polygons[r].alpha_element, n_v, 0.0);
            make_2d_array_double(&polygons[r].alpha, n_region, n_v, 0.0);
            make_2d_array_double(&polygons[r].beta, n_week, n_v, 0.0);
            
            make_2d_array_double(&polygons[r].X, n_week, n_cov, 0.0);
            make_3d_array_double(&polygons[r].ps, n_s, n_week, n_v, 0.5);
            
        }
    }
    fclose(file);
    
    file = fopen("data/south5_provincecode.txt", "r");
    if (file==NULL) printf("Can't open province code file!\n");
    else
    {
        rewind(file);
        for (r = 0; r < n_region; r++)
        {
            fscanf(file, "%d", &polygons[r].province);
            polygons[r].province--;
        }
    }
    
    /* Read neighborhood information */
    file = fopen("data/south5_nb.txt", "r");
    if (file==NULL) {
        printf("Can't open the neighbor info file!\n");
    }
    else
    {   rewind(file);
        for (r=0; r < n_region; r++)
        {
            fscanf(file, "%d", &polygons[r].n_nb);
            make_1d_array_int(&polygons[r].nb, polygons[r].n_nb,0);
            for (j = 0; j < polygons[r].n_nb; j++)
            {
                fscanf(file, "%d", &(polygons[r].nb[j]));
                polygons[r].nb[j]--;
            }
        }
    }
    fclose(file);
    
    /* Read eigenvalues and eigenvectors of the neighborhood matrix */
    
    file = fopen("data/south5_nb_eval_prov.txt", "r");
    if (file == NULL) printf("Can't open the nb_eval_prov file!\n");
    else
    {
        for (r=0; r < n_prov; r++)
        {
            fscanf(file, "%lf", &nb_eval_prov[r]);
        }
    }
    fclose(file);
    
    file = fopen("data/south5_nb_evec_prov.txt", "r");
    if (file == NULL) printf("Can't open the nb_evec_prov file!\n");
    else
    {
        for (r=0; r < n_prov; r++)
        {
            for (j=0; j < n_prov; j++)
            {
                fscanf(file, "%lf", &nb_evec_prov[j][r]);
            }
        }
    }
    fclose(file);
    
    /* Read Yv */
    file = fopen("data/2009_south5_Yv_realdata.txt", "r");
    if (file == NULL)
    {
        printf("Can't open Yv file!\n");
    }
    else
    {
        for (t=0; t < n_week; t++)
        {
            for (s=0; s < n_s; s++)
            {
                for (r=0; r < n_region; r++)
                {
                    fscanf(file, "%d", &(polygons[r].Yv[s][t]));
                }
            }
        }
    }
    fclose(file);
    
    /* Read Z */
    file = fopen("data/2009_south5_Z_realdata.txt", "r");
    if (file == NULL)
    {
        printf("Can't open Z file!\n");
    }
    else
    {
        for (v=0; v < n_v; v++)
        {
            for (t=0; t < n_week; t++)
            {
                for (s=0; s < n_s; s++)
                {
                    for (r=0; r < n_region; r++)
                    {
                        fscanf(file, "%d", &(polygons[r].Z[s][t][v]));
                    }
                }
            }
        }
    }
    fclose(file);
    
    file = fopen("data/2009_south5_Z_agg_realdata.txt", "r");
    if (file == NULL)
    {
        printf("Can't open Z file!\n");
    }
    else
    {
        for (v=0; v < n_v; v++)
        {
            for (t=0; t < n_week; t++)
            {
                for (s=0; s < n_s; s++)
                {
                    for (r=0; r < n_region; r++)
                    {
                        fscanf(file, "%lf", &(polygons[r].Z_agg[s][t][v]));
                    }
                }
            }
        }
    }
    fclose(file);
    
    /* Read Covariates */
    file = fopen("data/covariates/2009_south5_ws_std.txt", "r");
    if (file == NULL) printf("Can't open wind speed file!\n");
    else
    {
        for (t=0; t < n_week; t++)
        {
            for (r=0; r < n_region; r++)
            {
                fscanf(file, "%lf", &(polygons[r].X[t][0]));
                polygons[r].X[t][4] = polygons[r].X[t][0]*polygons[r].X[t][0];
                polygons[r].X[t][8] = polygons[r].X[t][4]*polygons[r].X[t][0];
            }
        }
    }
    fclose(file);

    file = fopen("data/covariates/2009_south5_temp_std.txt", "r");
    if (file == NULL) printf("Can't open temperature file!\n");
    else
    {
        for (t=0; t < n_week; t++)
        {
            for (r=0; r < n_region; r++)
            {
                fscanf(file, "%lf", &(polygons[r].X[t][1]));
                polygons[r].X[t][5] = polygons[r].X[t][1]*polygons[r].X[t][1];
                polygons[r].X[t][9] = polygons[r].X[t][5]*polygons[r].X[t][1];
            }
        }
    }
    fclose(file);
    
    file = fopen("data/covariates/2009_south5_rh_std.txt", "r");
    if (file == NULL) printf("Can't open relative humidity file!\n");
    else
    {
        for (t=0; t < n_week; t++)
        {
            for (r=0; r < n_region; r++)
            {
                fscanf(file, "%lf", &(polygons[r].X[t][2]));
                polygons[r].X[t][6] = polygons[r].X[t][2]*polygons[r].X[t][2];
                polygons[r].X[t][10] = polygons[r].X[t][6]*polygons[r].X[t][2];
            }
        }
    }
    fclose(file);
    
    file = fopen("data/covariates/2009_south5_log_pop_density_std.txt", "r");
    if (file == NULL) printf("Can't open the population file!\n");
    else
    {
        for (r=0; r < n_region; r++)
        {
            fscanf(file, "%lf", &(polygons[r].X[0][3]));
            polygons[r].X[0][7] = polygons[r].X[0][3]*polygons[r].X[0][3];
            polygons[r].X[0][11] = polygons[r].X[0][7]*polygons[r].X[0][3];
            for (t=2; t < n_week; t++)
            {
                polygons[r].X[t][3] = polygons[r].X[0][3];
                polygons[r].X[t][7] = polygons[r].X[0][7];
                polygons[r].X[t][11] = polygons[r].X[0][11];
            }
        }
    }
    fclose(file);
    
    /* Calculate design matrix */
    initialize_matrix(&design_1);
    initialize_matrix(&design_0);
    initialize_matrix(&inv_design_1);
    initialize_matrix(&inv_design_0);
    for (v=0; v < n_v; v++)
    {
        initialize_matrix(&(inv_design_E[v]));
        initialize_matrix(&(inv_design_S[v]));
    }
    
    inflate_matrix(&design_1, n_cov+1, n_cov+1, 0.0);
    inflate_matrix(&design_0, n_cov, n_cov, 0.0);
    inflate_matrix(&inv_design_1, n_cov+1, n_cov+1, 0.0);
    inflate_matrix(&inv_design_0, n_cov, n_cov, 0.0);
    for (v=0; v < n_v; v++)
    {
        inflate_matrix(&(inv_design_E[v]), n_cov, n_cov, 0.0);
        inflate_matrix(&(inv_design_S[v]), n_cov, n_cov, 0.0);
    }
    
    for (r = 0; r < n_cov; r++)
    {
        for (j = r; j < n_cov; j++)
        {
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_t=0; index_t < n_week; index_t++)
                {
                    design_1.data[r][j] = design_1.data[r][j] + polygons[index_r].X[index_t][r] * polygons[index_r].X[index_t][j];
                }
            }
            design_0.data[r][j] = design_1.data[r][j];
            if (j != r)
            {
                design_1.data[j][r] = design_1.data[r][j];
                design_0.data[j][r] = design_0.data[r][j];
            }
        }
    }
    
    for (r=0; r < n_cov; r++)
    {
        for (index_r=0; index_r < n_region; index_r++)
        {
            for (index_t=0; index_t < n_week; index_t++)
            {
                design_1.data[n_cov][r] = design_1.data[n_cov][r] + polygons[index_r].X[index_t][r];
            }
        }
        design_1.data[r][n_cov] = design_1.data[n_cov][r];
    }
    
    design_1.data[n_cov][n_cov] = (double)(n_region * n_week);
    
    inverse(&design_1, &inv_design_1);
    inverse(&design_0, &inv_design_0);

    initialize_matrix(&design_basis);
    initialize_matrix(&inv_design_basis);
    inflate_matrix(&design_basis, n_knots+3, n_knots+3,0.0);
    inflate_matrix(&inv_design_basis, n_knots+3, n_knots+3, 0.0);
    
    for (v=0; v < n_v; v++)
    {
        initialize_matrix(&(inv_design_beta[v]));
        inflate_matrix(&(inv_design_beta[v]), n_knots+3, n_knots+3, 0.0);
    }
    
    for (k=0; k < n_knots+3; k++)
    {
        for (j=0; j < n_knots+3; j++)
        {
            for (t=0; t < n_week; t++)
            {
                design_basis.data[k][j] = design_basis.data[k][j] + basis[k][t]*basis[j][t];
            }
        }
    }
    inverse(&design_basis, &inv_design_basis);
    
    
    /* Initialization */
    // eta_P0
    for (v=0; v < n_v; v++)
    {
        eta_P0[v] = rnorm(5.5, 0.2, &seed);
    }
    
    // ps
    for (r = 0; r < n_region; r++)
    {
        for (t = 0; t < n_week; t++)
        {
            for (v=0; v < n_v; v++)
            {
                polygons[r].ps[0][t][v] = inv_logit(eta_P0[v]);
                polygons[r].ps[1][t][v] = 1.0 - polygons[r].ps[0][t][v];
            }
        }
    }
    
    // alpha
    for (r=0; r< n_region; r++)
    {
        for (j=0; j < n_region; j++)
        {
            for (v=0; v < n_v; v++)
            {
                polygons[r].alpha[j][v] = polygons[r].alpha_element[v] + nu * polygons[j].alpha_element[v];
            }
        }
    }
    
    // Y
    for (r=0; r < n_region; r++)
    {
        for (s=0; s < n_s; s++)
        {
            for (t=0; t < n_week; t++)
            {
                sum_Z = 0.0;
                for (v=0; v < n_v; v++)
                {
                    sum_Z = sum_Z + polygons[r].Z_agg[s][t][v];
                }
                polygons[r].Y[s][t][n_v-1] = polygons[r].Yv[s][t];
                for (v=0; v < n_v-1; v++)
                {
                    polygons[r].Y[s][t][v] = floor(polygons[r].Z_agg[s][t][v] + (polygons[r].Yv[s][t]-sum_Z)/3);
                    polygons[r].Y[s][t][n_v-1] = polygons[r].Y[s][t][n_v-1] - polygons[r].Y[s][t][v];
                }
            }
        }
    }
    
    // Ys
    for (r=0; r < n_region; r++)
    {
        for (t=0; t < n_week; t++)
        {
            for (v=0; v < n_v; v++)
            {
                polygons[r].Ys[t][v] = 0;
                for (s=0; s < n_s; s++)
                {
                    polygons[r].Ys[t][v] = polygons[r].Ys[t][v] + polygons[r].Y[s][t][v];
                }
            }
        }
    }
    
    // beta_fix
    for (v=0; v < n_v; v++)
    {
        for (k=0; k < n_knots+3; k++)
        {
            beta_fix[k][v] = rnorm(0.0, 0.2, &seed);
        }
    }

    // beta
    for (r=0; r < n_region; r++)
    {
        for (t=0; t < n_week; t++)
        {
            for (v=0; v < n_v; v++)
            {
                for (k=0; k < n_knots+3; k++)
                {
                    polygons[r].beta[t][v] = polygons[r].beta[t][v] + beta_fix[k][v] * basis[k][t];
                }
            }
        }
    }
    
    // transmission rates
    for (v=0; v < n_v; v++)
    {
        gamma0[v] = rgamma(10.0, 0.1, &seed);
        gamma1[v] = gamma0[v]/2.0;
        gamma2[v] = gamma0[v]/10.0;
    }

    // Adaption
    // initialize tuning parameter
    for (v=0; v < n_v; v++)
    {
        tuning_gamma0[v] = 0.1;
        tuning_gamma1[v] = 0.01;
        tuning_gamma2[v] = 0.01;
        tuning_eta_E[v] = 1.0;
        tuning_eta_S[v] = 0.1;
        tuning_beta_fixs[v] = 0.1;
        tuning_alphas[v] = 0.001;
    }
    
    printf("Adaption begins!\n");
    
    while (adapt_status == 0 && count_batch < max_batch)
    {
        count_batch = count_batch + 1;
        printf("Adaption continues... %d\n", count_batch);
        adapt_status = 1;
        for (i=0; i < n_batch; i++)
        {
            // Spatial random effects: alpha_element
            // province-level spatial effects
            for (v=0; v < n_v; v++)
            {
                for (index_r=0; index_r < n_prov; index_r++)
                {
                    tuned_eval_alpha_prov[index_r] = nb_eval_prov[index_r]/tuning_alphas[v];
                }
                
                for (index_r=0; index_r < n_prov; index_r++)
                {
                    old_values_prov[index_r] = alpha_element_prov[index_r][v];
                }
                
                rgmrf(new_values_prov, n_prov, n_comp, old_values_prov, tuned_eval_alpha_prov, nb_evec_prov, &seed);
                
                logacprate = -1.0 * alphas_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,sigma2_alpha, eta_E, eta_S, c);
                
                for (index_r=0; index_r < n_region; index_r++)
                {
                    polygons[index_r].alpha_element[v] = new_values_prov[polygons[index_r].province];
                }
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_j=0; index_j < n_region; index_j++)
                    {
                        polygons[index_r].alpha[index_j][v] = polygons[index_r].alpha_element[v] + nu * polygons[index_j].alpha_element[v];
                    }
                }
                for (index_r=0; index_r < n_prov; index_r++)
                {
                    alpha_element_prov[index_r][v] = new_values_prov[index_r];
                }
                n_accept_alphas[v] = n_accept_alphas[v] + 1;
                
                logacprate = logacprate + alphas_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,sigma2_alpha, eta_E, eta_S, c);
                
                if (runiform(&seed) > exp(logacprate))
                {
                    for (index_r=0; index_r < n_region; index_r++)
                    {
                        polygons[index_r].alpha_element[v] = old_values_prov[polygons[index_r].province];
                    }
                    for (index_r=0; index_r < n_region; index_r++)
                    {
                        for (index_j=0; index_j < n_region; index_j++)
                        {
                            polygons[index_r].alpha[index_j][v] = polygons[index_r].alpha_element[v] + nu * polygons[index_j].alpha_element[v];
                        }
                    }
                    for (index_r=0; index_r < n_prov; index_r++)
                    {
                        alpha_element_prov[index_r][v] = old_values_prov[index_r];
                    }
                    n_accept_alphas[v] = n_accept_alphas[v] - 1;
                }
            }

            // Variance of spatial random effects
            for (v=0; v < n_v; v++)
            {
                shape_sigma2_alpha = (n_region - n_comp - 2.0)/2 + prior_shape_alpha[v];
                rate_sigma2_alpha = 0;
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_nb=0; index_nb < polygons[index_r].n_nb; index_nb++)
                    {
                        rate_sigma2_alpha = rate_sigma2_alpha + pow(polygons[index_r].alpha_element[v]-polygons[polygons[index_r].nb[index_nb]].alpha_element[v],2);
                    }
                }
                rate_sigma2_alpha = rate_sigma2_alpha/4.0; // Above double loop count each pair of neighbors twice , plus the 2 in the denominator, divide the above results by 4
                rate_sigma2_alpha += prior_rate_alpha[v];
                
                sigma2_alpha[v] = 1/rgamma(shape_sigma2_alpha, 1.0/rate_sigma2_alpha, &seed);
                
            }
            
            // Temporal fix effects, updata beta_fix, then calculate beta; normal prior.
            //block sampling
            for (v=0; v < n_v; v++)
            {
                for (k=0; k < n_knots+3; k++)
                {
                    for (j=0; j < n_knots+3; j++)
                    {
                        inv_design_beta[v].data[k][j] = tuning_beta_fixs[v] * inv_design_basis.data[k][j];
                    }
                }
                
                for (k=0; k < n_knots+3; k++)
                {
                    old_beta_fixs[k] = beta_fix[k][v];
                }
                
                rmultinorm(new_beta_fixs, n_knots+3, old_beta_fixs, &(inv_design_beta[v]), &seed);
                
                for (k=0; k < n_knots+3; k++)
                {
                    diff_beta_fixs[k] = new_beta_fixs[k] - old_beta_fixs[k];
                }
                
                logacprate = -1.0 * beta_fixs_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,beta_fix,beta_fix_mean,beta_fix_var, eta_E, eta_S, c);
                
                for (k=0; k < n_knots+3; k++)
                {
                    beta_fix[k][v] = new_beta_fixs[k];
                }
                
                for (index_t=0; index_t < n_week; index_t++)
                {
                    for (index_r=0; index_r < n_region; index_r++)
                    {
                        for (k=0; k < n_knots+3; k++)
                        {
                            polygons[index_r].beta[index_t][v] = polygons[index_r].beta[index_t][v] + diff_beta_fixs[k]*basis[k][index_t];
                        }
                    }
                }
                logacprate = logacprate + beta_fixs_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,beta_fix,beta_fix_mean,beta_fix_var, eta_E, eta_S, c);
                
                n_accept_beta_fixs[v] = n_accept_beta_fixs[v] + 1;
                
                if (runiform(&seed) > exp(logacprate))
                {
                    n_accept_beta_fixs[v] = n_accept_beta_fixs[v] - 1;
                    
                    for (k=0; k < n_knots+3; k++)
                    {
                        beta_fix[k][v] = old_beta_fixs[k];
                    }
                    
                    for (index_t=0; index_t < n_week; index_t++)
                    {
                        for (index_r=0; index_r < n_region; index_r++)
                        {
                            for (k=0; k < n_knots+3; k++)
                            {
                                polygons[index_r].beta[index_t][v] = polygons[index_r].beta[index_t][v] - diff_beta_fixs[k]*basis[k][index_t];
                            }
                        }
                    }
                }
            }

            // Transmission rates, one by one.
            for (v=0; v < n_v; v++)
            {
                //update gamma0
                old_gamma0 = gamma0[v];
                new_gamma0 = rnorm(old_gamma0, tuning_gamma0[v], &seed);
                if (new_gamma0 > gamma_lb)
                {
                    logacprate = -1 * gamma0_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape0, rate0, eta_E, eta_S, c);
                    gamma0[v] = new_gamma0;
                    logacprate = logacprate + gamma0_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape0, rate0, eta_E, eta_S, c);
                    logacprate = logacprate + log((1-pnorm(0,old_gamma0, tuning_gamma0[v]))/(1-pnorm(0,new_gamma0,tuning_gamma0[v])));
                    
                    n_accept_gamma0[v] = n_accept_gamma0[v] + 1;
                    
                    if(runiform(&seed) > exp(logacprate))
                    {
                        n_accept_gamma0[v] = n_accept_gamma0[v] - 1;
                        
                        gamma0[v] = old_gamma0;
                    }
                }
                
                
                //update gamma1
                old_gamma1 = gamma1[v];
                new_gamma1 = rnorm(old_gamma1, tuning_gamma1[v], &seed);
                if (new_gamma1 > gamma_lb)
                {
                    logacprate = -1 * gamma1_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape1, rate1, eta_E, eta_S, c);
                    gamma1[v] = new_gamma1;
                    logacprate = logacprate + gamma1_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape1, rate1, eta_E, eta_S, c);
                    logacprate = logacprate + log((1-pnorm(0,old_gamma1, tuning_gamma1[v]))/(1-pnorm(0,new_gamma1,tuning_gamma1[v])));
                    
                    n_accept_gamma1[v] = n_accept_gamma1[v] + 1;
                    
                    if(runiform(&seed) > exp(logacprate))
                    {
                        n_accept_gamma1[v] = n_accept_gamma1[v] - 1;
                        
                        gamma1[v] = old_gamma1;
                    }
                }
                
                
                //update gamma2
                old_gamma2 = gamma2[v];
                new_gamma2 = rnorm(old_gamma2, tuning_gamma2[v], &seed);
                if (new_gamma2 > gamma_lb)
                {
                    logacprate = -1 * gamma2_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape2, rate2, eta_E, eta_S, c);
                    gamma2[v] = new_gamma2;
                    logacprate = logacprate + gamma2_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape2, rate2, eta_E, eta_S, c);
                    logacprate = logacprate + log((1-pnorm(0,old_gamma2, tuning_gamma2[v]))/(1-pnorm(0,new_gamma2,tuning_gamma2[v])));
                    
                    n_accept_gamma2[v] = n_accept_gamma2[v] + 1;
                    
                    if(runiform(&seed) > exp(logacprate))
                    {
                        n_accept_gamma2[v] = n_accept_gamma2[v] - 1;
                        gamma2[v] = old_gamma2;
                    }
                    
                }
            }
            

            // Y
            // block sampling; multinomial proposal
            for (t=2; t < n_week-2; t++)
            {
                for (s=0; s < n_s; s++)
                {
                    for (r=0; r < n_region; r++)
                    {
                        if (polygons[r].Yv[s][t] == 0) continue;
                        
                        sum_Z_int = 0;
                        for (v = 0; v < n_v; v++) sum_Z_int += ceil(polygons[r].Z_agg[s][t][v]);
                        if (polygons[r].Yv[s][t] == sum_Z_int)
                        {
                            for (v = 0; v < n_v; v++) polygons[r].Y[s][t][v] = ceil(polygons[r].Z_agg[s][t][v]);
                            continue;
                        }
                        
                        sum_pmf_Y = 0.0;
                        for (v = 0; v < n_v; v++)
                        {
                            pmf_Y[v] = polygons[r].ps[s][t][v] * lambda_fun(polygons, r, t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
                            sum_pmf_Y += pmf_Y[v];
                        }
                        for (v = 0; v < n_v; v++) pmf_Y[v] = pmf_Y[v] / sum_pmf_Y;
                        
                        for (v = 0; v < n_v; v++) old_Y_block[v] = polygons[r].Y[s][t][v];
                        
                        rmultinomial(n_v, polygons[r].Yv[s][t] - sum_Z_int, pmf_Y, new_Y_block, &seed);
                        for (v = 0; v < n_v; v++)
                        {
                            new_Y_block[v] += ceil(polygons[r].Z_agg[s][t][v]);
                        }
                        
                        logacprate = 0.0;
                        for (v = 0; v < n_v; v++) logacprate += logpmfY_fun(old_Y_block[v],polygons,r,s,t,v,gamma0,gamma1,gamma2,eta_E,eta_S,c) - (old_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v])) * log(pmf_Y[v]) + lgamma(old_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v]) + 1);
                        
                        logacprate = -1.0*logacprate;
                        
                        for (v = 0; v < n_v; v++) logacprate += logpmfY_fun(new_Y_block[v],polygons,r,s,t,v,gamma0,gamma1,gamma2,eta_E,eta_S,c) - (new_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v])) * log(pmf_Y[v]) + lgamma(new_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v]) + 1);
                        
                        if (runiform(&seed) < exp(logacprate))
                        {
                            for (v = 0; v < n_v; v++)
                            {
                                polygons[r].Y[s][t][v] = new_Y_block[v];
                                polygons[r].Ys[t][v] = polygons[r].Ys[t][v] - old_Y_block[v] + new_Y_block[v];
                            }
                        }
                    }
                }
            }
            
            // pathogenicity
            for (v=0; v < n_v; v++)
            {
                para1_beta = 1.0;
                para2_beta = 1.0;
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_t=0; index_t < n_week; index_t++)
                    {
                        para1_beta += polygons[index_r].Y[0][index_t][v];
                        para2_beta += polygons[index_r].Y[1][index_t][v];
                    }
                }
                
                new_ps0 = rbeta(para1_beta, para2_beta, &seed);
                
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_t=0; index_t < n_week; index_t++)
                    {
                        polygons[index_r].ps[0][index_t][v] = new_ps0;
                        polygons[index_r].ps[1][index_t][v] = 1.0 - new_ps0;
                    }
                }
            }

            // covariates effects on susceptibility
            for (v=0; v < n_v; v++)
            {
                for (r=0; r < n_cov; r++)
                {
                    for (j=0; j < n_cov; j++)
                    {
                        inv_design_E[v].data[r][j] = tuning_eta_E[v] * inv_design_0.data[r][j];
                        inv_design_S[v].data[r][j] = tuning_eta_S[v] * inv_design_0.data[r][j];
                    }
                }
            }
            
            for (v=0; v < n_v; v++)
            {
                for (index_j=0; index_j < n_cov; index_j++)
                {
                    old_eta_S[index_j] = eta_S[index_j][v];
                }
                
                rmultinorm(new_eta_S, n_cov, old_eta_S, &inv_design_S[v], &seed);
                
                logacprate = 0.0;
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_t=2; index_t < n_week; index_t++)
                    {
                        lambda_old = lambda_fun(polygons, index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
                        for (index_j=0; index_j < n_cov; index_j++) eta_S[index_j][v] = new_eta_S[index_j];
                        lambda_new = lambda_fun(polygons, index_r, index_t, v, gamma0, gamma1,gamma2, eta_E, eta_S, c);
                        for (index_j=0; index_j < n_cov; index_j++) eta_S[index_j][v] = old_eta_S[index_j];
                        logacprate = logacprate + lambda_old - lambda_new + polygons[index_r].Ys[index_t][v] * log(lambda_new/lambda_old);
                    }
                }
                
                if (runiform(&seed) < exp(logacprate))
                {
                    n_accept_eta_S[v] = n_accept_eta_S[v] + 1;
                    
                    for (index_j=0; index_j < n_cov; index_j++)
                    {
                        eta_S[index_j][v] = new_eta_S[index_j];
                    }
                }
            }
            
            // covariates effects on environment
            for (v=0; v < n_v; v++)
            {
                for (index_j=0; index_j < n_cov; index_j++)
                {
                    old_eta_E[index_j] = eta_E[index_j][v];
                }
                
                rmultinorm(new_eta_E, n_cov, old_eta_E, &inv_design_E[v], &seed);
                
                logacprate = -1.0 * eta_E_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,eta_E,eta_S,eta_E_mean,eta_E_var,c);
                
                for (index_j=0; index_j < n_cov; index_j++)
                {
                    eta_E[index_j][v] = new_eta_E[index_j];
                }
                
                logacprate = logacprate + eta_E_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,eta_E,eta_S,eta_E_mean,eta_E_var,c);
                
                n_accept_eta_E[v] = n_accept_eta_E[v] + 1;
                
                if (runiform(&seed) > exp(logacprate))
                {
                    n_accept_eta_E[v] = n_accept_eta_E[v] - 1;
                    for (index_j=0; index_j < n_cov; index_j++)
                    {
                        eta_E[index_j][v] = old_eta_E[index_j];
                    }
                }
            }
        
        }
    
        // adjust tuning parameters
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_gamma0[v]) / (double)n_batch;
            if (accept_rate < single_lower)
            {
                tuning_gamma0[v] = tuning_gamma0[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > single_upper)
            {
                tuning_gamma0[v] = tuning_gamma0[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            n_accept_gamma0[v] = 0;
        }
        
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_gamma1[v]) / (double)n_batch;
            if (accept_rate < single_lower)
            {
                tuning_gamma1[v] = tuning_gamma1[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > single_upper)
            {
                tuning_gamma1[v] = tuning_gamma1[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            n_accept_gamma1[v] = 0;
        }

        for (v=0; v < n_v; v++)
        {
            
            accept_rate = ((double)n_accept_gamma2[v]) / (double)n_batch;
            if (accept_rate < single_lower)
            {
                tuning_gamma2[v] = tuning_gamma2[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > single_upper)
            {
                tuning_gamma2[v] = tuning_gamma2[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            n_accept_gamma2[v] = 0;
        }
        
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_alphas[v]) / (double)n_batch;
            if (accept_rate < group_lower)
            {
                tuning_alphas[v] = tuning_alphas[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > group_upper)
            {
                tuning_alphas[v] = tuning_alphas[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            
            n_accept_alphas[v] = 0;
            
        }
        
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_beta_fixs[v]) / (double)n_batch;
            if (accept_rate < group_lower)
            {
                tuning_beta_fixs[v] = tuning_beta_fixs[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > group_upper)
            {
                tuning_beta_fixs[v] = tuning_beta_fixs[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            
            n_accept_beta_fixs[v] = 0;
        }
        
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_eta_E[v]) / (double)n_batch;
            if (accept_rate < group_lower)
            {
                tuning_eta_E[v] = tuning_eta_E[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > group_upper)
            {
                tuning_eta_E[v] = tuning_eta_E[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            
            n_accept_eta_E[v] = 0;
        }
        
        for (v=0; v < n_v; v++)
        {
            accept_rate = ((double)n_accept_eta_S[v]) / (double)n_batch;
            if (accept_rate < group_lower)
            {
                tuning_eta_S[v] = tuning_eta_S[v] * (1.0 - tuning_adjust);
                adapt_status = 0;
            }
            else if (accept_rate > group_upper)
            {
                tuning_eta_S[v] = tuning_eta_S[v] * (1.0 + tuning_adjust);
                adapt_status = 0;
            }
            
            n_accept_eta_S[v] = 0;
        }
    }
    
    if (count_batch == max_batch) printf("Adaption aborted!\n");
    else printf("Adaption done!\n");
    
    // Gibbs Sampler
    for (i=0; i < n_sim + n_burn; i++)
    {
        if ((i+1)%1000==0) printf("Iteration %d\n", i+1);
        // Spatial random effects: alpha_element
        // province-level spatial effects
        for (v=0; v < n_v; v++)
        {
            for (index_r=0; index_r < n_prov; index_r++)
            {
                tuned_eval_alpha_prov[index_r] = nb_eval_prov[index_r]/tuning_alphas[v];
            }
            
            for (index_r=0; index_r < n_prov; index_r++)
            {
                old_values_prov[index_r] = alpha_element_prov[index_r][v];
            }
            
            rgmrf(new_values_prov, n_prov, n_comp, old_values_prov, tuned_eval_alpha_prov, nb_evec_prov, &seed);
            
            logacprate = -1.0 * alphas_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,sigma2_alpha, eta_E, eta_S, c);
            
            for (index_r=0; index_r < n_region; index_r++)
            {
                polygons[index_r].alpha_element[v] = new_values_prov[polygons[index_r].province];
            }
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_j=0; index_j < n_region; index_j++)
                {
                    polygons[index_r].alpha[index_j][v] = polygons[index_r].alpha_element[v] + nu * polygons[index_j].alpha_element[v];
                }
            }
            for (index_r=0; index_r < n_prov; index_r++)
            {
                alpha_element_prov[index_r][v] = new_values_prov[index_r];
            }
            n_accept_alphas[v] = n_accept_alphas[v] + 1;
            
            logacprate = logacprate + alphas_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,sigma2_alpha, eta_E, eta_S, c);
            
            if (runiform(&seed) > exp(logacprate))
            {
                for (index_r=0; index_r < n_region; index_r++)
                {
                    polygons[index_r].alpha_element[v] = old_values_prov[polygons[index_r].province];
                }
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (index_j=0; index_j < n_region; index_j++)
                    {
                        polygons[index_r].alpha[index_j][v] = polygons[index_r].alpha_element[v] + nu * polygons[index_j].alpha_element[v];
                    }
                }
                for (index_r=0; index_r < n_prov; index_r++)
                {
                    alpha_element_prov[index_r][v] = old_values_prov[index_r];
                }
                n_accept_alphas[v] = n_accept_alphas[v] - 1;
            }
        }
        
        // Variance of spatial random effects
        for (v=0; v < n_v; v++)
        {
            shape_sigma2_alpha = (n_region - n_comp - 2.0)/2 + prior_shape_alpha[v];
            rate_sigma2_alpha = 0;
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_nb=0; index_nb < polygons[index_r].n_nb; index_nb++)
                {
                    rate_sigma2_alpha = rate_sigma2_alpha + pow(polygons[index_r].alpha_element[v]-polygons[polygons[index_r].nb[index_nb]].alpha_element[v],2);
                }
            }
            rate_sigma2_alpha = rate_sigma2_alpha/4.0; // Above double loop count each pair of neighbors twice , plus the 2 in the denominator, divide the above results by 4
            rate_sigma2_alpha += prior_rate_alpha[v];
            sigma2_alpha[v] = 1/rgamma(shape_sigma2_alpha, 1.0/rate_sigma2_alpha, &seed);
            
        }
        
        // Temporal fix effects, updata beta_fix, then calculate beta; normal prior.
        //block sampling
        for (v=0; v < n_v; v++)
        {
            for (k=0; k < n_knots+3; k++)
            {
                for (j=0; j < n_knots+3; j++)
                {
                    inv_design_beta[v].data[k][j] = tuning_beta_fixs[v] * inv_design_basis.data[k][j];
                }
            }
            
            for (k=0; k < n_knots+3; k++)
            {
                old_beta_fixs[k] = beta_fix[k][v];
            }
            
            rmultinorm(new_beta_fixs, n_knots+3, old_beta_fixs, &(inv_design_beta[v]), &seed);
            
            for (k=0; k < n_knots+3; k++)
            {
                diff_beta_fixs[k] = new_beta_fixs[k] - old_beta_fixs[k];
            }
            
            logacprate = -1.0 * beta_fixs_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,beta_fix,beta_fix_mean,beta_fix_var, eta_E, eta_S, c);
            
            for (k=0; k < n_knots+3; k++)
            {
                beta_fix[k][v] = new_beta_fixs[k];
            }
            
            for (index_t=0; index_t < n_week; index_t++)
            {
                for (index_r=0; index_r < n_region; index_r++)
                {
                    for (k=0; k < n_knots+3; k++)
                    {
                        polygons[index_r].beta[index_t][v] = polygons[index_r].beta[index_t][v] + diff_beta_fixs[k]*basis[k][index_t];
                    }
                }
            }
            logacprate = logacprate + beta_fixs_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,beta_fix,beta_fix_mean,beta_fix_var, eta_E, eta_S, c);
            
            n_accept_beta_fixs[v] = n_accept_beta_fixs[v] + 1;
            
            if (runiform(&seed) > exp(logacprate))
            {
                n_accept_beta_fixs[v] = n_accept_beta_fixs[v] - 1;
                
                for (k=0; k < n_knots+3; k++)
                {
                    beta_fix[k][v] = old_beta_fixs[k];
                }
                
                for (index_t=0; index_t < n_week; index_t++)
                {
                    for (index_r=0; index_r < n_region; index_r++)
                    {
                        for (k=0; k < n_knots+3; k++)
                        {
                            polygons[index_r].beta[index_t][v] = polygons[index_r].beta[index_t][v] - diff_beta_fixs[k]*basis[k][index_t];
                        }
                    }
                }
            }
        }
        
        // Transmission rates
        // one by one
         for (v=0; v < n_v; v++)
         {
             //update gamma0
             old_gamma0 = gamma0[v];
             new_gamma0 = rnorm(old_gamma0, tuning_gamma0[v], &seed);
             if (new_gamma0 > gamma_lb)
             {
                 n_accept_gamma0[v] = n_accept_gamma0[v] + 1;
                 
                 logacprate = -1 * gamma0_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape0, rate0, eta_E, eta_S, c);
                 gamma0[v] = new_gamma0;
                 logacprate = logacprate + gamma0_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape0, rate0, eta_E, eta_S, c);
                 logacprate = logacprate + log((1-pnorm(0,old_gamma0, tuning_gamma0[v]))/(1-pnorm(0,new_gamma0,tuning_gamma0[v])));
         
                 if(runiform(&seed) > exp(logacprate))
                 {
                     gamma0[v] = old_gamma0;
                     n_accept_gamma0[v] = n_accept_gamma0[v] - 1;
                 }
             }

             //update gamma1
             old_gamma1 = gamma1[v];
             new_gamma1 = rnorm(old_gamma1, tuning_gamma1[v], &seed);
             if (new_gamma1 > gamma_lb)
             {
                 n_accept_gamma1[v] = n_accept_gamma1[v] + 1;
                 logacprate = -1 * gamma1_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape1, rate1, eta_E, eta_S, c);
                 gamma1[v] = new_gamma1;
                 logacprate = logacprate + gamma1_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape1, rate1, eta_E, eta_S, c);
                 logacprate = logacprate + log((1-pnorm(0,old_gamma1, tuning_gamma1[v]))/(1-pnorm(0,new_gamma1,tuning_gamma1[v])));
         
                 if(runiform(&seed) > exp(logacprate))
                 {
                     gamma1[v] = old_gamma1;
                     n_accept_gamma1[v] = n_accept_gamma1[v] - 1;
                 }
             }

             //update gamma2
             old_gamma2 = gamma2[v];
             new_gamma2 = rnorm(old_gamma2, tuning_gamma2[v], &seed);
             if (new_gamma2 > gamma_lb)
             {
                 n_accept_gamma2[v] = n_accept_gamma2[v] + 1;
                 logacprate = -1 * gamma2_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape2, rate2, eta_E, eta_S, c);
                 gamma2[v] = new_gamma2;
                 logacprate = logacprate + gamma2_logdensity_fun(polygons, v, gamma0, gamma1, gamma2, shape2, rate2, eta_E, eta_S, c);
                 logacprate = logacprate + log((1-pnorm(0,old_gamma2, tuning_gamma2[v]))/(1-pnorm(0,new_gamma2,tuning_gamma2[v])));
         
                 if(runiform(&seed) > exp(logacprate))
                 {
                     gamma2[v] = old_gamma2;
                     n_accept_gamma2[v] = n_accept_gamma2[v] - 1;
                 }
             }
         }
        
        // Y
        // block sampling; multinomial proposal
        for (t=2; t < n_week-2; t++)
        {
            for (s=0; s < n_s; s++)
            {
                for (r=0; r < n_region; r++)
                {
                    if (polygons[r].Yv[s][t] == 0) continue;
                    
                    sum_Z_int = 0;
                    for (v = 0; v < n_v; v++) sum_Z_int += ceil(polygons[r].Z_agg[s][t][v]);
                    if (polygons[r].Yv[s][t] == sum_Z_int)
                    {
                        for (v = 0; v < n_v; v++) polygons[r].Y[s][t][v] = ceil(polygons[r].Z_agg[s][t][v]);
                        continue;
                    }
                    
                    sum_pmf_Y = 0.0;
                    for (v = 0; v < n_v; v++)
                    {
                        pmf_Y[v] = polygons[r].ps[s][t][v] * lambda_fun(polygons, r, t, v, gamma0, gamma1, gamma2, eta_E, eta_S, c);
                        sum_pmf_Y += pmf_Y[v];
                    }
                    for (v = 0; v < n_v; v++) pmf_Y[v] = pmf_Y[v] / sum_pmf_Y;
                        
		    for (v = 0; v < n_v; v++) old_Y_block[v] = polygons[r].Y[s][t][v];
                   
		    rmultinomial(n_v, polygons[r].Yv[s][t] - sum_Z_int, pmf_Y, new_Y_block, &seed);
                    for (v = 0; v < n_v; v++)
                    {
                        new_Y_block[v] += ceil(polygons[r].Z_agg[s][t][v]);
                    }
                    
                    logacprate = 0.0;
                    for (v = 0; v < n_v; v++) logacprate += logpmfY_fun(old_Y_block[v],polygons,r,s,t,v,gamma0,gamma1,gamma2,eta_E,eta_S,c) - (old_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v])) * log(pmf_Y[v]) + lgamma(old_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v]) + 1);
                    
                    logacprate = -1.0*logacprate;
                    
                    for (v = 0; v < n_v; v++) logacprate += logpmfY_fun(new_Y_block[v],polygons,r,s,t,v,gamma0,gamma1,gamma2,eta_E,eta_S,c) - (new_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v])) * log(pmf_Y[v]) + lgamma(new_Y_block[v] - ceil(polygons[r].Z_agg[s][t][v]) + 1);
                    
                    if (runiform(&seed) < exp(logacprate))
                    {
                        for (v = 0; v < n_v; v++)
                        {
                            polygons[r].Y[s][t][v] = new_Y_block[v];
                            polygons[r].Ys[t][v] = polygons[r].Ys[t][v] - old_Y_block[v] + new_Y_block[v];
                        }
                    }
                }
            }
        }

        // pathogenicity
        // no covariates in pathogenicity part
        for (v=0; v < n_v; v++)
        {
            para1_beta = 1.0;
            para2_beta = 1.0;
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_t=0; index_t < n_week; index_t++)
                {
                    para1_beta += polygons[index_r].Y[0][index_t][v];
                    para2_beta += polygons[index_r].Y[1][index_t][v];
                }
            }
            
            new_ps0 = rbeta(para1_beta, para2_beta, &seed);
            
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_t=0; index_t < n_week; index_t++)
                {
                    polygons[index_r].ps[0][index_t][v] = new_ps0;
                    polygons[index_r].ps[1][index_t][v] = 1.0 - new_ps0;
                }
            }
        }


        // covariates effects on susceptibility
        for (v=0; v < n_v; v++)
        {
            for (index_j=0; index_j < n_cov; index_j++)
            {
                old_eta_S[index_j] = eta_S[index_j][v];
            }
            
            rmultinorm(new_eta_S, n_cov, old_eta_S, &inv_design_S[v], &seed);
            
            logacprate = 0.0;
            for (index_r=0; index_r < n_region; index_r++)
            {
                for (index_t=2; index_t < n_week; index_t++)
                {
                    lambda_old = lambda_fun(polygons, index_r, index_t, v, gamma0, gamma1, gamma2, eta_E,eta_S,c);
                    for (index_j=0; index_j < n_cov; index_j++) eta_S[index_j][v] = new_eta_S[index_j];
                    lambda_new = lambda_fun(polygons, index_r, index_t, v, gamma0, gamma1, gamma2, eta_E, eta_S,c);
                    for (index_j=0; index_j < n_cov; index_j++) eta_S[index_j][v] = old_eta_S[index_j];
                    logacprate = logacprate + lambda_old - lambda_new + polygons[index_r].Ys[index_t][v] * log(lambda_new/lambda_old);
                }
            }
            
            if (runiform(&seed) < exp(logacprate))
            {
                n_accept_eta_S[v] = n_accept_eta_S[v] + 1;
                
                for (index_j=0; index_j < n_cov; index_j++)
                {
                    eta_S[index_j][v] = new_eta_S[index_j];
                }
            }
        }


        // covariates effects on environment
        for (v=0; v < n_v; v++)
        {
            for (index_j=0; index_j < n_cov; index_j++)
            {
                old_eta_E[index_j] = eta_E[index_j][v];
            }
            
            rmultinorm(new_eta_E, n_cov, old_eta_E, &inv_design_E[v], &seed);
            
            logacprate = -1.0 * eta_E_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,eta_E,eta_S,eta_E_mean,eta_E_var,c);
            
            for (index_j=0; index_j < n_cov; index_j++)
            {
                eta_E[index_j][v] = new_eta_E[index_j];
            }
            
            logacprate = logacprate + eta_E_logdensity_fun(polygons,v,gamma0,gamma1,gamma2,eta_E,eta_S,eta_E_mean,eta_E_var,c);
            
            n_accept_eta_E[v] = n_accept_eta_E[v] + 1;
            
            if (runiform(&seed) > exp(logacprate))
            {
                n_accept_eta_E[v] = n_accept_eta_E[v] - 1;
                
                for (index_j=0; index_j < n_cov; index_j++)
                {
                    eta_E[index_j][v] = old_eta_E[index_j];
                }
            }
        }
        
        //save to file, use scan in R to read the files for plots
        if (i >= n_burn && (i+1)%thin==0)
        {
            sprintf(filename, "real_data_analysis/results/Y_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open Y output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    for (t=0; t < n_week; t++)
                    {
                        for (s=0; s < n_s; s++)
                        {
                            for (r=0; r < n_region; r++)
                            {
                                fprintf(file, "%d ", polygons[r].Y[s][t][v]);
                                if (r == n_region-1) fprintf(file, "\n");
                            }
                        }
                    }
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/gamma0_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open gamma0 output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    fprintf(file, "%f ", gamma0[v]);
                    if (v == n_v - 1) fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/gamma1_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open gamma1 output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    fprintf(file, "%f ", gamma1[v]);
                    if (v == n_v - 1) fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/gamma2_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open gamma2 output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    fprintf(file, "%f ", gamma2[v]);
                    if (v == n_v - 1) fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/alpha_element_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open alpha_element output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    for (r=0; r < n_region; r++)
                    {
                        fprintf(file, "%f ", polygons[r].alpha_element[v]);
                    }
                    fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/sigma2_alpha_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open sigma2_alpha output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    fprintf(file, "%f ", sigma2_alpha[v]);
                    if (v == n_v - 1) fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/beta_fix_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL)
            {
                printf("Can't open beta_fix output file!\n");
            }
            else
            {
                for (v=0; v < n_v; v++)
                {
                    for (k=0; k < n_knots+3; k++)
                    {
                        fprintf(file, "%f ", beta_fix[k][v]);
                    }
                    fprintf(file, "\n");
                }
            }
            fclose(file);
            
            sprintf(filename, "real_data_analysis/results/eta_P0_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL) printf("Can't open eta_P0 output file!\n");
            else
            {
                for (v=0; v < n_v; v++)
                {
                    fprintf(file, "%f ", logit(polygons[0].ps[0][0][v]));
                }
                fprintf(file, "\n");
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/eta_S_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL) printf("Can't open eta_S output file!\n");
            else
            {
                for (v=0; v < n_v; v++)
                {
                    for (j=0; j < n_cov; j++)
                    {
                        fprintf(file, "%f ", eta_S[j][v]);
                    }
                    fprintf(file, "\n");
                }
            }
            fclose(file);

            sprintf(filename, "real_data_analysis/results/eta_E_output_realdata_%d.txt", index_chain);
            file = fopen(filename,"a");
            if (file == NULL) printf("Cant't open eta_E output file!\n");
            else
            {
                for (v=0; v < n_v; v++)
                {
                    for (j=0; j < n_cov; j++)
                    {
                        fprintf(file, "%f ", eta_E[j][v]);
                    }
                    fprintf(file,"\n");
                }
            }
            fclose(file);
        }
    }
    
    // acceptance rate output
    sprintf(filename, "real_data_analysis/results/acceptance_rate_output_realdata_%d.txt", index_chain);
    file = fopen(filename,"a");
    fprintf(file, "gamma0\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_gamma0[v]/(n_sim + n_burn));
    
    fprintf(file, "gamma1\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_gamma1[v]/(n_sim + n_burn));
    
    fprintf(file, "gamma2\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_gamma2[v]/(n_sim + n_burn));
    
    fprintf(file, "eta_E\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_eta_E[v]/(n_sim + n_burn));
    
    fprintf(file, "eta_S\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_eta_S[v]/(n_sim + n_burn));

    fprintf(file, "alpha\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_alphas[v]/(n_sim + n_burn));

    fprintf(file, "beta_fix\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, (double)n_accept_beta_fixs[v]/(n_sim + n_burn));
    
    fclose(file);
    
    sprintf(filename, "real_data_analysis/results/tuning_para_output_realdata_%d.txt", index_chain);
    file = fopen(filename,"a");
    fprintf(file, "gamma0\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_gamma0[v]);
    
    fprintf(file, "gamma1\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_gamma1[v]);
    
    fprintf(file, "gamma2\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_gamma2[v]);
    
    fprintf(file, "eta_E\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_eta_E[v]);
    
    fprintf(file, "eta_S\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_eta_S[v]);
    
    fprintf(file, "alpha\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_alphas[v]);
    
    fprintf(file, "beta_fix\n");
    for (v=0; v < n_v; v++) fprintf(file, "v=%d %f\n", v, tuning_beta_fixs[v]);
    
    fclose(file);
    
    for(r=0; r < n_region; r++)
    {
        free(polygons[r].nb);
        free_2d_array_int(polygons[r].Yv);
        free_2d_array_int(polygons[r].Ys);
        free_3d_array_int(polygons[r].Z);
        free_3d_array_double(polygons[r].Z_agg);
        free_3d_array_int(polygons[r].Y);
        
        free(polygons[r].alpha_element);
        free_2d_array_double(polygons[r].alpha);
        free_2d_array_double(polygons[r].beta);
        
        free_2d_array_double(polygons[r].X);
        free_3d_array_double(polygons[r].ps);
        
    }
    
    free(polygons);
    free(gamma0);
    free(gamma1);
    free(gamma2);
    free_2d_array_double(alpha_element_prov);
    free(old_values_prov);
    free(new_values_prov);
    free(nb_eval_prov);
    free(tuned_eval_alpha_prov);
    free(sigma2_alpha);
    free_2d_array_double(basis);
    free_2d_array_double(nb_evec_prov);
    free_2d_array_double(beta_fix);
    free_2d_array_double(beta_fix_mean);
    free_2d_array_double(beta_fix_var);
    
    free_2d_array_double(eta_E);
    free_2d_array_double(eta_S);
    free(eta_P0);
    
    deflate_matrix(&design_0);
    deflate_matrix(&design_1);
    for (v=0; v < n_v; v++)
    {
        deflate_matrix(&(inv_design_E[v]));
        deflate_matrix(&(inv_design_S[v]));
    }
    deflate_matrix(&design_basis);
    deflate_matrix(&inv_design_basis);
    
    t2 = time(NULL);
    printf("\n use time: %f seconds\n", difftime(t2, t1));
    
    return 0;

}



