//
//  datastructure.h
//  
//
//  Created by Xueying Tang on 11/13/13.
//
//

/* data struture of a prefecture */

typedef struct {
    int id; // id of prefectures
    int code; // zonecode of prefectures
    
    int n_nb; // number of neighbors
    int *nb; // neighbors
    
    int **Yv; // data Y_rt^{+s}
    int ***Z; // data lab test cases
    
    int ***Y;
    int **Ys; // determined by Y
    
    double *alpha_element; // h2h spatial random effects
    double **alpha; // determined by alpha_element
    double **beta; // h2h temporal effects; determined by beta_fix (parameters) and beta_random
    
    double **X; // covariates
    
    double ***ps; // probability for mild cases
    
}POLYGON;