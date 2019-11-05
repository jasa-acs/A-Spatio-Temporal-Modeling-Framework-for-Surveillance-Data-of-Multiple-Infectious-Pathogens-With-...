//
//  datastructure_agg.h
//  
//
//  Created by Tang Xueying  on 11/13/13.
//
//

/*struture of a prefecture*/
typedef struct {
    int id; // id of prefectures
    int code; // zonecode of prefectures
    int province;
    
    int n_nb; // number of neighbors
    int *nb; // neighbors
    
    int **Yv;
    int ***Z;
    double ***Z_agg;
    int ***Y;
    int **Ys;
    
    double *alpha_element;
    double **alpha;
    double **beta;
    
    double **X;
    
    double ***ps;
    
}POLYGON;