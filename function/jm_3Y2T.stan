/*
This code was adapted from the codes in the following paper, 
Brilleman, S., Crowther, M., Moreno-Betancur, M., Buros Novik, J., & Wolfe, R. (2018). Joint longitudinal and time-to-event models via Stan. Proceedings of StanCon 2018, 1-18.

The original code joint models two continuous longitudinal outcomes and one event outcome using stan. 
In the code below, we extended to 3 longitudinal outcomes (Normal, Normal, and Bernulli distributions) and one main event outcome subject to a dependent censoring outcome.
*/
  
functions {
  /*
  * Evaluate the linear predictor for the glmer submodel
  * @param X Design matrix (include intercept) for fe; 
  * @param Z Design matrix for re, for a single grouping factor
  * @param Z_id Group indexing for Z
  * @param beta Vector of population level parameters
  * @param bMat Matrix of group level params (random effects)
  * @param shift Num of columns to shift in bMat to link with the m-th biomarker 
  * @return A vector containing the linear predictor for the glmer submodel
  */  
  vector evaluate_eta(matrix X, vector[] Z, int[] Z_id, 
                      vector beta, matrix bMat, int shift) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num FE params for specific Y
    int p = size(Z);    // num RE params for specific Y, size() returns the nrow()
    vector[N] eta;
    
    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);
    
    for (k in 1:p)
      for (n in 1:N)
        eta[n] = eta[n] + (bMat[Z_id[n], k + shift]) * Z[k,n];
    return eta;
  } 

  /*
  * Get the indices corresponding to the lower tri of a square matrix
  * @param dim The number of rows in the square matrix
  * @return A vector of indices
  */ 
  int[] lower_tri_indices(int dim) {
    int indices[dim + choose(dim, 2)];
    int mark = 1;
    for (r in 1:dim) {
      for (c in r:dim) {
        indices[mark] = (r - 1) * dim + c;
        mark = mark + 1;
      }
    }
    return indices;
  } 
}


data {
  //-----------------------------------------------------------//
  //----- Longitudinal submodels
  
  int<lower=2,upper=3> K;       // number of long. submodels 
  
  // population level dimensions
  int<lower=0> y_N[K];          // num observations
  int<lower=0> y_K[K];          // num predictors

  // population level data
  vector[y_N[1]] y1;            // response vectors
  vector[y_N[2]] y2;
  int<lower=0, upper=1> y3[y_N[3]];   
  matrix[y_N[1],y_K[1]] y1_X;   // fe design matrix
  matrix[y_N[2],y_K[2]] y2_X; 
  matrix[y_N[3],y_K[3]] y3_X; 
  
  // group level dimensions
  int<lower=0> b_N;             // num groups
  int<lower=0> b_K;             // total num params
  int<lower=0> b_KM[K];         // num params in each submodel

  // group level data
  vector[y_N[1]] y1_Z[b_KM[1]]; // re design matrix
  vector[y_N[2]] y2_Z[b_KM[2]];
  vector[y_N[3]] y3_Z[b_KM[3]];
  int<lower=0> y1_Z_id[y_N[1]]; // group ids (long id)
  int<lower=0> y2_Z_id[y_N[2]];
  int<lower=0> y3_Z_id[y_N[3]];
  
  //-----------------------------------------------------------//
  //----- Event submodel for main event
  
  // data for calculating event submodel linear predictor in GK quadrature
  int<lower=0> e_K;             // num. of predictors in event submodel
  int<lower=0> a_K;             // num. of association parameters
  int<lower=0> Npat;            // num. individuals 
  int<lower=0> Nevents;         // num. events (ie. not censored) 
  int<lower=0> qnodes;          // num. of nodes for GK quadrature 
  int<lower=0> Npat_times_qnodes; 
  int<lower=0> nrow_e_Xq;       // num. rows in event submodel predictor matrix
  vector[nrow_e_Xq] e_times;    // event times and quadrature points 
  matrix[nrow_e_Xq,e_K] e_Xq;   // predictor matrix (event submodel) 
  int<lower=0> basehaz_df;      // df for B-splines baseline hazard 
  matrix[nrow_e_Xq,basehaz_df] basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[Npat_times_qnodes] qwts;         // GK quadrature weights with (b-a)/2 scaling

  // data for calculating long. submodel linear predictor in GK quadrature
  int<lower=0> nrow_y_Xq[K]; // num. rows in long. predictor matrix at quadpoints 
  matrix[nrow_y_Xq[1],y_K[1]] y1_Xq;   // fe design matrix at quadpoints
  matrix[nrow_y_Xq[2],y_K[2]] y2_Xq; 
  matrix[nrow_y_Xq[3],y_K[3]] y3_Xq;
  vector[nrow_y_Xq[1]] y1_Zq[b_KM[1]]; // re design matrix at quadpoints
  vector[nrow_y_Xq[2]] y2_Zq[b_KM[2]];
  vector[nrow_y_Xq[3]] y3_Zq[b_KM[3]];
  int<lower=0> y1_Zq_id[nrow_y_Xq[1]]; // group indexing for re design matrix
  int<lower=0> y2_Zq_id[nrow_y_Xq[2]]; 
  int<lower=0> y3_Zq_id[nrow_y_Xq[3]];
  
  //-----------------------------------------------------------//
  //----- Event submodel for dependent censoring
  
  // data for calculating event submodel linear predictor in GK quadrature
  int<lower=0> e_K2;             // num. of predictors in event submodel
  int<lower=0> a_K2;             // num. of association parameters
  int<lower=0> Npat2;            // num. individuals 
  int<lower=0> Nevents2;         // num. events (ie. not censored) 
  int<lower=0> qnodes2;          // num. of nodes for GK quadrature 
  int<lower=0> Npat_times_qnodes2; 
  int<lower=0> nrow_e_Xq2;       // num. rows in event submodel predictor matrix
  vector[nrow_e_Xq2] e_times2;   // event times and quadrature points
  matrix[nrow_e_Xq2,e_K2] e_Xq2; // predictor matrix (event submodel)
  int<lower=0> basehaz_df2;      // df for B-splines baseline hazard
  matrix[nrow_e_Xq2,basehaz_df2] basehaz_X2; // design matrix (basis terms) for baseline hazard
  vector[Npat_times_qnodes2] qwts2;          // GK quadrature weights with (b-a)/2 scaling

  // data for calculating long. submodel linear predictor in GK quadrature
  int<lower=0> nrow_y_Xq2[K]; // num. rows in long. predictor matrix at quadpoints
  matrix[nrow_y_Xq2[1],y_K[1]] y1_Xq2;   // fe design matrix at quadpoints
  matrix[nrow_y_Xq2[2],y_K[2]] y2_Xq2; 
  matrix[nrow_y_Xq2[3],y_K[3]] y3_Xq2;
  vector[nrow_y_Xq2[1]] y1_Zq2[b_KM[1]]; // re design matrix at quadpoints
  vector[nrow_y_Xq2[2]] y2_Zq2[b_KM[2]];
  vector[nrow_y_Xq2[3]] y3_Zq2[b_KM[3]];
  int<lower=0> y1_Zq_id2[nrow_y_Xq2[1]]; // group indexing for re design matrix
  int<lower=0> y2_Zq_id2[nrow_y_Xq2[2]]; 
  int<lower=0> y3_Zq_id2[nrow_y_Xq2[3]];
}


transformed data {
  // indexing used to extract lower tri of RE covariance matrix
  int b_cov_idx[b_K + choose(b_K, 2)];
  if (b_K > 0) b_cov_idx = lower_tri_indices(b_K);
}


parameters {
  vector[y_K[1]] y1_beta;     // coefs in long. submodels
  vector[y_K[2]] y2_beta;     
  vector[y_K[3]] y3_beta;     
  real<lower=0> y1_aux;       // residual error SDs 
  real<lower=0> y2_aux;       
							  
  vector[e_K] e_beta;         // coefs in event submodel (log hazard ratios)
  vector[a_K] a_beta;         // assoc params (log hazard ratios)
  vector[basehaz_df] e_aux;   // coefs for baseline hazard  

  vector[e_K2] e_beta2;       // coefs in event submodel (log hazard ratios)
  vector[a_K2] a_beta2;       // assoc params (log hazard ratios)
  vector[basehaz_df2] e_aux2; // coefs for baseline hazard    

  // group level params     
  matrix[b_N, b_K] b_mat;               // group level params 
  cholesky_factor_cov[b_K] b_cholesky;  // cholesky factor of cov matrix
}


model {
  //-----------------------------------------------------------//
  //----- Log-likelihood for longitudinal submodels
  {
    // declare linear predictors
    vector[y_N[1]] y1_eta; 
    vector[y_N[2]] y2_eta;
	vector[y_N[3]] y3_eta;

    // evaluate linear predictor for each long. submodel
    y1_eta = evaluate_eta(y1_X, y1_Z, y1_Z_id, y1_beta, b_mat, 0);
    y2_eta = evaluate_eta(y2_X, y2_Z, y2_Z_id, y2_beta, b_mat, b_KM[1]);
	y3_eta = inv_logit(evaluate_eta(y3_X, y3_Z, y3_Z_id, y3_beta, b_mat, b_KM[1] + b_KM[2]));
    
    // increment the target with the log-lik
    target += normal_lpdf(y1 | y1_eta, y1_aux);
    target += normal_lpdf(y2 | y2_eta, y2_aux);
	target += bernoulli_lpmf(y3 | y3_eta);
  }
  
  //-----------------------------------------------------------//
  //----- Log-likelihood for main event submodel
  {
    vector[nrow_y_Xq[1]] y1_eta_q; 
    vector[nrow_y_Xq[2]] y2_eta_q;
	vector[nrow_y_Xq[3]] y3_eta_q;
    vector[nrow_e_Xq] e_eta_q; 
    vector[nrow_e_Xq] log_basehaz;  // log baseline hazard AT event time and quadrature points
    vector[nrow_e_Xq] log_haz_q;    // log hazard AT event time and quadrature points
    vector[Nevents] log_haz_etimes; // log hazard AT the event time only
    vector[Npat_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points
    
    // Event submodel: linear predictor at event time and quadrature points
    e_eta_q = e_Xq * e_beta;
    
    // Long. submodel: linear predictor at event time and quadrature points
    y1_eta_q = evaluate_eta(y1_Xq, y1_Zq, y1_Zq_id, y1_beta, b_mat, 0);
    y2_eta_q = evaluate_eta(y2_Xq, y2_Zq, y2_Zq_id, y2_beta, b_mat, b_KM[1]);
	y3_eta_q = evaluate_eta(y3_Xq, y3_Zq, y3_Zq_id, y3_beta, b_mat, b_KM[1] + b_KM[2]);
  
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    e_eta_q = e_eta_q + a_beta[1] * y1_eta_q + a_beta[2] * y2_eta_q + a_beta[3] * y3_eta_q;
    
    // Log baseline hazard at event time and quadrature points
    log_basehaz = basehaz_X * e_aux; 
    
    // Log hazard at event time and quadrature points
    log_haz_q = log_basehaz + e_eta_q;
  
    // Log hazard at event times only
    log_haz_etimes = head(log_haz_q, Nevents);
  
    // Log hazard at quadrature points only
    log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);
 
    // Log likelihood for event submodel
    target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));  
  }    
  
  //-----------------------------------------------------------//
  //----- Log-likelihood for dependent censoring submodel 
  {
    vector[nrow_y_Xq2[1]] y1_eta_q2; 
    vector[nrow_y_Xq2[2]] y2_eta_q2;
	vector[nrow_y_Xq2[3]] y3_eta_q2;
    vector[nrow_e_Xq2] e_eta_q2; 
    vector[nrow_e_Xq2] log_basehaz2;  // log baseline hazard AT event time and quadrature points
    vector[nrow_e_Xq2] log_haz_q2;    // log hazard AT event time and quadrature points
    vector[Nevents2] log_haz_etimes2; // log hazard AT the event time only
    vector[Npat_times_qnodes2] log_haz_qtimes2; // log hazard AT the quadrature points
    
    // Event submodel: linear predictor at event time and quadrature points
    e_eta_q2 = e_Xq2 * e_beta2;
    
    // Long. submodel: linear predictor at event time and quadrature points
    y1_eta_q2 = evaluate_eta(y1_Xq2, y1_Zq2, y1_Zq_id2, y1_beta, b_mat, 0);
    y2_eta_q2 = evaluate_eta(y2_Xq2, y2_Zq2, y2_Zq_id2, y2_beta, b_mat, b_KM[1]);
	y3_eta_q2 = evaluate_eta(y3_Xq2, y3_Zq2, y3_Zq_id2, y3_beta, b_mat, b_KM[1] + b_KM[2]);
  
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    e_eta_q2 = e_eta_q2 + a_beta2[1] * y1_eta_q2 + a_beta2[2] * y2_eta_q2 + a_beta2[3] * y3_eta_q2;
    
    // Log baseline hazard at event time and quadrature points
    log_basehaz2 = basehaz_X2 * e_aux2; 
    
    // Log hazard at event time and quadrature points
    log_haz_q2 = log_basehaz2 + e_eta_q2;
  
    // Log hazard at event times only
    log_haz_etimes2 = head(log_haz_q2, Nevents2);
  
    // Log hazard at quadrature points only
    log_haz_qtimes2 = tail(log_haz_q2, Npat_times_qnodes2);
 
    // Log likelihood for event submodel
    target += sum(log_haz_etimes2) - dot_product(qwts2, exp(log_haz_qtimes2));  
  }    

  //-----------------------------------------------------------//
  // ----- Log-likelihood for random effects
  for(n in 1:b_N)
	target += multi_normal_cholesky_lpdf(b_mat[n, 1:b_K] | rep_vector(0, b_K), b_cholesky); 
}
