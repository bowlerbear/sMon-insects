data {
  
  //number of sites
  int<lower = 1> n_site;
  
  // survey-level detection covariates - multiple surveys at each site
  int<lower = 1> total_surveys;
  int<lower = 1> m_p;// number of covariates
  matrix[total_surveys, m_p] X_p;// detection covariate matrix
  
  // survey level information 
  int<lower = 1, upper = n_site> site[total_surveys];
  int<lower = 0, upper = 1> y[total_surveys];//survey-level observation - presence/absence
  int<lower = 0, upper = total_surveys> start_idx[n_site];
  int<lower = 0, upper = total_surveys> end_idx[n_site];
  
  // spline data for occupancy model for site variation
  int Y[n_site];  // site-level observations - presence/absence
  // data for splines
  int Ks;  // number of linear effects
  matrix[n_site, Ks] Xs;  // design matrix for the linear effects
  // data for spline t2(x,y)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[n_site, knots_1[1]] Zs_1_1;
  matrix[n_site, knots_1[2]] Zs_1_2;
  matrix[n_site, knots_1[3]] Zs_1_3;
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[n_site];
  
  // number of surveys at each site
  int<lower = 0> n_survey[n_site];
  
  //data for predictions
  int<lower = 1> complete_N;
  // data for splines
  int complete_Ks;  // number of linear effects
  matrix[complete_N, complete_Ks] complete_Xs;  // design matrix for the linear effects
  // data for spline t2(x,y,yearIndex)
  int complete_nb_1;  // number of bases
  int complete_knots_1[complete_nb_1];  // number of knots
  // basis function matrices
  matrix[complete_N, complete_knots_1[1]] complete_Zs_1_1;
  matrix[complete_N, complete_knots_1[2]] complete_Zs_1_2;
  matrix[complete_N, complete_knots_1[3]] complete_Zs_1_3;
}
parameters {
  
  //detection covariates
  vector[m_p] beta_p;
  
  //spline parameters
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // spline coefficients
  // parameters for spline t2(x,y,yearIndex)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  vector[knots_1[2]] zs_1_2;
  vector[knots_1[3]] zs_1_3;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  real<lower=0> sds_1_2;  // standard deviations of spline coefficients
  real<lower=0> sds_1_3;  // standard deviations of spline coefficients
}
transformed parameters {
  
  // detection linear predictor
  vector[total_surveys]logit_p = X_p * beta_p;
  
  //spline transformed parameters
  vector[knots_1[1]] s_1_1;
  vector[knots_1[2]] s_1_2;
  vector[knots_1[3]] s_1_3;
  s_1_1 = sds_1_1 * zs_1_1;
  s_1_2 = sds_1_2 * zs_1_2;
  s_1_3 = sds_1_3 * zs_1_3;
}
model {
  
  //occupancy model linear predictor
  vector[n_site] logit_psi = Intercept + rep_vector(0.0, n_site) + Xs * bs + Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3;
  
  //transform to log scale
  vector[n_site] log_psi = log_inv_logit(logit_psi);
  vector[n_site] log1m_psi = log1m_inv_logit(logit_psi);
  
  //priors for detection model
  beta_p ~ normal(0, 1);
  
  //priors for spline occupancy model
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_1_2);
  target += std_normal_lpdf(zs_1_3);
  
  //likelihood

  for (i in 1:n_site) {
      if (any_seen[i]) {
        // if species is seen at the site at least once - site is occupied
        target += log_psi[i] 
        + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] | 
                                 logit_p[start_idx[i]:end_idx[i]]);
      } else {
        // if species never seen at the site - site may or may not be occupied
        target += log_sum_exp(
          log_psi[i] + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] |
                                              logit_p[start_idx[i]:end_idx[i]]), 
          log1m_psi[i]
        );
      }
    }
}

// to change eventually
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
  
  //psi on original scale
  vector[complete_N] psi = inv_logit(Intercept + rep_vector(0.0, complete_N) + complete_Xs * bs + complete_Zs_1_1 * s_1_1 + complete_Zs_1_2 * s_1_2 + complete_Zs_1_3 * s_1_3);
  
}
