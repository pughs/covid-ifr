data {
  int num_loc;
  int num_country;
  
  // death data
  int length_dstar; 
  int D_star[length_dstar]; // # recorded deaths for each age bin/loc
  vector[length_dstar] N; // population for each age bin/loc for deaths
  
  // seroprevalence data
  int length_rstar; 
  int R_star[length_rstar]; // # tested positive each age bin/loc
  int n[length_rstar]; // # tested each age bin/loc
  
  // sensitivity and specificity 
  int num_test;
  int<lower=0> sens_n[num_test]; // number of positive controls tested
  int<lower=0> sens_x[num_test]; // number of positive controls that tested positive
  int<lower=0> spec_n[num_test]; // number of negative controls tested
  int<lower=0> spec_x[num_test]; // number of negative controls that tested negative
  
  // numeric integration
  int n_int; // number of integration bins in [0,100]
  real h; // the width of the integration bins
  matrix[n_int, num_loc] f_expanded; // population age distribution for ages on the integration grid
  vector[length_rstar] f_normalizing_constants_rstar;
  vector[length_dstar] f_normalizing_constants_dstar;
  
  // Covariates
  int ncolX;
  int ncolZ;
  matrix[n_int, ncolX] X; 
  matrix[n_int, ncolZ] Z;
  
  // The rows of X that correspond to death age bins
  int min_match_X_to_death[length_dstar]; 
  int max_match_X_to_death[length_dstar];
  
  // The rows of Z that correspond to the seroprevalence age bins
  int min_match_Z_to_sero[length_rstar];
  int max_match_Z_to_sero[length_rstar];
  
  // Indicators
  int study_match_pi[num_loc]; // match test used to pi_A
  int study_match_dstar[length_dstar]; // which location
  int study_match_rstar[length_rstar]; // which location
  
  // 
  int num_mult_loc;
  int num_not_mult_loc;
  int num_country_mult_loc;
  int mult_loc_ind[num_mult_loc];
  int not_mult_loc_ind[num_not_mult_loc];
  int country_inds[num_mult_loc]; 
  
  vector[ncolZ-1] gamma_world; // mean of seroprevalence prior
  vector<lower=0>[ncolZ-1] tau; // standard deviation of seroprevalence prior 
}

parameters {
  // Prevalence model
  matrix[ncolZ, num_loc] gamma_raw; // prevalence covariates for each location
  
  // Test characteristics
  vector<lower=0, upper=1>[num_test] sens; // sensitivity for each test
  vector<lower=0, upper=1>[num_test] spec; // specificity for each test
  
  // IFR
  matrix[ncolX, num_loc] beta_raw; // death covariates for each location
  vector[ncolX] beta_world;
  vector[num_country_mult_loc] beta_country;
  vector<lower=0>[ncolX] sigma; // sd for betas
  real<lower=0> sigma_country; // sd for country intercepts
}

transformed parameters {
  matrix[ncolZ, num_loc] gamma;
  matrix[ncolX, num_loc] beta;
  
  gamma[1,] = -1 + 1.5*gamma_raw[1,];
  for (i in 2:ncolZ) {
    gamma[i,] = gamma_world[i-1] + tau[i-1]*gamma_raw[i,];
  }
  
  beta[1,mult_loc_ind] = to_row_vector(beta_world[1] + beta_country[country_inds]) + to_row_vector(sigma[1]*beta_raw[1, mult_loc_ind]);
  beta[1,not_mult_loc_ind] = to_row_vector(beta_world[1] + sqrt(sigma[1]^2 + sigma_country^2)*beta_raw[1, not_mult_loc_ind]);
  for (i in 2:ncolX) { 
    beta[i,] = beta_world[i] + sigma[i]*beta_raw[i,];
  }
}

model {
  // Integration
  matrix[n_int, num_loc] p_a; // test positivity
  matrix[n_int, num_loc] ifr_a;
  matrix[n_int, num_loc] pi_a;
  matrix[n_int, num_loc] p_a_int;
  matrix[n_int, num_loc] lambda_a_int;
  vector[length_rstar] p_A;
  vector[length_dstar] lambda_A;
  
  // Seropositivity //
  pi_a = inv_logit(Z * gamma);
  
  // Positivity //
  for (i in 1:num_loc) {
    p_a[,i] = sens[study_match_pi[i]] * pi_a[,i] + 
        (1 - spec[study_match_pi[i]]) * (1-pi_a[,i]); 
  }
        
  // IFR //
  ifr_a = exp(X*beta);
        
  // Integrate positivity
  p_a_int = p_a .* f_expanded;
  for (i in 1:length_rstar) {
    p_A[i] = f_normalizing_constants_rstar[i]*(
                (0.5*(p_a_int[min_match_Z_to_sero[i], study_match_rstar[i]]+
                     p_a_int[max_match_Z_to_sero[i], study_match_rstar[i]])+
                sum(p_a_int[min_match_Z_to_sero[i]:(max_match_Z_to_sero[i]-1), study_match_rstar[i]]))
            );
  }
  
  // Integrate to lambda 
  lambda_a_int = ifr_a .* pi_a .* f_expanded;
  for (i in 1:length_dstar) {
    lambda_A[i] = f_normalizing_constants_dstar[i]*(
        0.5*(lambda_a_int[min_match_X_to_death[i], study_match_dstar[i]] + 
              lambda_a_int[max_match_X_to_death[i], study_match_dstar[i]]) + 
        sum(lambda_a_int[min_match_X_to_death[i]:(max_match_X_to_death[i]-1), study_match_dstar[i]])
      );
  }
  
  
  
  // Model //
  D_star ~ poisson(N .* lambda_A );
  R_star ~ binomial(n, p_A); 
  
  // Priors //
    
  // Prevalence
  for (i in 1:ncolZ) {
    gamma_raw[i,] ~ std_normal();
  }
  
  // IFR 
  for (i in 1:(ncolX)) { 
    beta_raw[i,] ~ std_normal();
  }
  
  beta_world ~ normal(0, 5); 
  sigma ~ normal(0, 2); 
  
  beta_country ~ normal(0, sigma_country); 
  sigma_country ~ normal(0,2);
  
  // Test characteristics
  sens_x ~ binomial(sens_n, sens);
  spec_x ~ binomial(spec_n, spec);
  
  sens ~ beta(10,1);
  spec ~ beta(50, 1); 
}
