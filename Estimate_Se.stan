data {
	
	int N;						// Number of observations
	int J;						// Number of individuals
	
	int participant_ID_index[N];
	
	real<lower=0.0> days_since_seroconv[N];
	
	real response[N];
	
	int hosp_by_ind[J];
	
	int N_rep;
	
	int K;
	int days_for_Se[K];
	
	real cutoff;
	
}
parameters {
	
	real lambda;
	
	vector[2] beta_severity;
	
	vector[J] intercept_raw;
	
	real<lower=0.0> sigma;
	
	real<lower=0.0> tau[2];
	
}
transformed parameters {
	
	vector[J] intercept;
	
	vector[N] mu;
	
	// Re-parametrization for efficiency
	for(j in 1:J) {
		
		intercept[j] = beta_severity[hosp_by_ind[j]] + (tau[hosp_by_ind[j]] * intercept_raw[j]);
		
	}
	
	for(n in 1:N) {
		
		mu[n] = intercept[participant_ID_index[n]] - (lambda * days_since_seroconv[n]);
		
	}
	
}
model {
	
	// Likelihood
	response ~ normal(mu, sigma);
	
	// Prior
	intercept_raw ~ normal(0.0, 1.0);
	
}
generated quantities {
	
	matrix<lower=0.0, upper=1.0>[2,K] sensitivity_by_severity_time;	// Sensitivity: 2 classes of severity, K days
	
	real ratio_sensitivity_by_time[K];								// Sensitivity ratio: K days
	
	// Calculate time-varying sensitivity by severity class
	for(k in 1:K) {
		
		vector[N_rep] intercept_new_tau1;
		vector[N_rep] intercept_new_tau2;
		
		vector[N_rep] intercept_new_sigma;
		
		vector[N_rep] Y_tilde_1;
		vector[N_rep] Y_tilde_2;
		
		vector[N_rep] Y_tilde_binary_1;
		vector[N_rep] Y_tilde_binary_2;
		
		// Simulate new data
		for(x in 1:N_rep) {
			
			intercept_new_tau1[x] = normal_rng(0.0, tau[1]);
			intercept_new_tau2[x] = normal_rng(0.0, tau[2]);
			intercept_new_sigma[x] = normal_rng(0.0, sigma);
			
			Y_tilde_1[x] = beta_severity[1] + intercept_new_tau1[x] + intercept_new_sigma[x] - (lambda * (days_for_Se[k] * 1.0));
			Y_tilde_2[x] = beta_severity[2] + intercept_new_tau2[x] + intercept_new_sigma[x] - (lambda * (days_for_Se[k] * 1.0));
			
			Y_tilde_binary_1[x] = Y_tilde_1[x] >= cutoff ? 1.0 : 0.0;
			Y_tilde_binary_2[x] = Y_tilde_2[x] >= cutoff ? 1.0 : 0.0;
			
		}
		
		// Calculate sensitivity
		sensitivity_by_severity_time[1,k] = sum(Y_tilde_binary_1) / (1.0 * N_rep);
		sensitivity_by_severity_time[2,k] = sum(Y_tilde_binary_2) / (1.0 * N_rep);
		
		ratio_sensitivity_by_time[k] = sensitivity_by_severity_time[1,k] / sensitivity_by_severity_time[2,k];
	}
	
}
