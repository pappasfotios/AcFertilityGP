data {
  int<lower=1> N;
  int<lower=1> N_animals;

  int<lower=0> n_eyed[N];
  int<lower=0> n_eggs[N];

  int<lower=1, upper=N_animals> dam_id[N];
  int<lower=1, upper=N_animals> sire_id[N];

  // fixed effects
  int<lower=1> K_year;
  int<lower=1> K_boost;
  int<lower=1, upper=K_year>  year_id[N];
  int<lower=1, upper=K_boost> boost_id[N];

  // permanent envs
  int<lower=1> n_fpe;
  int<lower=1> n_mpe;

  int<lower=1, upper=n_fpe> fpe_id[N];
  int<lower=1, upper=n_mpe> mpe_id[N];

  // seperate residual terms
  int<lower=1> n_e_proxy;
  int<lower=1, upper=n_e_proxy> e_proxy_id[N];

  // additive genetic relationships matrix
  matrix[N_animals, N_animals] A;
}

transformed data {
  matrix[N_animals, N_animals] A_jit = A;
  for (i in 1:N_animals) A_jit[i,i] += 1e-9;
  matrix[N_animals, N_animals] L_A = cholesky_decompose(A_jit);
}

parameters {
  // intercepts
  real alpha_female;
  real alpha_capacity;

  // fixed effects
  vector[K_year]  year_female_raw;
  vector[K_year]  year_capacity_raw;
  vector[K_boost] boost_capacity_raw;

  // ## Bivariate ##
  // Z_a ~ iid N(0,1); then a = L_A * Z_a * L_Sigma_a' (N_animals x 2)
  matrix[N_animals, 2] Z_a;
  cholesky_factor_corr[2] L_Ga;   // genetic correlation (Cholesky)
  vector<lower=0>[2] tau_a;       // sqrt(var_A) per trait

  // PEs
  vector[n_fpe] z_fpe;
  vector[n_mpe] z_mpe;
  real<lower=0> sigma_fpe;
  real<lower=0> sigma_mpe;

  // residual + covariance structure
  matrix[2, n_e_proxy] Z_e;
  cholesky_factor_corr[2] L_Omega_e;
  vector<lower=0>[2] tau_e;

  // beta–binomial overdispersion
  real<lower=0> conc;
}

transformed parameters {
  // fixed effects
  vector[K_year]  year_female    = year_female_raw  - mean(year_female_raw);
  vector[K_year]  year_capacity  = year_capacity_raw- mean(year_capacity_raw);
  vector[K_boost] boost_capacity = boost_capacity_raw - mean(boost_capacity_raw);

  // == Additive genetic effects ==
  matrix[2,2] L_Sigma_a = diag_pre_multiply(tau_a, L_Ga); // lower-Cholesky of Ga
  matrix[N_animals, 2] a_animal = L_A * Z_a * L_Sigma_a'; // N_animals × 2
  // 1 = female-liability BV, 2 = capacity BV
  vector[N_animals] a_female   = a_animal[,1];
  vector[N_animals] a_capacity = a_animal[,2];

  // PEs
  vector[n_fpe] a_fpe = sigma_fpe * z_fpe;
  vector[n_mpe] a_mpe = sigma_mpe * z_mpe;

  // residuals
  matrix[2, n_e_proxy] E = diag_pre_multiply(tau_e, L_Omega_e) * Z_e;
  vector[n_e_proxy] e_f = (E[1])';
  vector[n_e_proxy] e_c = (E[2])';

  // linear predictors and corresponding values
  vector[N] female_lin;
  vector[N] capacity_lin;
  vector<lower=0, upper=1>[N] FemaleFertility;
  vector<lower=0>[N]          Capacity;
  vector<lower=0, upper=1>[N] p;

  for (i in 1:N) {
    female_lin[i] =
      alpha_female +
      year_female[ year_id[i] ] +
      a_female[   dam_id[i] ] +
      a_fpe[      fpe_id[i] ] +
      e_f[        e_proxy_id[i] ];

    capacity_lin[i] =
      alpha_capacity +
      year_capacity[ year_id[i] ] +
      boost_capacity[ boost_id[i] ] +
      a_capacity[ sire_id[i] ] +
      a_mpe[      mpe_id[i] ] +
      e_c[        e_proxy_id[i] ];

    FemaleFertility[i] = Phi_approx(female_lin[i]);      // fast probit link
    Capacity[i]        = log1p_exp(capacity_lin[i]);     // softplus (ln(1 + exp(x))) to ensure positivity
    p[i] = FemaleFertility[i] * (1 - exp(-Capacity[i] / n_eggs[i]));  // formula
  }
}

model {
  // priors
  alpha_female   ~ normal(0, 0.7);
  alpha_capacity ~ normal(0, 0.7);

  year_female_raw    ~ normal(0, 0.7);
  year_capacity_raw  ~ normal(0, 0.7);
  boost_capacity_raw ~ normal(0, 0.7);

  // additive genetic + correlation
  tau_a[1] ~ normal(0, 1);   // female-liability sd
  tau_a[2] ~ normal(0, 0.3); // capacity sd (scale)
  L_Ga     ~ lkj_corr_cholesky(2);
  to_vector(Z_a) ~ std_normal();

  // PEs
  sigma_fpe ~ normal(0, 1);
  sigma_mpe ~ normal(0, 1);
  z_fpe ~ std_normal();
  z_mpe ~ std_normal();

  // residuals
  tau_e     ~ normal(0, 1);
  L_Omega_e ~ lkj_corr_cholesky(2);
  to_vector(Z_e) ~ std_normal();

  // beta–binomial likelhood
  conc ~ gamma(2, 0.5);
  for (i in 1:N) {
    real alpha = p[i] * conc;
    real beta  = (1 - p[i]) * conc;
    n_eyed[i] ~ beta_binomial(n_eggs[i], alpha, beta);
  }
}

generated quantities {
  vector[N] log_lik;
  int n_eyed_rep[N];

  // == Genetic (co-)variances on liability scales ==
  matrix[2,2] Omega_a = multiply_lower_tri_self_transpose(L_Ga); // genetic correlation matrix
  matrix[2,2] Ga      = multiply_lower_tri_self_transpose(diag_pre_multiply(tau_a, L_Ga)); // genetic covariance matrix
  real rho_A          = Omega_a[1,2]; // genetic correlation (female vs capacity)
  real CovA_fc        = Ga[1,2];

  // Trait-specific variance components
  real Va_female   = Ga[1,1];
  real Va_capacity = Ga[2,2];

  real Vpe_female  = square(sigma_fpe);
  real Vpe_capacity= square(sigma_mpe);

  real Vres_female = square(tau_e[1]); // proxy/resid on liability scale
  real Vres_capacity = square(tau_e[2]);

  // Latent-scale heritabilities
  real Vp_female   = Va_female + Vpe_female + Vres_female + 1; // +1 for probit residual
  real Vp_capacity = Va_capacity + Vpe_capacity + Vres_capacity;

  real h2_female   = Va_female   / Vp_female;
  real h2_capacity = Va_capacity / Vp_capacity;

  // residual correlation
  matrix[2,2] Omega_e = multiply_lower_tri_self_transpose(L_Omega_e);
  real rho_e = Omega_e[1,2];

  // Pointwise log-likelihood and posterior checks
  for (i in 1:N) {
    real alpha = p[i] * conc;
    real beta  = (1 - p[i]) * conc;
    log_lik[i] = beta_binomial_lpmf(n_eyed[i] | n_eggs[i], alpha, beta);
    n_eyed_rep[i] = beta_binomial_rng(n_eggs[i], alpha, beta);
  }
}
