data {
  int<lower=1> N;
  int<lower=1> N_animals;

  // counts
  int<lower=0> n_eyed[N];
  int<lower=1> n_eggs[N];

  // IDs
  int<lower=1, upper=N_animals> dam_id[N];
  int<lower=1, upper=N_animals> sire_id[N];

  // fixed effects
  int<lower=1> K_year;
  int<lower=1> K_boost;
  int<lower=1, upper=K_year>  year_id[N];
  int<lower=1, upper=K_boost> boost_id[N];

  // ejaculate effect
  int<lower=1> n_ejac;
  int<lower=1, upper=n_ejac> ejaculate_id[N];

  // per-observation proxy
  int<lower=1> n_e_proxy;
  int<lower=1, upper=n_e_proxy> e_proxy_id[N];

  // additive relationship matrix
  matrix[N_animals, N_animals] A;
}

transformed data {
  matrix[N_animals, N_animals] A_jit = A;
  for (i in 1:N_animals) A_jit[i,i] += 1e-9;
  matrix[N_animals, N_animals] L_A = cholesky_decompose(A_jit);
}

parameters {
  // intercepts
  real alpha_female;     // probit
  real alpha_capacity;   // log(capacity)

  // fixed effects
  vector[K_year]  year_female_raw;
  vector[K_year]  year_capacity_raw;
  vector[K_boost] boost_capacity_raw;

  // additive genetic
  matrix[N_animals, 2] Z_a;
  vector<lower=0>[2]   tau_a;
  cholesky_factor_corr[2] L_Omega_a;        // genetic correlation

  // random effects
  vector[n_ejac]    a_ejac;
  vector[n_e_proxy] e_c;
  vector[n_e_proxy] e_f;

  // their scales
  real<lower=0> sigma_ejac;
  real<lower=0> sigma_ec;
  real<lower=0> sigma_ef;
}

transformed parameters {
  // centering
  vector[K_year]  year_female    = year_female_raw   - mean(year_female_raw);
  vector[K_year]  year_capacity  = year_capacity_raw - mean(year_capacity_raw);
  vector[K_boost] boost_capacity = boost_capacity_raw - mean(boost_capacity_raw);

  // G = L_Sigma_a * L_Sigma_a'
  matrix[2,2] L_Sigma_a = diag_pre_multiply(tau_a, L_Omega_a);

  // 1=female, 2=capacity
  matrix[N_animals, 2] a_mat = L_A * Z_a * L_Sigma_a';
  vector[N_animals] a_female   = col(a_mat, 1);
  vector[N_animals] a_capacity = col(a_mat, 2);

  // linear predictors
  vector[N] female_lin;
  vector[N] capacity_lin;
  vector<lower=0, upper=1>[N] FemaleFertility;
  vector<lower=0>[N]          r;    // lambda rate
  vector<lower=0, upper=1>[N] p;

  for (i in 1:N) {
    female_lin[i] =
      alpha_female +
      year_female[ year_id[i] ] +
      a_female[   dam_id[i] ] +
      e_f[        e_proxy_id[i] ];

    capacity_lin[i] =
      alpha_capacity +
      year_capacity[ year_id[i] ] +
      boost_capacity[ boost_id[i] ] +
      a_capacity[ sire_id[i] ] +
      a_ejac[     ejaculate_id[i] ] +
      e_c[        e_proxy_id[i] ];

    // links
    FemaleFertility[i] = Phi_approx(female_lin[i]);

    // re-center on log-r
    {
      real cap_eta = capacity_lin[i] - log(n_eggs[i]);  // log r
      r[i] = exp(cap_eta);
    }

    p[i] = FemaleFertility[i] * (1 - exp(-r[i]));
  }
}
  
model {
  // intercepts and fixed effects
  alpha_female   ~ normal(0, 1);
  alpha_capacity ~ normal(9, 0.5);
  year_female_raw    ~ normal(0, 0.5);
  year_capacity_raw  ~ normal(0, 0.5);
  boost_capacity_raw ~ normal(0, 0.5);

  // additive genetics
  tau_a[1] ~ normal(0.30, 0.15);   // female
  tau_a[2] ~ normal(0.30, 0.15);   // capacity
  L_Omega_a ~ lkj_corr_cholesky(6);
  to_vector(Z_a) ~ std_normal();

  // random effects
  a_ejac ~ normal(0, sigma_ejac);
  e_c    ~ normal(0, sigma_ec);
  e_f    ~ normal(0, sigma_ef);

  sigma_ejac ~ normal(0.30, 0.15);
  sigma_ec   ~ normal(0.20, 0.10);
  sigma_ef   ~ normal(0.25, 0.15);

  // BINOMIAL likelihood
  for (i in 1:N)
    n_eyed[i] ~ binomial(n_eggs[i], p[i]);
}

generated quantities {
  vector[N] log_lik;
  int n_eyed_rep[N];

  // genetic (co)variance matrix and correlation
  matrix[2,2] Sigma_a =
    quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_a), tau_a);
  real Va_female      = Sigma_a[1,1];
  real Va_capacity    = Sigma_a[2,2];
  real CovA_f_c       = Sigma_a[1,2];
  real rhoA_f_c       = CovA_f_c / sqrt(Va_female * Va_capacity);

  // non-genetic variance components (liability scales)
  real Vejac_capacity = square(sigma_ejac);
  real Vres_capacity  = square(sigma_ec);
  real Vres_female    = square(sigma_ef);

  real Vp_female      = Va_female + Vres_female + 1;           // probit residual = 1
  real h2_female      = Va_female / Vp_female;

  // male capacity
  real Vp_capacity    = Va_capacity + Vejac_capacity + Vres_capacity;
  real h2_capacity    = Va_capacity / Vp_capacity;

  // log-lik & PPC (binomial)
  for (i in 1:N) {
    log_lik[i]    = binomial_lpmf(n_eyed[i] | n_eggs[i], p[i]);
    n_eyed_rep[i] = binomial_rng(n_eggs[i], p[i]);
  }
}
