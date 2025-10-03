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
  real alpha_female;            // probit
  real alpha_male;              // probit

  // fixed effects
  vector[K_year]  year_female_raw;
  vector[K_year]  year_male_raw;
  vector[K_boost] boost_male_raw;

  // additive genetics (female, male)
  matrix[N_animals, 2] Z_a;
  vector<lower=0>[2]   tau_a;
  cholesky_factor_corr[2] L_Omega_a;        // correlation

  // RE
  vector[n_ejac]    a_ejac; 
  vector[n_e_proxy] e_m;                    // male resdual
  vector[n_e_proxy] e_f;                    // female residual

  real<lower=0> sigma_ejac;
  real<lower=0> sigma_em;
  real<lower=0> sigma_ef;
}

transformed parameters {
  // centerin
  vector[K_year]  year_female = year_female_raw - mean(year_female_raw);
  vector[K_year]  year_male   = year_male_raw   - mean(year_male_raw);
  vector[K_boost] boost_male  = boost_male_raw  - mean(boost_male_raw);

  // genetic Cholesky
  matrix[2,2] L_Sigma_a = diag_pre_multiply(tau_a, L_Omega_a);

  // additive genetic effects
  matrix[N_animals, 2] a_mat = L_A * Z_a * L_Sigma_a';
  vector[N_animals] a_female = col(a_mat, 1);
  vector[N_animals] a_male   = col(a_mat, 2);

  // linear predictors
  vector[N] female_lin;
  vector[N] male_lin;
  vector<lower=0, upper=1>[N] FemaleFertility;
  vector<lower=0, upper=1>[N] MaleFertility;
  vector<lower=0, upper=1>[N] p;

  for (i in 1:N) {
    female_lin[i] =
      alpha_female +
      year_female[ year_id[i] ] +
      a_female[   dam_id[i] ] +
      e_f[        e_proxy_id[i] ];

    male_lin[i] =
      alpha_male +
      year_male[ year_id[i] ] +
      boost_male[ boost_id[i] ] +
      a_male[    sire_id[i] ] +
      a_ejac[    ejaculate_id[i] ] +
      e_m[       e_proxy_id[i] ];

    FemaleFertility[i] = Phi_approx(female_lin[i]);
    MaleFertility[i]   = Phi_approx(male_lin[i]);
    p[i]               = FemaleFertility[i] * MaleFertility[i];
  }
}

model {
  // intercepts
  alpha_female ~ normal(0, 1);
  alpha_male   ~ normal(0, 1);

  // fixed effects
  year_female_raw ~ normal(0, 0.5);
  year_male_raw   ~ normal(0, 0.5);
  boost_male_raw  ~ normal(0, 0.5);

  // genetics
  tau_a       ~ normal(0.30, 0.15);
  L_Omega_a   ~ lkj_corr_cholesky(6);
  to_vector(Z_a) ~ std_normal();

  // random effects
  a_ejac ~ normal(0, sigma_ejac);
  e_m    ~ normal(0, sigma_em);
  e_f    ~ normal(0, sigma_ef);

  sigma_ejac ~ normal(0.25, 0.15);
  sigma_em   ~ normal(0.25, 0.15);
  sigma_ef   ~ normal(0.25, 0.15);

  // binomial liklihood
  for (i in 1:N)
    n_eyed[i] ~ binomial(n_eggs[i], p[i]);
}

generated quantities {
  vector[N] log_lik;
  int n_eyed_rep[N];

  // genetic correlation
  matrix[2,2] Sigma_a =
    quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_a), tau_a);
  real Va_female   = Sigma_a[1,1];
  real Va_male     = Sigma_a[2,2];
  real CovA_f_m    = Sigma_a[1,2];
  real rhoA_f_m    = CovA_f_m / sqrt(Va_female * Va_male);

  // non-gen VCs
  real Vejac_male  = square(sigma_ejac);
  real Vres_male   = square(sigma_em);
  real Vres_female = square(sigma_ef);

  // heritabilities on probit scales
  real Vp_female      = Va_female + Vres_female + 1;
  real h2_female      = Va_female / Vp_female;

  real Vp_male        = Va_male + Vejac_male + Vres_male + 1;
  real h2_male        = Va_male / Vp_male;

  // log-lik & PPC
  for (i in 1:N) {
    log_lik[i]    = binomial_lpmf(n_eyed[i] | n_eggs[i], p[i]);
    n_eyed_rep[i] = binomial_rng(n_eggs[i], p[i]);
  }
}
