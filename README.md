# AcFertilityGP
This project aims to estimate genetic parameters for observed and latent fertility traits of Arctic (superior) charr. It includes a multi-trait model to analyze sperm phenotypes together with growth ([Sperm_traits](https://github.com/pappasfotios/AcFertilityGP/tree/main/Sperm_traits)), and a second multivariate model to evaluate female and male gamete outputs and respective body sizes ([SexLimited](https://github.com/pappasfotios/AcFertilityGP/tree/main/SexLimited)). Additionally, a sophisticated Bayesian framework is presented to account for distinct female and male genetic contributions to overall reproductive success outcomes ([FertRate_Models](https://github.com/pappasfotios/AcFertilityGP/tree/main/FertRate_Models)).


![Untitled (11)](https://github.com/user-attachments/assets/f0f92cf3-4865-4af9-a1d6-276b492bd190)

## Example usage of stan model (model 2):

```
library(rstan)

stan_data <- readRDS("stan_data.rds")

init = function() list(
  alpha_capacity = 8.5,   # ~log(5000)
  alpha_female   = 0,
  tau_a = c(0.2, 0.2),
  sigma_ejac = 0.2, sigma_ec = 0.2, sigma_ef = 0.2,
  L_Omega_a = diag(2)
)


options(mc.cores = 8)
rstan_options(auto_write = TRUE)

fit <- stan(
  file = "FM2.stan",
  data = stan_data,
  chains = 4,
  iter = 10000,
  warmup = 4000,
  init = init,
  thin = 1,
  cores = 4,
  seed = 321,
  control = list(adapt_delta=0.99, max_treedepth=15)
)

```

<img width="1875" height="2520" alt="Obs_vs_latent_EBVs (1)" src="https://github.com/user-attachments/assets/209f0206-30d7-4426-9f49-e8ef58df6d8c" />
