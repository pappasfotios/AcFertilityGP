# AcFertilityGP
Genetic parameters for fertility traits of Arctic (superior) charr


## Part 1: data from cross success over multiple genarations (full-records, egg count and eyed embryos count)
```
library(brms)

fertility_model <- brm(
  bf(perc_eye ~ FemaleFertility * MaleFertility,
     FemaleFertility ~ 1 + (1 | gr(Dam, cov = A_full)) + Year + (1 | fpe) + (1 | e_proxy),
     MaleFertility ~ 1 + boost + (1 | gr(Sire, cov = A_full)) + Year + (1 | mpe) + (1 | e_proxy),
     nl = TRUE
  ),

  family = zero_one_inflated_beta(),

  data = eggs_corrected,
  data2 = list(A_full = Ap),

  prior = c(
    prior(normal(0.5, 0.3), class = "b", nlpar = "FemaleFertility", lb = 0, ub = 1),
    prior(normal(0.5, 0.3), class = "b", nlpar = "MaleFertility", lb = 0, ub = 1),
    
    prior(normal(0, 1), class = "sd", nlpar = "FemaleFertility", lb = 0),
    prior(normal(0, 1), class = "sd", nlpar = "MaleFertility", lb = 0)
  ),
  cores = 3, chains = 1, iter = 500000, warmup = 200000, thin = 150,
  control = list(adapt_delta = 0.98, max_treedepth = 15)
)
```
**OR**
```
fertility_model <- brm(
  bf(
    perc_eye ~ FemaleFertility * (1 - exp(-Capacity / n_eggs)),
    
    FemaleFertility ~ 1 + (1 | gr(Dam, cov = A_full)) + Year + (1 | fpe) + (1 | e_proxy),   
    Capacity ~ 1 + boost + (1 | gr(Sire, cov = A_full)) + Year + (1 | mpe) + (1 | e_proxy), 
    #Survival ~ 1 + s(Inb) + (1 | e_proxy),
    
    nl = TRUE
  ),
  
  family = zero_one_inflated_beta(),
  
  data = eggs_corrected,
  data2 = list(A_full = Ap),
  
  prior = c(
    prior(normal(0.5, 0.3), class = "b", nlpar = "FemaleFertility", lb = 0, ub = 1),
    prior(normal(0, 1), class = "sd", nlpar = "FemaleFertility", lb = 0),
    
    # Capacity prior
    prior(normal(12000, 6000), class = "b", nlpar = "Capacity", lb = 0),
    prior(normal(0, 4000), class = "sd", nlpar = "Capacity", lb = 0)
  ),
  
  cores = 3, chains = 1, iter = 500000, warmup = 200000, thin = 150,
  control = list(adapt_delta = 0.98, max_treedepth = 15)
)
```

![Untitled (11)](https://github.com/user-attachments/assets/f0f92cf3-4865-4af9-a1d6-276b492bd190)


**Current struggles in preliminary analysis of cross data**

**struggle 1:** missing parents. When filtering out totally invalid entries the data set consists of 624 crosses where at least one of the parents is valid (correspond to pedigree entries). When strict, the data-set consists of 546 crosses where both breeders are valid.

**possible approaches I have tried:** 1) work with 546 clean records. 2) work with the 624 data-points by assigning phantoms (cannot assess the effect of this). 
**final solution:** I took option 3 XD. Only kept phantom dams. My line of thought was that I can probably afford given the (for the most part) 2 dams to 1 sire breeding scheme.


**struggle 2:** priors for fertility parameters. Truncated normal seems promissing but also other options such as beta and normal (inverse logit link) have been tested.

**struggle 3:** year as fixed or random effect. 8 years, including 2 where the boost solution was used (strong intercept for male fertrility in all models). Treating it as fixed will probably reduce complexity and overparameterization.

**struggle 4:** for 2025, number of spawned eggs is missing for some of the failed families.

**solution 4:** use annual mean or "impute" from e.g. female length, weight, Kf, Year. Should not be very critical either way for 0 observations.

<br>
<br>

## Part 2: sperm analysis records (363 from 2020 and 720 from 2024). Selected variables to analyze are log(concentration + 1), curvilinear velocity and straightness of path.
```
sperm_trait ~ year(fixed) + day_of_sampling(covariate) + animal + residual
```
![image](https://github.com/user-attachments/assets/fe3bdef3-1418-4da1-a529-ad0ba95e100e)


**Current struggles in preliminary analysis of CASA data**

**struggle 1:** --While single-step yields more or less expected variance components for sperm traits,
VCE using only additive relationships from pedigree yields very high heritability estimates (~0.6-0.95) and correlations in the case of multitrait analysis including both years or only 2024. 
In a univariate fashion that pattern breaks when including data from both years but persists in the case of 2024.--

**solution 1:** Year should be coded as cross alpha and not cross numer.


**struggle 2:** CASA values from 2020 seem a bit odd... there are no 0s at all for VCL but a lot of decimal values close to it. From last 2 years I know that 0s are in fact not uncommon at all... Cross-check.

**struggle 3:** Martin mentioned that but have not given it a lot of thought. The boost solution introducices a scale problem or dispersion as well?

**struggle 4:** sperm velocities for samples with concentration bellow 200 were not measured. Should I assign 0 or have them missing?

CP: Better missing.
FP: Ok, I agree.

**struggle - fun fact 5:** Presence of replicates. Apparently we have recaptured a few males more than once, likely because the stuff accidentally returned them to the sampling tank. I should probably only keep the first occurance (?)

CP: Probably the first measurement. Would be interesting to see whether recordings vary a lot.

FP: Would repeated measures enhance power? Shall we include an autoregressive term in this case?

CP: You can try if you want. But I don't will help much in our case.

![image](https://github.com/user-attachments/assets/8db04171-6d9c-470f-8377-07fb8c20e851)

![Obs_vs_latent_EBVs](https://github.com/user-attachments/assets/327f4641-de9f-4883-a70a-00af1c6f942e)


