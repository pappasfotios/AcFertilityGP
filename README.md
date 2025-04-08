# AcFertilityGP
Genetic parameters for fertility traits of Arctic (superior) charr


## Part 1: data from cross success over multiple genarations (full-records, egg count and eyed embryos count)
```
library(brms)

fertility_model <- brm(
  # Response model for the proportion of eyed eggs (Id added as residual term)
  bf(perc_eye ~ FemaleFertility * MaleFertility,
     FemaleFertility ~ 1 + (1|gr(Dam, cov = A_full)) + (1|Year) + (1|Id),
     MaleFertility ~ 1 + (1|gr(Sire, cov = A_full)) + Year + (1|Id),
     # could also add an embryo survival coefficient modelled as a non-liner function of inbreeding
     nl = TRUE
  ),
  
  family = zero_one_inflated_beta(),
  
  data = eggs_corrected,
  data2 = list(A_full = A_full),
  
  prior = c(
    prior(normal(0.5, 0.3), class = "b", nlpar = "FemaleFertility", lb = 0, ub = 1),
    prior(normal(0.5, 0.3), class = "b", nlpar = "MaleFertility", lb = 0, ub = 1),
    prior(exponential(1), class = "sd", nlpar = "FemaleFertility"),
    prior(exponential(1), class = "sd", nlpar = "MaleFertility")
  ),
  
  cores = 6, chains = 3, iter = 60000, warmup = 20000,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
```
**OR**
```
fertility_model <- brm(
  bf(
    perc_eye ~ FemaleFertility * (1 - exp(-Capacity / n_eggs)),
    
    FemaleFertility ~ 1 + (1|gr(Dam, cov = A_full)) + (1|Year) + (1|Id),
    Capacity ~ 1 + boost + (1|gr(Sire, cov = A_full)) + (1|Year) + (1|Id),
    
    nl = TRUE
  ),
  
  family = zero_one_inflated_beta(),
  
  data = eggs_corrected,
  data2 = list(A_full = A_full),
  
  prior = c(
    prior(normal(0.5, 0.3), class = "b", nlpar = "FemaleFertility", lb = 0, ub = 1),
    prior(exponential(1), class = "sd", nlpar = "FemaleFertility"),
    
    prior(normal(2000, 1000), class = "b", nlpar = "Capacity", lb = 0),
    prior(exponential(0.01), class = "sd", nlpar = "Capacity")  # weak prior to allow variability
  ),
  
  cores = 6, chains = 3, iter = 60000, warmup = 20000,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
```

![Untitled (11)](https://github.com/user-attachments/assets/f0f92cf3-4865-4af9-a1d6-276b492bd190)


**Current struggles in preliminary analysis of cross data**

**struggle 1:** missing parents. When filtering out totally invalid entries the data set consists of 624 crosses where at least one of the parents is valid (correspond to pedigree entries). When strict, the data-set consists of 546 crosses where both breeders are valid.

**possible approaches I have tried:** 1) work with 546 clean records. 2) work with the 624 data-points by assigning phantoms (cannot assess the effect of this). 


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
