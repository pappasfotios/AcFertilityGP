# AcFertilityGP
Genetic parameters for fertility traits of Arctic (superior) charr


**Data structure**

-546 cross records (full-records, egg count and eyed embryos count)

-1083 sperm analysis records (363 from 2020 and 720 from 2024) : selected variables to analyze are log(concentration), curvilinear velocity and straightness (only for 2024).


**current struggles in preliminary analysis**

While single-step yields more or less expected variance components for sperm traits,
VCE using only additive relationships from pedigree yields very high heritability estimates (~0.6-0.95) and correlations in the case of multitrait analysis including both years or only 2024. 
In a univariate fashion that pattern breaks when including data from both years but persists in the case of 2024.

Example output:
```
h2  - Function: g_3_3_1_1/(g_3_3_1_1+r_1_1)
  Mean:   0.69330    
  Sample Mean:   0.69119    
  Sample SD:   0.51365E-01
  
h2  - Function: g_3_3_2_2/(g_3_3_2_2+r_2_2)
  Mean:   0.83571    
  Sample Mean:   0.83410    
  Sample SD:   0.36757E-01
  
h2  - Function: g_3_3_3_3/(g_3_3_3_3+r_3_3)
  Mean:   0.96017    
  Sample Mean:   0.96008    
  Sample SD:   0.13189E-01
  
rg12  - Function: g_3_3_1_2/((g_3_3_1_1*g_3_3_2_2)^0.5)
  Mean:   0.91496    
  Sample Mean:   0.91613    
  Sample SD:   0.24162E-01
  
rg13  - Function: g_3_3_1_3/((g_3_3_1_1*g_3_3_3_3)^0.5)
  Mean:   0.93371    
  Sample Mean:   0.93478    
  Sample SD:   0.19944E-01
  
rg23  - Function: g_3_3_2_3/((g_3_3_2_2*g_3_3_3_3)^0.5)
  Mean:   0.96322    
  Sample Mean:   0.96393    
  Sample SD:   0.10933E-01
```





