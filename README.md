# Spline-based AFT Model
R code for implementing the method in article “Pang M, Platt R, Schuster T, Abrahamowicz M. Spline-based Accelerated Failure Time Model. Under review at *Statistics in Medicine* 2020".

## Content
### Description
A semi-parametric accelerated failure time (AFT) model that models the baseline hazard using regression B-splines.

The code has been written using R with the following version information:<br/>
- R version 3.6.3 (2016-06-21)<br/> 
- Platform x86_64-apple-darwin15.6.0 (64-bit)<br/> 
- Using R packages:<br/> 
  - survival version 3.1-12
  - splines version v3.6.2
  
#### Code to implement the spline-based AFT model:
##### `SplineAFT.R`
This program includes all necessary functions to provide estimates of:
- covariate effects (time ratios)
- Hazard function and survival curve conditional on an arbitrary covariate pattern

The program is called by the program `Example.R` and `Diagnostics.R`. 

#### Code to run the spline-based AFT model:
##### `Example.R`
The program uses the dataset `colon` in the R `survival` package, and generates the results save in  `colon_splineaft.rda`.

#### Code to check the AFT and PH assumptions:
##### `Diagnostics.R`
The program loads `colon_splineaft.rda`.

For questions or comments about the code please contact Menglan Pang (menglan.pang at mail.mcgill.ca).
