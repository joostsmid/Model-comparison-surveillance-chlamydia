# Comparison of mathematical models developed to estimate C. trachomatis prevalence and incidence using the available public health surveillance data about tests and diagnoses at national or sub-national level

This project contains the code for two models different models ("Lewis"<sup>1,2</sup> and "Ali"<sup>3</sup>) for the estimation of prevalence and incidence of C. trachomatis in England between 2000-2015 using English surveillance data about testing and diagnosis (used in [2]), and the code to compare the results of these models.
The original model code can be found on other github repositories: (https://github.com/joanna-lewis/ct_trends) for the "Lewis" model and (https://github.com/leftygray/Chlamydia_incidence_model) for the "Ali" model.
Some modifications have been made to the original model code. In particular, the "Lewis" model, originally written in Python, was put into **R** language. The "Ali" model, originaly developed to describe prevalence and incidence of C. trachomatis in Australia in 6 age cohorts, was modified to describe prevalence and incidence of C. trachomatis in England in 2 age cohorts, instead.
The data used come from Chandra _et al._ (2017)<sup>4</sup> providing maximum and minimum estimates for the number of tests and diagnoses in England from 2000-2012.

### Code organization ###
All the model scripts, input data and outputs are stored in the main directory and 3 sub-directories. The main directory also contains files to run the code producing prevalence and incidence figures estimated by the two models (in one plot).

_Main directory files_

The following two R scripts run both models. They generate plots describing prevalence and incidence between 2000-2015 estimated by the two models. These plots are stored in the `figures/` directory.

- `compare_Lewis_Ali_maxdata.R:` Uses the data on maximum number of tests and diagnoses.

- `compare_Lewis_Ali_mindata.R:` Uses the data on minimum number of tests and diagnoses.

_Lewis_

Contains all the specific R functions and data used in the "Lewis" model. This folder contains the file `run_Lewis.R` which is used for the model comparison. Alternatively, one can run the `LewisWhite_Lancet.Rmd` file which produced a pdf file with additional figures for the "Lewis" model only. Most of these figures were originally produced for the original research article.

_Ali_

Contains all the specific R functions and data used in the "Ali" model. Not all functions should be ran. I produced the Bayesian model outputs for the English data already as the files `Ali/output/posterior_maxdata.R` and `Ali/output/posterior_mindata.R`. One can also produce these outputs oneself by uncommenting lines 25-27 in `compare_Lewis_Ali_maxdata.R` and `compare_Lewis_Ali_mindata.R`.

_figures_

Plots describing prevalence and incidence in England between 2000-2015 estimated by the two models are stored in the `figures/` directory.

### References ###

1.	Lewis, J. et al. (2017) [Estimating Local Chlamydia Incidence and Prevalence Using Surveillance Data](http://dx.doi:10.1097/Ede.0000000000000655). Epidemiology 28, 492-502
2.	Lewis, J. et al. (2018) [Changes in chlamydia prevalence and duration of infection estimated from testing and diagnosis rates in England: a model-based analysis using surveillance data](http://dx.doi:10.1016/S2468-2667(18)30071-9) Lancet Public Health 3, E271-E278
3.	Ali, H. et al. (2015) [A new approach to estimating trends in chlamydia incidence](http://dx.doi:10.1136/sextrans-2014-051631) Sexually Transmitted Infections 91, 513-519
4.	Chandra, N. L. et al. (2017) [Filling in the gaps: estimating numbers of chlamydia tests and diagnoses by age group and sex before and during the implementation of the English National Screening Programme, 2000 to 2012](http://dx.doi:10.2807/1560-7917.ES.2017.22.5.30453) Euro Surveill 22
