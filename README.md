# Heterogeneity in subnational mortality in the context of the COVID-19 pandemic: The case of Belgian districts
Benjamin-Samuel Schlüter, Bruno Masquelier, Carlo Giovanni Camarda, January 2022

Paper submitted to Archives of Public Health (BMC) scientific journal. Pre-print [link](https://assets.researchsquare.com/files/rs-1160099/v1_covered.pdf?c=1639757516)



## Purpose

This file and all other files in the repository should enable you to reproduce our results, plots, tables, ... from our joint paper. Our paper assesses how the mortality of the Belgian districts has been impacted in the context of the Covid19  pandemic. We assess how the heterogeneity in mortality levels within Belgium has evolved in 2020 relatively to previous years. In addition, we propose a methodology to measure the shock on mortality at the district level, accounting for different sources of uncertainty.

## Directory structure
There should be a single parent directory on your file system, where you have to set your R working directory and which will be divided into sub-directories named:

/code

&nbsp;&nbsp;&nbsp;&nbsp; /functions

&nbsp;&nbsp;&nbsp;&nbsp;/stan

/data

&nbsp;&nbsp;&nbsp;&nbsp;/raw

&nbsp;&nbsp;&nbsp;&nbsp;/tidy

&nbsp;&nbsp;&nbsp;&nbsp;/estimates

The repository available on Github where this README stands is the repository containing these folders. Make sure to exactly reproduce this structure on your laptop before running the R scripts. 

## Data

The analyses use Belgian register data and death certificates. Legally, I cannot share the data publicly. I am however evaluating the possibility to put a random sample of the anonymized data. You also have open data freely available [here](https://statbel.fgov.be/en/open-data?category=50) but the age groups available are much wider than in this study.

The data set used for these analyses is called `df_ageyearsexadmins_extrapol2020exp` located in the ./data/tidy folder. It has been created with the script `extract_data_covidsubBE.R` which links the register data with death certificates and manipulates the data to obtain a data frame of death count and exposure by age, sex, year region, province and district.  For the year 2021, we use population data freely available [here](https://statbel.fgov.be/en/open-data?category=23).

## R scripts

Script `smr_subnational.R` (within /code) in combination with script `smr_subnational.stan` (inside /code/stan), allows to assess the subnational heterogeneity in mortality over the period 2015-2020. The metrics used are standardized mortality ratios. Figures about this part of the analysis are produced by the script.

Script `e0_nat_projection_LC.R` (within /code) estimates the Lee-Carter model at the national level in a Bayesian context using `lee_carter.stan` (Stan code written by Monica Alexander available on her Github). It stores the national posterior life expectancy in 2020 that will be used to construct subnational counterfactual scenario in 2020.

Script `e0_subnat_projection` (within /code) in combination with `diff_e0_subnat_trend.stan` measures the shock on mortality at the district level during the year 2020. As a shock measure, we use the life expectancy at birth. The developed methodology tries to assess, as much as possible, for all uncertainty sources.

Finally, `maps_creation.R` creates maps with outputs stored while running `e0_subnat_projection`.  

