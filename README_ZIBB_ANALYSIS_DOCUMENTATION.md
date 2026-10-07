# Squamata ZIBB-PGLMM Analysis Documentation: 

## Working analysis record

This document summarizes the current status of the Squamata comparative oncology analysis.

It is intended to document:

- where the project started
- why a Zero-Inflated Beta-Binomial model was added
- how the new model is structured
- what models have been run
- the current results
- model diagnostics
- limitations
- unresolved questions
- reasonable next steps

### Like I stated above, this is just a working record for needed documentation. Interpretation to come.

---

<img width="1080" height="1080" alt="image" src="https://github.com/user-attachments/assets/3107e973-08df-4a1b-9d51-932adc564492" />

## Table of Contents

1. [Project background](#1-project-background)
2. [Why add a ZIBB model?](#2-why-add-a-zibb-model)
3. [Current response definitions](#3-current-response-definitions)
4. [Main analysis script](#4-main-analysis-script)
5. [Input files](#5-input-files)
6. [Data cleaning](#6-data-cleaning)
7. [Phylogenetic structure](#7-phylogenetic-structure)
8. [Current life-history predictors](#8-current-life-history-predictors)
9. [Why beta-binomial?](#9-why-beta-binomial)
10. [Why zero inflation?](#10-why-zero-inflation)
11. [Current model structure](#11-current-model-structure)
12. [Bayesian sampling settings](#12-bayesian-sampling-settings)
13. [Current priors](#13-current-priors)
14. [Models originally defined](#14-models-originally-defined)
15. [Minimum sample-size safeguard](#15-minimum-sample-size-safeguard)
16. [Multivariable sample sizes](#16-multivariable-sample-sizes)
17. [Successfully fitted models](#17-successfully-fitted-models)
18. [Current summary results](#18-current-summary-results)
19. [Neoplasia interpretation](#19-neoplasia-interpretation)
20. [Malignancy results](#20-malignancy-results)
21. [Malignancy interpretation](#21-malignancy-interpretation)
22. [Overall current biological result](#22-overall-current-biological-result)
23. [Relationship to the previous PGLS results](#23-relationship-to-the-previous-pgls-results)
24. [Bayesian interpretation](#24-bayesian-interpretation)
25. [MCMC convergence](#25-mcmc-convergence)
26. [Posterior predictive checks](#26-posterior-predictive-checks)
27. [Understanding the zero-count plots](#27-understanding-the-zero-count-plots)
28. [This is not a convergence problem](#28-this-is-not-a-convergence-problem)
29. [Detailed example: Neoplasia + longevity](#29-detailed-example-neoplasia--longevity)
30. [Neoplasia + longevity main effect](#30-neoplasia--longevity-main-effect)
31. [Zero-inflation intercept](#31-zero-inflation-intercept)
32. [Sampling effort and zero inflation](#32-sampling-effort-and-zero-inflation)
33. [Beta-binomial dispersion](#33-beta-binomial-dispersion)
34. [Phylogenetic random effect](#34-phylogenetic-random-effect)
35. [Important question raised by the diagnostics](#35-important-question-raised-by-the-diagnostics)
36. [Current model should not yet be treated as final](#36-current-model-should-not-yet-be-treated-as-final)
37. [Recommended next sensitivity analysis](#37-recommended-next-sensitivity-analysis)
38. [Why start with Neoplasia + longevity?](#38-why-start-with-neoplasia--longevity)
39. [Model comparison caution](#39-model-comparison-caution)
40. [Multivariable modeling limitation](#40-multivariable-modeling-limitation)
41. [Should the N = 18 model eventually be run?](#41-should-the-n--18-model-eventually-be-run)
42. [Current interpretation compared with manuscript claims](#42-current-interpretation-compared-with-manuscript-claims)
43. [Current strongest statement](#43-current-strongest-statement)
44. [Statements that are NOT currently justified](#44-statements-that-are-not-currently-justified)
45. [Current output directory](#45-current-output-directory)
46. [Primary summary file](#46-primary-summary-file)
47. [Model-specific datasets](#47-model-specific-datasets)
48. [Model summaries](#48-model-summaries)
49. [Diagnostics files](#49-diagnostics-files)
50. [Effect plots](#50-effect-plots)
51. [Posterior predictive plots](#51-posterior-predictive-plots)
52. [Zero-count posterior predictive plots](#52-zero-count-posterior-predictive-plots)
53. [Other output files](#53-other-output-files)
54. [Fitted model objects](#54-fitted-model-objects)
55. [Reproducibility](#55-reproducibility)
56. [GitHub organization](#56-github-organization)
57. [General workflow](#57-general-workflow)
58. [Where the project started](#58-where-the-project-started)
59. [Where the project is now](#59-where-the-project-is-now)
60. [Where to next?](#60-where-to-next)
61. [Questions still unresolved](#61-questions-still-unresolved)
62. [Current conclusions](#62-current-conclusions)
63. [Short summary — Gucci but not as Gucci as one could be](#63-short-summary--gucci-but-not-as-gucci-as-one-could-be)
64. [Current status](#64-current-status)
# 1. Project background

The Squamata project examines cancer prevalence across lizards and snakes in relation to life-history traits.

The existing project primarily used **Phylogenetic Generalized Least Squares (PGLS)** models to examine associations between cancer prevalence and traits such as:

- maximum longevity
- gestation length
- adult body mass
- litter size
- age at maturity
- other reproductive and life-history variables

The two primary cancer outcomes are:

- **Neoplasia**
- **Malignancy**

The existing analysis also includes:

- phylogenetic visualization
- Pagel's lambda
- univariate PGLS models
- multivariable PGLS models
- ordinary linear models
- Variance Inflation Factor (VIF) checks
- interaction testing
- stepwise AIC model selection

The original analyses generally modeled already-calculated prevalence values.

For example:

```r
NeoplasiaPrevalence ~ predictor
```

or:

```r
MalignancyPrevalence ~ predictor
```

---

# 2. Why add a ZIBB model?

The cancer data are fundamentally **count data**.

For each species, the dataset contains a number of tumor cases observed among a known number of necropsy records.

For example:

```text
3 tumor cases out of 40 records
```

This produces a prevalence of:

```text
3 / 40 = 0.075
```

However, modeling prevalence alone does not directly preserve the amount of sampling behind each prevalence estimate.

For example:

```text
2 / 20 = 0.10
```

and:

```text
20 / 200 = 0.10
```

have the same prevalence, but the second estimate is based on considerably more information.

For this reason, the new analysis models the **number of tumor cases out of the number of available records**.

The selected framework is a:

# Zero-Inflated Beta-Binomial Phylogenetic Generalized Linear Mixed Model

or:

# ZIBB-PGLMM

This framework allows the analysis to account for:

1. Different sampling effort among species
2. Tumor counts rather than prevalence alone
3. Overdispersion
4. Excess zero counts
5. Shared evolutionary history
6. Life-history predictors

---

# 3. Current response definitions

The current ZIBB analysis uses:

```r
NeoplasiaCases = NeoplasiaWithDenominators
MalignancyCases = Malignant
Trials = RecordsWithDenominators
```

Therefore, the response is modeled using:

```r
cases | trials(Trials)
```

rather than:

```r
NeoplasiaPrevalence
```

or:

```r
MalignancyPrevalence
```

alone.

---

# 4. Main analysis script

The new analysis is contained in:

```text
6. ZIBB_PGLMM_Squamata.R
```

The existing PGLS scripts were INTENTIONALLY LEFT UNCHANGED.

This ZIBB analysis was added as a separate workflow so that the original analysis remains available for comparison for backtracking and because I literally didnt have access to the repo xD

---

# 5. Input files

The current ZIBB pipeline uses:

```text
squamata_data.csv
min20Fixed516.nwk
```

The complete source dataset is also available as:

```text
min20-2022.05.16.csv
```
see -> https://aacrjournals.org/cancerdiscovery/article/15/1/227/750844/Cancer-Prevalence-across-VertebratesCancer-across 

The script prefers the existing Squamata-only dataset.

If that file is unavailable, the complete dataset can be filtered to Squamata.

---

# 6. Data cleaning

Negative numerical values used as missing-value codes are converted to `NA`.

Species names are standardized by:

- trimming whitespace
- replacing spaces with underscores

For example:

```text
Pantherophis guttatus
```

becomes:

```text
Pantherophis_guttatus
```

The analysis also checks for:

- duplicated species
- missing species names
- species absent from the phylogeny
- tumor case counts larger than denominators
- non-integer tumor counts
- non-integer denominators
- missing predictor values

Species missing from the phylogenetic tree are documented separately.

---

# 7. Phylogenetic structure

Species are matched to:

```text
min20Fixed516.nwk
```

For each fitted model the phylogenetic tree is restricted to the species actually included in that model.

The phylogenetic covariance/correlation matrix is generated with:

```r
A_model <- ape::vcv(tree_model, corr = TRUE)
```

The ZIBB model then includes:

```r
(1 | gr(Species, cov = A_model))
```

This allows related species to share covariance due to common evolutionary history.

Closely related species are therefore not automatically treated as completely independent observations.

---

# 8. Current life-history predictors

The current pipeline examines:

- Maximum longevity
- Gestation length
- Adult body mass
- Litter size

The relevant standardized variables are:

```text
longevity_s
gestation_s
mass_s
litter_size_s
```

Positive life-history variables are log10 transformed before standardization.

For example:

```r
log_longevity = log10(max_longevity_months)
longevity_s = scale(log_longevity)
```

and:

```r
log_mass = log10(adult_weight_g)
mass_s = scale(log_mass)
```

Standardization means that coefficients describe approximately a **one-standard-deviation change in the log-transformed predictor**.

This is important for interpretation.

For example, an odds ratio of:

```text
1.53
```

does **not** mean that one additional month of longevity increases the odds by 53%.

It refers to approximately a one-standard-deviation increase in standardized log longevity.

---

# 9. Why beta-binomial?

A standard binomial model assumes that variability is determined entirely by:

- the underlying probability
- the number of trials

Biological comparative data often contain more variability than this assumption allows.

Species may differ for reasons not completely represented by the measured predictors, including:

- genetics
- environment
- husbandry
- diagnostic practices
- institution
- demographic structure
- unmeasured life-history traits

A **beta-binomial** distribution allows additional variation beyond a simple binomial model.

This is commonly described as **overdispersion**.

---

# 10. Why zero inflation?

Some species have zero observed tumor cases.

A standard count model may generate zeros naturally, but some datasets can contain more zeros than expected.

A zero-inflated model adds an additional process that allows observations to enter an extra-zero state.

The current model allows this additional zero process to depend on sampling effort.

Sampling effort is represented as:

```r
log_trials_s = scale(log(Trials))
```

The zero-inflation model is:

```r
zi ~ log_trials_s
```

The purpose is to examine whether the probability of belonging to the additional-zero component changes with the amount of available sampling.

---

# 11. Current model structure

A typical current model is:

```r
cases | trials(Trials) ~ longevity_s +
  (1 | gr(Species, cov = A_model))

zi ~ log_trials_s
```

using:

```r
family = zero_inflated_beta_binomial()
```

The models are fitted in `brms`.

---

# 12. Bayesian sampling settings

The current pipeline uses:

```r
chains = 4
iter = 4000
warmup = 2000
```

with:

```r
adapt_delta = 0.99
max_treedepth = 15
```

and a fixed random seed for reproducibility.

---

# 13. Current priors

The analysis uses regularizing priors for:

- intercepts
- regression coefficients
- phylogenetic random-effect standard deviation
- zero-inflation parameters
- beta-binomial dispersion

The general goal is to avoid unnecessarily extreme parameter estimates while still allowing the data to inform the posterior.

---

# 14. Models originally defined

The pipeline was designed to fit the following models.

## Univariate models

```text
longevity_only
gestation_only
mass_only
litter_size_only
```

## Multivariable models

```text
longevity + gestation
longevity + gestation + mass
full life-history model
```

Each model was intended to be evaluated separately for:

```text
Neoplasia
Malignancy
```

---

# 15. Minimum sample-size safeguard

The script currently contains:

```r
MIN_SPECIES = 20
```

Models containing fewer than 20 usable species are automatically skipped 
The current modeling as of Oct 3rd, 2026 is as shown and further explained in the following sections.


<img width="627" height="227" alt="Screenshot 2026-10-03 at 9 11 34 PM" src="https://github.com/user-attachments/assets/9c59d25e-e0d4-49f1-8e0e-065930df7bd4" />


This was included because the ZIBB-PGLMM is relatively parameter-rich.

The model may contain:

- one or more regression coefficients
- a phylogenetic random effect
- beta-binomial dispersion
- a zero-inflation intercept
- a zero-inflation sampling-effort coefficient

Trying to estimate all of these components with very small sample sizes can produce unstable or poorly identified models.

---

# 16. Multivariable sample sizes

The attempted multivariable models produced the following usable sample sizes:

| Model | Number of species |
|---|---:|
| Longevity + gestation | 18 |
| Longevity + gestation + mass | 14 |
| Full life-history model | 12 |

Because all were below:

```text
MIN_SPECIES = 20
```

these models were automatically skipped.

This occurred for both:

- Neoplasia
- Malignancy

These are not model failures.

They are intentional safeguards.

---

# 17. Successfully fitted models

Eight univariate ZIBB-PGLMMs were successfully fitted.

## Neoplasia

- Longevity
- Gestation
- Adult body mass
- Litter size

## Malignancy

- Longevity
- Gestation
- Adult body mass
- Litter size

---

<img width="1313" height="171" alt="Screenshot 2026-10-03 at 8 37 28 PM" src="https://github.com/user-attachments/assets/aa6aa6fb-fc43-4a68-bc63-67ffffc6a35c" />


# 18. Current summary results

## Neoplasia

| Predictor | N | Median beta | 95% CrI | Median OR | 95% OR CrI | Clear effect? |
|---|---:|---:|---|---:|---|---|
| Longevity | 49 | 0.423 | 0.058 to 0.833 | 1.53 | 1.06 to 2.30 | Yes |
| Gestation | 20 | 0.072 | -0.411 to 0.574 | 1.07 | 0.66 to 1.78 | No |
| Body mass | 42 | 0.152 | -0.236 to 0.503 | 1.16 | 0.79 to 1.65 | No |
| Litter size | 35 | -0.118 | -0.502 to 0.237 | 0.89 | 0.61 to 1.27 | No |

---

# 19. Neoplasia interpretation

## Longevity

The longevity model used:

```text
N = 49 species
```

Posterior estimate:

```text
Median beta = 0.423
95% credible interval = 0.058 to 0.833
```

Odds ratio:

```text
Median OR = 1.53
95% credible interval = 1.06 to 2.30
```

The 95% credible interval does not include zero on the coefficient scale.

Current interpretation:

> Longer-lived Squamata species show higher neoplasia odds in the current ZIBB phylogenetic model.

Because longevity was standardized after log transformation, the odds ratio corresponds approximately to a one-standard-deviation increase in log longevity.

---

## Gestation

The gestation model used:

```text
N = 20 species
```

Posterior estimate:

```text
Median beta = 0.072
95% credible interval = -0.411 to 0.574
```

Odds ratio:

```text
Median OR = 1.07
95% credible interval = 0.66 to 1.78
```

The credible interval includes zero.

Current interpretation:

> There is currently no clear evidence for an association between gestation length and neoplasia in the ZIBB model.

The relatively small sample size should also be considered when interpreting the wide interval.

---

## Adult body mass

The body-mass model used:

```text
N = 42 species
```

Posterior estimate:

```text
Median beta = 0.152
95% credible interval = -0.236 to 0.503
```

Odds ratio:

```text
Median OR = 1.16
95% credible interval = 0.79 to 1.65
```

Current interpretation:

> There is currently no clear association between adult body mass and neoplasia in the univariate ZIBB model.

---

## Litter size

The litter-size model used:

```text
N = 35 species
```

Posterior estimate:

```text
Median beta = -0.118
95% credible interval = -0.502 to 0.237
```

Odds ratio:

```text
Median OR = 0.89
95% credible interval = 0.61 to 1.27
```

Current interpretation:

> There is currently no clear association between litter size and neoplasia.

---

# 20. Malignancy results

| Predictor | N | Median beta | 95% CrI | Median OR | 95% OR CrI | Clear effect? |
|---|---:|---:|---|---:|---|---|
| Longevity | 49 | 0.367 | 0.003 to 0.789 | 1.44 | ~1.00 to 2.20 | Borderline positive |
| Gestation | 20 | 0.054 | -0.445 to 0.552 | 1.05 | 0.64 to 1.74 | No |
| Body mass | 42 | 0.039 | -0.394 to 0.467 | 1.04 | 0.67 to 1.59 | No |
| Litter size | 35 | -0.202 | -0.654 to 0.188 | 0.82 | 0.52 to 1.21 | No |

---

# 21. Malignancy interpretation

## Longevity

The longevity model used:

```text
N = 49 species
```

Posterior estimate:

```text
Median beta = 0.367
95% credible interval = 0.003 to 0.789
```

Odds ratio:

```text
Median OR = 1.44
```

The lower bound is only slightly above zero on the coefficient scale.

Current interpretation:

> Longer-lived Squamata species appear to have higher malignancy odds, although the evidence is weaker and more borderline than the neoplasia result.

This result should therefore be described more cautiously.

---

## Gestation

Posterior estimate:

```text
Median beta = 0.054
95% credible interval = -0.445 to 0.552
```

Current interpretation:

> No clear association between gestation and malignancy is currently supported.

---

## Adult body mass

Posterior estimate:

```text
Median beta = 0.039
95% credible interval = -0.394 to 0.467
```

Current interpretation:

> No clear association between adult body mass and malignancy is currently supported.

---

## Litter size

Posterior estimate:

```text
Median beta = -0.202
95% credible interval = -0.654 to 0.188
```

Current interpretation:

> No clear association between litter size and malignancy is currently supported.

---

# 22. Overall current biological result

At the present stage, **maximum longevity is the most consistent predictor**.

Longevity shows a positive association with:

- Neoplasia
- Malignancy

The neoplasia result is more clearly supported.

The malignancy result is positive but comparatively borderline.

The other tested predictors currently have 95% credible intervals that include zero:

- Gestation
- Adult body mass
- Litter size

Therefore, the present univariate ZIBB results support longevity more consistently than the other examined life-history traits.

---

# 23. Relationship to the previous PGLS results

Previous analyses using PGLS suggested associations involving life-history traits including longevity and gestation.

The ZIBB analysis does not reproduce every previous result.

The most important current difference is:

```text
Longevity remains supported.
Gestation does not currently show a clear ZIBB effect.
```

This does not automatically mean the previous PGLS analysis was incorrect.

The two approaches ask related questions using different assumptions.

The ZIBB analysis:

- models tumor counts directly
- includes denominators directly
- allows beta-binomial overdispersion
- contains a zero-inflation component
- includes phylogenetic covariance

The PGLS analysis modeled calculated prevalence as a continuous response.

Differences between the results should therefore be examined rather than treated as errors.

---

# 24. Bayesian interpretation

These models are Bayesian.

Therefore, the primary inference should be based on:

- posterior estimates
- 95% credible intervals
- posterior probabilities
- posterior predictive checks

The results should not be described using conventional frequentist language such as:

```text
p < 0.05
```

unless referring specifically to the older PGLS analyses.

For the ZIBB analysis, language such as the following is preferable:

> The 95% credible interval excluded zero.

or:

> The posterior supports a positive association.

or:

> The credible interval included zero, so the effect is not clearly supported.

---

# 25. MCMC convergence

The currently fitted models showed very good sampling diagnostics.

Maximum Rhat values were approximately:

```text
1.002
```

and the fitted models had:

```text
0 divergent transitions
```

Rhat values close to 1 indicate that the MCMC chains mixed and converged well.

The absence of divergent transitions is also encouraging.

Therefore:

> There is no obvious MCMC convergence problem in the currently fitted models.

However, good convergence does not automatically mean that the chosen probability model perfectly reproduces the biological data.

Model fit must also be evaluated.

---

# 26. Posterior predictive checks

The pipeline produces posterior predictive checks for each model.

These include:

```text
ppcheck_*.png
```

and:

```text
zero_check_*.png
```

The zero-count diagnostic specifically compares:

- the observed proportion of species with zero tumor cases
- the zero proportions expected in simulated datasets generated from the fitted posterior

---

# 27. Understanding the zero-count plots

In the zero-count posterior predictive plots:

```text
Red dashed line = observed zero-case fraction
Blue distribution = posterior-predicted zero-case fraction
```

Across several models, the red line occurs toward the left side of the posterior predictive distribution.

This means that:

> The current ZIBB models often generate somewhat more zero-case species than are actually present in the observed data.

This pattern is visible for both:

- Neoplasia
- Malignancy

across multiple life-history predictors.

---

# 28. This is not a convergence problem

The zero-count discrepancy should not be confused with MCMC failure.

The models showed:

```text
Rhat ≈ 1.00
0 divergences
```

Therefore, the sampler appears to have successfully estimated the specified model.

The posterior predictive diagnostic instead asks a different question:

> Does the fitted statistical model generate datasets that resemble the observed dataset?

The zero-count plots suggest that the model's treatment of zeros deserves additional examination.

---

# 29. Detailed example: Neoplasia + longevity

The full summary for the neoplasia-longevity model showed:

```text
Family: zero_inflated_beta_binomial

Formula:
cases | trials(Trials) ~ longevity_s +
  (1 | gr(Species, cov = A_model))

zi ~ log_trials_s
```

Number of observations:

```text
49
```

---

# 30. Neoplasia + longevity main effect

The model estimated:

```text
longevity_s = 0.43
95% credible interval = 0.06 to 0.83
```

This supports the positive longevity association seen in the summary table.

---

# 31. Zero-inflation intercept

For the same model:

```text
zi_Intercept = -2.93
95% credible interval = -4.63 to -1.65
```

Because the zero-inflation component uses a logit link, the intercept can be converted to a probability using:

```r
plogis(-2.93)
```

which is approximately:

```text
0.05
```

Therefore, at average standardized sampling effort, the model estimates approximately:

```text
5% extra structural-zero probability
```

The exact posterior distribution should be used for formal reporting, but this provides an intuitive interpretation.

---

# 32. Sampling effort and zero inflation

The model estimated:

```text
zi_log_trials_s = 0.16
95% credible interval = -1.09 to 1.34
```

The credible interval crosses zero broadly.

Current interpretation:

> There is no clear evidence in this model that sampling effort strongly predicts the additional zero-inflation process.

---

# 33. Beta-binomial dispersion

For the neoplasia-longevity model:

```text
phi = 7.54
95% credible interval = 4.25 to 12.05
```

The beta-binomial component allows extra species-to-species variation beyond an ordinary binomial model.

The current results indicate that this overdispersion component is an important part of the fitted distribution.

---

# 34. Phylogenetic random effect

For the same model, the estimated phylogenetic random-effect standard deviation was approximately:

```text
0.25
```

with a broad credible interval.

This suggests that some residual phylogenetic structure may remain after accounting for longevity.

The phylogenetic effect should be interpreted cautiously because the interval is relatively broad.

---

# 35. Important question raised by the diagnostics

The current ZIBB models:

- converge well
- estimate only modest explicit zero inflation in at least the longevity-neoplasia model
- still tend to generate somewhat more zero-case species than observed

This raises an important statistical question:

> Is the additional zero-inflation process actually necessary, or does the beta-binomial distribution already account adequately for the observed zeros?

This has not yet been formally tested.

---

# 36. Current model should not yet be treated as final

The ZIBB analysis was selected because the data contain:

- heterogeneous denominators
- overdispersion
- zero tumor counts

However, the presence of zeros alone does not prove that a zero-inflated model is necessary.

A beta-binomial model can naturally generate zeros.

Therefore:

> The current ZIBB model is a reasonable working model, but the need for explicit zero inflation should be evaluated before calling it the final model.

---

# 37. Recommended next sensitivity analysis

Useful next analysis could be to compare current ZIBB vs beta-binomial but given the amount of zeros and what Maley is asking for it can probably be forgone. Beta-Binomial was previously done for another comp-onc ACE project and we eventually stuck with ZIBB per Dr. Maleys guidance. However if it were to be done see the below:

## Model A — current ZIBB

```r
family = zero_inflated_beta_binomial()
```

against:

## Model B — beta-binomial without zero inflation

```r
family = beta_binomial()
```

The two models should be fit using:

- the exact same species
- the exact same tumor outcome
- the exact same predictor
- the exact same phylogeny
- the exact same covariance matrix
- comparable priors

The first comparison should probably use:

```text
Neoplasia + longevity
```

because this model currently has:

```text
N = 49
```

and shows the clearest life-history association.

---

# 38. Why start with Neoplasia + longevity?

The neoplasia-longevity model is currently the strongest test case because:

- it uses 49 species
- it converged well
- it had zero divergent transitions
- the longevity credible interval excluded zero
- the zero-inflation component was modest
- posterior predictive zero checks raised an interpretable model-fit question

IF a beta-binomial model produces similar longevity estimates while reproducing the observed zero frequency better, this would suggest that explicit zero inflation may not be necessary but once again per guidance beta-binomial can probably be forgone.

---

# 39. Model comparison caution

Formal model comparison should ONLY be performed when competing models use the **same observations**.

For example:

```text
longevity-only model: N = 49
gestation-only model: N = 20
```

These models should NOT automatically be ranked against each other using LOO or another predictive criterion because they were fit to different datasets.

A valid comparison should use the same species in both models.

This is especially important for future multivariable analysis!!

---

# 40. Multivariable modeling limitation

The largest current limitation (Oct 3rd 2026) for multivariable modeling is missing life-history data.

Current complete-case sample sizes fall rapidly:

```text
Longevity alone: 49
Gestation alone: 20
Longevity + gestation: 18
Longevity + gestation + mass: 14
Full life-history: 12
```

This means that adding predictors changes both:

- the statistical model
- the species being analyzed

Therefore, an apparent change in a coefficient could result from either:

1. controlling for another biological variable, or
2. analyzing a different subset of species

This should be kept in mind when future multivariable models are considered! 

---

# 41. Should the N = 18 model eventually be run?

Maybe, Possibly, Possibly so

The:

```text
longevity + gestation
```

model has:

```text
N = 18
```

This MAY eventually be worth fitting as an **exploratory or sensitivity analysis**.

However, it should NOT AUTOMATICALLY be treated as the primary model because:

- the sample size is small
- the model is relatively complex
- it contains multiple predictors
- it contains phylogenetic covariance
- it contains beta-binomial dispersion
- it contains a zero-inflation process

The current conservative choice is to keep it excluded from the main analysis UNTIL the simpler model structure is better understood.

---

# 42. Current interpretation compared with manuscript claims

The existing manuscript emphasized relationships involving longevity and gestation.

The current ZIBB results suggest:

```text
Longevity: supported
Gestation: not clearly supported
Body mass: not clearly supported
Litter size: not clearly supported
```

Therefore, manuscript language regarding gestation should be revisited once the ZIBB analysis is finalized. 
The current ZIBB results DO NOT clearly support a gestation effect, although the gestation models 
are based on a relatively small number of species and additional sensitivity analyses may still be informative.

The manuscript should NOT be revised until the remaining ZIBB diagnostics and planned sensitivity checks are completed. Regular Beta-Bi may be used as a sensitivity check if wanted.


---

# 43. Current strongest statement

The clearest statement supported by the current ZIBB results is:

> Maximum longevity is positively associated with neoplasia prevalence across the analyzed Squamata species after accounting for sampling effort, overdispersion, excess-zero structure, and phylogenetic relatedness.

A more cautious corresponding malignancy statement is:

> Maximum longevity also shows a positive association with malignancy, although the posterior support is weaker than for neoplasia.

---

# 44. Statements that are NOT currently justified

At the present stage, the analysis should not claim that:

```text
Zero inflation is definitely necessary.
```

or:

```text
Gestation is a robust predictor under the new model.
```

or:

```text
Body mass has no biological effect.
```

or:

```text
Litter size has no biological effect.
```

or:

```text
The ZIBB model proves longevity causes cancer.
```

The analysis identifies statistical associations, not causation.

A credible interval including zero also does not prove that the biological effect is exactly zero.

---

# 45. Current output directory

The new results are stored in:

```text
ZIBB_PGLMM_Squamata/
```

---

# 46. Primary summary file

The primary coefficient summary is:

```text
ZIBB_PGLMM_Squamata_summary.csv
```

https://github.com/WillWeston94/squamates/blob/main/ZIBB_PGLMM_Squamata/ZIBB_PGLMM_Squamata_summary.csv

<img width="1313" height="171" alt="Screenshot 2026-10-03 at 8 37 28 PM" src="https://github.com/user-attachments/assets/a41ef5e7-b593-42b2-ad64-922cccbdc35d" />


This contains information including:

- response
- model
- predictor
- sample size
- zero fraction
- median number of trials
- posterior coefficient interval
- odds ratio interval
- posterior probability of positive effect
- posterior probability of negative effect
- whether the 95% credible interval excludes zero
- beta-binomial `phi`
- maximum Rhat
- number of divergences

---

# 47. Model-specific datasets

Files such as:

```text
data_used_Neoplasia_longevity_only.csv
```

contain the exact species and variables used in each fitted model.

These are important because different predictors have different amounts of missing data.

Saving model-specific datasets makes it possible to identify exactly which species contributed to each analysis.

---

# 48. Model summaries

Files such as:

```text
summary_Neoplasia_longevity_only.txt
```

contain the full `brms` model summaries.

These include:

- regression coefficients
- credible intervals
- phylogenetic random-effect parameters
- zero-inflation coefficients
- beta-binomial dispersion
- Rhat
- effective sample size

---

# 49. Diagnostics files

Files such as:

```text
diagnostics_Neoplasia_longevity_only.csv
```

contain posterior diagnostic information for model parameters.

These can be used to inspect:

- Rhat
- bulk effective sample size
- tail effective sample size

---

# 50. Effect plots

Files such as:

```text
plot_Neoplasia_longevity_only_longevity_s.png
```

show the relationship between a predictor and the fitted outcome.

The plots include:

- observed data
- posterior prediction
- 95% Bayesian credible interval

Other standardized predictors, when present, are held at their mean.

Phylogenetic random effects are excluded from population-level prediction curves.

---

# 51. Posterior predictive plots

Files such as:

```text
ppcheck_Neoplasia_longevity_only.png
```

provide general posterior predictive diagnostics.

These should be inspected along with numerical convergence diagnostics.

---

# 52. Zero-count posterior predictive plots

Files such as:

```text
zero_check_Neoplasia_longevity_only.png
```

specifically evaluate whether the model reproduces the observed number of zero-case species.

These plots were added because zero inflation is one of the main reasons for considering the ZIBB family.

---

# 53. Other output files

The directory also contains:

```text
data_coverage.csv
tumor_count_audit.csv
species_missing_from_tree.csv
sessionInfo.txt
```

These provide additional reproducibility and data-quality information.

---

# 54. Fitted model objects

The fitted `brms` model objects are saved locally as:

```text
fit_*.rds
```

These objects can be large.

They are currently excluded from GitHub using:

```text
.gitignore
```

with a rule similar to:

```text
ZIBB_PGLMM_Squamata/*.rds
```

The models can be regenerated from the analysis script.

---

# 55. Reproducibility

The pipeline saves:

```text
sessionInfo.txt
```

This records the R environment and package versions used for the analysis.

The script also uses a fixed random seed.

Together with the input files, script, and model-specific datasets, this should make the analysis substantially easier to reproduce.

---

# 56. GitHub organization

The ZIBB analysis was developed as a separate addition to the existing Squamata project.

The intended repository structure is approximately:

```text
squamates/
│
├── 6. ZIBB_PGLMM_Squamata.R
├── README_ZIBB_ANALYSIS_DOCUMENTATION.md
├── squamata_data.csv
├── min20Fixed516.nwk
│
├── ZIBB_PGLMM_Squamata/
│   ├── ZIBB_PGLMM_Squamata_summary.csv
│   ├── data_coverage.csv
│   ├── tumor_count_audit.csv
│   ├── species_missing_from_tree.csv
│   ├── data_used_*.csv
│   ├── diagnostics_*.csv
│   ├── summary_*.txt
│   ├── plot_*.png
│   ├── ppcheck_*.png
│   ├── zero_check_*.png
│   └── sessionInfo.txt
│
└── existing PGLS scripts and project files
```

The existing scripts were not overwritten.

---

# 57. General workflow

The current overall analysis workflow can be summarized as:

```text
Original prevalence-based PGLS analysis
                ↓
Review existing data and scripts
                ↓
Identify tumor counts and denominators
                ↓
Build ZIBB-PGLMM pipeline
                ↓
Match species to phylogeny
                ↓
Fit univariate ZIBB models
                ↓
Check convergence
                ↓
Generate posterior effect estimates
                ↓
Run posterior predictive checks
                ↓
Inspect zero-count behavior
                ↓
Question whether explicit zero inflation is necessary
                ↓
Compare ZIBB against beta-binomial (Probably not needed according to what Dr. Maley wants but still useful to learn if wanted!)
                ↓
Reassess multivariable modeling
                ↓
Update manuscript interpretation
```

---

# 58. Where the project started

The project began with:

- prevalence-based PGLS
- Pagel's lambda
- phylogenetic tree visualization
- ordinary regression
- VIF analysis
- stepwise AIC selection
- manuscript conclusions focused particularly on longevity and gestation

These analyses provided the initial biological direction.

---

# 59. Where the project is now

The project now also has a functioning Bayesian phylogenetic count-modeling pipeline.

The current ZIBB analysis:

- models actual tumor case counts
- incorporates the number of available records
- allows overdispersion
- includes a zero-inflation component
- incorporates phylogenetic covariance
- standardizes life-history predictors
- saves model-specific datasets
- exports effect estimates
- exports odds ratios
- generates posterior predictive checks
- checks zero-count behavior
- checks Rhat
- checks effective sample size
- checks divergences
- saves reproducibility information

Eight univariate models have successfully completed.

---

# 60. Where to next? 

<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/109f0edd-943c-4349-a84d-723bea8c30f2" />

Potential next steps include:

1. Start with:
   ```text
   Neoplasia + longevity
   ```

2. Repeat the comparison for:
   ```text
   Malignancy + longevity
   ```

3. Examine whether:
   ```text
   zi ~ log_trials_s
   ```
   is actually useful.

4. Compare that formulation with:
   ```text
   zi ~ 1
   ```

5. Examine additional posterior predictive diagnostics.

6. Use formal predictive comparison only when models contain identical species.

7. Consider an exploratory:
   ```text
   longevity + gestation
   ```
   model with N = 18.

8. Avoid making the N = 12 full model a primary inferential model unless additional data become available.

9. Review whether additional life-history data can be obtained for currently incomplete species.

10. Compare the final count-model results directly with the existing PGLS findings.

11. Review manuscript wording after model selection is finalized.

---

# 61. Questions still unresolved

The following questions remain open:


### Does sampling effort explain extra zeros?

There is currently little evidence for this in the neoplasia-longevity model.

### Is longevity robust?

Currently, yes, particularly for neoplasia.

### Is the malignancy longevity effect robust?

It is positive, but weaker and more borderline.

### Is gestation robust?

Not in current (Oct 3rd 2026) univariate ZIBB analysis.

### Should very small multivariable models be fitted?

Possibly as sensitivity analyses, but not automatically as primary models.

---


# 62. Current conclusions

The current ZIBB-PGLMM pipeline is functioning well and producing stable posterior sampling.

The most consistent current biological finding is:

> Longer-lived Squamata species tend to show greater neoplasia prevalence and, more weakly, greater malignancy prevalence.

The currently tested univariate ZIBB models do **not** show similarly clear effects for:

- gestation
- adult body mass
- litter size

The current models also converged well, with:

- Rhat values near 1.00
- 0 divergent transitions in the fitted models

Posterior predictive checks suggest that some ZIBB models may generate somewhat more zero-case species than are observed in the actual data.

This does **not** indicate that the ZIBB models failed or that zero inflation should automatically be removed.

Instead, it suggests that the zero-inflation component should be examined more closely before the analysis is considered fully finalized.

Therefore:

> The current (Oct 3rd, 2026) ZIBB results are informative and biologically interpretable, but additional sensitivity and model-checking work may still improve the final analysis.

The primary analysis framework will remain the **phylogenetic Zero-Inflated Beta-Binomial model**.

Potential follow-up work includes:

- examining whether `zi ~ log_trials_s` is the most appropriate zero-inflation structure
- comparing it with a simpler `zi ~ 1` formulation
- examining posterior predictive fit in greater detail
- running selected lower-sample multivariable models as exploratory analyses
- optionally fitting a non-zero-inflated beta-binomial model as a sensitivity check

A plain beta-binomial model would therefore be used as a **robustness comparison**, not as a replacement for the primary ZIBB framework.

---

# 63. Short summary ( Gucci but not as Gucci as one could be )

## Where it started

The project began with prevalence-based PGLS analyses examining relationships between cancer prevalence and life-history traits across Squamata.

Those analyses provided the original biological direction, particularly around traits such as longevity and gestation.

## Where the group is now

The project now has a working Bayesian phylogenetic Zero-Inflated Beta-Binomial pipeline using:

- tumor case counts
- necropsy denominators
- beta-binomial overdispersion
- zero inflation
- phylogenetic relatedness
- standardized life-history predictors
- posterior diagnostics
- posterior predictive checks

Eight univariate ZIBB-PGLMMs have been successfully fitted.

The strongest and most consistent current signal is:

> Maximum longevity is positively associated with neoplasia and, more weakly, malignancy.

Gestation, adult body mass, and litter size do not currently show similarly clear effects in the univariate ZIBB models.

The sampler itself is behaving well, with good convergence and no divergent transitions in the fitted models.

## Where one could go?

The next stage is not to replace the ZIBB framework, but to **refine and stress-test it because Trust but Verify**.

Possible next analyses include:

1. examining the current zero-inflation structure more closely

2. comparing:

```text
zi ~ log_trials_s
```

with:

```text
zi ~ 1
```

3. optionally fitting a standard phylogenetic beta-binomial model as a sensitivity analysis

4. examining whether the main longevity result remains stable across reasonable model specifications

5. selectively running the smaller multivariable models as exploratory analyses

6. comparing the finalized ZIBB results with the original PGLS results

7. using the finalized ZIBB interpretation to guide manuscript revisions

The goal is to determine...

> How robust are the biological conclusions within a ZIBB-centered analysis framework?

---

# 64. Current status

**Analysis status:** Ongoing

**Primary modeling framework:** Phylo ZIBB GLMM

**ZIBB pipeline:** Shi Works

**Successfully fitted models:** 8 univariate models

**MCMC convergence:** Good

**Divergences:** 0 in currently fitted models

**Strongest current predictor:** Maximum longevity

**Strongest current outcome association:** Neoplasia + longevity

**Malignancy + longevity:** Positive, but weaker than the neoplasia result

**Gestation:** No clear effect in the current univariate ZIBB models

**Adult body mass:** No clear effect in the current univariate ZIBB models

**Litter size:** No clear effect in the current univariate ZIBB models

**Zero-count posterior checks:** Some models predict somewhat more zero-case species than observed

**Beta-binomial without zero inflation:** Optional sensitivity analysis, BUT not the primary framework

**Multivariable analysis:** Currently limited by missing life-history data and small complete-case sample sizes

**Skipped model fits:** 6 multivariable model instances were intentionally skipped because `N < 20`

**The above 6 small-sample models:** Can still be run later as exploratory or sensitivity analyses if biologically justified

**Manuscript interpretation:** Should remain provisional (but its close) until the ZIBB diagnostics, sensitivity checks, and selected follow-up analyses are complete

And as always,

Trust, but verify.

(•_•)

( •_•)>⌐■-■

(⌐■_■)
