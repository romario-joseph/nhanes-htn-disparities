# Racial Disparities in Hypertension Among U.S. Adults
### A design-based cross-sectional analysis of NHANES 2017-2018

**Author:** Romario Joseph, MPH (BU SPH, Epidemiology & Biostatistics)
**Stack:** R 4.x, `survey`, `srvyr`, tidyverse, gtsummary, Quarto
**Data:** National Health and Nutrition Examination Survey (NHANES) 2017-2018. Public-use files released by the U.S. CDC / NCHS.

I built this project to do one thing properly: estimate the Black-White gap in hypertension prevalence, awareness, treatment, and control from NHANES, and do it with the survey-design machinery that the data actually require, not the unweighted shortcut that most early-career analyses settle for.

## Epidemiological Objective

Hypertension is still the single largest attributable risk factor for cardiovascular mortality in the United States, and the burden is not distributed evenly. Non-Hispanic Black adults experience earlier onset, higher prevalence, and lower treatment-control rates than non-Hispanic White adults. That gap is on the *Healthy People 2030* priority list for a reason, and any honest equity claim has to come from an analysis that respects the NHANES sampling design.

The three questions I am asking, in order:

1. Prevalence. What is the design-adjusted prevalence of hypertension (JNC-8 / ACC-AHA thresholds) across OMB race/ethnicity categories, and which gaps remain after adjustment for age, sex, education, and insurance status?
2. Cascade of care. Among adults with hypertension, what fraction are aware of their diagnosis, on antihypertensive pharmacotherapy, and meeting blood-pressure control targets, stratified by race/ethnicity?
3. Inferential validity. Are the observed gaps robust to the complex survey design of NHANES (stratification, clustering, unequal probability of selection)?

The third question is the one that decides whether the first two mean anything, so I built the whole repo around the design-based estimation framework.

## Methodological Framework

NHANES is not a simple random sample. It is a multistage, stratified, clustered probability sample with deliberate oversampling of selected demographic subgroups. Any analysis that ignores that structure produces biased point estimates and badly miscalibrated standard errors. I therefore implemented the analysis end-to-end in the Lumley `survey` package, which is the R reference implementation of the Horvitz-Thompson estimator and Taylor-series linearization for complex surveys.

**Survey design object.** I run every analysis on a single `survey.design2` object declared with:

- Primary Sampling Unit: `SDMVPSU` (the masked NHANES PSU identifier)
- Stratum: `SDMVSTRA` (the masked design stratum)
- Probability weight: `WTMEC2YR` (the two-year MEC examination weight; `WTINT2YR` for interview-only outcomes)
- Nesting flag: `nest = TRUE` so PSUs are interpreted as nested within strata, exactly as the NCHS documentation requires

Those four arguments are the mathematical machinery that converts the NHANES sample into an unbiased estimator of the U.S. non-institutionalized civilian population.

**Point estimators.** I estimate prevalence, awareness, treatment, and control proportions with `survey::svymean()` and `survey::svyciprop()`. For the cascade-of-care outcomes I do subdomain estimation via `survey::svyby()` with the full design object subset, never by physically filtering the data, because dropping rows turns a valid subdomain estimator into an invalid one. That distinction is exactly the kind of detail a methods reviewer will look for.

**Variance estimation.** Standard errors come from Taylor-series linearization, which is the default in `survey` and the NCHS-recommended approach for NHANES. I use the logit transformation in `svyciprop(... method = "logit")` so the 95% intervals on bounded proportions stay inside (0, 1) and do not pathologize in small subgroups.

**Adjusted comparisons.** I fit race/ethnicity contrasts adjusted for age, sex, education, and insurance with `survey::svyglm()` using a quasi-binomial family and a logit link. Adjusted prevalence ratios come from marginal standardization (predicted probabilities averaged across the empirical covariate distribution), which sidesteps the odds-ratio inflation you get when the outcome is common.

**Degrees of freedom.** All Wald-type test statistics use `degf(design)` (PSUs minus strata) rather than n minus 1, which matters because NHANES has a modest number of PSUs and the wrong df makes the p-values lie.

**Subgroup minimums.** Per NCHS analytic guidelines, I suppress point estimates when the unweighted denominator falls below 30 or the relative standard error exceeds 30%, and the suppression rule is documented in the output tables so a reader can see exactly which cells were dropped and why.

The point of being this explicit is that the framework is mathematical (Horvitz-Thompson estimation, design-based variance, subdomain inference), not just a black-box call to the `survey` package.

## Data Architecture

The pipeline ingests NHANES public-use `.XPT` files directly from the NCHS file server, harmonizes them into a single analytic frame, and emits both the headline tables and the diagnostics I need to defend the design specification.

**Stage 1, ingest (`R/01_load_nhanes.R`).** Pulls the relevant 2017-2018 components: `DEMO_J` (demographics, weights, PSU, stratum), `BPX_J` (blood pressure examination), `BPQ_J` (blood pressure questionnaire, for awareness and medication), and the auxiliary files for education (`DMDEDUC2`), insurance (`HIQ011`), and body measures (`BMX_J`). Files come in via `haven::read_xpt()` and join on `SEQN`.

**Stage 2, cleaning and harmonization (`R/02_clean.R`).**

*Cohort restriction.* I restrict to adults aged 18 and over with a non-missing MEC two-year weight (`WTMEC2YR > 0`). Participants who were interviewed but not examined are dropped before the design declaration. They have a MEC weight of zero by construction, so this is correct, but I document it so the reviewer does not have to chase the question.

*Blood pressure construction.* Following the NHANES convention, systolic and diastolic blood pressure are the mean of the second and third valid readings (`BPXSY2/BPXSY3`, `BPXDI2/BPXDI3`). The first reading is discarded to mitigate white-coat effect.

*Hypertension phenotype.* A participant is hypertensive if (a) mean SBP is at or above 130 mmHg, OR (b) mean DBP is at or above 80 mmHg, OR (c) the participant self-reports current use of antihypertensive medication (`BPQ050A == 1`). Sensitivity analyses re-run the headline tables at the JNC-7 threshold (140/90) and report both, so the conclusions do not hinge on the threshold.

*Awareness, treatment, control.* Defined per NHANES analytic conventions: awareness is "ever told by a provider"; treatment is current antihypertensive medication; control is on treatment AND mean BP under 130/80.

*Race/ethnicity recoding.* I collapse `RIDRETH3` into five mutually exclusive OMB-aligned categories: Non-Hispanic White, Non-Hispanic Black, Hispanic (Mexican-American plus Other Hispanic), Non-Hispanic Asian, and Other / Multiracial. I keep the Asian category because NHANES has oversampled this group since 2011-2012, and its inclusion is exactly what `WTMEC2YR` is calibrated to support.

*Missing data.* I handle item-level missingness with complete-case analysis on the analytic variables and log the unweighted retention rate at each step to `outputs/cohort_flow.csv`. NHANES non-response is already baked into `WTMEC2YR` by the NCHS weighting protocol, so I intentionally avoid additional imputation on the design variables.

**Stage 3, design declaration (`R/03_design.R`).** Builds the single `survey.design2` object that downstream scripts consume. I declare it exactly once so the design specification is not re-implemented inconsistently across analyses.

**Stage 4, estimation (`R/04_estimate.R`).** Produces (i) overall and race-stratified prevalence with logit-method 95% CIs, (ii) the four-stage cascade-of-care tables, (iii) adjusted prevalence ratios via marginal standardization from `svyglm`, and (iv) a diagnostic table comparing weighted and unweighted prevalence so the magnitude of the design adjustment is visible to a reader.

**Stage 5, reporting (`reports/nhanes_disparities.qmd`).** Quarto report with the publication-ready figures and tables.

```
nhanes-htn-disparities/
├── README.md
├── R/
│   ├── 01_load_nhanes.R         # ingest DEMO_J, BPX_J, BPQ_J, BMX_J
│   ├── 02_clean.R               # cohort restriction, BP phenotype, race recoding
│   ├── 03_design.R              # survey.design2 with PSU, stratum, WTMEC2YR
│   └── 04_estimate.R            # svyciprop, svyby, svyglm + marginal standardization
├── data/
│   └── README.md                # NHANES file URLs and SHA hashes
├── reports/
│   └── nhanes_disparities.qmd
├── outputs/
│   ├── tables/
│   └── figures/
└── LICENSE
```

## Reproduce

```r
install.packages(c("survey","srvyr","haven","tidyverse","gtsummary","quarto"))
source("R/01_load_nhanes.R")   # downloads NHANES 2017-2018 .XPT files
source("R/02_clean.R")         # cohort + BP phenotype + race recoding
source("R/03_design.R")        # declares survey.design2 with PSU, stratum, WTMEC2YR
source("R/04_estimate.R")      # design-adjusted prevalence + cascade + svyglm
quarto::quarto_render("reports/nhanes_disparities.qmd")
```

## Why this matters

A Black-White hypertension gap estimated from unweighted NHANES data is not a population estimate. It is a sample statistic that confuses the oversampling design with a population signal. Implementing the design-based framework above is what licenses any equity claim made from these data. The same `survey`-package pipeline carries over directly to other complex surveys (BRFSS, MEPS, the Belgian Health Interview Survey), where the same PSU / stratum / weight machinery applies.

## Disclaimer

NHANES public-use files contain no direct identifiers and are released by NCHS for unrestricted secondary analysis. Results in this repository should not be interpreted as a substitute for the NCHS Data Briefs.

## Contact

Romario Joseph, rjoseph3@bu.edu, [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
