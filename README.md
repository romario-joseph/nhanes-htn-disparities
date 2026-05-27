# Racial Disparities in Hypertension Among U.S. Adults
### A Design-Based Cross-Sectional Analysis of NHANES 2017–2018

**Author:** Romario Joseph, MPH · BU SPH (Epidemiology & Biostatistics)
**Stack:** R 4.x · `survey` · `srvyr` · tidyverse · gtsummary · Quarto
**Data:** National Health and Nutrition Examination Survey (NHANES) 2017–2018 — public-use files released by the U.S. CDC / NCHS.

---

## Epidemiological Objective

Hypertension remains the single largest attributable risk factor for cardiovascular mortality in the United States, and its burden is not distributed equally. Non-Hispanic Black adults experience earlier onset, higher prevalence, and lower treatment-control rates than Non-Hispanic White adults — a gradient that the U.S. Department of Health and Human Services has formally targeted as a *Healthy People 2030* priority. This repository estimates **race- and ethnicity-stratified prevalence, awareness, treatment, and control of hypertension among U.S. adults (≥ 18 years)** using the most recent pre-pandemic NHANES cycle (2017–2018), and asks three interlocking equity questions:

1. **Prevalence.** What is the design-adjusted prevalence of hypertension (JNC-8 / ACC-AHA thresholds) across OMB race/ethnicity categories, and which gaps remain after adjustment for age, sex, education, and insurance status?
2. **Cascade of care.** Among adults with hypertension, what fraction are aware of their diagnosis, on antihypertensive pharmacotherapy, and meeting blood-pressure control targets — stratified by race/ethnicity?
3. **Inferential validity.** Are the observed gaps robust to the **complex survey design** of NHANES, including stratification, clustering, and unequal probability of selection?

Because the third question dictates the inferential validity of the first two, the entire repository is structured to make the design-based statistical adjustment the headline methodological contribution.

---

## Methodological Framework

NHANES is **not** a simple random sample. It is a multistage, stratified, clustered probability sample with deliberate oversampling of selected demographic subgroups. Any analysis that ignores this structure will produce biased point estimates and dramatically miscalibrated standard errors. This repository therefore implements a fully **design-based estimation framework** using the Lumley `survey` package, which is the R-language reference implementation of the Horvitz–Thompson estimator and Taylor-series linearization for complex surveys.

- **Survey design object.** All analyses are conducted on a single `survey.design2` object declared with:
  - **Primary Sampling Unit:** `SDMVPSU` (the masked NHANES PSU identifier),
  - **Stratum:** `SDMVSTRA` (the masked design stratum),
  - **Probability weight:** `WTMEC2YR` (the two-year MEC examination weight; `WTINT2YR` is used for interview-only outcomes),
  - **Nesting flag:** `nest = TRUE` so PSUs are interpreted as nested within strata, as the NCHS documentation requires.

  These four arguments are what mathematically convert the NHANES sample into an unbiased estimator of the U.S. non-institutionalized civilian population.

- **Point estimators.** Prevalence, awareness, treatment, and control proportions are estimated with `survey::svymean()` / `svyciprop()`. For the cascade-of-care outcomes, subdomain estimation is performed via `survey::svyby()` *with the full design object subset*, not by physically filtering the data — the distinction matters because dropping rows changes the variance estimator from a valid subdomain estimator to an invalid one.
- **Variance estimation.** Standard errors are obtained by **Taylor-series linearization**, the default in `survey` and the NCHS-recommended approach for NHANES. 95% confidence intervals use the logit transformation (`svyciprop(… method = "logit")`) so that intervals on bounded proportions remain inside [0, 1] and do not pathologize in small subgroups.
- **Adjusted comparisons.** Race/ethnicity contrasts adjusted for age, sex, education, and insurance are fit with `survey::svyglm()` using a quasi-binomial family and a logit link. Adjusted prevalence ratios are computed by marginal standardization (predicted probabilities averaged across the empirical covariate distribution), avoiding the well-documented odds-ratio inflation that occurs when outcomes are common.
- **Degrees of freedom.** All Wald-type test statistics use `degf(design)` (i.e., PSUs minus strata) rather than `n − 1`, which is essential for valid hypothesis testing in surveys with a modest number of PSUs.
- **Subgroup minimums.** Per NCHS analytic guidelines, point estimates are suppressed when the unweighted denominator falls below 30 or the relative standard error exceeds 30%, with the suppression rule documented in the output tables.

This framework explicitly signals fluency in the **mathematical theory of complex survey designs** — Horvitz–Thompson estimation, design-based variance, and subdomain inference — not merely in the use of the `survey` package as a black box.

---

## Data Architecture

The pipeline ingests NHANES public-use `.XPT` files directly from the NCHS file server, harmonizes them into a single analytic frame, and emits both the headline tables and the diagnostic outputs needed to defend the design specification.

**Stage 1 — Ingest (`R/01_load_nhanes.R`).** Pulls the relevant 2017–2018 components: `DEMO_J` (demographics, weights, PSU, stratum), `BPX_J` (blood pressure examination), `BPQ_J` (blood pressure questionnaire — awareness and medication), and the auxiliary files for education (`DMDEDUC2`), insurance (`HIQ011`), and body measures (`BMX_J`). Files are read with `haven::read_xpt()` and joined on `SEQN`.

**Stage 2 — Cleaning & harmonization (`R/02_clean.R`).**
- **Cohort restriction.** Analysis is restricted to adults aged ≥ 18 with a non-missing MEC two-year weight (`WTMEC2YR > 0`). Participants who were interviewed but not examined are dropped *before* design declaration; this is correct because their MEC weight is zero by construction, but documenting it eliminates an obvious reviewer question.
- **Blood pressure construction.** Per the JNC-8 / ACC-AHA convention used by NHANES, systolic and diastolic blood pressure are computed as the mean of the second and third valid readings (`BPXSY2/BPXSY3`, `BPXDI2/BPXDI3`), discarding the first measurement to mitigate white-coat effect.
- **Hypertension phenotype.** A participant is coded as hypertensive if (a) mean SBP ≥ 130 mmHg, OR (b) mean DBP ≥ 80 mmHg, OR (c) the participant self-reports current use of antihypertensive medication (`BPQ050A == 1`). Sensitivity analyses re-run the headline tables at the JNC-7 threshold (140/90) and report both, so the conclusions are not threshold-dependent.
- **Awareness, treatment, control.** Defined per the NHANES analytic conventions: awareness = ever told by a provider; treatment = currently on antihypertensive medication; control = on treatment AND mean BP < 130/80.
- **Race/ethnicity recoding.** `RIDRETH3` is collapsed into five mutually exclusive OMB-aligned categories: Non-Hispanic White, Non-Hispanic Black, Hispanic (Mexican-American + Other Hispanic combined), Non-Hispanic Asian, and Other / Multiracial. The Asian category is preserved because NHANES has oversampled this group since 2011–2012, and its inclusion is precisely what `WTMEC2YR` is calibrated to support.
- **Missing data.** Item-level missingness is handled via complete-case analysis on the analytic variables, with the unweighted retention rate at each step logged to `outputs/cohort_flow.csv`. Because NHANES non-response is already incorporated into `WTMEC2YR` by the NCHS weighting protocol, additional imputation on the design variables is intentionally avoided.

**Stage 3 — Design declaration (`R/03_design.R`).** Builds the single `survey.design2` object that downstream scripts consume. Declared exactly once so the design specification is not re-implemented inconsistently across analyses.

**Stage 4 — Estimation (`R/04_estimate.R`).** Produces (i) overall and race-stratified prevalence with logit-method 95% CIs; (ii) the four-stage cascade-of-care tables; (iii) adjusted prevalence ratios via marginal standardization from `svyglm`; and (iv) a diagnostic table comparing weighted vs unweighted prevalence so the magnitude of the design adjustment is transparent.

**Stage 5 — Reporting (`reports/nhanes_disparities.qmd`).** Quarto report rendering the publication-ready figures and tables.

```
nhanes-htn-disparities/
├── README.md
├── R/
│   ├── 01_load_nhanes.R         ← ingest DEMO_J, BPX_J, BPQ_J, BMX_J
│   ├── 02_clean.R               ← cohort restriction, BP phenotype, race recoding
│   ├── 03_design.R              ← survey.design2 with PSU, stratum, WTMEC2YR
│   └── 04_estimate.R            ← svyciprop, svyby, svyglm + marginal standardization
├── data/
│   └── README.md                ← NHANES file URLs and SHA hashes
├── reports/
│   └── nhanes_disparities.qmd
├── outputs/
│   ├── tables/
│   └── figures/
└── LICENSE
```

---

## Reproduce

```r
install.packages(c("survey","srvyr","haven","tidyverse","gtsummary","quarto"))
source("R/01_load_nhanes.R")   # downloads NHANES 2017–2018 .XPT files
source("R/02_clean.R")         # cohort + BP phenotype + race recoding
source("R/03_design.R")        # declares survey.design2 with PSU, stratum, WTMEC2YR
source("R/04_estimate.R")      # design-adjusted prevalence + cascade + svyglm
quarto::quarto_render("reports/nhanes_disparities.qmd")
```

---

## Why this matters

A Black–White hypertension gap estimated from unweighted NHANES data is not a population estimate — it is a sample statistic that confuses the oversampling design with a population signal. Implementing the design-based framework above is what licenses any equity claim made from these data, and the same `survey`-package pipeline carries over directly to other complex surveys (BRFSS, MEPS, the Belgian Health Interview Survey) where the same PSU / stratum / weight machinery applies.

---

## Disclaimer

NHANES public-use files contain no direct identifiers and are released by NCHS for unrestricted secondary analysis. Results in this repository should not be interpreted as a substitute for the NCHS Data Briefs.

## Contact

Romario Joseph · rjoseph3@bu.edu · [LinkedIn](https://www.linkedin.com/in/romariojosephpublichealth/)
