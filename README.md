# Latent class analysis in R — a worked blueprint

A complete, runnable template for a latent class analysis (LCA) of attitudinal
survey data: from raw items through measurement checks and class enumeration to
robustness checks and publication-ready figures.

It was developed in the course of applied survey research and generalised
here: the structure, the modelling decisions and the robustness checks are
the ones used in practice, with indicators renamed neutrally
(`attitude_*`, `policy_*`, `norms_*`) so the workflow can be adapted to any
attitudinal item battery.

No survey data is included: the notebook runs on a synthetic dataset built
with a known four-class structure, so every method below can be checked
against a ground truth before you point it at real data.

**[→ Read the rendered notebook](lca_blueprint.md)**

## What it covers

1. **Recoding.** Non-response codes to `NA`, reversal of negatively worded
   items with a guard against unrecoded codes, collapsing to three categories.
   All rules live in one function so observed and imputed data are treated
   identically.
2. **Measurement checks.** Cronbach's alpha and largest pairwise correlation
   per item block, a one-factor CFA for the summated scale, and response
   category shares to justify collapsing.
3. **Class enumeration.** `poLCA` fitted across K = 2–6 for eight candidate
   indicator sets, compared on BIC, AIC, relative entropy (Ramaswamy et al.,
   1993), and minimum class size. Respondents with item non-response are kept
   (`na.rm = FALSE`) rather than dropped listwise.
4. **Local independence.** Within-class polychoric correlations as a heuristic
   check that the class structure explains the item associations.
5. **Robustness.**
   - *Survey weights:* row replication for pseudo-maximum-likelihood estimates
     (Patterson, Dayton & Graubard, 2002), with information criteria rescaled
     to the original N — otherwise the inflated likelihood makes BIC favour
     ever more classes.
   - *Missing data:* multiple imputation with `mice` (proportional-odds
     models), indicators rebuilt with the same rules, agreement measured after
     matching class labels.
   - *Indicator set:* re-estimation with an additional item.
   - Class labels are arbitrary across solutions, so they are matched with the
     Hungarian algorithm (`clue::solve_LSAP`) before agreement is computed.
6. **Figures.** Conditional item-response probabilities, weighted covariate
   composition by class, socio-demographic profiles, weighted age
   distributions.

## Files

```
lca_blueprint.Rmd      the analysis notebook
lca_blueprint.md       rendered output with figures
R/functions.R          recoding, measurement, LCA and plotting helpers
R/simulate_data.R      synthetic data generator (known 4-class structure)
```

## Running it

```r
install.packages(c("poLCA", "dplyr", "tidyr", "purrr", "ggplot2", "haven",
                   "psych", "lavaan", "mice", "clue", "here", "scales",
                   "tibble", "rmarkdown"))
rmarkdown::render("lca_blueprint.Rmd")
```

A full run takes roughly 15–20 minutes, mostly in the random starts and the
robustness checks. For a quick look:

```r
rmarkdown::render("lca_blueprint.Rmd",
                  params = list(nrep = 2, run_robustness = FALSE))
```

## Using it on your own data

Place a `.sav` file at the path in `params$data_file` with variables named
`attitude_*`, `policy_*`, `norms_*`, plus `weight`, `age`, `sex`,
`education_raw` and `vote`. Otherwise, adjust the names in the
`recoding-rules` and `export` chunks — those two chunks hold all the
project-specific decisions, and everything downstream follows from them.

The notebook assumes the final weighted model is estimated in LatentGOLD and
falls back to the `poLCA` solution when LatentGOLD output is absent.

## References

Linzer, D. A., & Lewis, J. B. (2011). poLCA: An R package for polytomous
variable latent class analysis. *Journal of Statistical Software*, 42(10), 1–29.

Patterson, B. H., Dayton, C. M., & Graubard, B. I. (2002). Latent class
analysis of complex sample survey data: Application to dietary data.
*Journal of the American Statistical Association*, 97(459), 721–741.

Ramaswamy, V., DeSarbo, W. S., Reibstein, D. J., & Robinson, W. T. (1993). An
empirical pooling approach for estimating marketing mix elasticities with PIMS
data. *Marketing Science*, 12(1), 103–124.

## Author

Julia Gorny · [LinkedIn](https://linkedin.com/in/julia-gorny)
