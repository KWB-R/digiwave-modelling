# DigiWave — Predictive UV-dose control for wastewater disinfection (public release)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21882515.svg)](https://doi.org/10.5281/zenodo.21882515)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reduce the energy of the effluent **UV disinfection** stage by dosing to demand:
predict the E. coli concentration entering the UV unit from plant sensor signals,
combine it with a UV **dose–response** relation, back-calculate the **minimum dose**
to meet the discharge limit, and quantify the **energy saving** vs. a constant dose.

## Data anonymization

The plant **sensor features** (all `*_roll3h` predictors) are anonymized: each is
affine-transformed with per-variable coefficients that are **withheld**, so their
absolute values, units, and scale cannot be recovered from this release. The fitted
models are transformed with the same map, so predictions reproduce as in the
original study. **E. coli measurements, timestamps, and the UV dose–response are
real and unchanged.**

Raw sensor data are available from the authors on reasonable request.

## Reproducing the results

```
Rscript R/_run_all.R      # runs stages 03–06, regenerates every figure in plots/
```

Stages **01 (clean raw plant data)** and **02 (train models)** require the
confidential raw sensor data, so they are shipped as **documented code only** and
are not executed by `_run_all.R`. Stages 03–06 run from the sanitized artifacts.

## Pipeline

| Stage | Does | Runnable here |
|-------|------|:---:|
| `01_clean_plant_data.R` | Clean plant CSV → 3 h rolling-mean features | code-only |
| `02_train_models.R` | Build training frame, QRF variable selection, fit models | code-only |
| `03_dose_response.R` | Fit UV dose → LRV model (`stanDR`) | ✅ |
| `04_dose_schedule.R` | Reconstruct applied UV-dose schedule | ✅ |
| `05_predict_validate.R` | Validate predictions vs. lab/AquaBio, overlay dose | ✅ |
| `06_energy_savings.R` | Recommended dose + energy saving | ✅ |

## The three candidate models

Fit in `02` by Quantile-Random-Forest importance-ranked forward selection (inlet
temperature `KA_Inf_T` is excluded — not available at scoring time). Of 19 rolled
sensor predictors:

| Model | # predictors | Colour | Selection rule |
|---|:---:|---|---|
| `bestR2` | 11 | blue | prefix length maximizing cross-val R² |
| `best7` | 6 | purple | top-7 importance prefix |
| `best4of7` | 3 | green | high marginal-gain subset |

Comparison (`plots/model_comparison.csv`): bestR2 is best (RMSE 0.25, R² 0.51,
Pearson 0.80); best7 (RMSE 0.30, R² 0.28); best4of7 is too sparse (R² −0.23).

## Dependencies

R packages: `tidyverse, data.table, zoo, quantregForest, Metrics, rstanarm,
lubridate, gridExtra`.

## Citation

Seis, W. (2026). *DigiWave — Predictive UV-dose control for wastewater
disinfection.* Kompetenzzentrum Wasser Berlin. https://doi.org/10.5281/zenodo.21882515

The DOI above always resolves to the latest release. Version-specific DOIs are
listed under "Versions" on the [Zenodo record](https://doi.org/10.5281/zenodo.21882515).
