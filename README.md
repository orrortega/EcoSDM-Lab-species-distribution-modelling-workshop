# EcoSDM Lab

**An interactive R Shiny workshop for learning species distribution modelling (SDM).**

EcoSDM Lab helps students explore how environmental predictors relate to species occurrence, fit and interpret three common modelling approaches, and compare predictions on held-out data. The app uses a **fictional grassland butterfly** and a simulated landscape by default, but also accepts users' own presence–absence CSV files.

Developed by **Dr Oscar Rodriguez de Rivera**, Senior Lecturer in Statistics, University of Exeter, for teaching species distribution modelling and ecological restoration.

## What you can do

The app has five tabs:

| Tab | Features |
| --- | --- |
| **1. Data & upload** | Explore survey data, environmental variables and a correlation matrix; use the built-in example, simulate a new species, or upload CSVs. |
| **2. GLM** | Fit a binomial generalised linear model, select predictors, examine coefficients and predicted occurrence maps. |
| **3. GAM** | Fit a binomial generalised additive model with a choice of linear and cubic-spline terms; inspect effective degrees of freedom and partial effects. |
| **4. Random Forest** | Fit a probability forest, adjust the number of trees and leaf size, and examine permutation variable importance. |
| **5. Compare & diagnose** | Compare held-out performance, ROC curves, calibration, residuals, predicted maps, model disagreement and candidate restoration sites. |

Key teaching features include:

- **Two validation strategies:** stratified random holdout and a single eastern spatial holdout.
- **Model performance measures:** AUC, Brier score, log loss, sensitivity and specificity.
- **A simulated-species designer:** change environmental preferences, occurrence prevalence, survey sample size and sampling bias near roads.
- **Known simulation truth:** reveal the true simulated occurrence probabilities to compare against fitted predictions.
- **Data export:** download example CSVs, the active survey records and model predictions.
- **Equations and explanations:** mathematical expressions and interpretation notes accompany each modelling approach.

## Model fitting and shared hosting

Model fitting is **button-controlled** to reduce unnecessary computation during workshops:

- Each model fits once when first needed so the starting example is available.
- After changing predictors or model settings, click **Fit GLM**, **Fit GAM** or **Fit Random Forest** to update that model.
- Use **Fit all three models** in *Compare & diagnose* after changing the dataset or the validation split.
- Display-only controls, such as changing the plotted predictor or classification threshold, do **not** refit models.
- If the selected settings differ from those used for a fit, the app prompts you to refit rather than silently presenting outdated results.

## Project structure

```text
EcoSDM-Lab/
├── app.R                          # Shiny UI, model functions and server
├── data/
│   ├── simulated_landscape.csv    # Fictional prediction landscape
│   └── simulated_observations.csv # Fictional presence–absence surveys
├── .gitignore                     # Files excluded from deployment
└──README.md                       # This file
```

The included example contains **3,640 simulated landscape cells and 750 simulated surveys**. These represent no real species or geographic location.

## Using your own data

Select **My own CSV files** under *Data & upload*.

1. Upload a **survey CSV** containing one row per survey, a binary `0`/`1` presence–absence column, numeric X and Y coordinates, and numeric environmental predictors.
2. Optionally upload a **prediction-location CSV** containing the same coordinate and predictor columns. No presence–absence response is needed for this second file.
3. Choose the response and coordinate fields in the app, select predictors and fit the models.

For example, survey data might look like:

```csv
presence,x_km,y_km,temperature,rainfall,grassland
1,10,25,13.4,980,72
0,13,27,15.1,860,31
1,18,32,12.7,1100,66
```

These are **illustrative rows, not a complete valid dataset**. The workshop app requires at least **40 complete survey records**, including at least **8 presences and 8 absences**, and limits uploads to **10,000 survey records** and **25,000 prediction locations** after cleaning. It supports at most **35 numeric predictors**, with at least five distinct observed values per predictor. CSV upload size is limited to 20 MB.

Without a separate prediction-location file, predictions are shown **at surveyed locations only**, not across unsampled areas.


## Interpretation and limitations

EcoSDM Lab is designed for **teaching**, not as a validated operational conservation decision system.

- The simulation assumes **perfect detection**; in real monitoring, non-detection does not necessarily mean absence.
- The models use **presence–absence** data, not presence-only records.
- Predictor significance is not the same as ecological importance; the p-value filtering option is exploratory and should not be used for confirmatory inference.
- The GAM allows a mix of linear and spline terms. An estimated smooth with EDF close to 1 is approximately linear.
- Differences among model predictions are **not** formal uncertainty intervals.
- Random and spatial holdouts illustrate different validation questions but do not replace a full external validation study.
- A high predicted occurrence probability alone is not a sufficient basis for choosing a restoration site.

## Online application

The deployed version of the app is available here:

[Launch EcoSDM Lab](https://orderortega.shinyapps.io/EcoSDM_Lab_SDM_worshop/)

## Source code

The source code is available on GitHub:

[GitHub repository](https://github.com/orrortega/EcoSDM-Lab-species-distribution-modelling-workshop)

## Author and contact

**Dr Oscar Rodriguez de Rivera Ortega**  
Senior Lecturer in Statistics, University of Exeter  
[Email](mailto:o.m.rodriguez-de-rivera-ortega@exeter.ac.uk)

