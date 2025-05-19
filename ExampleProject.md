# ExampleProject

2025-05-19

## Install or update the **GLMMPreprocessor** package

The first step guarantees that you are working with the most recent development build of the package from GitHub.  If an older version is already installed it is removed, and the latest commit is re‑installed.

```r
if ("GLMMPreprocessor" %in% rownames(installed.packages())) {
  remove.packages("GLMMPreprocessor")
}

devtools::install_github("alexChernikov1/GLMMPreprocessor", force = TRUE, upgrade = "never")
```

## Load the analysis toolkit

Below we load **GLMMPreprocessor** together with a set of tidy‑verse, modelling, and visualisation packages that will be used throughout the workflow.

```r
library(GLMMPreprocessor)
library(usdm)
library(FactoMineR)
library(factoextra)
library(rpart.plot)
library(glmnet)
library(randomForest)
library(glmnet)
library(caret)
library(tidyr)
library(dplyr)
library(nlme)
library(lubridate)
library(lme4)
library(car)
library(lmerTest)
library(ggplot2)
library(DescTools)
library(mice)
library(httr)
library(jsonlite)
library(dplyr)
library(vcd)  
library(pheatmap)
library(dplyr)
library(ggplot2)
library(lubridate)
library(httr)
library(jsonlite)
library(FactoMineR)
library(factoextra)
library(glmnet)
library(dplyr)
library(randomForest)
library(pdp)
library(R6)
library(lme4)
library(dplyr)
library(tibble)
library(foreach)
library(doParallel)
```

## Bring the example data into the session

We use the full U.S. domestic‐flights data set bundled with the package.  The chunk reads the `.RData` file and leaves us with a 4‑million‑row tibble.  A preview of the structure is shown automatically when the chunk runs.

```r
load("~/GLMMPreprocessor_data/airtraffic.RData")
```

tibble $4,045,608 × 17$ (S3: tbl\_df/tbl/data.frame)

YEAR : num $1:4045608$ 2015 2015 2015 2015 2015 …

MONTH : num $1:4045608$ 1 1 1 1 1 1 1 1 1 1 …

DAY : num $1:4045608$ 1 1 1 1 1 1 1 1 1 1 …

DAY\_OF\_WEEK : Factor w/ 7 levels “1”,“2”,“3”,“4”,..: 4 4 4 4 4 4 4 4 4 4 …

AIRLINE : Factor w/ 14 levels “AA”,“AS”,“B6”,..: 2 1 12 1 2 4 9 12 1 4 …

TAIL\_NUMBER : Factor w/ 4888 levels “7819A”,“7820L”,..: 1618 1552 422 1512 2127 1137 2760 2405 1557 3929 …

ORIGIN\_AIRPORT : Factor w/ 76 levels “ABQ”,“ANC”,“ATL”,..: 2 36 67 36 66 67 35 36 67 35 …

DESTINATION\_AIRPORT : Factor w/ 76 levels “ABQ”,“ANC”,“ATL”,..: 66 54 15 44 2 46 46 15 21 3 …

SCHEDULED\_DEPARTURE : num $1:4045608$ 5 10 20 20 25 25 25 30 30 30 …

SCHEDULED\_TIME : num $1:4045608$ 205 280 286 285 235 217 181 273 195 221 …

DISTANCE : num $1:4045608$ 1448 2330 2296 2342 1448 …

ORIGIN\_LATITUDE : num $1:4045608$ 61.2 33.9 37.6 33.9 47.4 …

ORIGIN\_LONGITUDE : num $1:4045608$ -150 -118 -122 -118 -122 …

DESTINATION\_LATITUDE : num $1:4045608$ 47.4 26.7 35.2 25.8 61.2 …

DESTINATION\_LONGITUDE: num $1:4045608$ -122.3 -80.1 -80.9 -80.3 -150 …

ARRIVAL\_DELAY\_BINOM : int $1:4045608$ 0 0 0 0 0 0 0 0 0 0 …

ROUTE : Factor w/ 2735 levels “ABQ\_ATL”,“ABQ\_BWI”,..: 38 1375 2420 1366 2356 2443 1308 1342 2425 1271 …

## 1 · Continuous‑predictor screening

Here we remove numeric predictors that are redundant or strongly collinear.  The helper `remove_cont_multicollinearity()` (part of **GLMMPreprocessor**)

* drops user‑specified columns,\*
* prunes features highly correlated with the **target**,\*
* clusters the remaining variables, keeping one representative per cluster,\*
* and provides an optional VIF check.\*
  The heat‑map makes it easy to spot any remaining high pairwise correlations.

```r
cont_result <- invisible(remove_cont_multicollinearity(
  data                  = airtraffic,
  target                = "ARRIVAL_DELAY_BINOM",
  target_cor_threshold  = 0.7,
  target_corr_last      = TRUE,
  cor_threshold         = 0.7,
  vif_threshold         = 5,
  verbose               = TRUE,
  keep_cols             = c(),
  drop_cols             = c("DAY"),
  draw_corr             = TRUE
))
```

Dropped user‑specified columns: DAY

Dropping 'DISTANCE' corr> 0.7 with 'SCHEDULED\_TIME'.

Feature Group summary (chosen / others):

Feature Group 1: MONTH

Feature Group 2: SCHEDULED\_DEPARTURE

Feature Group 3: SCHEDULED\_TIME  /  DISTANCE

Feature Group 4: ORIGIN\_LATITUDE

Feature Group 5: ORIGIN\_LONGITUDE

Feature Group 6: DESTINATION\_LATITUDE

Feature Group 7: DESTINATION\_LONGITUDE

variable cluster removed

MONTH MONTH 1 FALSE

SCHEDULED\_DEPARTURE SCHEDULED\_DEPARTURE 2 FALSE

SCHEDULED\_TIME SCHEDULED\_TIME 3 FALSE

DISTANCE DISTANCE 3 TRUE

ORIGIN\_LATITUDE ORIGIN\_LATITUDE 4 FALSE

ORIGIN\_LONGITUDE ORIGIN\_LONGITUDE 5 FALSE

DESTINATION\_LATITUDE DESTINATION\_LATITUDE 6 FALSE

DESTINATION\_LONGITUDE DESTINATION\_LONGITUDE 7 FALSE

![](ExampleProject_files/figure-gfm/unnamed-chunk-3-1.png)

Before: 9 numeric cols; After: 7 retained.
Retained columns:
“MONTH” “SCHEDULED\_DEPARTURE” “SCHEDULED\_TIME”
“ORIGIN\_LATITUDE” “ORIGIN\_LONGITUDE” “DESTINATION\_LATITUDE”
“DESTINATION\_LONGITUDE”

## 2 · Categorical‑predictor screening

Next we tackle high‑cardinality factors and overlapping categorical information using `remove_factor_multicollinearity()`.  The procedure clusters factors with similar level patterns (measured by Cramer’s V) and keeps a single representative of each cluster.  Optionally, association with the target can be used as a last pruning step.

```r
cat_result <- invisible(remove_factor_multicollinearity(
  df                       = airtraffic,
  target_col               = "ARRIVAL_DELAY_BINOM",
  prune_target_assoc       = 0.7,
  prune_target_assoc_last  = TRUE,
  keep_cols                = c(),
  drop_cols                = c("DAY"),
  k                        = 4,
  verbose                  = TRUE
))
```

Dropped user‑specified columns: DAY

Pruning numeric predictors: YEAR, MONTH, SCHEDULED\_DEPARTURE, SCHEDULED\_TIME, DISTANCE, ORIGIN\_LATITUDE, ORIGIN\_LONGITUDE, DESTINATION\_LATITUDE, DESTINATION\_LONGITUDE

![](ExampleProject_files/figure-gfm/unnamed-chunk-4-1.png)

— Cluster Membership and NA Counts —
Feature Group summary (chosen / others):
Feature Group 1: DAY\_OF\_WEEK
Feature Group 2: AIRLINE / TAIL\_NUMBER
Feature Group 3: ORIGIN\_AIRPORT / ROUTE
Feature Group 4: DESTINATION\_AIRPORT

![](ExampleProject_files/figure-gfm/unnamed-chunk-4-2.png)

Before: 6 factor cols; After: 4 retained.
Retained columns: “DAY\_OF\_WEEK” “AIRLINE” “ORIGIN\_AIRPORT” “DESTINATION\_AIRPORT”

## 3 · Assemble the modelling data set

We gather the retained numeric and factor predictors, bring back the binary **target**, and drop the geographical latitude/longitude fields that will not be used in this particular example.

```r
# Step 3: Extract selected column names
numeric_cols     <- setdiff(names(cont_result$pruned_data), "ARRIVAL_DELAY_BINOM")
categorical_cols <- names(cat_result$pruned_data)

# Step 4: Combine selected columns

df_cleaned <- airtraffic[, c("ARRIVAL_DELAY_BINOM", numeric_cols, categorical_cols), drop = FALSE]
```

```r
str(df_cleaned)
```

tibble $4,045,608 × 12$ (S3: tbl\_df/tbl/data.frame)

ARRIVAL\_DELAY\_BINOM : int $1:4045608$ 0 0 0 0 0 0 0 0 0 0 …

MONTH : num $1:4045608$ 1 1 1 1 1 1 1 1 1 1 …

SCHEDULED\_DEPARTURE : num $1:4045608$ 5 10 20 20 25 25 25 30 30 30 …

SCHEDULED\_TIME : num $1:4045608$ 205 280 286 285 235 217 181 273 195 221 …

ORIGIN\_LATITUDE : num $1:4045608$ 61.2 33.9 37.6 33.9 47.4 …

ORIGIN\_LONGITUDE : num $1:4045608$ -150 -118 -122 -118 -122 …

DESTINATION\_LATITUDE : num $1:4045608$ 47.4 26.7 35.2 25.8 61.2 …

DESTINATION\_LONGITUDE: num $1:4045608$ -122.3 -80.1 -80.9 -80.3 -150 …

DAY\_OF\_WEEK : Factor w/ 7 levels “1”,“2”,“3”,“4”,..: 4 4 4 4 4 4 4 4 4 4 …

AIRLINE : Factor w/ 14 levels “AA”,“AS”,“B6”,..: 2 1 12 1 2 4 9 12 1 4 …

ORIGIN\_AIRPORT : Factor w/ 76 levels “ABQ”,“ANC”,“ATL”,..: 2 36 67 36 66 67 35 36 67 35 …

DESTINATION\_AIRPORT : Factor w/ 76 levels “ABQ”,“ANC”,“ATL”,..: 66 54 15 44 2 46 46 15 21 3 …

```r
to_drop <- c("ORIGIN_LATITUDE","ORIGIN_LONGITUDE","DESTINATION_LATITUDE","DESTINATION_LONGITUDE")

df_final <- df_cleaned %>% select(-all_of(to_drop))
```

## 4 · Configure the GLMM pre‑processor

We define the grouping variables that will have random intercepts, set the train/test split, and instantiate **GLMMPreprocessor** for a binary outcome.

```r
random_effects  <- c("MONTH", "AIRLINE","ORIGIN_AIRPORT","DESTINATION_AIRPORT")
no_interactions <- c()
TRAIN_PCT       <- 0.005
RAND_SEED       <- 663

glmm <- GLMMPreprocessor$new(target = "ARRIVAL_DELAY_BINOM")
```

## 5 · Prepare the training / test split

The `preprocess()` method scales numeric predictors, encodes factors, and (optionally) balances the training set via undersampling so that both classes are represented equally.

```r
prep_out <- glmm$preprocess(
  df_final,
  train_pct         = TRAIN_PCT,
  undersample_train = TRUE,
  random_state      = RAND_SEED
)
```

After undersampling: train rows = 7558 (balanced 0/1)

## 6 · Stepwise BIC‑guided model search

`stepwise_bic_glmm()` explores polynomial terms, interactions, and random‑effect structures, greedily adding or dropping fixed effects according to the Bayesian Information Criterion.  The search stops when no candidate improves (or worsens within 50 %) of the best BIC so far.

```r
fit_result <- glmm$stepwise_bic_glmm(
  X                 = prep_out$X,
  y                 = prep_out$y,
  random_effects    = random_effects,
  no_interactions   = no_interactions,
  poly              = 3,
  bic_tol           = 0.5,
  max_add_steps     = 1000,
  max_drop_steps    = 1000,
  cores             = max(1, parallel::detectCores() - 1),
  re_var_threshold  = 1e-6,
  int_var_threshold = 1e-6
)
```

Fixed‑effect search space: 40 terms

(+ / – output truncated for brevity)

## 7 · Inspect the final model

Below is the summary of the selected GLMM, including random‑effect standard deviations and fixed‑effect coefficients.

```r
fit_result$fit
```

Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss‑Hermite Quadrature, nAGQ = 0) \[glmerMod]

(... output truncated ...)

## 8 · Visualise marginal relationships

The convenience wrapper `plot_marginal()` draws the population‑average predicted probability curve for a single focal predictor, holding all other variables at typical values.

```r
glmm$plot_marginal(
  model_fit = fit_result$fit,
  df        = glmm$train_df,
  focal_var = "SCHEDULED_DEPARTURE"
)
```

![](ExampleProject_files/figure-gfm/unnamed-chunk-11-1.png)

## 9 · Explore two‑way interactions

Finally, `plot_marginal_interaction()` shows how the predicted probability surface changes across the joint range of two predictors—in this case, departure time and longitude.

```r
glmm$plot_marginal_interaction(
  model_fit = fit_result$fit,
  df        = glmm$train_df,
  vars      = c("SCHEDULED_DEPARTURE", "ORIGIN_LONGITUDE")
)
```

![](ExampleProject_files/figure-gfm/unnamed-chunk-12-1.png)

