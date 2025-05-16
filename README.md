Multicollinearity Tools
=======================

This package provides robust tools for identifying and removing multicollinearity in predictors, both numeric and categorical. It includes VIF-based filtering, correlation pruning, and hierarchical clustering for categorical variables based on Cramer's V.

Features
--------

- Removes multicollinear numeric features using correlation and VIF filtering
- Groups highly related features into clusters and selects representative features
- Handles NA values and ties intelligently
- Categorical variable reduction using Cramer's V and hierarchical clustering
- Dendrogram visualizations before and after selection

Usage
-----

### Numeric Feature Filtering

    clean <- remove_cont_multicollinearity(
      data,                      # data.frame containing numeric predictors and target
      target,                    # name of the target column
      target_cor_threshold = 0.7,# threshold for correlation with target to drop a feature
      cor_threshold        = 0.7,# threshold for correlation between features to drop one
      vif_threshold        = 5,  # VIF threshold for further filtering
      verbose              = TRUE,   # whether to print messages
      keep_cols            = character(), # columns to always keep
      drop_cols            = character(), # columns to drop before filtering
      draw_corr            = FALSE  # whether to draw a heatmap
    )

### Categorical Feature Filtering

    reduced_factors <- factor_remove_collinearity(
      df,                        # full data frame with factor columns
      target_col     = "TargetColumn",  # name of the target variable
      drop_cols      = "date",   # columns to exclude from analysis
      keep_cols      = character(), # preferred factor columns to retain
      k              = 5         # number of clusters to form
    )


GLMMPreprocessor  
================

**GLMMPreprocessor** is an R6-based utility designed for preprocessing and modeling binary response data using Generalized Linear Mixed Models (GLMMs). It supports automated feature engineering, stepwise BIC-guided model selection, and marginal effect visualization. It is especially useful in contexts where fixed and random effects must be distinguished, polynomial and interaction terms are explored, and interpretability is key.

Core Capabilities  
------------------

- Binary response preprocessing (factor or numeric 0/1)
- Train/test split with optional undersampling to balance classes
- Automatic handling of date and factor variables
- Numeric standardization with saved scaling parameters
- Polynomial expansion (e.g., squared, cubic terms)
- Optional interaction terms with exclusion rules
- Stepwise forward and backward feature selection via BIC
- Mixed model fitting using `glmer()` or fallback to `glm()` if needed
- Support for random intercepts via `(1 | group)` syntax
- Detection and removal of negligible variance random effects
- Fixed-effects marginal prediction plots (1D and 2D)

Typical Usage  
-------------

```r
library(GLMMPreprocessor)

# Set modeling controls
random_effects  <- c("date_month", "Service_Range_Desc_Household")
no_interactions <- c("Unemployment_Rate", "Interest_Rate")
TRAIN_PCT       <- 0.20
RAND_SEED       <- 663

# Prepare preprocessor
glmm <- GLMMPreprocessor$new(target = "LTR_Binary")

# Apply to a single data frame
prep_out <- glmm$preprocess(
  df,
  train_pct         = TRAIN_PCT,
  undersample_train = TRUE,
  random_state      = RAND_SEED
)

# Run stepwise BIC-guided GLMM search
fit_result <- glmm$stepwise_bic_glmm(
  X                 = prep_out$X,
  y                 = prep_out$y,
  random_effects    = random_effects,
  no_interactions   = no_interactions,
  poly              = 3,
  bic_tol           = 1,
  max_add_steps     = 10000,
  max_drop_steps    = 10000,
  cores             = max(1, parallel::detectCores() - 1),
  re_var_threshold  = 1e-6,
  int_var_threshold = 1e-6
)

# Visualize marginal effect of a predictor
glmm$plot_marginal(
  model_fit = fit_result$fit,
  df        = glmm$train_df,
  focal_var = "Unemployment_Rate"
)

# Joint marginal effect (e.g., cubic schedule × distance)
glmm$plot_marginal_interaction(
  model_fit = fit_result$fit,
  df        = glmm$train_df,
  vars      = c("SCHEDULED_DEPARTURE", "DISTANCE")
)
```

Function Reference  
------------------

- `GLMMPreprocessor$new(target)` — initializes the preprocessor with a binary target
- `$preprocess(df, train_pct, undersample_train, random_state)` — preprocesses the data, splits into train/test, and scales numerics
- `$stepwise_bic_glmm(X, y, ...)` — performs stepwise BIC model selection with polynomial and interaction terms
- `$plot_marginal(model_fit, df, focal_var)` — 1D marginal effect plot
- `$plot_marginal_interaction(model_fit, df, vars)` — 2D marginal effect heatmap for two variables

Plotting Tools
--------------

- `plot_marginal()` – Visualize the marginal effect of a single variable
- `plot_marginal_interaction()` – Visualize the joint marginal effects of two numeric predictors using a heatmap with contours

Dependencies
------------

The package uses the following key packages:

R6, lme4, ggplot2, foreach, doParallel, tibble, car, mice, jsonlite, httr, randomForest, caret, pdp, lmerTest, FactoMineR, factoextra, vcd, pheatmap, rpart.plot, dplyr, tidyr, lubridate

Make sure these are installed:

    install.packages(c(
      "R6", "lme4", "lmerTest", "ggplot2", "tibble", "foreach", "doParallel",
      "car", "mice", "jsonlite", "httr", "randomForest", "caret", "pdp",
      "FactoMineR", "factoextra", "vcd", "pheatmap", "rpart.plot", "dplyr", "tidyr", "lubridate"
    ))

Contributing
------------

Contributions are welcome. Please open an issue or submit a pull request if you'd like to suggest a feature or fix a bug.

License
-------

MIT © Alexander Chernikov

