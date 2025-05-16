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

```r
# Step 1: Filter numeric features
cont_result <- remove_cont_multicollinearity(
  data                  = df,             # full dataframe
  target                = "TargetColumn", # name of the target variable (must be 0/1 numeric or binary factor)
  target_cor_threshold  = 0.7,            # threshold to drop predictors highly correlated with target
  cor_threshold         = 0.7,            # threshold to drop one of any two highly correlated predictors
  vif_threshold         = 5,              # max acceptable variance inflation factor
  verbose               = TRUE,           # print messages about dropped variables
  keep_cols             = character(),    # predictors you want to force-keep
  drop_cols             = character(),    # predictors you want to drop before filtering
  draw_corr             = FALSE           # whether to plot the correlation heatmap
)

# Step 2: Filter categorical features
cat_result <- factor_remove_collinearity(
  df             = df,                   # original dataframe with factors
  target_col     = "TargetColumn",       # name of the target variable
  drop_cols      = "date",               # drop date before computing factor similarity
  keep_cols      = character(),          # factors you want to force-keep
  k              = 5                     # number of clusters to use in hierarchical clustering
)

# Step 3: Extract selected column names
numeric_cols     <- setdiff(names(cont_result$pruned_data), "TargetColumn") # selected numeric features
categorical_cols <- colnames(cat_result)                                    # selected factor features

# Step 4: Combine selected columns from original df
# Note: 'date' column is manually added back here
df_cleaned <- df[, c("TargetColumn", numeric_cols, categorical_cols, "date"), drop = FALSE]

```



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
random_effects  <- c("random_effect_1", "random_effect_2")   # grouping variables for random effects
no_interactions <- c("fixed_effect_1", "fixed_effect_2")     # fixed effects to exclude from interaction terms
TRAIN_PCT       <- 0.20                                       # percentage of data to use for training
RAND_SEED       <- 663                                        # seed for reproducibility

# Prepare preprocessor
glmm <- GLMMPreprocessor$new(target = "target_feature")      # initialize with binary target variable

# Apply to a single data frame
prep_out <- glmm$preprocess(
  df,
  train_pct         = TRAIN_PCT,        # fraction of rows to use in training set
  undersample_train = TRUE,             # whether to undersample majority class
  random_state      = RAND_SEED         # seed for reproducible sampling
)

# Run stepwise BIC-guided GLMM search
fit_result <- glmm$stepwise_bic_glmm(
  X                 = prep_out$X,               # preprocessed feature matrix
  y                 = prep_out$y,               # target vector
  random_effects    = random_effects,           # random effect variables
  no_interactions   = no_interactions,          # fixed effects to exclude from interactions
  poly              = 3,                        # include squared and cubic terms
  bic_tol           = 0.5,                      # allow up to 50% increase *from previous best BIC* when adding/removing terms
  max_add_steps     = 10000,                    # maximum terms to add during forward search
  max_drop_steps    = 10000,                    # maximum terms to drop during backward search
  cores             = max(1, parallel::detectCores() - 1),  # number of cores for parallel processing
  re_var_threshold  = 1e-6,                     # drop random effects with negligible variance
  int_var_threshold = 1e-6                      # exclude numeric features with low variance from interaction pool
)

# Visualize marginal effect of a predictor
glmm$plot_marginal(
  model_fit = fit_result$fit,                   # fitted model
  df        = glmm$train_df,                    # training data used for centering/scaling
  focal_var = "fixed_effect_3"                  # focal variable to visualize
)

# Joint marginal effect (e.g., cubic schedule × distance)
glmm$plot_marginal_interaction(
  model_fit = fit_result$fit,                   # fitted model
  df        = glmm$train_df,                    # training data
  vars      = c("fixed_effect_3", "fixed_effect_4")  # variables for interaction plot
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

