Multicollinearity Tools
=======================

This package provides robust tools for identifying and removing multicollinearity in predictors, both numeric and categorical. It includes VIF-based filtering, correlation pruning, and hierarchical clustering for categorical variables based on Cramer’s V. To see an implementation of the package: [ExampleProject](ExampleProject.md)

Features
--------

- Removes multicollinear **numeric** features using target‐correlation, pairwise correlation and VIF filtering  
- Groups highly related numeric features into clusters and selects representative features  
- Categorical variable reduction using Cramer’s V and hierarchical clustering  
- Intelligent handling of NA values, ties and single‐level factors  
- Dendrogram visualizations before and after selection  

Usage
-----

### Numeric feature reduction

```r
cont_result <- remove_cont_multicollinearity(
  data                  = df,             # full dataframe
  target                = "TargetColumn", # name of the target variable (0/1 or binary factor)
  target_cor_threshold  = 0.7,            # drop predictors with |corr| > 0.7 vs target
  target_corr_last      = TRUE,           # remove features correlated with target after VIF pruning
  cor_threshold         = 0.7,            # drop one of any two predictors with |corr| > 0.7
  vif_threshold         = 5,              # drop predictors until all VIFs < 5
  verbose               = TRUE,           # print diagnostic messages
  drop_cols             = character(),    # predictors to drop before filtering
  keep_cols             = character(),    # predictors to force-keep
  draw_corr             = FALSE           # plot final correlation heatmap
)

# cont_result is an invisible list with:
#  • pruned_data: data.frame(target, retained numeric predictors)
#  • cluster_df:   data.frame(variable, cluster, removed)
```

### Categorical feature reduction

```r
cat_result <- remove_factor_multicollinearity(
  df             = df,                   # original dataframe
  target_col     = "TargetColumn",       # name of the target (for filtering non-factors)
  prune_target_assoc        = 0.7,       # new: Cramer's V threshold vs target
  prune_target_assoc_last   = TRUE,      # new: run target–assoc prune after clustering?
  keep_cols      = character(),          # factors to force-keep
  drop_cols      = "date",               # columns to drop before clustering
  k              = 5,                    # number of clusters
  verbose        = TRUE                  # print diagnostic messages
)

# cat_result is an invisible list with:
#  • pruned_data: data.frame of retained factor predictors
#  • cluster_df:   data.frame(variable, cluster, NA_count, removed)
```

### Combine selected features

```r
numeric_cols     <- setdiff(names(cont_result$pruned_data), "TargetColumn")
categorical_cols <- colnames(cat_result$pruned_data)

df_cleaned <- df[, c("TargetColumn", numeric_cols, categorical_cols, "date"), drop = FALSE]
```

GLMMPreprocessor  
================

An R6‐based utility designed for preprocessing and modeling binary response data with Generalized Linear Mixed Models (GLMMs). Supports automated feature engineering, train/test splitting, stepwise BIC‐guided model selection, and marginal effect visualization.

Core Capabilities
------------------

- Binary response checks and recoding (factor → 0/1)  
- Train/test split with optional undersampling  
- Automatic date conversion and single-level factor removal  
- Numeric standardization (center & scale)  
- Polynomial expansion (up to cubic by default)  
- Optional interaction terms with exclusion rules  
- Stepwise forward/backward selection via BIC (with `bic_tol` from previous best)  
- Mixed‐model fitting with `glmer()` (or fallback `glm()`)  
- Random‐effect variance pruning  
- Marginal effect plots (1D and 2D)  

Typical Usage
-------------

```r
library(GLMMPreprocessor)

# 1) Initialize with binary target
glmm <- GLMMPreprocessor$new(target = "target_feature")

# 2) Preprocess: split, scale, (optional) undersample
prep_out <- glmm$preprocess(
  df,
  train_pct         = 0.20,
  undersample_train = TRUE,
  random_state      = 663
)
# returns a list:
#  • X_train, y_train, X_test, y_test, X, y

# 3) Stepwise BIC GLMM search
fit_result <- glmm$stepwise_bic_glmm(
  X                 = prep_out$X,
  y                 = prep_out$y,
  random_effects    = c("group1", "group2"),
  no_interactions   = c("cov1", "cov2"),
  poly              = 3,
  bic_tol           = 0.5,          # allow up to 50% increase from previous best BIC
  max_add_steps     = 10000,
  max_drop_steps    = 10000,
  cores             = parallel::detectCores() - 1,
  re_var_threshold  = 1e-6,
  int_var_threshold = 1e-6
)
# returns an invisible list:
#  • fit, formula, n_add, n_drop, final_BIC

# 4) 1D marginal effect plot
p1 <- glmm$plot_marginal(
  model_fit    = fit_result$fit,
  df           = glmm$train_df,
  focal_var    = "cov3"
)
# returns a ggplot2 object

# 5) 2D joint marginal effect heatmap
p2 <- glmm$plot_marginal_interaction(
  model_fit = fit_result$fit,
  df        = glmm$train_df,
  vars      = c("cov3", "cov4")
)
# returns a ggplot2 object
```

Function Reference
------------------

### `remove_cont_multicollinearity(data, target, …)`

Filters numeric predictors by target correlation, pairwise correlation, and VIF.  
**Returns**: invisible list with  
- `pruned_data` — data.frame(target, retained numeric predictors)  
- `cluster_df` — data.frame(variable, cluster, removed)

### `remove_factor_multicollinearity(df, target_col, …)`

Clusters factor predictors by Cramer’s V and selects one per cluster.  
**Returns**: invisible list with  
- `pruned_data` — data.frame of retained factor predictors  
- `cluster_df` — data.frame(variable, cluster, NA_count, removed)

### `GLMMPreprocessor$new(target)`

Creates a new preprocessor object for binary response modeling.

### `$preprocess(df, train_pct, undersample_train, random_state)`

Splits, optionally undersamples, drops single‐level factors, converts dates, and scales numerics.  
**Returns**: list with  
- `X_train`, `y_train`, `X_test`, `y_test` — training/testing features and target  
- `X`, `y` — training features and target (for model search)

### `$stepwise_bic_glmm(X, y, random_effects, …)`

Performs forward/backward stepwise GLMM selection using BIC.  
**Returns** (invisible list) with  
- `fit` — final fitted model object  
- `formula` — model formula as string  
- `n_add`, `n_drop` — number of terms added/dropped  
- `final_BIC` — BIC of the selected model

### `$plot_marginal(model_fit, df, focal_var, …)`

Generates a 1D marginal effect plot (fixed effects only).  
**Returns**: ggplot2 object

### `$plot_marginal_interaction(model_fit, df, vars, …)`

Generates a 2D joint marginal effect heatmap.  
**Returns**: ggplot2 object

Dependencies
------------

Requires the following R packages:

```r
install.packages(c(
  "R6", "lme4", "lmerTest", "ggplot2", "tibble",
  "foreach", "doParallel", "car", "mice", "jsonlite",
  "httr", "randomForest", "caret", "pdp", "FactoMineR",
  "factoextra", "vcd", "pheatmap", "rpart.plot",
  "dplyr", "tidyr", "lubridate", "glmnet", "DescTools"
))
```

Contributing
------------

Contributions are welcome. Please open an issue or submit a pull request to suggest features or fixes.

License
-------

MIT © Alexander Chernikov
