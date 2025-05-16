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

GLMMPreprocessor is an R6-based package for preprocessing data and performing stepwise BIC-guided feature selection for Generalized Linear Mixed Models (GLMMs). It includes tools for data splitting, scaling, polynomial and interaction term modeling, and marginal effect visualization.

Features
--------

- Preprocessing of binary response data for GLMMs
- Train/test split with optional class balancing
- Scaling of numeric variables
- Handles date and factor columns automatically
- Stepwise BIC feature selection using `glmer` or `glm`
- Support for polynomial terms and interactions
- Visualizations for 1D and 2D marginal effects

Quick Start
-----------

    library(GLMMPreprocessor)

    # Initialize preprocessor
    glmm <- GLMMPreprocessor$new(target = "ARRIVAL_DELAY_BINOM")

    # Preprocess data
    out <- glmm$preprocess(my_data)

    # Run stepwise BIC-guided GLMM search
    fit_result <- glmm$stepwise_bic_glmm(
      X = out$X, y = out$y,
      random_effects = c("AIRLINE", "ROUTE")
    )

    # Plot marginal effect of a focal predictor
    glmm$plot_marginal(
      model_fit = fit_result$fit,
      df = glmm$train_df,
      focal_var = "SCHEDULED_DEPARTURE"
    )

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

