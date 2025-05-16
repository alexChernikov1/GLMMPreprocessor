GLMMPreprocessor
=================

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

Installation
------------

To install the development version directly from GitHub:

    # install.packages("remotes") if not already installed
    remotes::install_github("YOUR_USERNAME/GLMMPreprocessor")

Replace "YOUR_USERNAME" with your actual GitHub username.

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
