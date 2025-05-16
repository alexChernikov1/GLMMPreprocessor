#' @importFrom R6 R6Class
#' @importFrom lme4 glmer glmerControl VarCorr
#' @importFrom stats glm BIC predict as.formula
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#' @importFrom ggplot2 scale_fill_viridis_c geom_raster geom_contour
#' @importFrom tibble tibble
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom car vif
#' @importFrom lmerTest lmer
#' @importFrom DescTools CramerV
#' @importFrom mice complete mice
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#' @importFrom vcd assocstats
#' @importFrom pheatmap pheatmap
#' @importFrom pdp partial
#' @importFrom caret train
#' @importFrom lubridate ymd
#' @importFrom dplyr select mutate arrange filter group_by summarise
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_biplot
#' @importFrom glmnet glmnet
#' @importFrom randomForest randomForest
#' @importFrom rpart.plot prp
NULL


`%||%` <- function(x, y) if (is.null(x)) y else x

GLMMPreprocessor <- R6Class(
  "GLMMPreprocessor",
  lock_objects = TRUE,
  public = list(
    target        = NULL,
    scaler_center = NULL,
    scaler_scale  = NULL,
    train_df      = NULL,
    test_df       = NULL,

    initialize = function(target = "ARRIVAL_DELAY_BINOM") self$target <- target,

    preprocess = function(df,
                          train_pct         = 0.05,
                          undersample_train = FALSE,
                          random_state      = 663) {

      ## ---- response checks --------------------------------------------------
      if (!self$target %in% names(df))
        stop("Target column ‘", self$target, "’ not found in data.")

      ycol <- df[[self$target]]
      if (is.factor(ycol)) {
        if (nlevels(ycol) != 2) stop("Target factor must have exactly 2 levels.")
        df[[self$target]] <- as.numeric(ycol == levels(ycol)[2L])
      } else {
        bad <- setdiff(unique(ycol), c(0, 1, NA))
        if (length(bad))
          stop("Target must be coded 0/1 (found values: ",
               paste(head(bad, 5L), collapse = ", "), ").")
      }

      ## ---- drop single‑level factors ---------------------------------------
      one_lvl <- names(df)[vapply(df, \(x) is.factor(x) && nlevels(x) < 2, logical(1))]
      if (length(one_lvl)) {
        message("Dropping single‑level factor(s): ", paste(one_lvl, collapse = ", "))
        df[one_lvl] <- NULL
      }

      ## ---- date → numeric ---------------------------------------------------
      date_cols <- names(df)[vapply(df, inherits, logical(1),
                                    what = c("Date", "POSIXct", "POSIXt"))]
      if (length(date_cols)) {
        for (dc in date_cols) df[[dc]] <- as.numeric(as.Date(df[[dc]]))
        message("Converted Date columns to numeric: ", paste(date_cols, collapse = ", "))
      }

      ## ---- train / test split ----------------------------------------------
      set.seed(random_state)
      idx0 <- which(df[[self$target]] == 0)
      idx1 <- which(df[[self$target]] == 1)
      train_idx <- c(sample(idx0, ceiling(length(idx0) * train_pct)),
                     sample(idx1, ceiling(length(idx1) * train_pct)))

      train_df <- df[train_idx, ]
      test_df  <- df[-train_idx, ]

      ## ---- OPTIONAL undersampling (train only) ------------------------------
      if (undersample_train) {
        tbl <- table(train_df[[self$target]])
        maj_class <- as.integer(names(tbl)[which.max(tbl)])
        min_class <- as.integer(names(tbl)[which.min(tbl)])

        maj_idx <- which(train_df[[self$target]] == maj_class)
        min_idx <- which(train_df[[self$target]] == min_class)

        set.seed(random_state)
        train_df <- train_df[c(min_idx, sample(maj_idx, length(min_idx))), , drop = FALSE]

        message("After undersampling: train rows = ", nrow(train_df),
                " (balanced 0/1)")
      } else {
        message("train rows: ", nrow(train_df), "; test rows: ", nrow(test_df))
      }

      ## ---- scaling ----------------------------------------------------------
      num_cols <- names(train_df)[vapply(train_df, is.numeric, logical(1)) &
                                    names(train_df) != self$target]

      self$scaler_center <- vapply(train_df[num_cols], mean, numeric(1), na.rm = TRUE)
      self$scaler_scale  <- vapply(train_df[num_cols], sd,   numeric(1), na.rm = TRUE)
      self$scaler_scale[self$scaler_scale == 0] <- 1

      scale_df <- function(d) {
        d[num_cols] <- sweep(
          sweep(d[num_cols], 2, self$scaler_center, `-`),
          2, self$scaler_scale, `/`
        )
        d
      }

      self$train_df <- train_df
      self$test_df  <- test_df

      list(
        X_train = scale_df(select(train_df, -all_of(self$target))),
        y_train = train_df[[self$target]],
        X_test  = scale_df(select(test_df , -all_of(self$target))),
        y_test  = test_df [[self$target]],
        X       = scale_df(select(train_df, -all_of(self$target))),
        y       = train_df[[self$target]]
      )
    },

    # ------------------------------------------------------------------------
    #  Step‑wise BIC GLMM search
    #     * new: any time a polynomial term I(x^k) is present,
    #            the *base* x is silently added to the formula so that
    #            it will live in fit@frame → avoids “object '<x>' not found”
    # ------------------------------------------------------------------------
    # ──────────────────────────────────────────────────────────────────────────
    stepwise_bic_glmm = function(
    X, y,
    random_effects,
    no_interactions   = NULL,
    poly              = 3,
    bic_tol           = 0.02,
    max_add_steps     = 200,
    max_drop_steps    = 200,
    cores             = max(1, parallel::detectCores() - 1),
    re_var_threshold  = 1e-6,
    int_var_threshold = 1e-6
    ) {
      target_col <- self$target

      safe <- function(nm) if (make.names(nm) == nm && !grepl(":", nm, fixed = TRUE))
        nm else paste0("`", nm, "`")
      sTerm <- function(t) {
        if (grepl("^I\\(", t))                       return(t)                 # KEEP ‘I(x^k)’
        if (grepl(":", t, fixed = TRUE))
          return(paste(vapply(strsplit(t, ":", fixed = TRUE)[[1]], safe, ""),
                       collapse = ":"))
        safe(t)
      }

      poly_base <- function(term) sub("^I\\(([^\\^]+)\\^.+\\)$", "\\1", term)

      fit_ic <- function(fx_terms, re_terms) {
        poly_terms   <- grep("^I\\(.+\\^\\d+\\)$", fx_terms, value = TRUE)
        base_needed  <- setdiff(vapply(poly_terms, poly_base, character(1)), fx_terms)
        fx_all       <- c(fx_terms, base_needed)

        rhs <- if (length(fx_all))
          paste(vapply(fx_all, sTerm, ""), collapse = " + ") else "1"
        if (length(re_terms))
          rhs <- paste(rhs,
                       paste0("(1|", vapply(re_terms, safe, ""), ")",
                              collapse = " + "),
                       sep = " + ")
        fmla <- as.formula(paste0(target_col, " ~ ", rhs))

        d <- X
        d[[target_col]] <- y
        d[re_terms]     <- lapply(d[re_terms], as.factor)

        fit <- suppressWarnings(tryCatch(
          if (length(re_terms))
            glmer(fmla, d, family = binomial(), nAGQ = 0,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl   = list(maxfun = 1e5)))
          else
            glm(fmla, d, family = binomial()),
          error = \(e) NULL))

        if (!is.null(fit) && length(fx_all) > 0) {
          ph <- try(predict(fit, type = "response"), silent = TRUE)
          if (inherits(ph, "try-error") || length(unique(ph > 0.5)) < 2)
            fit <- NULL
        }
        if (is.null(fit) || !is.finite(logLik(fit)))
          return(list(bic = Inf, fit = NULL))
        list(bic = BIC(fit), fit = fit)
      }

      ## ---- candidate fixed‑effect pool (as in previous version) -------------
      blocked   <- unique(random_effects)
      cont_cols <- names(X)[vapply(X, is.numeric, logical(1)) & !(names(X) %in% blocked)]
      cat_cols  <- setdiff(names(X), c(cont_cols, blocked))

      low_cont  <- cont_cols[sapply(cont_cols,
                                    \(v) var(X[[v]], na.rm = TRUE) < int_var_threshold)]
      cont_cols <- setdiff(cont_cols, low_cont)

      base_terms <- c(cont_cols, cat_cols)
      poly_terms <- if (poly > 1L)
        unlist(lapply(cont_cols,
                      \(cn) paste0("I(", cn, "^", 2:poly, ")")),
               FALSE)
      else character(0)

      allowed_int <- setdiff(base_terms, unique(no_interactions %||% character(0)))
      interaction_terms <- if (length(allowed_int) >= 2L) {
        pairs <- combn(allowed_int, 2, simplify = FALSE)
        unlist(lapply(pairs,
                      \(z) if (any(z %in% low_cont)) NULL else paste0(z[1], ":", z[2])),
               FALSE)
      } else character(0)

      candidate_fixed <- unique(c(base_terms, poly_terms, interaction_terms))
      message("Fixed‑effect search space: ", length(candidate_fixed), " terms")

      ## ---- forward + backward search (identical search logic) --------------
      current_fx <- character(0)
      current_re <- random_effects
      best       <- fit_ic(current_fx, current_re)
      if (is.infinite(best$bic)) {
        message("Intercept‑only mixed model failed → retrying without random effects …")
        current_re <- character(0); best <- fit_ic(current_fx, current_re)
      }
      if (is.infinite(best$bic)) stop("Even intercept‑only GLM failed – check response.")

      current_bic  <- best$bic
      current_fit  <- best$fit
      best_fx  <- current_fx; best_re <- current_re; best_bic <- current_bic; best_fit <- current_fit

      cl <- makeCluster(cores); registerDoParallel(cl); on.exit(stopCluster(cl))

      add_steps <- 0
      repeat {
        add_steps <- add_steps + 1
        if (add_steps > max_add_steps) break
        pool <- setdiff(candidate_fixed, current_fx); if (!length(pool)) break

        bic_tbl <- foreach(term = pool, .combine = rbind, .packages = "lme4") %dopar% {
          out <- fit_ic(c(current_fx, term), current_re)
          data.frame(term = term, bic = out$bic, stringsAsFactors = FALSE)
        }
        cand <- arrange(bic_tbl, bic) |> slice(1)

        if (cand$bic < best_bic) {
          best_bic <- cand$bic; best_fx <- c(current_fx, cand$term)
          best_re  <- current_re;       best_fit <- fit_ic(best_fx, best_re)$fit
        }

        if (cand$bic <= best_bic * (1 + bic_tol)) {
          message(" +  ", cand$term, "  (BIC = ", round(cand$bic, 2), ")")
          current_fx  <- c(current_fx, cand$term)
          current_bic <- cand$bic
          current_fit <- fit_ic(current_fx, current_re)$fit
        } else break
      }

      drop_steps <- 0
      repeat {
        drop_steps <- drop_steps + 1
        if (drop_steps > max_drop_steps || !length(current_fx)) break

        bic_tbl <- foreach(term = current_fx, .combine = rbind, .packages = "lme4") %dopar% {
          out <- fit_ic(setdiff(current_fx, term), current_re)
          data.frame(term = term, bic = out$bic, stringsAsFactors = FALSE)
        }
        cand <- arrange(bic_tbl, bic) |> slice(1)

        if (cand$bic < best_bic) {
          best_bic <- cand$bic; best_fx <- setdiff(current_fx, cand$term)
          best_re  <- current_re;       best_fit <- fit_ic(best_fx, best_re)$fit
        }

        if (cand$bic <= best_bic * (1 + bic_tol)) {
          message(" -  ", cand$term, "  (BIC = ", round(cand$bic, 2), ")")
          current_fx  <- setdiff(current_fx, cand$term)
          current_bic <- cand$bic
          current_fit <- fit_ic(current_fx, current_re)$fit
        } else break
      }

      repeat {
        if (!length(current_re)) break
        vc <- as.data.frame(VarCorr(current_fit))
        tiny <- vc$grp[vc$vcov < re_var_threshold & vc$var1 == "(Intercept)"]
        tiny <- intersect(tiny, current_re); if (!length(tiny)) break
        message("Dropping near‑zero RE: ", paste(tiny, collapse = ", "))
        current_re <- setdiff(current_re, tiny)
        tmp <- fit_ic(current_fx, current_re)
        current_fit <- tmp$fit; current_bic <- tmp$bic
        if (current_bic < best_bic) { best_bic <- current_bic; best_fx <- current_fx
        best_re  <- current_re;  best_fit <- current_fit }
      }

      rhs <- if (length(best_fx))
        paste(vapply(best_fx, sTerm, ""), collapse = " + ") else "1"
      if (length(best_re))
        rhs <- paste(rhs,
                     paste0("(1|", vapply(best_re, safe, ""), ")",
                            collapse = " + "),
                     sep = " + ")

      invisible(list(
        fit       = best_fit,
        formula   = paste0(target_col, " ~ ", rhs),
        n_add     = add_steps  - 1,
        n_drop    = drop_steps - 1,
        final_BIC = best_bic
      ))
    },



    plot_marginal = function(model_fit,
                             df,
                             focal_var     = "SCHEDULED_DEPARTURE",
                             n_points      = 200,
                             hold_distance = "mean") {
      if (!inherits(model_fit, "glmerMod") && !inherits(model_fit, "glm"))
        stop("`model_fit` must be a fitted model returned by stepwise_bic_glmm()")

      # grab scaling constants used in preprocess()
      c0 <- self$scaler_center[focal_var]
      s0 <- self$scaler_scale [focal_var]

      # if the interaction SCHEDULED_DEPARTURE:DISTANCE exists, we also need distance
      if ("DISTANCE" %in% names(self$scaler_center)) {
        cD <- self$scaler_center["DISTANCE"]
        sD <- self$scaler_scale ["DISTANCE"]
        dist_val <- switch(
          hold_distance,
          mean  = mean(df$DISTANCE, na.rm = TRUE),
          median = median(df$DISTANCE, na.rm = TRUE),
          stop("`hold_distance` must be 'mean' or 'median'")
        )
        dist_scaled <- (dist_val - cD) / sD
      }

      # range for the focal predictor on its ORIGINAL scale
      x_seq <- seq(min(df[[focal_var]]), max(df[[focal_var]]), length.out = n_points)

      # get factor levels from the model frame
      mf   <- model_fit@frame
      fac_defs <- lapply(mf, function(col)
        if (is.factor(col))
          levels(col))

      # construct new data grid
      newdat <- tibble::tibble(focal_orig = x_seq)
      newdat[[focal_var]] <- (x_seq - c0) / s0
      # add numeric covariates that the model needs but are not the focal var
      for (num in names(self$scaler_center)) {
        if (num == focal_var)
          next            # skip the one we’re varying
        if (!num %in% names(newdat)) {
          # not already added (e.g., DISTANCE)
          c_i <- self$scaler_center[num]
          s_i <- self$scaler_scale [num]
          val <- mean(df[[num]], na.rm = TRUE)   # or median(df[[num]])
          newdat[[num]] <- (val - c_i) / s_i     # scaled baseline
        }
      }

      # add distance if needed
      if (exists("dist_scaled", inherits = FALSE))
        newdat$DISTANCE <- dist_scaled

      # add one baseline value for every *factor* that appears in the model
      for (nm in names(fac_defs)) {
        levs <- fac_defs[[nm]]
        if (!is.null(levs) &&
            length(levs) > 0) {
          # <- skip numerics
          newdat[[nm]] <- factor(levs[1], levels = levs)   # length-1 ->recycled
        }
      }


      # predict FE probability
      newdat$prob <- predict(model_fit,
                             newdata = newdat,
                             type    = "response",
                             re.form = NA)  # population average (drops RE)

      # build ggplot
      p <- ggplot2::ggplot(newdat, ggplot2::aes(x = focal_orig, y = prob)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::labs(
          x = paste(focal_var, "(original scale)"),
          y = "P(delay = 1)  (fixed effects only)",
          title = paste("Marginal effect of", focal_var)
        ) +
        ggplot2::theme_minimal()

      return(p)
    },
    ## 2D marginal plot
    plot_marginal_interaction = function(model_fit,
                                         df,
                                         vars        = c("SCHEDULED_DEPARTURE", "DISTANCE"),
                                         n_grid      = 60,
                                         hold_method = c(mean = mean, median = median)) {
      stopifnot(length(vars) == 2)
      v1 <- vars[1]
      v2 <- vars[2]

      if (!inherits(model_fit, "glmerMod") &&
          !inherits(model_fit, "glm"))
        stop("model_fit must be the object returned by stepwise_bic_glmm()")

      ## helpers
      scale_one <- function(x, c_, s_)
        (x - c_) / s_

      # scaling constants
      c1 <- self$scaler_center[v1]
      s1 <- self$scaler_scale[v1]
      c2 <- self$scaler_center[v2]
      s2 <- self$scaler_scale[v2]

      # grid on original scale
      x1 <- seq(min(df[[v1]], na.rm = TRUE), max(df[[v1]], na.rm = TRUE), length.out = n_grid)
      x2 <- seq(min(df[[v2]], na.rm = TRUE), max(df[[v2]], na.rm = TRUE), length.out = n_grid)

      grid <- expand.grid(v1_orig = x1,
                          v2_orig = x2,
                          KEEP.OUT.ATTRS = FALSE)

      # numeric scaled variables
      grid[[v1]] <- scale_one(grid$v1_orig, c1, s1)
      grid[[v2]] <- scale_one(grid$v2_orig, c2, s2)

      # add numeric predictors at their mean
      for (num in names(self$scaler_center)) {
        if (num %in% vars)
          next
        if (!num %in% names(grid)) {
          c_ <- self$scaler_center[num]
          s_ <- self$scaler_scale[num]
          val <- hold_method$mean(df[[num]], na.rm = TRUE)   # default = mean
          grid[[num]] <- scale_one(val, c_, s_)
        }
      }

      # factor predictors at first level
      mf <- model_fit@frame
      fac_defs <- lapply(mf, function(col)
        if (is.factor(col))
          levels(col))
      for (nm in names(fac_defs)) {
        levs <- fac_defs[[nm]]
        if (!is.null(levs) &&
            length(levs) > 0 && !nm %in% names(grid)) {
          grid[[nm]] <- factor(levs[1], levels = levs)
        }
      }

      # predict fixed-effect probability
      grid$prob <- predict(model_fit,
                           newdata = grid,
                           type    = "response",
                           re.form = NA)

      # plot
      ggplot(grid, aes(x = v1_orig, y = v2_orig, z = prob)) +
        geom_raster(aes(fill = prob), interpolate = TRUE) +
        geom_contour(colour = "white", alpha = 0.6) +
        scale_fill_viridis_c(name = "P(delay = 1)") +
        labs(
          x = paste(v1, "(original)"),
          y = paste(v2, "(original)"),
          title = paste("Joint marginal effect:", v1, "×", v2)
        ) +
        theme_minimal()
    }

  )
)
