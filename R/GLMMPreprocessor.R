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

#' Null coalescing operator `%||%`
#'
#' Returns `y` if `x` is `NULL`, otherwise returns `x`.
#'
#' @param x An object that may be NULL.
#' @param y A default to return if `x` is NULL.
#'
#' @return `x` if not NULL, otherwise `y`.
#' @export
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @export
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

#' @export
remove_cont_multicollinearity <- function(
    data,
    target,
    target_cor_threshold = 0.7,
    cor_threshold        = 0.70,
    vif_threshold        = 5,
    verbose              = TRUE,
    drop_cols            = character(),
    keep_cols            = character(),
    draw_corr            = FALSE
) {
  # --------------------------------------------------------
  # 0) Helper functions
  # --------------------------------------------------------
  uf_init <- function(items) {
    parent <- setNames(seq_along(items), items)
    rank   <- setNames(rep(0, length(items)), items)
    list(parent = parent, rank = rank)
  }
  uf_find <- function(x, uf) {
    if (uf$parent[x] != x) uf$parent[x] <- uf_find(uf$parent[x], uf)
    uf$parent[x]
  }
  uf_union <- function(a, b, uf) {
    rootA <- uf_find(a, uf); rootB <- uf_find(b, uf)
    if (rootA != rootB) {
      if      (uf$rank[rootA] < uf$rank[rootB]) uf$parent[rootA] <- rootB
      else if (uf$rank[rootA] > uf$rank[rootB]) uf$parent[rootB] <- rootA
      else {
        uf$parent[rootB] <- rootA
        uf$rank[rootA]  <- uf$rank[rootA] + 1
      }
    }
    uf
  }
  na_count <- function(x) sum(is.na(x))

  # --------------------------------------------------------
  # 1) Validate & prepare data
  # --------------------------------------------------------
  if (!target %in% names(data)) {
    stop("Target column '", target, "' not found in data.")
  }
  if (length(drop_cols) > 0) {
    data <- data[, setdiff(names(data), drop_cols), drop = FALSE]
    if (verbose) message("Dropped user-specified columns: ", paste(drop_cols, collapse = ", "))
  }
  #--------------------------------------------------
  # 1.5) Warn for columns that are neither numeric nor factor
  #--------------------------------------------------
  not_num_or_factor <- names(data)[sapply(data, function(x) !is.numeric(x) && !is.factor(x))]
  all_ignored <- union(drop_cols, target)
  not_num_or_factor <- setdiff(not_num_or_factor, all_ignored)
  if (length(not_num_or_factor) > 0 && verbose) {
    warning(paste0("Variables that are neither numeric nor factor: ", paste(not_num_or_factor, collapse = ", "), "\n"))
  }

  target_vec <- data[[target]]
  df_preds   <- data[, setdiff(names(data), target), drop = FALSE]

  # -- New: Remove non-numeric predictors --
  non_numeric_cols <- names(df_preds)[!sapply(df_preds, is.numeric)]
  if (length(non_numeric_cols) > 0 && verbose) {
    warning(paste0("Pruning non-numeric predictors: ", paste(non_numeric_cols, collapse = ", "), "\n"))
  }
  df_preds <- df_preds[, sapply(df_preds, is.numeric), drop = FALSE]


  # keep only numeric with variation
  df_preds   <- df_preds[, sapply(df_preds, is.numeric), drop = FALSE]
  df_preds   <- df_preds[, sapply(df_preds, function(x) {
    !all(is.na(x)) && length(unique(na.omit(x))) > 1
  }), drop = FALSE]

  # --------------------------------------------------------
  # 2) Target correlation filtering
  # --------------------------------------------------------
  if (ncol(df_preds) > 0) {
    tcorrs <- sapply(df_preds, function(x) cor(x, target_vec, use = "complete.obs"))
    drop_target <- names(tcorrs)[
      abs(tcorrs) > target_cor_threshold &
        !names(tcorrs) %in% keep_cols &
        !is.na(tcorrs)
    ]
    if (length(drop_target) > 0) {
      if (verbose) message(
        "Dropping ", length(drop_target),
        " predictors with |corr| > ", target_cor_threshold,
        " vs target: ", paste(drop_target, collapse = ", ")
      )
      df_preds <- df_preds[, setdiff(names(df_preds), drop_target), drop = FALSE]
    }
  }
  if (ncol(df_preds) < 2) {
    if (verbose) message("Not enough predictors after target filtering.")
    return(list(pruned_data = data.frame(setNames(list(target_vec), target)),
                cluster_df  = NULL))
  }

  all_cols       <- colnames(df_preds)
  uf             <- uf_init(all_cols)
  removed_status <- setNames(rep(FALSE, length(all_cols)), all_cols)

  # --------------------------------------------------------
  # 3) Correlation filtering among predictors
  # --------------------------------------------------------
  corr_mat <- suppressWarnings(cor(df_preds, use = "pairwise.complete.obs"))
  high_cor <- which(abs(corr_mat) > cor_threshold & upper.tri(corr_mat), arr.ind = TRUE)
  while (nrow(high_cor) > 0) {
    i      <- high_cor[1,1]; j <- high_cor[1,2]
    c1     <- colnames(corr_mat)[i]; c2 <- colnames(corr_mat)[j]
    na1    <- na_count(df_preds[[c1]]); na2 <- na_count(df_preds[[c2]])
    in1    <- c1 %in% keep_cols; in2 <- c2 %in% keep_cols
    if      (in1 && !in2) { drop_col <- c2; keep_col <- c1 }
    else if (!in1 && in2) { drop_col <- c1; keep_col <- c2 }
    else if (na1 > na2)    { drop_col <- c1; keep_col <- c2 }
    else if (na2 > na1)    { drop_col <- c2; keep_col <- c1 }
    else {
      drop_col <- c2; keep_col <- c1
      if (in1 && in2) warning(
        "Tie for keep_cols in '", c1, "' & '", c2,
        "'; dropping '", drop_col, "'."
      )
    }
    uf <- uf_union(drop_col, keep_col, uf)
    removed_status[drop_col] <- TRUE
    if (verbose) message(
      "Dropping '", drop_col, "' corr> ", cor_threshold,
      " with '", keep_col, "'."
    )
    df_preds[[drop_col]] <- NULL
    if (ncol(df_preds) <= 1) break
    corr_mat  <- suppressWarnings(cor(df_preds, use = "pairwise.complete.obs"))
    high_cor  <- which(abs(corr_mat) > cor_threshold & upper.tri(corr_mat), arr.ind = TRUE)
  }

  # --------------------------------------------------------
  # 4) VIF filtering
  # --------------------------------------------------------
  if (ncol(df_preds) > 1) {
    set.seed(123); fake_y <- rnorm(nrow(df_preds))
    repeat {
      form     <- as.formula(paste("fake_y ~", paste(names(df_preds), collapse = " + ")))
      fit_vif  <- lm(form, data = df_preds)
      all_vifs <- car::vif(fit_vif)
      max_vif  <- max(all_vifs)
      if (max_vif < vif_threshold) break
      worst_col  <- names(which.max(all_vifs))
      other_cols <- setdiff(names(df_preds), worst_col)
      if (length(other_cols) == 0) {
        removed_status[worst_col] <- TRUE
        if (verbose) message("Dropping '", worst_col, "' (VIF=", round(max_vif,2), ") no partner.")
        df_preds[[worst_col]] <- NULL; break
      }
      w_corrs     <- sapply(other_cols, function(cc)
        cor(df_preds[[worst_col]], df_preds[[cc]], use = "pairwise.complete.obs"))
      partner_col <- names(which.max(abs(w_corrs)))
      uf <- uf_union(worst_col, partner_col, uf)
      removed_status[worst_col] <- TRUE
      if (verbose) message(
        "Dropping '", worst_col, "' (VIF=", round(max_vif,2),
        ") corr with '", partner_col, "'."
      )
      df_preds[[worst_col]] <- NULL
      if (ncol(df_preds) < 2) break
    }
  }

  # --------------------------------------------------------
  # 5) Cluster summary & final output
  # --------------------------------------------------------
  final_cols  <- names(uf$parent)
  final_roots <- sapply(final_cols, uf_find, uf = uf)
  clusters    <- as.integer(factor(final_roots))
  df_clusters <- data.frame(
    variable = final_cols,
    cluster  = clusters,
    removed  = removed_status[final_cols],
    stringsAsFactors = FALSE
  )

  if (verbose) {
    # print per‐cluster chosen vs others
    cluster_summary <- df_clusters %>%
      group_by(cluster) %>%
      summarise(
        chosen = variable[!removed][1],
        others = paste(variable[removed], collapse = ", "),
        .groups = "drop"
      )
    message("\nFeature Group summary (chosen / others):")
    apply(cluster_summary, 1, function(row) {
      msg <- paste0("  Feature Group ", row["cluster"], ": ",
                    row["chosen"],
                    if (nzchar(row["others"])) paste0("  /  ", row["others"]) else "")
      message(msg)
    })
    # also show full table
    print(df_clusters)
  }

  pruned_data <- cbind(
    setNames(data.frame(target_vec), target),
    df_preds
  )

  if (draw_corr) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      warning("Install 'pheatmap' to draw heatmap.")
    } else {
      cm <- cor(pruned_data, use = "pairwise.complete.obs")
      pheatmap::pheatmap(
        cm,
        color           = colorRampPalette(c("blue", "white", "red"))(50),
        main            = "Correlation Heatmap",
        display_numbers = TRUE,
        angle_col       = 90
      )
    }
  }


  #--------------------------------------------------
  # 6) Final summary & return
  #--------------------------------------------------
  if (verbose) {
    original_cols <- setdiff(names(data), c(target, drop_cols))
    numeric_cols <- names(data[, original_cols])[sapply(data[, original_cols], is.numeric)]
    cat("\nBefore:", length(numeric_cols), "numeric cols;",
        "After:", ncol(df_preds), "retained.\n")
    cat("Retained columns:\n")
    print(colnames(df_preds))
  }

  invisible(list(pruned_data = pruned_data, cluster_df = df_clusters))
}

#' @export
remove_factor_multicollinearity <- function(
    df,
    target_col     = "Likelihood_to_Recommend__Relationship",
    drop_cols      = "date",
    keep_cols      = character(),
    k              = 5,
    verbose        = TRUE
) {
  #--------------------------------------------------
  # 0) Prep: require packages
  #--------------------------------------------------
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package.")
  }
  if (!requireNamespace("vcd", quietly = TRUE)) {
    stop("Please install the 'vcd' package (for assocstats).")
  }

  #--------------------------------------------------
  # 0.1) Drop user-specified columns
  #--------------------------------------------------
  if (length(drop_cols) > 0 && verbose) {
    message("Dropped user-specified columns: ", paste(drop_cols, collapse = ", "))
  }
  df <- df[, setdiff(names(df), drop_cols), drop = FALSE]

  #--------------------------------------------------
  # 0.5) Prune numeric predictors (we only want factors)
  #--------------------------------------------------
  numeric_preds <- setdiff(
    names(df)[sapply(df, is.numeric)],
    target_col
  )
  if (length(numeric_preds) > 0 && verbose) {
    message(
      "Pruning numeric predictors: ",
      paste(numeric_preds, collapse = ", ")
    )
  }
  df <- df[, setdiff(names(df), numeric_preds), drop = FALSE]

  #--------------------------------------------------
  # 1) Create factor-only data (exclude numeric target)
  #--------------------------------------------------
  model_data  <- df %>% dplyr::select(all_of(target_col), where(~ !is.numeric(.)))
  factor_data <- model_data %>% dplyr::select(-all_of(target_col))

  #--------------------------------------------------
  # 1.5) Prune any remaining non-factor predictors
  #--------------------------------------------------
  non_factor_cols <- names(factor_data)[!sapply(factor_data, is.factor)]
  if (length(non_factor_cols) > 0 && verbose) {
    message(
      "Pruning non-factor predictors: ",
      paste(non_factor_cols, collapse = ", ")
    )
  }
  factor_data <- factor_data[, sapply(factor_data, is.factor), drop = FALSE]

  #--------------------------------------------------
  # 2) Clean factors: char → NA → factor
  #--------------------------------------------------
  factor_data <- factor_data %>%
    dplyr::mutate(across(where(is.factor),    as.character)) %>%
    dplyr::mutate(across(everything(),        ~ na_if(., "NA"))) %>%
    dplyr::mutate(across(everything(),        ~ na_if(., ""))) %>%
    dplyr::mutate(across(everything(),        as.factor))

  #--------------------------------------------------
  # 3) Compute pairwise Cramer's V
  #--------------------------------------------------
  cramer_v <- function(x, y) {
    tbl <- table(x, y)
    vcd::assocstats(tbl)$cramer
  }
  vars       <- colnames(factor_data)
  n_vars     <- length(vars)
  cramer_mat <- matrix(0, n_vars, n_vars, dimnames = list(vars, vars))
  for (i in seq_len(n_vars - 1)) {
    for (j in seq(i + 1, n_vars)) {
      v <- cramer_v(factor_data[[i]], factor_data[[j]])
      cramer_mat[i, j] <- v
      cramer_mat[j, i] <- v
    }
  }

  #--------------------------------------------------
  # 4) Dendrogram BEFORE grouping
  #--------------------------------------------------
  d_mat <- as.dist(1 - cramer_mat)
  hc    <- hclust(d_mat, method = "complete")
  plot(hc, main = "Dendrogram BEFORE Factor Grouping", xlab = "", sub = "")

  #--------------------------------------------------
  # 5) Cut into k clusters
  #--------------------------------------------------
  clusters <- cutree(hc, k = k)

  #--------------------------------------------------
  # 6) Show cluster membership + NA counts
  #--------------------------------------------------
  df_clusters <- data.frame(
    variable = names(clusters),
    cluster  = as.character(clusters),
    stringsAsFactors = FALSE
  )
  df_clusters$NA_count <- vapply(
    df_clusters$variable,
    function(v) sum(is.na(factor_data[[v]])), integer(1)
  )
  df_clusters <- df_clusters %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(NA_count, .by_group = TRUE)

  if (verbose) {
    cat("\n--- Cluster Membership and NA Counts ---\n")
    print(df_clusters, n = Inf)
  }

  #--------------------------------------------------
  # 7) Select one “best” per cluster
  #--------------------------------------------------
  selected_cols <- vapply(
    unique(clusters),
    FUN.VALUE = character(1),
    function(cl) {
      members   <- names(clusters[clusters == cl])
      na_counts <- setNames(
        vapply(members, function(m) sum(is.na(factor_data[[m]])), integer(1)),
        members
      )
      kc_in_cl <- intersect(members, keep_cols)
      if (length(kc_in_cl) > 0) {
        return(kc_in_cl[which.min(na_counts[kc_in_cl])][1])
      }
      names(which.min(na_counts))[1]
    }
  )
  factor_data_reduced <- factor_data[, unname(selected_cols), drop = FALSE]

  #--------------------------------------------------
  # 7.5) Feature-Group summary
  #--------------------------------------------------
  if (verbose) {
    df_clusters$removed <- !df_clusters$variable %in% unname(selected_cols)
    cat("\nFeature Group summary (chosen / others):\n")
    for (cl in sort(unique(df_clusters$cluster))) {
      in_cl   <- df_clusters$variable[df_clusters$cluster == cl]
      kept    <- in_cl[!df_clusters$removed[df_clusters$cluster == cl]]
      dropped <- setdiff(in_cl, kept)
      cat(sprintf(
        "  Feature Group %s: %s%s\n",
        cl,
        kept,
        if (length(dropped)) paste0("  /  ", paste(dropped, collapse = ", ")) else ""
      ))
    }
  }

  #--------------------------------------------------
  # 8) Dendrogram AFTER grouping
  #--------------------------------------------------
  if (ncol(factor_data_reduced) > 1) {
    vars2 <- colnames(factor_data_reduced)
    cm2   <- matrix(0, length(vars2), length(vars2), dimnames = list(vars2, vars2))
    for (i in seq_len(length(vars2) - 1)) {
      for (j in seq(i + 1, length(vars2))) {
        v2 <- cramer_v(factor_data_reduced[[i]], factor_data_reduced[[j]])
        cm2[i, j] <- v2
        cm2[j, i] <- v2
      }
    }
    hc2 <- hclust(as.dist(1 - cm2), method = "complete")
    plot(hc2, main = "Dendrogram AFTER Factor Grouping Selection", xlab = "", sub = "")
  } else {
    message("\nOnly one factor retained; no second dendrogram.")
  }

  #--------------------------------------------------
  # 9) Final summary & return
  #--------------------------------------------------
  if (verbose) {
    cat(sprintf(
      "\nBefore: %d factor cols; After: %d retained.\n",
      ncol(factor_data), ncol(factor_data_reduced)
    ))
    cat("Retained columns:\n")
    print(colnames(factor_data_reduced))
  }

  invisible(list(pruned_data = factor_data_reduced, cluster_df = df_clusters))
}






