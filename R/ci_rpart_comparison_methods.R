# Custom rpart concentration-index tree methods extracted from
# R/RDC_analysis_CI_Tree_Trials.Rmd and adapted for the cleaned report data.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

ci_rpart_weighted_ci <- function(y, rank_var, wt) {
  keep <- is.finite(y) & is.finite(rank_var) & is.finite(wt) & wt > 0
  y <- y[keep]
  rank_var <- rank_var[keep]
  wt <- wt[keep]

  if (length(y) < 2L || sum(wt) <= 0) {
    return(0)
  }

  ord <- order(rank_var)
  y <- y[ord]
  wt <- wt[ord]
  W <- sum(wt)
  R <- (cumsum(wt) - 0.5 * wt) / W
  mu <- sum(wt * y) / W

  if (!is.finite(mu) || abs(mu) < .Machine$double.eps) {
    return(0)
  }

  out <- 2 / mu * sum(wt * (y - mu) * (R - 0.5)) / W
  if (is.finite(out)) out else 0
}

make_simple_ci_rpart_method <- function() {
  ci_init <- function(y, offset, parms = NULL, wt) {
    if (!is.matrix(y) || ncol(y) != 2L) {
      stop("Response must be cbind(outcome, rank_var).", call. = FALSE)
    }

    sfun <- function(yval, dev, wt, ylevel = NULL, digits = getOption("digits"), ...) {
      if (is.matrix(yval)) {
        paste0(
          "mean=", format(signif(yval[, 1], digits)),
          ", CI=", format(signif(yval[, 2], digits)),
          ", n=", round(wt)
        )
      } else {
        paste0(
          "mean=", format(signif(yval[1], digits)),
          ", CI=", format(signif(yval[2], digits)),
          ", n=", round(wt)
        )
      }
    }

    tfun <- function(yval, dev = NULL, wt = NULL, ylevel = NULL,
                     digits = getOption("digits") - 3L, n = NULL,
                     use.n = FALSE, ...) {
      mu <- if (is.matrix(yval)) yval[, 1] else yval[1]
      ci <- if (is.matrix(yval)) yval[, 2] else yval[2]
      lab <- paste0(format(signif(mu, digits)), "\nCI=", format(signif(ci, digits)))
      if (isTRUE(use.n) && !is.null(n)) {
        lab <- paste0(lab, "\nn=", n)
      }
      lab
    }

    environment(sfun) <- .GlobalEnv
    environment(tfun) <- .GlobalEnv

    list(
      y = y,
      parms = parms,
      numresp = 2L,
      numy = 2L,
      summary = sfun,
      text = tfun
    )
  }

  ci_eval <- function(y, wt, parms = NULL) {
    outcome <- y[, 1]
    rnk <- y[, 2]
    ci <- ci_rpart_weighted_ci(outcome, rnk, wt)
    W <- sum(wt)

    list(
      label = c(mean_y = sum(wt * outcome) / W, CI = ci),
      deviance = abs(ci) * W
    )
  }

  ci_split <- function(y, wt, x, parms = NULL, continuous) {
    outcome <- y[, 1]
    rnk <- y[, 2]
    n <- length(outcome)
    parent_dev <- abs(ci_rpart_weighted_ci(outcome, rnk, wt)) * sum(wt)

    if (continuous) {
      goodness <- numeric(n - 1L)
      direction <- numeric(n - 1L)

      for (i in seq_len(n - 1L)) {
        if (x[i] >= x[i + 1L]) {
          next
        }

        left <- seq_len(i)
        right <- (i + 1L):n
        wL <- sum(wt[left])
        wR <- sum(wt[right])

        if (wL == 0 || wR == 0) {
          next
        }

        dL <- abs(ci_rpart_weighted_ci(outcome[left], rnk[left], wt[left])) * wL
        dR <- abs(ci_rpart_weighted_ci(outcome[right], rnk[right], wt[right])) * wR
        goodness[i] <- parent_dev - (dL + dR)
        direction[i] <- sign(sum(wt[right] * outcome[right]) / wR -
          sum(wt[left] * outcome[left]) / wL)

        if (direction[i] == 0) {
          direction[i] <- 1
        }
      }

      list(goodness = pmax(goodness, 0), direction = direction)
    } else {
      lvls <- sort(unique(x))
      means <- vapply(lvls, function(level) {
        idx <- x == level
        if (sum(wt[idx]) == 0) {
          return(0)
        }
        sum(wt[idx] * outcome[idx]) / sum(wt[idx])
      }, numeric(1))
      ord <- order(means)
      xo <- match(x, lvls[ord])
      goodness <- numeric(length(lvls) - 1L)

      for (k in seq_len(length(lvls) - 1L)) {
        left <- xo <= k
        right <- !left
        wL <- sum(wt[left])
        wR <- sum(wt[right])

        if (wL == 0 || wR == 0) {
          next
        }

        dL <- abs(ci_rpart_weighted_ci(outcome[left], rnk[left], wt[left])) * wL
        dR <- abs(ci_rpart_weighted_ci(outcome[right], rnk[right], wt[right])) * wR
        goodness[k] <- parent_dev - (dL + dR)
      }

      list(goodness = pmax(goodness, 0), direction = ord)
    }
  }

  list(init = ci_init, eval = ci_eval, split = ci_split)
}

make_ci_rpart_method <- function(
    loss = c("absolute", "squared", "positive", "negative"),
    rerank_within_node = TRUE,
    maxcat_enum = 12L,
    eps = 1e-12,
    allow_negative = FALSE) {
  loss <- match.arg(loss)
  loss_codes <- c(absolute = 1, squared = 2, positive = 3, negative = 4)

  p0 <- c(
    loss = unname(loss_codes[loss]),
    rerank = as.numeric(rerank_within_node),
    maxcat = as.numeric(maxcat_enum),
    eps = as.numeric(eps),
    allow_negative = as.numeric(allow_negative)
  )

  normalise_parms <- function(parms = NULL) {
    p <- p0

    if (!is.null(parms) && length(parms) > 0L) {
      if (is.list(parms)) {
        if (!is.null(parms$loss) && is.character(parms$loss)) {
          parms$loss <- unname(loss_codes[match.arg(parms$loss, names(loss_codes))])
        }
        if (!is.null(parms$rerank_within_node)) {
          parms$rerank <- as.numeric(parms$rerank_within_node)
          parms$rerank_within_node <- NULL
        }
        if (!is.null(parms$maxcat_enum)) {
          parms$maxcat <- as.numeric(parms$maxcat_enum)
          parms$maxcat_enum <- NULL
        }
        u <- unlist(parms, use.names = TRUE)
      } else {
        u <- parms
      }

      if (!is.null(names(u))) {
        keep <- intersect(names(u), names(p))
        if (length(keep)) {
          p[keep] <- suppressWarnings(as.numeric(u[keep]))
        }
      }
    }

    if (!is.finite(p["loss"]) || !(round(p["loss"]) %in% 1:4)) {
      p["loss"] <- p0["loss"]
    }
    p["loss"] <- round(p["loss"])
    p["rerank"] <- as.numeric(is.finite(p["rerank"]) && p["rerank"] > 0)
    if (!is.finite(p["maxcat"]) || p["maxcat"] < 2) {
      p["maxcat"] <- p0["maxcat"]
    }
    p["maxcat"] <- max(2, round(p["maxcat"]))
    if (!is.finite(p["eps"]) || p["eps"] <= 0) {
      p["eps"] <- p0["eps"]
    }
    p["allow_negative"] <- as.numeric(is.finite(p["allow_negative"]) && p["allow_negative"] > 0)
    p
  }

  weighted_midrank01 <- function(ses, wt) {
    n <- length(ses)
    if (n == 0L) {
      return(numeric())
    }

    o <- order(ses)
    ses_o <- ses[o]
    wt_o <- wt[o]
    W <- sum(wt_o)

    if (!is.finite(W) || W <= 0) {
      return(rep(NA_real_, n))
    }

    starts <- c(1L, which(diff(ses_o) != 0) + 1L)
    ends <- c(starts[-1L] - 1L, n)
    cw <- c(0, cumsum(wt_o))
    r_o <- numeric(n)

    for (j in seq_along(starts)) {
      r_o[starts[j]:ends[j]] <- ((cw[starts[j]] + cw[ends[j] + 1L]) / 2) / W
    }

    r <- numeric(n)
    r[o] <- r_o
    r
  }

  concentration_index <- function(y, ses, wt, parms) {
    keep <- is.finite(y) & is.finite(ses) & is.finite(wt) & wt > 0
    y <- y[keep]
    ses <- ses[keep]
    wt <- wt[keep]

    if (length(y) <= 1L || sum(wt) <= 0) {
      return(0)
    }

    W <- sum(wt)
    mu <- sum(wt * y) / W

    if (!is.finite(mu) || abs(mu) < parms["eps"]) {
      return(0)
    }

    R <- if (parms["rerank"] > 0.5) weighted_midrank01(ses, wt) else ses
    Rbar <- sum(wt * R) / W
    ci <- 2 * sum(wt * (y - mu) * (R - Rbar)) / W / mu
    if (is.finite(ci)) ci else 0
  }

  loss_value <- function(ci, parms) {
    code <- as.integer(round(parms["loss"]))
    if (code == 1L) return(abs(ci))
    if (code == 2L) return(ci^2)
    if (code == 3L) return(max(ci, 0))
    if (code == 4L) return(max(-ci, 0))
    abs(ci)
  }

  node_risk <- function(y, wt, parms) {
    W <- sum(wt)
    if (!is.finite(W) || W <= 0) {
      return(0)
    }
    ci <- concentration_index(y[, 1], y[, 2], wt, parms)
    W * loss_value(ci, parms)
  }

  init_ci <- function(y, offset, parms = NULL, wt) {
    parms <- normalise_parms(parms)

    if (!is.matrix(y) || ncol(y) != 2L) {
      stop("Use a two-column response: cbind(outcome, ses_or_rank) ~ predictors.", call. = FALSE)
    }
    if (length(offset)) {
      warning("Offset ignored for this CI-based rpart method.", call. = FALSE)
    }

    y <- cbind(outcome = as.numeric(y[, 1]), ses = as.numeric(y[, 2]))
    if (any(!is.finite(y))) {
      stop("Outcome and SES/rank response columns must be finite.", call. = FALSE)
    }
    if (!as.logical(parms["allow_negative"]) && any(y[, 1] < 0)) {
      stop("Negative outcome values found; set allow_negative = TRUE only if intentional.", call. = FALSE)
    }

    fmt <- function(z, digits, nsmall = NULL) {
      digits <- max(1L, as.integer(digits))
      if (is.null(nsmall)) {
        format(signif(as.numeric(z), digits), trim = TRUE)
      } else {
        format(signif(as.numeric(z), digits), trim = TRUE, nsmall = nsmall)
      }
    }

    as_yval_matrix <- function(z) {
      if (is.null(dim(z))) matrix(z, nrow = 1L) else as.matrix(z)
    }

    pad3 <- function(z) {
      z <- as.numeric(z)
      c(z, rep(NA_real_, max(0L, 3L - length(z))))[1:3]
    }

    sfun <- function(yval, dev, wt, ylevel = NULL, digits = getOption("digits"),
                     nsmall = NULL, ...) {
      M <- as_yval_matrix(yval)
      out <- character(nrow(M))

      for (i in seq_len(nrow(M))) {
        yy <- pad3(M[i, ])
        dev_wt <- if (length(dev) >= i && length(wt) >= i && is.finite(wt[i]) && wt[i] > 0) {
          dev[i] / wt[i]
        } else {
          NA_real_
        }
        out[i] <- paste0(
          "mean=", fmt(yy[1], digits, nsmall),
          ", CI=", fmt(yy[2], digits, nsmall),
          ", loss=", fmt(yy[3], digits, nsmall),
          ", dev/wt=", fmt(dev_wt, digits, nsmall)
        )
      }
      out
    }

    pfun <- function(yval, ylevel = NULL, digits = getOption("digits"),
                     nsmall = NULL, ...) {
      M <- as_yval_matrix(yval)
      out <- character(nrow(M))

      for (i in seq_len(nrow(M))) {
        yy <- pad3(M[i, ])
        out[i] <- paste0(
          "mean=", fmt(yy[1], digits, nsmall),
          ", CI=", fmt(yy[2], digits, nsmall),
          ", loss=", fmt(yy[3], digits, nsmall)
        )
      }
      out
    }

    tfun <- function(yval, dev = NULL, wt = NULL, ylevel = NULL,
                     digits = getOption("digits") - 3L, n = NULL,
                     use.n = FALSE, nsmall = NULL, ...) {
      M <- as_yval_matrix(yval)
      out <- character(nrow(M))

      for (i in seq_len(nrow(M))) {
        yy <- pad3(M[i, ])
        out[i] <- paste0(
          "CI=", fmt(yy[2], digits, nsmall),
          "\nmean=", fmt(yy[1], digits, nsmall)
        )
        if (isTRUE(use.n) && length(n) >= i) {
          out[i] <- paste0(out[i], "\nn=", n[i])
        }
      }
      out
    }

    label_env <- environment()
    environment(sfun) <- label_env
    environment(pfun) <- label_env
    environment(tfun) <- label_env

    list(
      y = y,
      parms = parms,
      numresp = 3L,
      numy = 2L,
      summary = sfun,
      print = pfun,
      text = tfun
    )
  }

  eval_ci <- function(y, wt, parms) {
    parms <- normalise_parms(parms)
    W <- sum(wt)
    mean_y <- if (W > 0) sum(wt * y[, 1]) / W else NA_real_
    ci <- concentration_index(y[, 1], y[, 2], wt, parms)
    loss <- loss_value(ci, parms)
    list(label = c(mean_y = mean_y, CI = ci, loss = loss), deviance = W * loss)
  }

  split_ci <- function(y, wt, x, parms, continuous) {
    parms <- normalise_parms(parms)
    n <- nrow(y)
    parent_risk <- node_risk(y, wt, parms)

    if (!is.finite(parent_risk) || parent_risk <= parms["eps"] || n <= 1L) {
      if (continuous) {
        return(list(goodness = rep(0, max(0L, n - 1L)), direction = rep(-1, max(0L, n - 1L))))
      }
      cats <- sort(unique(x))
      return(list(goodness = rep(0, max(0L, length(cats) - 1L)), direction = cats))
    }

    split_gain <- function(left) {
      if (!any(left) || all(left)) {
        return(0)
      }

      left_risk <- node_risk(y[left, , drop = FALSE], wt[left], parms)
      right_risk <- node_risk(y[!left, , drop = FALSE], wt[!left], parms)
      max(0, (parent_risk - left_risk - right_risk) / parent_risk)
    }

    if (continuous) {
      goodness <- numeric(n - 1L)
      for (i in seq_len(n - 1L)) {
        if (x[i] < x[i + 1L]) {
          goodness[i] <- split_gain(seq_len(n) <= i)
        }
      }
      return(list(goodness = goodness, direction = rep(-1, n - 1L)))
    }

    cats <- sort(unique(x))
    m <- length(cats)
    if (m <= 1L) {
      return(list(goodness = numeric(), direction = cats))
    }

    goodness_for_order <- function(ord) {
      g <- numeric(m - 1L)
      for (j in seq_len(m - 1L)) {
        g[j] <- split_gain(x %in% ord[seq_len(j)])
      }
      g
    }

    if (m <= parms["maxcat"]) {
      best_g <- 0
      best_left <- cats[1L]
      max_mask <- 2^(m - 1L) - 2L

      for (mask in 0:max_mask) {
        bits <- as.logical(as.integer(intToBits(mask))[seq_len(m - 1L)])
        left_cats <- cats[c(TRUE, bits)]
        g <- split_gain(x %in% left_cats)
        if (g > best_g) {
          best_g <- g
          best_left <- left_cats
        }
      }

      ord <- c(best_left, setdiff(cats, best_left))
    } else {
      score <- vapply(cats, function(k) {
        idx <- x == k
        concentration_index(y[idx, 1], y[idx, 2], wt[idx], parms)
      }, numeric(1))
      mean_score <- vapply(cats, function(k) {
        idx <- x == k
        Wk <- sum(wt[idx])
        if (Wk > 0) sum(wt[idx] * y[idx, 1]) / Wk else mean(y[idx, 1])
      }, numeric(1))
      ord <- cats[order(score, mean_score)]
    }

    list(goodness = goodness_for_order(ord), direction = ord)
  }

  list(init = init_ci, eval = eval_ci, split = split_ci)
}

ci_rpart_formula <- function(outcome = "deadu5_num", rank = "wealth", predictors) {
  stats::as.formula(
    paste0(
      "cbind(", outcome, ", ", rank, ") ~ ",
      paste(predictors, collapse = " + ")
    )
  )
}

fit_ci_rpart_comparison <- function(
    data,
    predictors,
    outcome = "deadu5_num",
    rank = "wealth",
    weights = "sample_weight",
    simple_control = NULL,
    robust_control = NULL,
    robust_loss = "absolute",
    robust_rerank_within_node = TRUE,
    robust_maxcat_enum = 12L) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Install rpart to fit the comparison trees.", call. = FALSE)
  }

  data <- as.data.frame(data)
  needed <- unique(c(outcome, rank, weights, predictors))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  keep <- stats::complete.cases(data[, needed, drop = FALSE])
  data <- data[keep, , drop = FALSE]

  simple_control <- simple_control %||% rpart::rpart.control(
    minbucket = 30L,
    cp = 0.001,
    maxdepth = 6L,
    xval = 10L
  )
  robust_control <- robust_control %||% rpart::rpart.control(
    minsplit = 580L,
    minbucket = 280L,
    cp = 0.001,
    xval = 0L
  )

  formula <- ci_rpart_formula(outcome = outcome, rank = rank, predictors = predictors)
  data$.ci_rpart_weight <- data[[weights]]

  simple_fit <- rpart::rpart(
    formula = formula,
    data = data,
    weights = .ci_rpart_weight,
    method = make_simple_ci_rpart_method(),
    control = simple_control,
    model = TRUE,
    y = TRUE
  )

  robust_fit <- rpart::rpart(
    formula = formula,
    data = data,
    weights = .ci_rpart_weight,
    method = make_ci_rpart_method(
      loss = robust_loss,
      rerank_within_node = robust_rerank_within_node,
      maxcat_enum = robust_maxcat_enum
    ),
    control = robust_control,
    model = TRUE,
    y = TRUE
  )

  list(
    simple_abs = simple_fit,
    robust_abs = robust_fit
  )
}

summarise_ci_rpart_fit <- function(fit, method_name = deparse(substitute(fit))) {
  frame <- fit$frame
  yval <- as.data.frame(fit$frame$yval2)

  if (ncol(yval) == 2L) {
    names(yval) <- c("mean_y", "CI")
    yval$CI_loss <- abs(yval$CI)
  } else {
    names(yval)[seq_len(min(3L, ncol(yval)))] <- c("mean_y", "CI", "CI_loss")[seq_len(min(3L, ncol(yval)))]
  }

  out <- data.frame(
    method = method_name,
    node = rownames(frame),
    var = frame$var,
    n = frame$n,
    weight = frame$wt,
    deviance = frame$dev,
    yval,
    row.names = NULL
  )
  out$is_terminal <- out$var == "<leaf>"
  out
}

summarise_ci_rpart_comparison <- function(fits) {
  out <- do.call(rbind, Map(summarise_ci_rpart_fit, fits, names(fits)))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::as.data.table(out)
  } else {
    out
  }
}

plot_ci_rpart_fit <- function(fit, main = NULL, use_rpart_plot = TRUE, ...) {
  if (isTRUE(use_rpart_plot) && requireNamespace("rpart.plot", quietly = TRUE)) {
    rpart.plot::rpart.plot(fit, main = main, ...)
  } else {
    plot(fit, uniform = TRUE, margin = 0.1, main = main)
    text(fit, use.n = TRUE, cex = 0.8)
  }
}
