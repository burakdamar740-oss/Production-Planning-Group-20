# =============================================================================
# BOX-JENKINS (ARIMA) - SHINY MODULE
#
# Supports AR(p), MA(q), ARMA(p, q) and ARIMA(p, d, q) models, optionally with a
# seasonal part (P, D, Q)[m]. Orders can be chosen manually or searched
# automatically (information criteria via forecast::auto.arima(), or the lowest
# in-sample MAD / MSE / MAPE over a grid of candidate models).
#
# Defines pure helper functions and the Shiny module:
#   arima_ui(id)
#   arima_server(id, selected_series)
#
# `selected_series` must be a reactive that returns a numeric vector.
# Nothing is executed when this file is sourced.
#
# Requires the 'forecast' package (forecast::Arima, auto.arima, forecast).
#
# Error metrics (MAD, MSE, MAPE) are computed on the in-sample one-step-ahead
# fitted values. The first d + D * m periods are excluded because they are used
# up by differencing.
# =============================================================================


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------

# Number formatting that tolerates NA.
arima_format_number <- function(x, digits = 2) {
  if (length(x) == 0 || is.null(x) || is.na(x)) {
    "NA"
  } else {
    format(round(x, digits), nsmall = digits)
  }
}

# TRUE if x is a single whole number between low and high (inclusive).
arima_is_whole_in <- function(x, low, high) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x == floor(x) && x >= low && x <= high
}

# Checks that the 'forecast' package is installed.
arima_require_forecast <- function() {
  if (!requireNamespace("forecast", quietly = TRUE)) {
    stop("The 'forecast' package is required. Install it with install.packages(\"forecast\").")
  }
  invisible(TRUE)
}

# Checks that the series is a usable numeric vector.
arima_validate_series <- function(demand_vec) {
  if (!is.numeric(demand_vec)) {
    stop("The data must be numeric.")
  }
  if (any(!is.finite(demand_vec))) {
    stop("The data cannot contain missing or non-finite values.")
  }
  if (length(demand_vec) < 10) {
    stop("Box-Jenkins models need at least 10 valid observations.")
  }
  invisible(TRUE)
}

# Non-seasonal orders implied by the chosen model type.
#   ar -> (p, 0, 0), ma -> (0, 0, q), arma -> (p, 0, q), arima -> (p, d, q)
arima_effective_orders <- function(model_type, p_value, d_value, q_value) {
  switch(
    model_type,
    ar    = c(p_value, 0, 0),
    ma    = c(0, 0, q_value),
    arma  = c(p_value, 0, q_value),
    arima = c(p_value, d_value, q_value),
    stop("Unknown model type.")
  )
}

# Readable name of a fitted model, built from fit$arma:
# c(p, q, P, Q, m, d, D).
arima_model_label <- function(arma_vec) {
  p_value <- arma_vec[1]
  q_value <- arma_vec[2]
  P_value <- arma_vec[3]
  Q_value <- arma_vec[4]
  m_value <- arma_vec[5]
  d_value <- arma_vec[6]
  D_value <- arma_vec[7]
  
  has_seasonal <- (P_value + Q_value + D_value) > 0
  if (!has_seasonal && d_value == 0) {
    base_label <- if (p_value > 0 && q_value > 0) {
      sprintf("ARMA(%d,%d)", p_value, q_value)
    } else if (p_value > 0) {
      sprintf("AR(%d)", p_value)
    } else if (q_value > 0) {
      sprintf("MA(%d)", q_value)
    } else {
      "White noise"
    }
    return(base_label)
  }
  
  base_label <- sprintf("ARIMA(%d,%d,%d)", p_value, d_value, q_value)
  if (has_seasonal) {
    base_label <- paste0(base_label,
                         sprintf("(%d,%d,%d)[%d]", P_value, D_value, Q_value, m_value))
  }
  base_label
}

# Friendly label for a coefficient name used by arima() (ar1, ma1, sar1, ...).
arima_term_label <- function(term_name) {
  if (grepl("^ar[0-9]+$", term_name)) {
    return(paste0("AR(", sub("^ar", "", term_name), ")"))
  }
  if (grepl("^ma[0-9]+$", term_name)) {
    return(paste0("MA(", sub("^ma", "", term_name), ")"))
  }
  if (grepl("^sar[0-9]+$", term_name)) {
    return(paste0("Seasonal AR(", sub("^sar", "", term_name), ")"))
  }
  if (grepl("^sma[0-9]+$", term_name)) {
    return(paste0("Seasonal MA(", sub("^sma", "", term_name), ")"))
  }
  if (identical(term_name, "intercept")) {
    return("Intercept (mean)")
  }
  if (identical(term_name, "drift")) {
    return("Drift")
  }
  term_name
}


# ----------------------------------------------------------------------------
# Error metrics
# ----------------------------------------------------------------------------

# MAD, MSE and MAPE of the fitted values.
# - The first `burn_in` periods are excluded (used up by differencing).
# - Only periods where both actual and fitted values are available are used.
# - Periods with zero actual demand are excluded from MAPE only
#   (division by zero), and counted in `n_used_mape`.
arima_error_metrics <- function(actual_vec, fitted_vec, burn_in = 0) {
  if (length(actual_vec) != length(fitted_vec)) {
    stop("Actual and fitted vectors must have the same length.")
  }
  
  valid_idx <- !is.na(actual_vec) & !is.na(fitted_vec) & (seq_along(actual_vec) > burn_in)
  if (!any(valid_idx)) {
    return(list(MAD = NA_real_, MSE = NA_real_, MAPE = NA_real_,
                n_used = 0L, n_used_mape = 0L))
  }
  
  actual_valid <- actual_vec[valid_idx]
  fitted_valid <- fitted_vec[valid_idx]
  error_vec    <- actual_valid - fitted_valid
  
  nonzero_idx <- actual_valid != 0
  mape_value <- if (any(nonzero_idx)) {
    mean(abs(error_vec[nonzero_idx]) / abs(actual_valid[nonzero_idx])) * 100
  } else {
    NA_real_
  }
  
  list(
    MAD         = mean(abs(error_vec)),
    MSE         = mean(error_vec^2),
    MAPE        = mape_value,
    n_used      = as.integer(sum(valid_idx)),
    n_used_mape = as.integer(sum(nonzero_idx))
  )
}

# Metrics as a small display table.
arima_metrics_table <- function(metrics) {
  data.frame(
    Metric = c("MAD", "MSE", "MAPE (%)"),
    Value  = c(
      arima_format_number(metrics$MAD, 4),
      arima_format_number(metrics$MSE, 4),
      arima_format_number(metrics$MAPE, 4)
    ),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# Model estimation
# ----------------------------------------------------------------------------

# Fits an ARIMA model with the given orders.
#   order_vec          : c(p, d, q)
#   seasonal_order_vec : c(P, D, Q)
#   season_length      : m (use 1 for a non-seasonal model)
arima_fit_manual <- function(demand_vec, order_vec, seasonal_order_vec = c(0, 0, 0),
                             season_length = 1, include_constant = TRUE) {
  arima_require_forecast()
  arima_validate_series(demand_vec)
  
  ts_obj <- stats::ts(demand_vec, frequency = season_length)
  tryCatch(
    suppressWarnings(
      forecast::Arima(
        ts_obj,
        order            = order_vec,
        seasonal         = list(order = seasonal_order_vec, period = season_length),
        include.constant = include_constant,
        method           = "CSS-ML"
      )
    ),
    error = function(e) {
      stop(paste0("The model could not be estimated: ", conditionMessage(e),
                  " Try different orders."))
    }
  )
}

# Automatic search with forecast::auto.arima() using an information criterion
# ("AICc", "AIC" or "BIC").
arima_fit_auto_ic <- function(demand_vec, season_length, use_seasonal,
                              criterion = "AICc", max_p = 3, max_q = 3) {
  arima_require_forecast()
  arima_validate_series(demand_vec)
  
  ic_value <- switch(criterion, AICc = "aicc", AIC = "aic", BIC = "bic",
                     stop("Unknown information criterion."))
  ts_obj <- stats::ts(demand_vec, frequency = if (use_seasonal) season_length else 1)
  
  tryCatch(
    suppressWarnings(
      forecast::auto.arima(
        ts_obj,
        seasonal  = use_seasonal,
        ic        = ic_value,
        max.p     = max_p,
        max.q     = max_q,
        max.P     = 1,
        max.Q     = 1,
        stepwise  = length(demand_vec) > 100
      )
    ),
    error = function(e) {
      stop(paste0("The automatic search failed: ", conditionMessage(e)))
    }
  )
}

# Automatic search by error minimization: the differencing orders d and D are
# chosen with unit-root tests (ndiffs / nsdiffs), then every combination of
# p = 0..max_p, q = 0..max_q (and P, Q = 0..1 for seasonal models) is estimated.
# The model with the lowest in-sample value of the chosen criterion
# ("MAD", "MSE" or "MAPE") is kept.
# Note: in-sample errors always favour larger models, so information criteria
# are usually the safer choice.
arima_search_by_error <- function(demand_vec, season_length, use_seasonal,
                                  criterion = "MAD", max_p = 3, max_q = 3) {
  arima_require_forecast()
  arima_validate_series(demand_vec)
  if (!(criterion %in% c("MAD", "MSE", "MAPE"))) {
    stop("Criterion must be 'MAD', 'MSE' or 'MAPE'.")
  }
  
  frequency_value <- if (use_seasonal) season_length else 1
  ts_obj <- stats::ts(demand_vec, frequency = frequency_value)
  
  d_value <- tryCatch(min(forecast::ndiffs(ts_obj), 2), error = function(e) 0)
  D_value <- if (use_seasonal) {
    tryCatch(min(forecast::nsdiffs(ts_obj), 1), error = function(e) 0)
  } else {
    0
  }
  burn_in <- d_value + D_value * frequency_value
  
  candidate_grid <- expand.grid(
    p = 0:max_p,
    q = 0:max_q,
    P = if (use_seasonal) 0:1 else 0,
    Q = if (use_seasonal) 0:1 else 0
  )
  
  best_fit    <- NULL
  best_value  <- Inf
  table_rows  <- vector("list", nrow(candidate_grid))
  
  for (row_idx in seq_len(nrow(candidate_grid))) {
    candidate <- candidate_grid[row_idx, ]
    candidate_fit <- tryCatch(
      suppressWarnings(
        forecast::Arima(
          ts_obj,
          order            = c(candidate$p, d_value, candidate$q),
          seasonal         = list(order = c(candidate$P, D_value, candidate$Q),
                                  period = frequency_value),
          include.constant = TRUE,
          method           = "CSS-ML"
        )
      ),
      error = function(e) NULL
    )
    if (is.null(candidate_fit)) {
      next
    }
    
    metrics <- arima_error_metrics(demand_vec, as.numeric(stats::fitted(candidate_fit)), burn_in)
    criterion_value <- metrics[[criterion]]
    if (is.na(criterion_value) || !is.finite(criterion_value)) {
      next
    }
    
    table_rows[[row_idx]] <- data.frame(
      Model = arima_model_label(candidate_fit$arma),
      value = criterion_value,
      AIC   = candidate_fit$aic,
      stringsAsFactors = FALSE
    )
    if (criterion_value < best_value) {
      best_value <- criterion_value
      best_fit   <- candidate_fit
    }
  }
  
  if (is.null(best_fit)) {
    stop("None of the candidate models could be estimated. Try a different series or settings.")
  }
  
  search_table <- do.call(rbind, table_rows)
  search_table <- search_table[order(search_table$value), , drop = FALSE]
  search_table <- utils::head(search_table, 10)
  names(search_table)[2] <- criterion
  rownames(search_table) <- NULL
  
  list(fit = best_fit, table = search_table, d = d_value, D = D_value)
}


# ----------------------------------------------------------------------------
# Output helpers
# ----------------------------------------------------------------------------

# Coefficient table: estimates, standard errors, z values and p-values.
# Returns NULL if the model has no estimated coefficients.
arima_coef_table <- function(fit) {
  coefficient_vec <- stats::coef(fit)
  if (length(coefficient_vec) == 0) {
    return(NULL)
  }
  
  std_error <- tryCatch(sqrt(diag(fit$var.coef)), error = function(e) NULL)
  if (is.null(std_error) || length(std_error) != length(coefficient_vec)) {
    std_error <- rep(NA_real_, length(coefficient_vec))
  }
  z_value <- coefficient_vec / std_error
  p_value <- 2 * stats::pnorm(-abs(z_value))
  
  data.frame(
    "Term"       = vapply(names(coefficient_vec), arima_term_label, character(1)),
    "Estimate"   = as.numeric(coefficient_vec),
    "Std. Error" = as.numeric(std_error),
    "z value"    = as.numeric(z_value),
    "p-value"    = as.numeric(p_value),
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
}

# Ljung-Box test on the residuals. `fit_df` is the number of estimated ARMA
# parameters. Returns NULL if the test cannot be run.
arima_ljung_box <- function(residual_vec, fit_df) {
  residual_vec <- residual_vec[is.finite(residual_vec)]
  lag_value <- min(max(10, fit_df + 3), floor(length(residual_vec) / 2))
  if (length(residual_vec) < 5 || lag_value <= fit_df) {
    return(NULL)
  }
  test_result <- stats::Box.test(residual_vec, lag = lag_value,
                                 type = "Ljung-Box", fitdf = fit_df)
  list(
    lag       = lag_value,
    statistic = unname(test_result$statistic),
    p_value   = test_result$p.value
  )
}


# ----------------------------------------------------------------------------
# Plots (base R graphics)
# ----------------------------------------------------------------------------

# Sample ACF and PACF of the original series and, when differencing is used,
# of the differenced series.
arima_plot_acf_pacf <- function(demand_vec, d_value = 0, D_value = 0, season_length = 1) {
  has_differencing <- d_value > 0 || D_value > 0
  
  differenced_vec <- demand_vec
  if (D_value > 0 && season_length > 1) {
    differenced_vec <- diff(differenced_vec, lag = season_length, differences = D_value)
  }
  if (d_value > 0) {
    differenced_vec <- diff(differenced_vec, differences = d_value)
  }
  
  lag_limit <- function(x) {
    max(1, min(floor(length(x) / 2), max(20, 3 * season_length)))
  }
  
  old_par <- graphics::par(mfrow = c(if (has_differencing) 2 else 1, 2),
                           mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  
  stats::acf(demand_vec, lag.max = lag_limit(demand_vec),
             main = "ACF of the original series")
  stats::pacf(demand_vec, lag.max = lag_limit(demand_vec),
              main = "PACF of the original series")
  
  if (has_differencing) {
    differencing_text <- sprintf("d = %d, D = %d", d_value, D_value)
    if (length(differenced_vec) > 5) {
      stats::acf(differenced_vec, lag.max = lag_limit(differenced_vec),
                 main = paste0("ACF of the differenced series (", differencing_text, ")"))
      stats::pacf(differenced_vec, lag.max = lag_limit(differenced_vec),
                  main = paste0("PACF of the differenced series (", differencing_text, ")"))
    } else {
      for (panel_idx in 1:2) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Too few observations after differencing.")
      }
    }
  }
}

# Historical data, in-sample fitted values and out-of-sample forecasts with
# 80% and 95% prediction intervals.
arima_plot_forecast <- function(demand_vec, fitted_vec, forecast_mean,
                                lower_80, upper_80, lower_95, upper_95,
                                burn_in = 0, plot_title = "ARIMA Forecast") {
  data_length <- length(demand_vec)
  horizon     <- length(forecast_mean)
  future_idx  <- data_length + seq_len(horizon)
  
  fitted_plot <- fitted_vec
  if (burn_in > 0) {
    fitted_plot[seq_len(burn_in)] <- NA
  }
  
  y_range <- range(c(demand_vec, fitted_plot, lower_95, upper_95), na.rm = TRUE)
  
  plot(seq_len(data_length), demand_vec,
       type = "l", col = "blue", lwd = 2,
       xlim = c(1, data_length + horizon), ylim = y_range,
       xlab = "Period", ylab = "Value", main = plot_title)
  
  graphics::polygon(c(future_idx, rev(future_idx)), c(lower_95, rev(upper_95)),
                    col = grDevices::adjustcolor("darkgreen", alpha.f = 0.15), border = NA)
  graphics::polygon(c(future_idx, rev(future_idx)), c(lower_80, rev(upper_80)),
                    col = grDevices::adjustcolor("darkgreen", alpha.f = 0.30), border = NA)
  
  lines(seq_len(data_length), fitted_plot, col = "red", lty = 2, lwd = 1.5)
  lines(c(data_length, future_idx), c(demand_vec[data_length], forecast_mean),
        col = "darkgreen", lwd = 2)
  points(future_idx, forecast_mean, pch = 18, cex = 1.4, col = "darkgreen")
  
  legend("topleft",
         legend = c("Actual", "Fitted (in-sample)", "Forecast",
                    "80% interval", "95% interval"),
         col = c("blue", "red", "darkgreen",
                 grDevices::adjustcolor("darkgreen", alpha.f = 0.30),
                 grDevices::adjustcolor("darkgreen", alpha.f = 0.15)),
         lty = c(1, 2, 1, NA, NA), lwd = c(2, 1.5, 2, NA, NA),
         pch = c(NA, NA, 18, 15, 15), pt.cex = c(1, 1, 1.4, 2, 2), bty = "n")
}

# Residuals over time and their ACF (diagnostic checking).
arima_plot_residuals <- function(residual_vec, burn_in = 0) {
  residual_vec <- residual_vec[seq_along(residual_vec) > burn_in]
  
  old_par <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  
  plot(residual_vec, type = "h", col = "red",
       xlab = "Observation (after differencing)", ylab = "Residual",
       main = "Residuals")
  graphics::abline(h = 0, lty = 2, col = "gray40")
  stats::acf(residual_vec, lag.max = max(1, min(20, floor(length(residual_vec) / 2))),
             main = "ACF of the residuals")
}


# ----------------------------------------------------------------------------
# Shiny module: UI
# ----------------------------------------------------------------------------

# Small summary card (used by the server).
arima_make_card <- function(label, value) {
  shiny::column(
    4,
    shiny::div(
      class = "well well-sm text-center",
      shiny::tags$small(label),
      shiny::h3(value)
    )
  )
}

arima_ui <- function(id) {
  ns <- shiny::NS(id)
  
  js_equals <- function(input_name, value) {
    sprintf("input['%s'] == '%s'", ns(input_name), value)
  }
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::radioButtons(
        ns("mode"), "How should the model be chosen?",
        choices = c(
          "Manual (choose the orders)"          = "manual",
          "Auto (search for the best model)"    = "auto"
        ),
        selected = "manual"
      ),
      
      # --- Manual orders ----------------------------------------------------------
      shiny::conditionalPanel(
        condition = js_equals("mode", "manual"),
        shiny::selectInput(
          ns("model_type"), "Model type",
          choices = c(
            "AR(p)"         = "ar",
            "MA(q)"         = "ma",
            "ARMA(p, q)"    = "arma",
            "ARIMA(p, d, q)" = "arima"
          ),
          selected = "arima"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] != 'ma'", ns("model_type")),
          shiny::numericInput(ns("p_order"), "AR order p", value = 1, min = 0, max = 5, step = 1)
        ),
        shiny::conditionalPanel(
          condition = js_equals("model_type", "arima"),
          shiny::numericInput(ns("d_order"), "Differencing order d", value = 1, min = 0, max = 2, step = 1)
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] != 'ar'", ns("model_type")),
          shiny::numericInput(ns("q_order"), "MA order q", value = 1, min = 0, max = 5, step = 1)
        ),
        shiny::checkboxInput(
          ns("include_constant"),
          "Include a constant (mean if d = 0, drift if d = 1)",
          value = TRUE
        )
      ),
      
      # --- Automatic search -----------------------------------------------------------
      shiny::conditionalPanel(
        condition = js_equals("mode", "auto"),
        shiny::selectInput(
          ns("criterion"), "Select the best model by",
          choices = c(
            "AICc (forecast::auto.arima)" = "AICc",
            "AIC (forecast::auto.arima)"  = "AIC",
            "BIC (forecast::auto.arima)"  = "BIC",
            "Lowest in-sample MAD"        = "MAD",
            "Lowest in-sample MSE"        = "MSE",
            "Lowest in-sample MAPE"       = "MAPE"
          ),
          selected = "AICc"
        ),
        shiny::numericInput(ns("max_p"), "Maximum p to try", value = 3, min = 0, max = 5, step = 1),
        shiny::numericInput(ns("max_q"), "Maximum q to try", value = 3, min = 0, max = 5, step = 1),
        shiny::helpText(
          "The differencing orders d and D are chosen automatically with unit-root tests. ",
          "Information criteria (AICc, AIC, BIC) penalize model size; in-sample errors ",
          "do not, so they tend to favour larger models. Large search ranges can take a while."
        )
      ),
      
      shiny::hr(),
      
      # --- Seasonal part -----------------------------------------------------------------
      shiny::checkboxInput(ns("use_seasonal"), "Include a seasonal component", value = FALSE),
      shiny::conditionalPanel(
        condition = sprintf("input['%s']", ns("use_seasonal")),
        shiny::numericInput(
          ns("season_length"), "Seasonal period m (4 = quarterly, 12 = monthly, 7 = weekly)",
          value = 4, min = 2, step = 1
        ),
        shiny::conditionalPanel(
          condition = js_equals("mode", "manual"),
          shiny::numericInput(ns("sp_order"), "Seasonal AR order P", value = 0, min = 0, max = 2, step = 1),
          shiny::numericInput(ns("sd_order"), "Seasonal differencing order D", value = 1, min = 0, max = 1, step = 1),
          shiny::numericInput(ns("sq_order"), "Seasonal MA order Q", value = 1, min = 0, max = 2, step = 1)
        )
      ),
      
      shiny::hr(),
      
      shiny::numericInput(
        ns("horizon"), "Forecast horizon (periods ahead)",
        value = 6, min = 1, max = 120, step = 1
      ),
      shiny::helpText(
        "Error metrics use the in-sample fitted values and exclude the first ",
        "d + D * m periods that are consumed by differencing. ",
        "Results update automatically; the data comes from the upload panel at the top."
      )
    ),
    
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Summary",
          shiny::uiOutput(ns("summary")),
          shiny::h4("Error Metrics (in-sample fitted values)"),
          shiny::tableOutput(ns("metrics_table")),
          shiny::h4("Parameter Estimates"),
          shiny::uiOutput(ns("coef_ui")),
          shiny::h4("Model Output"),
          shiny::verbatimTextOutput(ns("model_text")),
          shiny::uiOutput(ns("search_ui"))
        ),
        shiny::tabPanel(
          "ACF / PACF",
          shiny::plotOutput(ns("acf_plot"), height = "560px"),
          shiny::helpText(
            "Use these plots to identify the orders: a PACF that cuts off after lag p suggests ",
            "an AR(p) part, an ACF that cuts off after lag q suggests an MA(q) part, and a ",
            "slowly decaying ACF suggests that differencing is needed. When d or D is larger ",
            "than 0, the plots of the differenced series are shown as well."
          )
        ),
        shiny::tabPanel(
          "Forecast",
          shiny::plotOutput(ns("forecast_plot"), height = "380px"),
          shiny::h4("Forecasts"),
          shiny::tableOutput(ns("forecast_table"))
        ),
        shiny::tabPanel(
          "Diagnostics",
          shiny::uiOutput(ns("ljung_ui")),
          shiny::plotOutput(ns("residual_plot"), height = "340px")
        ),
        shiny::tabPanel(
          "Detailed Results",
          shiny::tableOutput(ns("detail_table"))
        )
      )
    )
  )
}


# ----------------------------------------------------------------------------
# Shiny module: server
# ----------------------------------------------------------------------------

# selected_series: a reactive returning a numeric vector (supplied by app.R).
arima_server <- function(id, selected_series) {
  stopifnot(shiny::is.reactive(selected_series))
  
  shiny::moduleServer(id, function(input, output, session) {
    
    # Turns any problem (validation message, input error, estimation error) into
    # a message on the screen instead of crashing the session.
    safe_run <- function(expr_fn) {
      tryCatch(
        expr_fn(),
        error = function(e) {
          error_msg <- conditionMessage(e)
          if (!nzchar(error_msg)) {
            error_msg <- "Upload a data file and select a numeric column in the panel above."
          }
          list(ok = FALSE, message = error_msg)
        }
      )
    }
    
    warning_box <- function(outcome) {
      shiny::div(class = "alert alert-warning", outcome$message)
    }
    
    # ---- Data and common settings --------------------------------------------------
    prepared <- shiny::reactive({
      raw_series <- selected_series()
      shiny::req(!is.null(raw_series), length(raw_series) > 0)
      shiny::validate(shiny::need(is.numeric(raw_series), "The data must be numeric."))
      demand_vec <- raw_series[is.finite(raw_series)]
      shiny::validate(shiny::need(
        length(demand_vec) >= 10,
        "Box-Jenkins models need at least 10 valid observations."
      ))
      
      horizon <- input$horizon
      shiny::validate(shiny::need(
        arima_is_whole_in(horizon, 1, 120),
        "The forecast horizon must be a whole number between 1 and 120."
      ))
      
      use_seasonal  <- isTRUE(input$use_seasonal)
      season_length <- 1
      if (use_seasonal) {
        season_length <- input$season_length
        shiny::validate(shiny::need(
          is.numeric(season_length) && length(season_length) == 1 && !is.na(season_length) &&
            season_length == floor(season_length) && season_length >= 2,
          "The seasonal period must be a whole number of at least 2."
        ))
        shiny::validate(shiny::need(
          length(demand_vec) >= 2 * season_length,
          sprintf("A seasonal model needs at least 2 full seasons: %d observations for m = %d, but the data has only %d.",
                  2 * season_length, season_length, length(demand_vec))
        ))
      }
      
      list(
        y             = demand_vec,
        horizon       = horizon,
        use_seasonal  = use_seasonal,
        season_length = season_length
      )
    })
    
    # ---- Model estimation and all derived results ---------------------------------------
    arima_result <- shiny::reactive({
      safe_run(function() {
        prep <- prepared()
        y    <- prep$y
        
        search_table   <- NULL
        selection_text <- NULL
        
        if (identical(input$mode, "auto")) {
          max_p     <- input$max_p
          max_q     <- input$max_q
          criterion <- input$criterion
          shiny::validate(
            shiny::need(arima_is_whole_in(max_p, 0, 5), "Maximum p must be a whole number between 0 and 5."),
            shiny::need(arima_is_whole_in(max_q, 0, 5), "Maximum q must be a whole number between 0 and 5."),
            shiny::need(criterion %in% c("AICc", "AIC", "BIC", "MAD", "MSE", "MAPE"),
                        "Choose a model selection criterion.")
          )
          
          if (criterion %in% c("AICc", "AIC", "BIC")) {
            fit <- arima_fit_auto_ic(y, prep$season_length, prep$use_seasonal,
                                     criterion, max_p, max_q)
            selection_text <- paste0("Auto mode: model selected by ", criterion,
                                     " (forecast::auto.arima)")
          } else {
            search <- arima_search_by_error(y, prep$season_length, prep$use_seasonal,
                                            criterion, max_p, max_q)
            fit          <- search$fit
            search_table <- search$table
            selection_text <- paste0("Auto mode: model with the lowest in-sample ", criterion,
                                     " (d = ", search$d, ", D = ", search$D, " from unit-root tests)")
          }
        } else {
          model_type <- input$model_type
          p_order    <- input$p_order
          d_order    <- input$d_order
          q_order    <- input$q_order
          
          needs_p <- model_type %in% c("ar", "arma", "arima")
          needs_q <- model_type %in% c("ma", "arma", "arima")
          shiny::validate(
            shiny::need(!needs_p || arima_is_whole_in(p_order, 0, 5),
                        "The AR order p must be a whole number between 0 and 5."),
            shiny::need(model_type != "arima" || arima_is_whole_in(d_order, 0, 2),
                        "The differencing order d must be a whole number between 0 and 2."),
            shiny::need(!needs_q || arima_is_whole_in(q_order, 0, 5),
                        "The MA order q must be a whole number between 0 and 5."),
            shiny::need(model_type != "ar" || p_order >= 1, "An AR model needs p of at least 1."),
            shiny::need(model_type != "ma" || q_order >= 1, "An MA model needs q of at least 1."),
            shiny::need(model_type != "arma" || (p_order >= 1 && q_order >= 1),
                        "An ARMA model needs both p and q of at least 1.")
          )
          
          order_vec <- arima_effective_orders(model_type, p_order, d_order, q_order)
          
          seasonal_order_vec <- c(0, 0, 0)
          if (prep$use_seasonal) {
            sp_order <- input$sp_order
            sd_order <- input$sd_order
            sq_order <- input$sq_order
            shiny::validate(
              shiny::need(arima_is_whole_in(sp_order, 0, 2),
                          "The seasonal AR order P must be a whole number between 0 and 2."),
              shiny::need(arima_is_whole_in(sd_order, 0, 1),
                          "The seasonal differencing order D must be 0 or 1."),
              shiny::need(arima_is_whole_in(sq_order, 0, 2),
                          "The seasonal MA order Q must be a whole number between 0 and 2.")
            )
            seasonal_order_vec <- c(sp_order, sd_order, sq_order)
          }
          
          fit <- arima_fit_manual(y, order_vec, seasonal_order_vec,
                                  prep$season_length, isTRUE(input$include_constant))
          selection_text <- "Manual mode: orders chosen by the user"
        }
        
        arma_vec <- fit$arma
        orders <- list(p = arma_vec[1], q = arma_vec[2], P = arma_vec[3], Q = arma_vec[4],
                       m = arma_vec[5], d = arma_vec[6], D = arma_vec[7])
        burn_in <- orders$d + orders$D * orders$m
        
        fitted_vec   <- as.numeric(stats::fitted(fit))
        residual_vec <- as.numeric(stats::residuals(fit))
        metrics      <- arima_error_metrics(y, fitted_vec, burn_in)
        
        forecast_obj <- forecast::forecast(fit, h = prep$horizon, level = c(80, 95))
        
        list(
          ok             = TRUE,
          prep           = prep,
          fit            = fit,
          label          = arima_model_label(arma_vec),
          orders         = orders,
          burn_in        = burn_in,
          selection_text = selection_text,
          search_table   = search_table,
          fitted         = fitted_vec,
          residuals      = residual_vec,
          metrics        = metrics,
          coef_table     = arima_coef_table(fit),
          ljung          = arima_ljung_box(
            residual_vec[seq_along(residual_vec) > burn_in],
            sum(arma_vec[1:4])
          ),
          forecast = list(
            mean     = as.numeric(forecast_obj$mean),
            lower_80 = as.numeric(forecast_obj$lower[, 1]),
            upper_80 = as.numeric(forecast_obj$upper[, 1]),
            lower_95 = as.numeric(forecast_obj$lower[, 2]),
            upper_95 = as.numeric(forecast_obj$upper[, 2])
          )
        )
      })
    })
    
    # ---- Summary -----------------------------------------------------------------------------
    output$summary <- shiny::renderUI({
      outcome <- arima_result()
      if (!outcome$ok) {
        return(warning_box(outcome))
      }
      fit <- outcome$fit
      
      ljung_card_value <- if (is.null(outcome$ljung)) {
        "n/a"
      } else {
        arima_format_number(outcome$ljung$p_value, 3)
      }
      
      shiny::tagList(
        shiny::h4(outcome$selection_text),
        shiny::fluidRow(
          arima_make_card("Model", outcome$label),
          arima_make_card("AIC", arima_format_number(fit$aic)),
          arima_make_card("BIC", arima_format_number(fit$bic))
        ),
        shiny::fluidRow(
          arima_make_card("AICc", arima_format_number(fit$aicc)),
          arima_make_card("Log-likelihood", arima_format_number(fit$loglik)),
          arima_make_card("Residual variance (sigma^2)", arima_format_number(fit$sigma2, 4))
        ),
        shiny::fluidRow(
          arima_make_card("Ljung-Box p-value", ljung_card_value),
          arima_make_card("Periods used for errors", outcome$metrics$n_used),
          arima_make_card("Differencing periods excluded", outcome$burn_in)
        ),
        if (!is.na(outcome$metrics$n_used_mape) &&
            outcome$metrics$n_used_mape < outcome$metrics$n_used) {
          shiny::tags$p(shiny::em(
            "Note: periods with zero actual demand were excluded from MAPE."
          ))
        }
      )
    })
    
    output$metrics_table <- shiny::renderTable({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      arima_metrics_table(outcome$metrics)
    })
    
    output$coef_ui <- shiny::renderUI({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      if (is.null(outcome$coef_table)) {
        return(shiny::tags$p(shiny::em(
          "This model has no estimated coefficients (white noise without a constant)."
        )))
      }
      shiny::tableOutput(session$ns("coef_table"))
    })
    
    output$coef_table <- shiny::renderTable({
      outcome <- arima_result()
      shiny::req(outcome$ok, !is.null(outcome$coef_table))
      outcome$coef_table
    }, digits = 4, na = "-")
    
    output$model_text <- shiny::renderPrint({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      print(outcome$fit)
    })
    
    output$search_ui <- shiny::renderUI({
      outcome <- arima_result()
      if (!outcome$ok || is.null(outcome$search_table)) {
        return(NULL)
      }
      shiny::tagList(
        shiny::h4("Best Candidate Models"),
        shiny::tableOutput(session$ns("search_table"))
      )
    })
    
    output$search_table <- shiny::renderTable({
      outcome <- arima_result()
      shiny::req(outcome$ok, !is.null(outcome$search_table))
      outcome$search_table
    }, digits = 4)
    
    # ---- ACF / PACF -------------------------------------------------------------------------------
    output$acf_plot <- shiny::renderPlot({
      prep    <- prepared()
      outcome <- arima_result()
      d_value <- if (outcome$ok) outcome$orders$d else 0
      D_value <- if (outcome$ok) outcome$orders$D else 0
      arima_plot_acf_pacf(prep$y, d_value, D_value, prep$season_length)
    })
    
    # ---- Forecast ---------------------------------------------------------------------------------------
    output$forecast_plot <- shiny::renderPlot({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      fc <- outcome$forecast
      arima_plot_forecast(
        outcome$prep$y, outcome$fitted, fc$mean,
        fc$lower_80, fc$upper_80, fc$lower_95, fc$upper_95,
        outcome$burn_in,
        paste0(outcome$label, ": Actual, Fitted and Forecast")
      )
    })
    
    output$forecast_table <- shiny::renderTable({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      fc <- outcome$forecast
      data.frame(
        "Period"       = length(outcome$prep$y) + seq_along(fc$mean),
        "Forecast"     = fc$mean,
        "Lower 80%"    = fc$lower_80,
        "Upper 80%"    = fc$upper_80,
        "Lower 95%"    = fc$lower_95,
        "Upper 95%"    = fc$upper_95,
        check.names    = FALSE
      )
    }, digits = 2)
    
    # ---- Diagnostics --------------------------------------------------------------------------------------
    output$ljung_ui <- shiny::renderUI({
      outcome <- arima_result()
      if (!outcome$ok) {
        return(warning_box(outcome))
      }
      if (is.null(outcome$ljung)) {
        return(shiny::tags$p("The Ljung-Box test could not be run (too few residuals)."))
      }
      ljung <- outcome$ljung
      verdict <- if (ljung$p_value >= 0.05) {
        "No significant autocorrelation is left in the residuals (the model looks adequate)."
      } else {
        "The residuals still show significant autocorrelation (consider different orders)."
      }
      shiny::tagList(
        shiny::tags$p(
          shiny::strong("Ljung-Box test: "),
          sprintf("Q = %s, lag = %d, p-value = %s",
                  arima_format_number(ljung$statistic, 3), ljung$lag,
                  arima_format_number(ljung$p_value, 4))
        ),
        shiny::tags$p(verdict)
      )
    })
    
    output$residual_plot <- shiny::renderPlot({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      arima_plot_residuals(outcome$residuals, outcome$burn_in)
    })
    
    # ---- Detailed results ------------------------------------------------------------------------------------
    output$detail_table <- shiny::renderTable({
      outcome <- arima_result()
      shiny::req(outcome$ok)
      y          <- outcome$prep$y
      fitted_vec <- outcome$fitted
      if (outcome$burn_in > 0) {
        fitted_vec[seq_len(outcome$burn_in)] <- NA
      }
      data.frame(
        "Period" = seq_along(y),
        "Actual" = y,
        "Fitted" = fitted_vec,
        "Error"  = y - fitted_vec,
        check.names = FALSE
      )
    }, na = "-", digits = 2)
  })
}